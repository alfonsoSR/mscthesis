from dataclasses import dataclass
import typing
import types
from tudatpy.astro import time_representation as ttime
from tudatpy.dynamics.propagation_setup import (
    integrator as tigrs,
    propagator as tprops,
)
from tudatpy.estimation.observable_models_setup import (
    links as tlinks,
    light_time_corrections as tlight,
)
from tudatpy.estimation.observations_setup import ancillary_settings as tanc
import numpy as np
from pathlib import Path


ENUMERATIONS = (
    tigrs.CoefficientSets,
    tprops.TranslationalPropagatorType,
    tlinks.LinkEndType,
    tanc.FrequencyBands,
    tlight.LightTimeFailureHandling,
)


@dataclass
class SetupBase:

    _raw: dict[str, typing.Any]
    present: bool

    @classmethod
    def from_raw(cls, raw_config: dict[str, typing.Any]) -> "typing.Self":

        present = True
        if "present" in raw_config:
            present = raw_config["present"]

        return cls(raw_config, present)

    def __getattribute__(self, name: str) -> typing.Any:

        if (name not in SetupBase.__annotations__) and (
            name not in dir(SetupBase)
        ):
            raise AttributeError(f"Invalid argument {name}")

        return super().__getattribute__(name)

    def __getattr__(self, name: str) -> typing.Any:

        if name not in type(self).__annotations__:
            raise AttributeError(
                f"Invalid attribute {name} for {self.__class__.__name__}"
            )

        _raw = super().__getattribute__("_raw")
        if name not in _raw:
            raise AttributeError(f"Invalid attribute {name}")

        # Get type
        __type = type(self).__annotations__[name]
        if isinstance(__type, types.UnionType):

            valid = (len(__type.__args__) == 2) and (
                (__type.__args__[-1] is types.NoneType)
            )
            if not valid:
                raise ValueError(f"Invalid union type: {__type}: {name}")

            if _raw[name] is None:
                return None

            __type = __type.__args__[0]

        if _raw[name] is None:
            raise ValueError(f"{name} not set for {self.__class__.__name__}")

        if __type is ttime.Time and isinstance(_raw[name], str):
            return ttime.DateTime.from_iso_string(
                _raw[name]
            ).to_epoch_time_object()

        if __type in ENUMERATIONS and isinstance(_raw[name], str):
            return getattr(__type, _raw[name])

        if __type is np.ndarray and isinstance(_raw[name], list):
            return np.array(_raw[name])

        if isinstance(__type, types.GenericAlias):
            if isinstance(_raw[name], list):

                __item_type = __type.__args__[0]
                if hasattr(__item_type, "from_raw"):
                    return [
                        getattr(__item_type, "from_raw")(item)
                        for item in _raw[name]
                    ]
                else:
                    if __item_type == ttime.Time and isinstance(
                        _raw[name][0], str
                    ):
                        return [
                            ttime.DateTime.from_iso_string(
                                item
                            ).to_epoch_time_object()
                            for item in _raw[name]
                        ]
                    else:
                        return [__item_type(item) for item in _raw[name]]

            else:
                raise NotImplementedError(
                    f"Not implemented generic alias: {__type}"
                )

        return __type(_raw[name])


@dataclass
class SetupCollectionBase:

    @classmethod
    def from_raw(
        cls, raw_settings: dict[str, dict[str, typing.Any]]
    ) -> "typing.Self":

        arguments = {}
        for raw_name, raw_val in raw_settings.items():

            argname = raw_name.lower()
            if argname not in cls.__annotations__:
                raise AttributeError(
                    f"Invalid attribute {argname} for {cls.__name__}"
                )
            argtype = cls.__annotations__[argname]

            if isinstance(argtype, types.GenericAlias):

                # Identification flags
                is_dict = (len(argtype.__args__) == 2) and (
                    argtype.__args__[0] is str
                )
                is_list = len(argtype.__args__) == 1
                if is_dict:
                    arguments[argname] = {
                        key: getattr(argtype.__args__[-1], "from_raw")(val)
                        for key, val in raw_val.items()
                    }
                elif is_list:
                    arguments[argname] = [
                        getattr(argtype.__args__[0], "from_raw")(val)
                        for val in raw_val
                    ]
                else:
                    raise NotImplementedError(
                        f"Invalid generic alias: {argtype}"
                    )
            else:

                if hasattr(argtype, "from_raw"):
                    arguments[argname] = getattr(argtype, "from_raw")(raw_val)
                else:
                    arguments[argname] = argtype(raw_val)

        return cls(**arguments)
