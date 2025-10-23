from ..config import CaseSetup
from tudatpy.astro import time_representation as ttime


def get_propagation_start_epoch_from_config(config: "CaseSetup") -> ttime.Time:

    match config.propagation.integrator.general.starting_point:

        case "start":
            return config.time.initial_epoch

        case "end":
            return config.time.final_epoch

        case "middle":

            duration = config.time.final_epoch - config.time.initial_epoch
            return config.time.initial_epoch + duration / 2.0

        case "custom":

            custom_epoch = (
                config.propagation.integrator.general.custom_start_epoch
            )
            if custom_epoch is None:
                raise ValueError("Custom propagation epoch not specified")

            return custom_epoch

        case _:
            raise NotImplementedError("Invalid propagation start option")
