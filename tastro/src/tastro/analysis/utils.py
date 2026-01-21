from ..config import CaseSetup
from tudatpy.astro import time_representation as ttime
from pathlib import Path
import yaml


def get_propagaton_start_epoch_from_config_file(
    config_file: Path,
) -> ttime.Time:

    # Load raw contents of configuration file
    with config_file.open("rb") as buffer:
        raw_config = yaml.safe_load(buffer)

    # Extract information about the starting point of the propagation
    integrator_config = raw_config["Propagation"]["integrator"]["general"]
    starting_point = integrator_config["starting_point"]
    initial_epoch = ttime.DateTime.from_iso_string(
        raw_config["Time"]["initial_epoch"]
    ).to_epoch_time_object()
    final_epoch = ttime.DateTime.from_iso_string(
        raw_config["Time"]["final_epoch"]
    ).to_epoch_time_object()
    del raw_config

    # Find propagation start epoch based on configuration
    match starting_point:

        case "start":
            return initial_epoch
        case "end":
            return final_epoch
        case "middle":
            duration = final_epoch - initial_epoch
            return initial_epoch + (duration / 2.0)
        case "custom":
            if integrator_config["custom_start_epoch"] is None:
                raise ValueError("Custom propagation epoch not specified")
            return ttime.DateTime.from_iso_string(
                integrator_config["custom_start_epoch"]
            ).to_epoch_time_object()
        case _:
            raise NotImplementedError("Invalid propagation start option")


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
