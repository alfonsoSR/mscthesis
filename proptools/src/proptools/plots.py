import yaml
from pathlib import Path
from dataclasses import dataclass


class PlotSettings:

    def __init__(self, config_file: Path) -> None:

        config = yaml.safe_load(config_file.open())["Plotting"]

        self.ephemerides = config["ephemerides"]
        self.show = config["show"]
        self.save = config["save"]

        return None
