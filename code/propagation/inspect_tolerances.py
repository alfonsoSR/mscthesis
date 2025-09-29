import pickle
from nastro import types as nt, graphics as ng
import numpy as np
import argparse
from pathlib import Path
from proptools import config as pcon
import numpy as np
import pickle

parser = argparse.ArgumentParser()
parser.add_argument("config_file", help="Path to configuration file")


if __name__ == "__main__":

    # Load configuration
    config_path = Path(parser.parse_args().config_file).resolve()
    if not config_path.exists():
        raise FileNotFoundError("Configuration file does not exist")
    config = pcon.PropSettings(config_path)

    # Load data
    with (config_path.parent / "tolerance.pkl").open("rb") as buffer:
        data = pickle.load(buffer)

    epochs = data["epochs"]
    error = nt.CartesianState[nt.Double](*data["error"])
    orbit = nt.CartesianState[nt.Double](*data["orbit"])

    with ng.PlotCartesianState() as fig:
        fig.add_state(epochs, error)
