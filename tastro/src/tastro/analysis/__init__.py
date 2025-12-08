"""
This module contains the source code, and command line interface, of utilities
to perform analyses with tastro. Functionality includes running propagations,
estimations, and pre-fit analyses, as well as generating plots with their
outcomes
"""

from .utils import get_propagation_start_epoch_from_config

__all__ = ["get_propagation_start_epoch_from_config"]
