from tastro.config.general import CaseSetup
from ..config.environment.stations import StationSetup
from tudatpy.dynamics.environment_setup import (
    ground_station as tgss,
    convert_ground_station_state_between_itrf_frames,
)
from tudatpy.astro import element_conversion as telc
import numpy as np
from tudatpy.astro import time_representation as ttime
from ..core import SettingsGenerator
from ..io.observations.doppler import get_ground_station_reference_state_itrf


class StationSettings(SettingsGenerator[StationSetup]):

    def __motion_settings(
        self, cstate_itrf2014: np.ndarray
    ) -> list[tgss.GroundStationMotionSettings]:

        # Initialize empty container
        motion_settings: list[tgss.GroundStationMotionSettings] = []

        # Add linear motion settings
        if self.local.coordinates.linear_motion:
            motion_settings.append(
                tgss.LinearGroundStationMotionSettings(
                    linear_velocity=cstate_itrf2014[3:],
                    reference_epoch=self.local.coordinates.reference_epoch,
                )
            )

        # Add shape-deformation motion settings
        if self.local.coordinates.body_deformation:
            raise NotImplementedError(
                "Displacements based on body deformation not implemented"
            )

        return motion_settings

    def station_settings(self) -> tgss.GroundStationSettings:

        # Get reference state for station
        cstate = get_ground_station_reference_state_itrf(
            station=self.name, source=self.local.coordinates.available_position
        )

        # Transform reference state to ITRF2014
        cstate_itrf2014 = convert_ground_station_state_between_itrf_frames(
            ground_station_state=cstate,
            epoch=self.local.coordinates.reference_epoch,
            base_frame=self.local.coordinates.itrf_version,
            target_frame="ITRF2014",
        )

        # Define motion settings
        motion_settings = self.__motion_settings(cstate_itrf2014)

        # Return ground station settings
        return tgss.basic_station(
            station_name=self.name,
            station_nominal_position=cstate_itrf2014[:3],
            station_position_element_type=telc.PositionElementTypes.cartesian_position_type,
            station_motion_settings=motion_settings,
        )
