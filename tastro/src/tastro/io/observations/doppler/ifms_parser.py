from pathlib import Path
from ....logging import log
import traceback
from tudatpy import data as tdata
from tudatpy.astro import time_representation as ttime


class IFMSLoader:

    def __init__(self, ifms_file: Path) -> None:

        # Ensure IFMS file exists
        if not ifms_file.exists():
            log.fatal(f"Invalid IFMS file: {ifms_file}")
            log.fatal(traceback.extract_stack()[-2])
            exit(1)

        # Load raw data map
        self.data = tdata.read_ifms_file(
            str(ifms_file), apply_tropospheric_correction=False
        ).raw_datamap

        # Extract metadata from file name
        self.metadata = self._get_metadata_from_filename(ifms_file)

        return None

    @property
    def station_id(self) -> int:

        return int(self.metadata["station_id"])

    @staticmethod
    def _get_metadata_from_filename(source_file: Path) -> dict[str, str]:

        # Extract metadata from file name
        filename_components: dict[str, str] = {
            "mission_id": source_file.stem[0],
            "station_id": source_file.stem[1:3],
            "data_source_id": source_file.stem[3:7],
            "data_level": source_file.stem[7:10],
            "data_type": source_file.stem[11:14],
            "ref_epoch": source_file.stem[15:24],
            "version": source_file.stem[25:],
        }

        return filename_components

    def load_raw_observations(self) -> dict[ttime.Time, float]:

        return {
            ttime.Time(float(epoch)): float(value)
            for (epoch, value) in zip(
                self.data["tdb_seconds_since_j2000"],
                self.data["doppler_averaged_frequency_hz"],
            )
        }
