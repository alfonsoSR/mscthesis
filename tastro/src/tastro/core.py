from .config import CaseSetup


class SettingsGenerator[T]:

    def __init__(self, name: str, local_config: T, config: CaseSetup) -> None:

        self.name = name
        self.local = local_config
        self.config = config

        return None
