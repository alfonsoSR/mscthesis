import typing
from dataclasses import dataclass


if typing.TYPE_CHECKING:
    from .config import CaseSetup


class AutoDataclass(type):

    def __new__(mcs, name, bases, attrs):

        cls = super().__new__(mcs, name, bases, attrs)

        return dataclass(cls)  # type: ignore


class SettingsGenerator[T]:

    def __init__(self, name: str, local_config: T, config: "CaseSetup") -> None:

        self.name = name
        self.local = local_config
        self.config = config

        return None
