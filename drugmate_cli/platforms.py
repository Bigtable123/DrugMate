"""Metadata for the retained Python CLI tool and skill configuration."""

from collections import OrderedDict
from typing import NamedTuple


class PlatformInfo(NamedTuple):
    label: str
    default_toolset: str


PLATFORMS: OrderedDict[str, PlatformInfo] = OrderedDict([
    ("cli", PlatformInfo(label="Python CLI", default_toolset="drugmate-cli")),
])


def platform_label(key: str, default: str = "") -> str:
    info = PLATFORMS.get(key)
    return info.label if info is not None else default


def get_all_platforms() -> "OrderedDict[str, PlatformInfo]":
    return OrderedDict(PLATFORMS)
