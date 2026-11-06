"""Parse the Spectre-style ``.sp`` netlists and ``.sym`` matching groups
used by the ``mycase`` examples.

The parser intentionally supports a small, explicit subset:

* ``subckt NAME ...`` / ``ends``
* MOS:  ``M1 (d g s b) model w=2u l=1u multi=4 nf=10``
* CAP:  ``C1 (p n) model nr=6 l=1.5u w=1u ...``
* RES:  ``R1 (p n) model l=10u w=4u series=6 ...``

Parameters are converted to microns (for layout) using the usual SPICE
suffixes.  ``multi``/``m``/``MR`` are kept as the device multiplier; ``nf``
is kept as the number of fingers per multiplier.  A matching group line in
``.sym`` is a whitespace-separated list of fully-qualified device names.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re
from typing import Iterable

_UNITS = {
    "f": 1e-15,
    "p": 1e-12,
    "n": 1e-9,
    "u": 1e-6,
    "m": 1e-3,
    "k": 1e3,
    "meg": 1e6,
}

_SUBCKT_RE = re.compile(r"^\s*\.?subckt\s+(\w+)", re.IGNORECASE)
_ENDS_RE = re.compile(r"^\s*\.?ends\b", re.IGNORECASE)
_MOS_RE = re.compile(
    r"^\s*(M\w*)\s*\((.*?)\)\s*(\w+)\s*(.*)$", re.IGNORECASE
)
_PASSIVE_RE = re.compile(
    r"^\s*([CR]\w*)\s*\((.*?)\)\s*(\w+)\s*(.*)$", re.IGNORECASE
)
_PARAM_RE = re.compile(
    r"(\w+)\s*=\s*([0-9.+\-eE]+[a-zA-Z]*)", re.IGNORECASE
)


def parse_value(text: str) -> float:
    """Convert a SPICE value such as ``2u`` or ``320n`` to microns."""
    match = re.match(r"^\s*([0-9.+\-eE]+)([a-zA-Z]*)\s*$", text)
    if not match:
        raise ValueError(f"invalid SPICE value: {text!r}")
    number = float(match.group(1))
    unit = match.group(2).lower()
    if unit in _UNITS:
        return number * _UNITS[unit] / 1e-6
    # No suffix: treat the value as microns for layout parameters.
    return number


def _parse_params(text: str) -> dict[str, float]:
    params: dict[str, float] = {}
    for match in _PARAM_RE.finditer(text):
        key = match.group(1).lower()
        try:
            params[key] = parse_value(match.group(2))
        except ValueError:
            continue
    return params


@dataclass
class Device:
    name: str
    kind: str  # "mos" | "cap" | "res"
    model: str
    nodes: tuple[str, ...]
    params: dict[str, float] = field(default_factory=dict)
    subckt: str = ""

    # ---- multiplicity helpers -------------------------------------------
    @property
    def multi(self) -> int:
        for key in ("multi", "m", "mr"):
            if key in self.params:
                return int(self.params[key])
        return 1

    @property
    def nf(self) -> int:
        return int(self.params.get("nf", 1))

    @property
    def nr(self) -> int:
        return int(self.params.get("nr", 1))

    @property
    def series(self) -> int:
        return int(self.params.get("series", 1))

    # ---- device type / model attributes --------------------------------
    @property
    def nch(self) -> bool:
        return "nch" in self.model.lower()

    @property
    def attr(self) -> list[str]:
        model = self.model.lower()
        return [key for key in ("lvt", "hvt", "od", "ud", "na", "mac") if key in model]

    @property
    def unit_count(self) -> int:
        """Number of unit cells in the default (multi-level) convention."""
        return self.multi

    @property
    def unit_width(self) -> float:
        """Per-finger width in microns (Spectre convention)."""
        return float(self.params.get("w", 0.0))

    @property
    def unit_length(self) -> float:
        return float(self.params.get("l", 0.0))

    def __repr__(self) -> str:
        return (
            f"Device({self.name!r}, kind={self.kind!r}, model={self.model!r}, "
            f"w={self.unit_width:g}, l={self.unit_length:g}, "
            f"multi={self.multi}, nf={self.nf})"
        )


def parse_netlist(path: str | Path) -> list[Device]:
    """Parse a Spectre-style netlist and return all supported devices."""
    path = Path(path)
    devices: list[Device] = []
    current_subckt = ""

    for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.strip()
        if not line or line.startswith(("*", "//", "#")):
            continue

        subckt_match = _SUBCKT_RE.match(line)
        if subckt_match:
            current_subckt = subckt_match.group(1)
            continue

        if _ENDS_RE.match(line):
            current_subckt = ""
            continue

        mos_match = _MOS_RE.match(line)
        if mos_match:
            name, node_text, model, param_text = mos_match.groups()
            devices.append(
                Device(
                    name=name,
                    kind="mos",
                    model=model,
                    nodes=tuple(node_text.split()),
                    params=_parse_params(param_text),
                    subckt=current_subckt,
                )
            )
            continue

        passive_match = _PASSIVE_RE.match(line)
        if passive_match:
            name, node_text, model, param_text = passive_match.groups()
            kind = "cap" if name.upper().startswith("C") else "res"
            devices.append(
                Device(
                    name=name,
                    kind=kind,
                    model=model,
                    nodes=tuple(node_text.split()),
                    params=_parse_params(param_text),
                    subckt=current_subckt,
                )
            )

    return devices


def _common_prefix_ending_with_underscore(names: Iterable[str]) -> str:
    names = list(names)
    if not names:
        return ""
    first = names[0]
    candidates = [i for i, ch in enumerate(first) if ch == "_"]
    for i in reversed(candidates):
        prefix = first[: i + 1]
        if all(name.startswith(prefix) for name in names):
            return prefix
    return ""


def parse_sym_groups(path: str | Path) -> list[list[str]]:
    """Parse a ``.sym`` file into local device-name groups.

    ``OTA2_my_M0 OTA2_my_M7`` becomes ``["M0", "M7"]``.  Blank lines and
    comments are ignored.  A line with a single device is returned as a
    one-element group.
    """
    path = Path(path)
    lines: list[list[str]] = []
    all_names: list[str] = []

    for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        names = line.split()
        lines.append(names)
        all_names.extend(names)

    prefix = _common_prefix_ending_with_underscore(all_names)
    groups: list[list[str]] = []
    for names in lines:
        groups.append([name[len(prefix):] if name.startswith(prefix) else name for name in names])
    return groups


def load_case(
    netlist_path: str | Path,
    sym_path: str | Path,
    subckt: str | None = None,
) -> tuple[dict[str, Device], list[list[str]]]:
    """Convenience wrapper: return ``({name: Device}, groups)``."""
    devices = parse_netlist(netlist_path)
    if subckt is not None:
        devices = [d for d in devices if d.subckt == subckt]
    by_name = {d.name: d for d in devices}
    groups = parse_sym_groups(sym_path)
    return by_name, groups
