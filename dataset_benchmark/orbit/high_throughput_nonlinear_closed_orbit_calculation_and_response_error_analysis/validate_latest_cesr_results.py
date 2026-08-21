#!/usr/bin/env python3
"""Read-only consistency audit for the latest-CESR orbit-error workflows.

The validator deliberately does not import SciBmad/Julia or execute an
experiment.  It checks only checked-in TOML/CSV sidecars and result tables.
Smoke artifacts are useful evidence of a working pipeline, but they are not
silently promoted to production: absent production files are reported as
``MISSING_PRODUCTION``.  Use ``--require-production`` when that condition
should make the command fail (exit status 2).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
import tomllib
from collections import Counter
from pathlib import Path
from typing import Any


EXPECTED_RING = "latest_cesr"
LATEST_LATTICE_FRAGMENT = "latest_lattice/latest_cesr_scibmad_repaired.jl"
CLOSURE_TOLERANCE = 1.0e-8
PRODUCTION_RHO_TRIALS = 600

COUNT_ALIASES: dict[str, tuple[str, ...]] = {
    "control_count": ("control_count", "controls"),
    "detector_count": ("detector_count", "detectors"),
    "element_count": ("element_count", "elements"),
    "active_normal_sextupoles": (
        "active_normal_sextupoles",
        "active_normal_sextupole_count",
        "active_normal_sextupule_count",
        "sextupole_count",
    ),
}

LIST_ALIASES: dict[str, tuple[str, ...]] = {
    "control_count": ("control_names", "control_names_from_config"),
    "detector_count": ("detector_names",),
    "element_count": ("element_names",),
    "active_normal_sextupoles": (
        "active_normal_sextupole_names",
        "active_normal_sextupule_names",
        "normal_sextupole_names",
    ),
}

FORBIDDEN_HEADER_PATTERNS = (
    re.compile(r"(^|_)(hh|hv|vv)(_|$)"),
    re.compile(r"(^|_)(third|cubic)(_|$)"),
    re.compile(r"(^|_)order_?3(_|$)"),
    re.compile(r"(^|_)(block_share|share_block)(_|$)"),
)


def normalized_name(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def is_forbidden_output_name(value: object) -> bool:
    return any(pattern.search(normalized_name(value)) for pattern in FORBIDDEN_HEADER_PATTERNS)


def as_finite_number(value: object) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        number = float(value)
    else:
        text = str(value).strip()
        if not text:
            return None
        try:
            number = float(text)
        except (TypeError, ValueError):
            return None
    return number if math.isfinite(number) else math.nan


def as_count(value: object) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (list, tuple)):
        return len(value)
    number = as_finite_number(value)
    if number is None or math.isnan(number) or not number.is_integer():
        return None
    return int(number)


def bool_from_csv(value: object) -> bool | None:
    text = str(value).strip().lower()
    if text in {"true", "t", "yes", "y", "1"}:
        return True
    if text in {"false", "f", "no", "n", "0"}:
        return False
    return None


def path_from_value(value: object, project_root: Path) -> Path | None:
    if value is None or not str(value).strip():
        return None
    candidate = Path(str(value))
    if candidate.is_file():
        return candidate
    if not candidate.is_absolute():
        candidate = project_root / candidate
    return candidate


class LatestCesrAudit:
    """Accumulate a small, JSON-serializable audit report."""

    def __init__(self, project_root: Path) -> None:
        self.project_root = project_root
        self.records: list[dict[str, Any]] = []
        self.metadata: list[dict[str, Any]] = []

    def add(self, status: str, scope: str, message: str, path: Path | None = None) -> None:
        record: dict[str, Any] = {
            "status": status,
            "scope": scope,
            "message": message,
        }
        if path is not None:
            record["path"] = str(path)
        self.records.append(record)

    def _new_records(self, start: int) -> list[dict[str, Any]]:
        return self.records[start:]

    def _has_status(self, status: str, start: int = 0) -> bool:
        return any(record["status"] == status for record in self.records[start:])

    def read_toml(self, path: Path, scope: str) -> dict[str, Any] | None:
        if not path.is_file():
            self.add("FAIL", scope, "missing TOML metadata", path)
            return None
        try:
            with path.open("rb") as stream:
                value = tomllib.load(stream)
        except (OSError, tomllib.TOMLDecodeError) as exc:
            self.add("FAIL", scope, f"cannot parse TOML: {exc}", path)
            return None
        if not isinstance(value, dict):
            self.add("FAIL", scope, "TOML root is not a table", path)
            return None
        return value

    def read_json(self, path: Path, scope: str) -> dict[str, Any] | None:
        if not path.is_file():
            self.add("FAIL", scope, "missing JSON metadata", path)
            return None
        try:
            with path.open("r", encoding="utf-8") as stream:
                value = json.load(stream)
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            self.add("FAIL", scope, f"cannot parse JSON: {exc}", path)
            return None
        if not isinstance(value, dict):
            self.add("FAIL", scope, "JSON root is not an object", path)
            return None
        return value

    def read_rho_metadata(
        self,
        directory: Path,
        scope: str,
    ) -> tuple[dict[str, Any] | None, Path | None]:
        """Read either a single-run TOML sidecar or a unified merged JSON sidecar."""

        toml_path = directory / "rho_sweep_metadata.toml"
        if toml_path.is_file():
            return self.read_toml(toml_path, scope + "/metadata"), toml_path
        json_path = directory / "rho_sweep_metadata.json"
        if json_path.is_file():
            return self.read_json(json_path, scope + "/metadata"), json_path
        self.add("FAIL", scope + "/metadata", "missing rho_sweep_metadata.toml or rho_sweep_metadata.json", directory)
        return None, None

    def read_csv(self, path: Path, scope: str) -> tuple[list[str], list[dict[str, str]]] | None:
        if not path.is_file():
            self.add("FAIL", scope, "missing CSV output", path)
            return None
        try:
            with path.open("r", newline="", encoding="utf-8-sig") as stream:
                reader = csv.DictReader(stream)
                headers = list(reader.fieldnames or [])
                rows = list(reader)
        except (OSError, UnicodeError, csv.Error) as exc:
            self.add("FAIL", scope, f"cannot parse CSV: {exc}", path)
            return None
        if not headers:
            self.add("FAIL", scope, "CSV has no header", path)
            return None
        return headers, rows

    def check_csv(
        self,
        path: Path,
        scope: str,
        required_headers: tuple[str, ...] = (),
        allow_empty: bool = False,
    ) -> tuple[list[str], list[dict[str, str]]] | None:
        result = self.read_csv(path, scope)
        if result is None:
            return None
        headers, rows = result
        missing = [header for header in required_headers if header not in headers]
        if missing:
            self.add("FAIL", scope, f"missing required columns: {', '.join(missing)}", path)
        forbidden = [header for header in headers if is_forbidden_output_name(header)]
        if forbidden:
            self.add(
                "FAIL",
                scope,
                "paper-scope output contains forbidden block/order columns: "
                + ", ".join(forbidden),
                path,
            )
        if not rows and not allow_empty:
            self.add("FAIL", scope, "CSV has no data rows", path)
        nonfinite: list[str] = []
        for row_number, row in enumerate(rows, start=2):
            for header, value in row.items():
                if value is None or not str(value).strip():
                    continue
                number = as_finite_number(value)
                if number is not None and math.isnan(number):
                    nonfinite.append(f"{row_number}:{header}={value}")
                    if len(nonfinite) >= 5:
                        break
            if len(nonfinite) >= 5:
                break
        if nonfinite:
            self.add("FAIL", scope, "non-finite numeric cells: " + "; ".join(nonfinite), path)
        return headers, rows

    def metadata_counts(self, metadata: dict[str, Any]) -> dict[str, int]:
        counts: dict[str, int] = {}
        for canonical, aliases in COUNT_ALIASES.items():
            for alias in aliases:
                if alias not in metadata:
                    continue
                count = as_count(metadata[alias])
                if count is None:
                    self.add("FAIL", "metadata", f"invalid count value for {alias}")
                else:
                    counts[canonical] = count
                break
        return counts

    def check_metadata(
        self,
        metadata: dict[str, Any],
        path: Path,
        scope: str,
        *,
        require_version: bool = True,
    ) -> dict[str, int]:
        start = len(self.records)
        if metadata.get("ring_id") != EXPECTED_RING:
            self.add("FAIL", scope, f"ring_id is {metadata.get('ring_id')!r}, expected {EXPECTED_RING!r}", path)
        lattice = str(metadata.get("lattice_path", "")).replace("\\", "/").lower()
        if LATEST_LATTICE_FRAGMENT not in lattice:
            self.add("FAIL", scope, "lattice_path is not the repaired latest-CESR SciBmad lattice", path)
        if not str(metadata.get("format", "")).strip():
            self.add("FAIL", scope, "metadata is missing format", path)
        version = str(metadata.get("scibmad_version", "")).strip()
        if not version:
            engine = str(metadata.get("engine", "")).strip().lower()
            if engine == "scibmad":
                self.add("WARN", scope, "engine is SciBmad but scibmad_version is omitted; use companion metadata for the exact version", path)
            elif require_version:
                self.add("FAIL", scope, "metadata is missing scibmad_version", path)
        if "rf_on" in metadata and metadata["rf_on"] is not True:
            self.add("FAIL", scope, f"rf_on is {metadata['rf_on']!r}; latest-CESR result must be RF-on", path)
        elif "rf_on" not in metadata:
            self.add("WARN", scope, "metadata does not state rf_on", path)

        counts = self.metadata_counts(metadata)
        for canonical, count in counts.items():
            if count <= 0:
                self.add("FAIL", scope, f"{canonical} must be positive, got {count}", path)
            for list_key in LIST_ALIASES[canonical]:
                values = metadata.get(list_key)
                if isinstance(values, list) and len(values) != count:
                    self.add(
                        "FAIL",
                        scope,
                        f"{list_key} has {len(values)} entries but {canonical}={count}",
                        path,
                    )

        control_count = counts.get("control_count")
        horizontal = as_count(metadata.get("horizontal_control_count"))
        vertical = as_count(metadata.get("vertical_control_count"))
        if control_count is not None and horizontal is not None and vertical is not None:
            if horizontal + vertical != control_count:
                self.add("FAIL", scope, "horizontal_control_count + vertical_control_count does not equal control_count", path)

        detector_count = counts.get("detector_count")
        per_plane = as_count(metadata.get("detector_count_per_plane"))
        if detector_count is not None and per_plane is not None and detector_count != per_plane:
            self.add("FAIL", scope, "detector_count differs from detector_count_per_plane", path)
        observable_count = as_count(metadata.get("observable_count"))
        if detector_count is not None and observable_count is not None and observable_count != 2 * detector_count:
            self.add("FAIL", scope, "observable_count is not twice detector_count", path)
        state_dimension = as_count(metadata.get("state_dimension"))
        if state_dimension is not None and state_dimension <= 0:
            self.add("FAIL", scope, "state_dimension must be positive", path)

        forbidden_keys = [key for key in metadata if is_forbidden_output_name(key)]
        if forbidden_keys:
            self.add("FAIL", scope, "metadata contains forbidden output keys: " + ", ".join(forbidden_keys), path)

        self.metadata.append({"scope": scope, "path": path, "metadata": metadata, "counts": counts})
        if not self._has_status("FAIL", start):
            self.add("PASS", scope, "latest-ring provenance and dynamic metadata are valid", path)
        return counts

    def choose_output(
        self,
        directory: Path,
        canonical: str,
        legacy: str,
        scope: str,
    ) -> Path | None:
        canonical_path = directory / canonical
        if canonical_path.is_file():
            return canonical_path
        legacy_path = directory / legacy
        if legacy_path.is_file():
            self.add("WARN", scope, f"legacy output name {legacy!r} used; prefer {canonical!r}", legacy_path)
            return legacy_path
        self.add("FAIL", scope, f"missing {canonical} (and legacy fallback {legacy})", directory)
        return None

    def check_production_presence(self, directory: Path, scope: str, required: tuple[str, ...]) -> bool:
        missing = [name for name in required if not (directory / name).is_file()]
        if not missing:
            return True
        present = [name for name in required if (directory / name).is_file()]
        if not present:
            self.add(
                "MISSING_PRODUCTION",
                scope,
                "production output is absent; smoke artifacts are not promoted to production",
                directory,
            )
        else:
            self.add("FAIL", scope, "production output is incomplete; missing " + ", ".join(missing), directory)
        return False

    def check_rho_production(self, latest_directory: Path) -> None:
        """Accept either a single-run TOML root or the official merged JSON layout."""

        merged_directory = latest_directory / "merged"
        merged_required = (
            "rho_sweep_metadata.json",
            "rho_sweep_summary.csv",
            "rho_sweep_trial_errors.csv",
        )
        root_required = (
            "rho_sweep_metadata.toml",
            "rho_sweep_summary.csv",
            "rho_sweep_trial_errors.csv",
        )
        merged_complete = all((merged_directory / name).is_file() for name in merged_required)
        root_complete = all((latest_directory / name).is_file() for name in root_required)
        if merged_complete:
            self.validate_rho_directory(merged_directory, "rho/production/merged", production=True)
            return
        if root_complete:
            self.validate_rho_directory(latest_directory, "rho/production/root", production=True)
            return

        merged_present = [name for name in merged_required if (merged_directory / name).is_file()]
        root_present = [name for name in root_required if (latest_directory / name).is_file()]
        if merged_present:
            missing = [name for name in merged_required if name not in merged_present]
            self.add(
                "FAIL",
                "rho/production/merged",
                "merged production output is incomplete; missing " + ", ".join(missing),
                merged_directory,
            )
        elif root_present:
            missing = [name for name in root_required if name not in root_present]
            self.add(
                "FAIL",
                "rho/production/root",
                "single-run production output is incomplete; missing " + ", ".join(missing),
                latest_directory,
            )
        else:
            self.add(
                "MISSING_PRODUCTION",
                "rho/production",
                "production output is absent; checked both latest_cesr/merged (unified JSON) and latest_cesr root (single-run TOML)",
                latest_directory,
            )

    def validate_rho_directory(self, directory: Path, scope: str, *, production: bool) -> None:
        start = len(self.records)
        summary_path = directory / "rho_sweep_summary.csv"
        trial_path = directory / "rho_sweep_trial_errors.csv"
        metadata, metadata_path = self.read_rho_metadata(directory, scope)
        if metadata is not None:
            # metadata_path is guaranteed non-None whenever metadata loaded.
            assert metadata_path is not None
            counts = self.check_metadata(metadata, metadata_path, scope + "/metadata")
            if production:
                trials_key = next(
                    (
                        key
                        for key in (
                            "trials_per_rho_scenario",
                            "trials_per_positive_rho_scenario",
                            "trials_per_rho",
                        )
                        if key in metadata
                    ),
                    None,
                )
                trials = as_count(metadata.get(trials_key)) if trials_key else None
                if trials is None or trials < PRODUCTION_RHO_TRIALS:
                    self.add(
                        "FAIL",
                        scope,
                        f"production {trials_key or 'trials_per_rho_scenario'} is {trials!r}; expected at least {PRODUCTION_RHO_TRIALS}",
                        metadata_path,
                    )
            total_key = next(
                (key for key in ("total_samples", "total_unique_samples") if key in metadata),
                None,
            )
            total_samples = as_count(metadata.get(total_key)) if total_key else None
            converged_count = as_count(metadata.get("converged_count"))
            failed_count = as_count(metadata.get("failed_count"))
            if failed_count is not None and failed_count != 0:
                self.add("FAIL", scope, f"rho metadata reports failed_count={failed_count}", metadata_path)
            if total_samples is not None and converged_count is not None and total_samples != converged_count:
                self.add("FAIL", scope, "rho total_samples does not equal converged_count", metadata_path)
            detector_count = counts.get("detector_count")
            observable_count = as_count(metadata.get("observable_count"))
            if detector_count is not None and observable_count is not None and observable_count != 2 * detector_count:
                self.add("FAIL", scope, "rho observable count is inconsistent with the two detector planes", metadata_path)

        summary = self.check_csv(
            summary_path,
            scope + "/summary",
            ("scenario", "rho", "trials", "converged_trials", "max_closure_norm"),
        )
        trials = self.check_csv(
            trial_path,
            scope + "/trials",
            ("scenario", "rho", "converged", "closure_norm"),
        )
        if trials is not None:
            _, rows = trials
            unconverged = [str(row.get("sample_id", index + 1)) for index, row in enumerate(rows) if bool_from_csv(row.get("converged")) is False]
            if unconverged:
                self.add("FAIL", scope, "unconverged rho samples: " + ", ".join(unconverged[:8]), trial_path)
            expected = None
            if metadata:
                total_key = next(
                    (key for key in ("total_samples", "total_unique_samples") if key in metadata),
                    None,
                )
                expected = as_count(metadata.get(total_key)) if total_key else None
            if expected is not None and len(rows) != expected:
                self.add("FAIL", scope, f"trial CSV has {len(rows)} rows but metadata total_samples={expected}", trial_path)
            for row in rows:
                closure = as_finite_number(row.get("closure_norm"))
                if closure is not None and (math.isnan(closure) or closure < 0 or closure > CLOSURE_TOLERANCE):
                    self.add("FAIL", scope, f"rho closure_norm exceeds {CLOSURE_TOLERANCE:g}", trial_path)
                    break
        if summary is not None:
            _, rows = summary
            for row in rows:
                closure = as_finite_number(row.get("max_closure_norm"))
                if closure is not None and (math.isnan(closure) or closure < 0 or closure > CLOSURE_TOLERANCE):
                    self.add("FAIL", scope, f"rho summary max_closure_norm exceeds {CLOSURE_TOLERANCE:g}", summary_path)
                    break
        if not self._has_status("FAIL", start):
            self.add("PASS" if production else "SMOKE", scope, "rho sweep tables have finite values, closure, and latest provenance", directory)

    def validate_thick_directory(self, directory: Path, scope: str, plane: str, *, production: bool) -> None:
        start = len(self.records)
        metadata_path = directory / "metadata.toml"
        metadata = self.read_toml(metadata_path, scope + "/metadata")
        counts: dict[str, int] = {}
        if metadata is not None:
            counts = self.check_metadata(metadata, metadata_path, scope + "/metadata")

        reconstruction = self.check_csv(
            directory / "reconstruction_summary.csv",
            scope + "/reconstruction",
            ("trials", "output_plane", "output_coordinate", "elements", "active_normal_sextupoles", "detectors", "total_all_element_relative_closure"),
        )
        direction_closure = self.check_csv(
            directory / "direction_closure.csv",
            scope + "/direction_closure",
            ("trial", "all_element_total_relative_closure"),
        )
        element_path = self.choose_output(directory, "element_contribution_summary.csv", "thick_sextupole_contribution_summary.csv", scope + "/elements")
        element_direction_path = self.choose_output(directory, "element_direction_contributions.csv", "thick_sextupole_direction_contributions.csv", scope + "/element_directions")
        element_summary = self.check_csv(
            element_path,
            scope + "/elements",
            ("element_name", "element_type", "eta_total", "magnitude_total"),
        ) if element_path else None
        element_direction = self.check_csv(
            element_direction_path,
            scope + "/element_directions",
            ("trial", "element_name", "projection_numerator", "contribution_norm_m"),
        ) if element_direction_path else None
        family_summary = self.check_csv(
            directory / "family_contribution_summary.csv",
            scope + "/families",
            ("family", "element_count", "eta_total", "magnitude_total"),
        )
        family_direction = self.check_csv(
            directory / "family_direction_contributions.csv",
            scope + "/family_directions",
            ("trial", "family", "contribution_norm_m", "projection_numerator"),
        )

        if reconstruction is not None and metadata is not None:
            _, rows = reconstruction
            if len(rows) != 1:
                self.add("FAIL", scope, f"reconstruction_summary has {len(rows)} rows; expected one aggregate row", directory / "reconstruction_summary.csv")
            else:
                row = rows[0]
                expected_values = {
                    "trials": as_count(metadata.get("trials")),
                    "elements": counts.get("element_count"),
                    "active_normal_sextupoles": counts.get("active_normal_sextupoles"),
                    "detectors": counts.get("detector_count"),
                }
                for key, expected in expected_values.items():
                    actual = as_count(row.get(key))
                    if expected is not None and actual != expected:
                        self.add("FAIL", scope, f"reconstruction {key}={actual!r}, metadata says {expected!r}", directory / "reconstruction_summary.csv")
                if str(row.get("output_plane", "")).lower() != plane:
                    self.add("FAIL", scope, f"reconstruction output_plane is {row.get('output_plane')!r}, expected {plane!r}", directory / "reconstruction_summary.csv")
                for key, value in row.items():
                    number = as_finite_number(value)
                    if number is None:
                        continue
                    if "closure" in key and (math.isnan(number) or number < 0):
                        self.add("FAIL", scope, f"negative/non-finite closure metric {key}", directory / "reconstruction_summary.csv")
                    if "relative_closure" in key and number > CLOSURE_TOLERANCE:
                        self.add("FAIL", scope, f"{key} exceeds {CLOSURE_TOLERANCE:g}: {number:g}", directory / "reconstruction_summary.csv")

        trials_expected = as_count(metadata.get("trials")) if metadata else None
        elements_expected = counts.get("element_count")
        families_expected = len(metadata.get("source_families", [])) if metadata and isinstance(metadata.get("source_families"), list) else None
        if element_summary is not None and elements_expected is not None and len(element_summary[1]) != elements_expected:
            self.add("FAIL", scope, f"element summary has {len(element_summary[1])} rows but element_count={elements_expected}", element_path)
        if element_direction is not None and elements_expected is not None and trials_expected is not None:
            expected_rows = elements_expected * trials_expected
            if len(element_direction[1]) != expected_rows:
                self.add("FAIL", scope, f"element direction table has {len(element_direction[1])} rows but expected {expected_rows}", element_direction_path)
        if family_summary is not None and families_expected is not None and len(family_summary[1]) != families_expected:
            self.add("FAIL", scope, f"family summary has {len(family_summary[1])} rows but source_families has {families_expected}", directory / "family_contribution_summary.csv")
        if family_direction is not None and families_expected is not None and trials_expected is not None:
            expected_rows = families_expected * trials_expected
            if len(family_direction[1]) != expected_rows:
                self.add("FAIL", scope, f"family direction table has {len(family_direction[1])} rows but expected {expected_rows}", directory / "family_direction_contributions.csv")
        if direction_closure is not None and trials_expected is not None and len(direction_closure[1]) != trials_expected:
            self.add("FAIL", scope, f"direction_closure has {len(direction_closure[1])} rows but trials={trials_expected}", directory / "direction_closure.csv")
        if direction_closure is not None:
            for row in direction_closure[1]:
                for key, value in row.items():
                    if "closure" not in key:
                        continue
                    number = as_finite_number(value)
                    if number is None or math.isnan(number) or number < 0:
                        self.add("FAIL", scope, f"invalid closure metric {key}", directory / "direction_closure.csv")
                    elif "relative" in key and number > CLOSURE_TOLERANCE:
                        self.add("FAIL", scope, f"{key} exceeds {CLOSURE_TOLERANCE:g}: {number:g}", directory / "direction_closure.csv")

        if element_summary is not None and family_summary is not None:
            element_values = [as_finite_number(row.get("eta_total")) for row in element_summary[1]]
            family_values = [as_finite_number(row.get("eta_total")) for row in family_summary[1]]
            if any(value is None or math.isnan(value) for value in element_values + family_values):
                self.add("FAIL", scope, "eta_total contains a non-numeric or non-finite value", directory)
                element_eta = family_eta = math.nan
            else:
                element_eta = sum(value for value in element_values if value is not None)
                family_eta = sum(value for value in family_values if value is not None)
            if abs(element_eta - family_eta) > 1.0e-10 * max(1.0, abs(element_eta), abs(family_eta)):
                self.add("FAIL", scope, f"family eta_total sum {family_eta:g} does not close to element sum {element_eta:g}", directory / "family_contribution_summary.csv")
            if reconstruction is not None and reconstruction[1]:
                projection = as_finite_number(reconstruction[1][0].get("total_all_element_signed_projection"))
                if projection is not None and not math.isnan(projection) and abs(family_eta - projection) > 1.0e-10 * max(1.0, abs(projection)):
                    self.add("FAIL", scope, f"family eta_total sum {family_eta:g} does not close to aggregate signed projection {projection:g}", directory / "reconstruction_summary.csv")

        if not self._has_status("FAIL", start):
            self.add("PASS" if production else "SMOKE", scope, "complete-element tables, family sum, and vector-closure metrics are valid", directory)

    def validate_sextupole_directory(self, directory: Path, scope: str, *, production: bool) -> None:
        start = len(self.records)
        metadata_path = directory / "metadata.toml"
        metadata = self.read_toml(metadata_path, scope + "/metadata")
        counts: dict[str, int] = {}
        if metadata is not None:
            counts = self.check_metadata(metadata, metadata_path, scope + "/metadata")
        reconstruction = self.check_csv(
            directory / "reconstruction_summary.csv",
            scope + "/reconstruction",
            ("trials", "elements", "active_normal_sextupoles", "detectors", "all_element_relative_closure"),
        )
        closure = self.check_csv(
            directory / "direction_closure.csv",
            scope + "/direction_closure",
            ("trial", "all_element_relative_closure", "normal_sextupole_relative_closure"),
        )
        detectors = self.check_csv(directory / "detectors.csv", scope + "/detectors", ("detector_name", "s_m"))
        summary = self.check_csv(
            directory / "sextupole_contribution_summary.csv",
            scope + "/sextupoles",
            ("element_name", "eta_total", "magnitude_total"),
        )
        direction = self.check_csv(
            directory / "sextupole_direction_contributions.csv",
            scope + "/sextupole_directions",
            ("trial", "element_name", "total_projection_numerator", "total_contribution_norm_m"),
        )

        trials_expected = as_count(metadata.get("trials")) if metadata else None
        elements_expected = counts.get("element_count")
        active_expected = counts.get("active_normal_sextupoles")
        detectors_expected = counts.get("detector_count")
        if detectors is not None and detectors_expected is not None and len(detectors[1]) != detectors_expected:
            self.add("FAIL", scope, f"detector table has {len(detectors[1])} rows but detector_count={detectors_expected}", directory / "detectors.csv")
        if summary is not None and active_expected is not None and len(summary[1]) != active_expected:
            self.add("FAIL", scope, f"sextupole summary has {len(summary[1])} rows but active_normal_sextupoles={active_expected}", directory / "sextupole_contribution_summary.csv")
        if direction is not None and active_expected is not None and trials_expected is not None:
            expected_rows = active_expected * trials_expected
            if len(direction[1]) != expected_rows:
                self.add("FAIL", scope, f"sextupole direction table has {len(direction[1])} rows but expected {expected_rows}", directory / "sextupole_direction_contributions.csv")
        if closure is not None and trials_expected is not None and len(closure[1]) != trials_expected:
            self.add("FAIL", scope, f"direction_closure has {len(closure[1])} rows but trials={trials_expected}", directory / "direction_closure.csv")
        if reconstruction is not None:
            if len(reconstruction[1]) != 1:
                self.add("FAIL", scope, "reconstruction_summary must contain one aggregate row", directory / "reconstruction_summary.csv")
            else:
                row = reconstruction[1][0]
                for key, expected in {
                    "trials": trials_expected,
                    "elements": elements_expected,
                    "active_normal_sextupoles": active_expected,
                    "detectors": detectors_expected,
                }.items():
                    actual = as_count(row.get(key))
                    if expected is not None and actual != expected:
                        self.add("FAIL", scope, f"reconstruction {key}={actual!r}, metadata says {expected!r}", directory / "reconstruction_summary.csv")
                for key, value in row.items():
                    if "relative_closure" not in key:
                        continue
                    number = as_finite_number(value)
                    if number is None or math.isnan(number) or number < 0:
                        self.add("FAIL", scope, f"invalid closure metric {key}", directory / "reconstruction_summary.csv")
                    elif "all_element" in key and number > CLOSURE_TOLERANCE:
                        self.add("FAIL", scope, f"{key} exceeds {CLOSURE_TOLERANCE:g}: {number:g}", directory / "reconstruction_summary.csv")
        if closure is not None:
            for row in closure[1]:
                for key, value in row.items():
                    if "relative_closure" not in key:
                        continue
                    number = as_finite_number(value)
                    if number is None or math.isnan(number) or number < 0:
                        self.add("FAIL", scope, f"invalid closure metric {key}", directory / "direction_closure.csv")
                    elif "all_element" in key and number > CLOSURE_TOLERANCE:
                        self.add("FAIL", scope, f"{key} exceeds {CLOSURE_TOLERANCE:g}: {number:g}", directory / "direction_closure.csv")

        if not self._has_status("FAIL", start):
            self.add("PASS" if production else "SMOKE", scope, "normal-sextupole tables retain all-element vector closure and finite values", directory)

    def validate_beta_directory(self, directory: Path, scope: str, *, production: bool) -> None:
        start = len(self.records)
        required = (
            "metadata.toml",
            "nominal_optics_metadata.toml",
            "nominal_optics_points.csv",
            "direction_optics_metadata.toml",
            "direction_optics_points.csv",
            "direction_optics_tunes.csv",
            "direction_optics_variation_summary.csv",
            "direction_sextupole_transport_factors.csv",
            "sextupole_detector_optics.csv",
            "correlation_summary.csv",
            "direction_element_correlation_data.csv",
            "element_correlation_data.csv",
        )
        missing = [name for name in required if not (directory / name).is_file()]
        if missing:
            status = "FAIL" if production else "PARTIAL_SMOKE"
            self.add(status, scope, "beta-phase output is incomplete; missing " + ", ".join(missing), directory)
            return

        metadata_path = directory / "metadata.toml"
        nominal_path = directory / "nominal_optics_metadata.toml"
        direction_path = directory / "direction_optics_metadata.toml"
        metadata = self.read_toml(metadata_path, scope + "/metadata")
        nominal_metadata = self.read_toml(nominal_path, scope + "/nominal_metadata")
        direction_metadata = self.read_toml(direction_path, scope + "/direction_metadata")
        base_counts: dict[str, int] = {}
        nominal_counts: dict[str, int] = {}
        direction_counts: dict[str, int] = {}
        if metadata is not None:
            base_counts = self.check_metadata(metadata, metadata_path, scope + "/metadata", require_version=False)
        if nominal_metadata is not None:
            nominal_counts = self.check_metadata(nominal_metadata, nominal_path, scope + "/nominal_metadata")
        if direction_metadata is not None:
            direction_counts = self.check_metadata(direction_metadata, direction_path, scope + "/direction_metadata")

        active = base_counts.get("active_normal_sextupoles") or nominal_counts.get("active_normal_sextupoles") or direction_counts.get("active_normal_sextupoles")
        detectors = base_counts.get("detector_count") or nominal_counts.get("detector_count") or direction_counts.get("detector_count")
        directions = as_count(metadata.get("direction_count")) if metadata else None
        trials = as_count(direction_metadata.get("trials")) if direction_metadata else None
        if trials is None:
            trials = directions
        if directions is None:
            directions = trials
        if active is None or detectors is None or trials is None:
            self.add("FAIL", scope, "beta-phase metadata does not provide active sextupole, detector, and direction/trial counts", directory)
        if metadata is not None:
            rows_declared = as_count(metadata.get("direction_optics_rows"))
            if rows_declared is not None and active is not None and directions is not None and rows_declared != active * directions:
                self.add("FAIL", scope, f"direction_optics_rows={rows_declared} but active sextupoles × directions={active * directions}", metadata_path)

        nominal_points = self.check_csv(
            directory / "nominal_optics_points.csv",
            scope + "/nominal_points",
            ("point_type", "element_name", "beta_1_m", "beta_2_m", "phi_1_turn", "phi_2_turn"),
        )
        direction_points = self.check_csv(
            directory / "direction_optics_points.csv",
            scope + "/direction_points",
            ("trial", "point_type", "element_name", "beta_1_m", "beta_2_m", "phi_1_turn", "phi_2_turn"),
        )
        tunes = self.check_csv(directory / "direction_optics_tunes.csv", scope + "/direction_tunes", ("trial", "full_tune_1_turn", "full_tune_2_turn"))
        variation = self.check_csv(directory / "direction_optics_variation_summary.csv", scope + "/variation", ("quantity", "direction_rows"))
        transport = self.check_csv(directory / "direction_sextupole_transport_factors.csv", scope + "/transport", ("trial", "element_name", "beta_x_sext_m", "beta_y_sext_m"))
        detector_optics = self.check_csv(directory / "sextupole_detector_optics.csv", scope + "/detector_optics", ("sextupole_name", "detector_name", "beta_x_sext_m", "beta_y_sext_m"))
        correlation = self.check_csv(directory / "correlation_summary.csv", scope + "/correlation", ("plane", "predictor", "directions", "elements"))
        direction_correlation = self.check_csv(directory / "direction_element_correlation_data.csv", scope + "/direction_correlation", ("plane", "trial", "element_name", "actual_contribution_norm_m"))
        element_correlation = self.check_csv(directory / "element_correlation_data.csv", scope + "/element_correlation", ("plane", "element_name", "actual_rms_relative_magnitude"))

        if active is not None and detectors is not None:
            expected_points = active + detectors
            if nominal_points is not None:
                type_counts = Counter(row.get("point_type") for row in nominal_points[1])
                if len(nominal_points[1]) != expected_points or type_counts.get("detector", 0) != detectors or type_counts.get("sextupole_exit", 0) != active:
                    self.add("FAIL", scope, f"nominal optics points do not match detector+sextupole counts ({detectors}+{active})", directory / "nominal_optics_points.csv")
            if direction_points is not None and trials is not None:
                expected = trials * expected_points
                type_counts = Counter(row.get("point_type") for row in direction_points[1])
                if len(direction_points[1]) != expected or type_counts.get("detector", 0) != trials * detectors or type_counts.get("sextupole_exit", 0) != trials * active:
                    self.add("FAIL", scope, f"direction optics points do not match trials × (detectors+sextupoles)={expected}", directory / "direction_optics_points.csv")
            if detector_optics is not None and len(detector_optics[1]) != active * detectors:
                self.add("FAIL", scope, f"sextupole-detector optics has {len(detector_optics[1])} rows, expected {active * detectors}", directory / "sextupole_detector_optics.csv")
        if trials is not None:
            if tunes is not None and len(tunes[1]) != trials:
                self.add("FAIL", scope, f"direction tunes has {len(tunes[1])} rows but trials={trials}", directory / "direction_optics_tunes.csv")
            if active is not None:
                if transport is not None and len(transport[1]) != trials * active:
                    self.add("FAIL", scope, f"transport table has {len(transport[1])} rows but trials × active sextupoles={trials * active}", directory / "direction_sextupole_transport_factors.csv")
                if direction_correlation is not None and len(direction_correlation[1]) != trials * active * 2:
                    self.add("FAIL", scope, f"direction correlation table has {len(direction_correlation[1])} rows but trials × active × 2={trials * active * 2}", directory / "direction_element_correlation_data.csv")
            if element_correlation is not None and active is not None and len(element_correlation[1]) != active * 2:
                self.add("FAIL", scope, f"element correlation table has {len(element_correlation[1])} rows but active × 2={active * 2}", directory / "element_correlation_data.csv")
        if correlation is not None:
            planes = {row.get("plane") for row in correlation[1]}
            if not {"x", "y"}.issubset(planes):
                self.add("FAIL", scope, "correlation summary does not contain both x and y planes", directory / "correlation_summary.csv")

        for key in ("contributions_csv", "closure_csv"):
            if metadata is None or key not in metadata:
                self.add("FAIL", scope, f"beta metadata is missing {key}", metadata_path)
                continue
            referenced = path_from_value(metadata[key], self.project_root)
            if referenced is None or not referenced.is_file():
                self.add("FAIL", scope, f"beta metadata {key} does not resolve to an existing artifact", metadata_path)

        if not self._has_status("FAIL", start):
            self.add("PASS" if production else "SMOKE", scope, "nominal/direction optics and beta-phase correlation tables are complete and finite", directory)

    def check_cross_consistency(self) -> None:
        scope = "cross/latest_cesr"
        start = len(self.records)
        fields = ("control_count", "detector_count", "element_count", "active_normal_sextupoles")
        observed: dict[str, dict[int, list[str]]] = {field: {} for field in fields}
        versions: dict[str, list[str]] = {}
        for item in self.metadata:
            metadata = item["metadata"]
            counts = item["counts"]
            label = item["scope"]
            for field in fields:
                value = counts.get(field)
                if value is not None:
                    observed[field].setdefault(value, []).append(label)
            version = str(metadata.get("scibmad_version", "")).strip()
            if version:
                versions.setdefault(version, []).append(label)
        for field, values in observed.items():
            if len(values) > 1:
                detail = "; ".join(f"{value}: {', '.join(labels[:3])}" for value, labels in sorted(values.items()))
                self.add("FAIL", scope, f"{field} disagrees across latest workflows ({detail})")
        if len(versions) > 1:
            self.add("FAIL", scope, "scibmad_version disagrees across latest workflows: " + "; ".join(sorted(versions)))
        if not self._has_status("FAIL", start):
            compact = ", ".join(f"{field}={next(iter(values))}" for field, values in observed.items() if len(values) == 1)
            version = next(iter(versions)) if len(versions) == 1 else "not uniformly recorded"
            self.add("PASS", scope, f"dynamic registry agreement: {compact}; SciBmad {version}")

    def run(self) -> None:
        error_root = self.project_root / "dataset_benchmark" / "orbit" / "error_analysis"

        rho_root = error_root / "response_rho_sweep_600"
        rho_latest = rho_root / EXPECTED_RING
        rho_smokes = ("gtpsa_smoke", "smoke", "smoke_recheck", "smoke_runner")
        found_rho_smoke = False
        for name in rho_smokes:
            directory = rho_latest / name
            if directory.is_dir():
                found_rho_smoke = True
                self.validate_rho_directory(directory, f"rho/{name}", production=False)
        if not found_rho_smoke:
            self.add("FAIL", "rho/smoke", "no latest-CESR rho smoke directory found", rho_latest)
        self.check_rho_production(rho_latest)

        thick_root = error_root / "thick_element_sextupole_sourcing"
        for plane in ("horizontal", "vertical"):
            output_plane = "x" if plane == "horizontal" else "y"
            smoke = thick_root / f"{plane}_results" / f"{EXPECTED_RING}_smoke"
            production = thick_root / f"{plane}_results" / EXPECTED_RING
            if smoke.is_dir():
                self.validate_thick_directory(smoke, f"thick/{plane}/smoke", output_plane, production=False)
            else:
                self.add("FAIL", f"thick/{plane}/smoke", "latest-CESR smoke directory is absent", smoke)
            required = ("metadata.toml", "reconstruction_summary.csv", "element_contribution_summary.csv", "element_direction_contributions.csv", "family_contribution_summary.csv", "family_direction_contributions.csv")
            if self.check_production_presence(production, f"thick/{plane}/production", required):
                self.validate_thick_directory(production, f"thick/{plane}/production", output_plane, production=True)

        sext_root = error_root / "sextupole_detector_contributions"
        sext_smoke = sext_root / "smoke" / EXPECTED_RING
        sext_production = sext_root / "results" / EXPECTED_RING
        if sext_smoke.is_dir():
            self.validate_sextupole_directory(sext_smoke, "sextupole_attribution/smoke", production=False)
        else:
            self.add("FAIL", "sextupole_attribution/smoke", "latest-CESR smoke directory is absent", sext_smoke)
        sext_required = ("metadata.toml", "reconstruction_summary.csv", "direction_closure.csv", "detectors.csv", "sextupole_contribution_summary.csv", "sextupole_direction_contributions.csv")
        if self.check_production_presence(sext_production, "sextupole_attribution/production", sext_required):
            self.validate_sextupole_directory(sext_production, "sextupole_attribution/production", production=True)

        beta_root = error_root / "sextupole_beta_phase_correlation"
        beta_smoke = beta_root / "smoke" / EXPECTED_RING
        beta_production = beta_root / "results" / EXPECTED_RING
        if beta_smoke.is_dir():
            self.validate_beta_directory(beta_smoke, "beta_phase/smoke", production=False)
        else:
            self.add("FAIL", "beta_phase/smoke", "latest-CESR smoke directory is absent", beta_smoke)
        beta_required = ("metadata.toml", "nominal_optics_metadata.toml", "nominal_optics_points.csv", "direction_optics_metadata.toml", "direction_optics_points.csv", "direction_optics_tunes.csv", "correlation_summary.csv")
        if self.check_production_presence(beta_production, "beta_phase/production", beta_required):
            self.validate_beta_directory(beta_production, "beta_phase/production", production=True)

        self.check_cross_consistency()

    def overall_status(self) -> str:
        statuses = {record["status"] for record in self.records}
        if "FAIL" in statuses:
            return "FAIL"
        if "MISSING_PRODUCTION" in statuses:
            return "MISSING_PRODUCTION"
        if "PARTIAL_SMOKE" in statuses:
            return "PARTIAL_SMOKE"
        if "WARN" in statuses:
            return "WARN"
        if "SMOKE" in statuses:
            return "SMOKE"
        return "PASS"


def default_project_root() -> Path:
    # .../CESR Project/dataset_benchmark/orbit/<this directory>/script.py
    return Path(__file__).resolve().parents[3]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=default_project_root(), help="CESR Project directory (default: inferred from this file)")
    parser.add_argument("--require-production", action="store_true", help="return exit status 2 when any ring-scoped production output is missing")
    parser.add_argument("--json", action="store_true", dest="as_json", help="emit the report as JSON")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    project_root = args.project_root.resolve()
    audit = LatestCesrAudit(project_root)
    audit.run()
    overall = audit.overall_status()
    report = {
        "overall": overall,
        "project_root": str(project_root),
        "production_required": bool(args.require_production),
        "records": audit.records,
    }
    if args.as_json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(f"OVERALL: {overall}")
        for record in audit.records:
            path = f" [{record['path']}]" if "path" in record else ""
            print(f"[{record['status']}] {record['scope']}: {record['message']}{path}")
    if overall == "FAIL":
        return 1
    if args.require_production and overall == "MISSING_PRODUCTION":
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
