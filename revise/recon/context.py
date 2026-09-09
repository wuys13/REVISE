from __future__ import annotations

import copy
import logging
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, ClassVar, Dict, List, Optional

from anndata import AnnData

from revise.svc import SVC


@dataclass
class PipelineContext:
    STAGES: ClassVar[tuple[tuple[str, str], ...]] = (
        ("validate_inputs", "pipeline.validate_inputs"),
        ("global_anchoring", "strategy.global_anchoring"),
        ("local_refinement", "strategy.local_refinement"),
        ("finalize", "pipeline.finalize"),
        ("evaluate", "pipeline.evaluate"),
    )

    merged_config: Dict[str, Any]
    profile: Optional[str]
    runtime: Dict[str, Any]
    route_key: str
    run_dir: Path
    logger: logging.Logger
    engine_defaults_hash: Optional[str] = None
    authority_hash: Optional[str] = None
    algorithm_config_hash: Optional[str] = None
    effective_config_hash: Optional[str] = None
    dry_run: bool = False
    finalize_callback: Optional[Callable[["PipelineContext"], None]] = None
    application_config_metadata: Dict[str, Any] = field(default_factory=dict)
    benchmark_config_metadata: Dict[str, Any] = field(default_factory=dict)

    runner_config: Any = None
    runner: Any = None
    input_specs: Any = None
    st_adata: Optional[AnnData] = None
    sc_ref_adata: Optional[AnnData] = None
    real_st_adata: Optional[AnnData] = None
    svc: Optional[SVC] = None

    quality_metrics: Dict[str, Any] = field(default_factory=dict)
    artifacts: Dict[str, Any] = field(default_factory=dict)
    provenance: Dict[str, Any] = field(default_factory=dict)
    artifact_records: List[Dict[str, Any]] = field(default_factory=list)
    sr_allocation_records: List[Dict[str, Any]] = field(default_factory=list)
    input_identities: List[Dict[str, str]] = field(default_factory=list)
    software_versions: Dict[str, str] = field(default_factory=dict)
    pm_on_cell: Any = None

    run_status: str = field(init=False, default="running")
    run_started_at: str = field(init=False)
    run_ended_at: Optional[str] = field(init=False, default=None)
    run_duration_seconds: Optional[float] = field(init=False, default=None)
    run_error: Optional[Dict[str, str]] = field(init=False, default=None)
    stage_records: List[Dict[str, Any]] = field(init=False)
    _stage_started_monotonic: Dict[str, float] = field(
        init=False, default_factory=dict, repr=False
    )
    _provenance_callback: Optional[Callable[["PipelineContext"], None]] = field(
        init=False, default=None, repr=False
    )
    _publication_commit: Optional[Callable[[], None]] = field(
        init=False, default=None, repr=False
    )
    _publication_rollback: Optional[Callable[[], None]] = field(
        init=False, default=None, repr=False
    )
    _run_started_monotonic: float = field(init=False, repr=False)
    local_refinement_record: Dict[str, Any] = field(init=False)

    def __post_init__(self) -> None:
        self.run_started_at = self._timestamp()
        self._run_started_monotonic = time.monotonic()
        self.stage_records = [
            {
                "name": name,
                "owner": owner,
                "status": "pending",
                "started_at": None,
                "duration_seconds": None,
                "reason": None,
                "error": None,
            }
            for name, owner in self.STAGES
        ]
        refinement = self.merged_config.get("local_refinement")
        self.local_refinement_record = {
            "route": self.route_key,
            "applied": False,
            "strength": (
                float(refinement["strength"])
                if refinement is not None
                else None
            ),
        }

    @staticmethod
    def _timestamp() -> str:
        return datetime.now(timezone.utc).isoformat()

    @staticmethod
    def _error_record(error: BaseException) -> Dict[str, str]:
        return {"type": type(error).__name__, "message": str(error)}

    @property
    def run_record(self) -> Dict[str, Any]:
        return {
            "status": self.run_status,
            "dry_run": bool(self.dry_run),
            "started_at": self.run_started_at,
            "ended_at": self.run_ended_at,
            "duration_seconds": self.run_duration_seconds,
            "error": self.run_error,
        }

    def set_provenance_callback(
        self,
        callback: Optional[Callable[["PipelineContext"], None]],
        *,
        notify: bool = True,
    ) -> None:
        self._provenance_callback = callback
        if notify:
            self._notify_provenance()

    def _notify_provenance(self) -> None:
        if self._provenance_callback is not None:
            self._provenance_callback(self)

    def set_pending_publication(
        self,
        *,
        commit: Callable[[], None],
        rollback: Callable[[], None],
    ) -> None:
        if self._publication_commit is not None:
            raise RuntimeError("A result publication is already pending")
        self._publication_commit = commit
        self._publication_rollback = rollback

    def commit_pending_publication(self) -> None:
        commit = self._publication_commit
        self._publication_commit = None
        self._publication_rollback = None
        if commit is None:
            return
        try:
            commit()
        except Exception:
            self.logger.warning(
                "[framework] could not remove the previous result backup",
                exc_info=True,
            )

    def rollback_pending_publication(self) -> None:
        rollback = self._publication_rollback
        self._publication_commit = None
        self._publication_rollback = None
        if rollback is not None:
            rollback()

    @contextmanager
    def _durable_transition(self):
        snapshot = {
            "run_status": self.run_status,
            "run_ended_at": self.run_ended_at,
            "run_duration_seconds": self.run_duration_seconds,
            "run_error": copy.deepcopy(self.run_error),
            "stage_records": copy.deepcopy(self.stage_records),
            "stage_started_monotonic": dict(self._stage_started_monotonic),
            "artifact_records": copy.deepcopy(self.artifact_records),
            "sr_allocation_records": copy.deepcopy(self.sr_allocation_records),
            "local_refinement_record": copy.deepcopy(
                self.local_refinement_record
            ),
        }
        try:
            yield
            self._notify_provenance()
        except BaseException:
            self.run_status = snapshot["run_status"]
            self.run_ended_at = snapshot["run_ended_at"]
            self.run_duration_seconds = snapshot["run_duration_seconds"]
            self.run_error = snapshot["run_error"]
            self.stage_records = snapshot["stage_records"]
            self._stage_started_monotonic = snapshot[
                "stage_started_monotonic"
            ]
            self.artifact_records = snapshot["artifact_records"]
            self.sr_allocation_records = snapshot["sr_allocation_records"]
            self.local_refinement_record = snapshot["local_refinement_record"]
            raise

    def _stage_record(self, name: str) -> Dict[str, Any]:
        for record in self.stage_records:
            if record["name"] == name:
                return record
        raise ValueError(f"Unknown lifecycle stage {name!r}")

    def start_stage(self, name: str) -> None:
        record = self._stage_record(name)
        if self.run_status != "running":
            raise RuntimeError(
                f"Cannot start stage {name!r} after run became {self.run_status!r}"
            )
        if record["status"] != "pending":
            raise RuntimeError(
                f"Stage {name!r} cannot start from {record['status']!r}"
            )
        with self._durable_transition():
            record["status"] = "running"
            record["started_at"] = self._timestamp()
            self._stage_started_monotonic[name] = time.monotonic()

    def succeed_stage(self, name: str) -> None:
        record = self._stage_record(name)
        if record["status"] != "running":
            raise RuntimeError(
                f"Stage {name!r} cannot succeed from {record['status']!r}"
            )
        with self._durable_transition():
            record["status"] = "succeeded"
            record["duration_seconds"] = max(
                0.0, time.monotonic() - self._stage_started_monotonic.pop(name)
            )

    def skip_stage(self, name: str, reason: str) -> None:
        record = self._stage_record(name)
        if record["status"] != "pending":
            raise RuntimeError(
                f"Stage {name!r} cannot be skipped from {record['status']!r}"
            )
        with self._durable_transition():
            record["status"] = "skipped"
            record["reason"] = str(reason)

    def skip_pending_stages(self, reason: str) -> None:
        pending = [
            record for record in self.stage_records if record["status"] == "pending"
        ]
        if not pending:
            return
        with self._durable_transition():
            for record in pending:
                record["status"] = "skipped"
                record["reason"] = str(reason)

    def terminate_stage(
        self,
        name: str,
        error: BaseException,
    ) -> None:
        record = self._stage_record(name)
        if record["status"] != "running":
            raise RuntimeError(
                f"Stage {name!r} cannot terminate from {record['status']!r}"
            )

        with self._durable_transition():
            record["status"] = "failed"
            record["duration_seconds"] = max(
                0.0, time.monotonic() - self._stage_started_monotonic.pop(name)
            )
            record["error"] = self._error_record(error)
            current_index = self.stage_records.index(record)
            for later in self.stage_records[current_index + 1 :]:
                if later["status"] == "pending":
                    later["status"] = "skipped"
                    later["reason"] = "upstream_failure"

            self.run_status = "failed"
            self.run_ended_at = self._timestamp()
            self.run_duration_seconds = max(
                0.0, time.monotonic() - self._run_started_monotonic
            )
            self.run_error = self._error_record(error)

    def terminate_run(
        self,
        error: BaseException,
    ) -> None:
        if self.run_status != "running":
            self._notify_provenance()
            return
        running = next(
            (
                record
                for record in self.stage_records
                if record["status"] == "running"
            ),
            None,
        )
        if running is not None:
            self.terminate_stage(str(running["name"]), error)
            return

        with self._durable_transition():
            for record in self.stage_records:
                if record["status"] == "pending":
                    record["status"] = "skipped"
                    record["reason"] = "upstream_failure"
            self.run_status = "failed"
            self.run_ended_at = self._timestamp()
            self.run_duration_seconds = max(
                0.0, time.monotonic() - self._run_started_monotonic
            )
            self.run_error = self._error_record(error)

    def mark_run_succeeded(self) -> None:
        if self.run_status != "running":
            raise RuntimeError(f"Cannot succeed a {self.run_status!r} run")
        unfinished = [
            record["name"]
            for record in self.stage_records
            if record["status"] in {"pending", "running"}
        ]
        if unfinished:
            raise RuntimeError(f"Cannot succeed run with unfinished stages: {unfinished}")
        with self._durable_transition():
            self.run_status = "succeeded"
            self.run_ended_at = self._timestamp()
            self.run_duration_seconds = max(
                0.0, time.monotonic() - self._run_started_monotonic
            )

    def record_artifact(self, artifact: Dict[str, Any]) -> None:
        with self._durable_transition():
            self.artifact_records.append(dict(artifact))

    def record_local_refinement(self, applied: bool) -> None:
        if not isinstance(applied, bool):
            raise TypeError("applied must be a bool")
        if not applied or self.local_refinement_record["applied"]:
            return
        with self._durable_transition():
            self.local_refinement_record["applied"] = True

    def record_sr_allocation(self, evidence: Dict[str, Any]) -> None:
        with self._durable_transition():
            self.sr_allocation_records.append(copy.deepcopy(evidence))

    @property
    def io(self) -> Dict[str, Any]:
        return self.merged_config.get("io", {})

    @property
    def columns(self) -> Dict[str, Any]:
        return self.merged_config.get("columns", {})

    @property
    def compatibility_mode(self) -> bool:
        return bool(self.runtime["compatibility_mode"])

    @property
    def route(self) -> Dict[str, Any]:
        route = {
            "mode": self.runtime.get("mode"),
            "task": self.runtime.get("task"),
            "svc_kind": self.runtime.get("svc_kind"),
            "strategy": self.runtime.get("strategy"),
            "compatibility_mode": self.runtime.get("compatibility_mode"),
        }
        if self.runtime.get("mode") == "application":
            route["application_route"] = self.runtime.get("application_route")
            route["application_mode"] = self.runtime.get("application_mode")
        else:
            route["confounding"] = self.runtime.get("confounding")
        return route
