"""Structured progress reporting for long-running experiments.

Writes a JSON progress file that can be read by monitor.sh or Claude Code
to determine what phase an experiment is in, how far through, and ETA.

Writes are throttled (default 2 s) to avoid I/O overhead when callbacks
fire once per corpus text.  Phase transitions and completion/failure always
write immediately.  Writes are atomic via os.replace().
"""

from __future__ import annotations

import json
import os
import time
import traceback
from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path


class ProgressReporter:
    """Report experiment progress to a JSON file.

    Usage::

        with ProgressReporter(path, "08_statistical_framework", total_phases=5) as progress:
            progress.begin_phase("load_texts", "Loading primary texts")
            load_texts()

            progress.begin_phase("build_corpus", "Building corpus", total_steps=1779)
            for i, text in enumerate(corpus, 1):
                process(text)
                progress.step(i, f"Processing text {i}/1779")

            progress.complete()

    If an exception propagates out of the ``with`` block the status is
    automatically set to ``failed`` with the error message.
    """

    def __init__(
        self,
        path: str | Path,
        experiment: str,
        total_phases: int,
        *,
        throttle_s: float = 2.0,
    ) -> None:
        self._path = Path(path)
        self._experiment = experiment
        self._total_phases = total_phases
        self._throttle_s = throttle_s

        self._started_at = datetime.now(timezone.utc)
        self._status = "running"
        self._phase: str | None = None
        self._phase_description: str | None = None
        self._phase_number = 0
        self._total_steps: int | None = None
        self._current_step: int | None = None
        self._message: str | None = None
        self._phase_start: float | None = None

        self._last_write_time = 0.0
        self._completed = False

    # -- context manager ------------------------------------------------

    def __enter__(self) -> ProgressReporter:
        self._write(force=True)
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: object,
    ) -> None:
        if exc_type is not None:
            self._status = "failed"
            self._message = f"{exc_type.__name__}: {exc_val}"
            self._write(force=True)
        elif not self._completed:
            self.complete()

    # -- public API -----------------------------------------------------

    def begin_phase(
        self,
        name: str,
        description: str,
        total_steps: int | None = None,
    ) -> None:
        """Start a new phase.  Always writes to disk immediately."""
        self._phase_number += 1
        self._phase = name
        self._phase_description = description
        self._total_steps = total_steps
        self._current_step = None
        self._message = None
        self._phase_start = time.monotonic()
        self._write(force=True)

    def step(self, current: int, message: str = "") -> None:
        """Update step progress within the current phase.

        Throttled: writes at most every *throttle_s* seconds.
        """
        self._current_step = current
        if message:
            self._message = message
        if self._should_write():
            self._write()

    def step_callback(self) -> Callable[[int, int], None]:
        """Return a ``(current, total)`` callback compatible with
        ``LemmatiserPipeline`` and ``SharedDistinctiveAnalyser``.

        On first call, sets *total_steps* from the ``total`` argument.
        """
        def _callback(current: int, total: int) -> None:
            if self._total_steps is None:
                self._total_steps = total
            self.step(current, f"Processing text {current}/{total}")

        return _callback

    def complete(self) -> None:
        """Mark the experiment as completed.  Always writes immediately."""
        self._status = "completed"
        self._completed = True
        self._write(force=True)

    # -- internals ------------------------------------------------------

    def _should_write(self) -> bool:
        return (time.monotonic() - self._last_write_time) >= self._throttle_s

    def _build_payload(self) -> dict:
        now = datetime.now(timezone.utc)
        elapsed = (now - self._started_at).total_seconds()

        payload: dict = {
            "experiment": self._experiment,
            "status": self._status,
            "started_at": self._started_at.isoformat(),
            "updated_at": now.isoformat(),
            "elapsed_s": round(elapsed, 1),
            "phase": self._phase,
            "phase_description": self._phase_description,
            "phase_number": self._phase_number,
            "total_phases": self._total_phases,
        }

        if (
            self._current_step is not None
            and self._total_steps is not None
            and self._total_steps > 0
        ):
            pct = self._current_step / self._total_steps * 100
            payload["step"] = self._current_step
            payload["total_steps"] = self._total_steps
            payload["step_pct"] = round(pct, 1)

            if self._phase_start is not None:
                phase_elapsed = time.monotonic() - self._phase_start
                if phase_elapsed > 0 and self._current_step > 0:
                    rate = self._current_step / phase_elapsed
                    remaining = self._total_steps - self._current_step
                    eta = remaining / rate if rate > 0 else 0
                    payload["rate"] = round(rate, 2)
                    payload["eta_s"] = round(eta)

        if self._message is not None:
            payload["message"] = self._message

        return payload

    def _write(self, *, force: bool = False) -> None:
        if not force and not self._should_write():
            return

        payload = self._build_payload()
        self._path.parent.mkdir(parents=True, exist_ok=True)

        tmp_path = self._path.with_name(".progress.json.tmp")
        tmp_path.write_text(json.dumps(payload, indent=2) + "\n")
        os.replace(tmp_path, self._path)

        self._last_write_time = time.monotonic()
