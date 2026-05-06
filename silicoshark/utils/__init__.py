"""silicoshark utility modules.

Currently exports the ProgressReporter pattern used across the
research-project ecosystem (see `~/.claude/CLAUDE.md` §Experiment
monitoring). Copied verbatim from `~/repo/workflow/lib/progress.py`
into the package so silicoshark has no path-dependent imports.
"""
from .progress import ProgressReporter

__all__ = ['ProgressReporter']
