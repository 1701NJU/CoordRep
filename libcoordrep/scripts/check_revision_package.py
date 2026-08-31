#!/usr/bin/env python3
"""Backward-compatible entry point for the current revision checker."""

from check_jacs_revision_release import main


if __name__ == "__main__":
    raise SystemExit(main())
