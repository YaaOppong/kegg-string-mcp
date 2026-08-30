"""Entry point for the static browser build. `demo/build_pages.py` rewrites the
imports in ui.py and replay.py so they resolve flat under Pyodide."""

from app.ui import build

demo = build()
