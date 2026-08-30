"""Hugging Face Space entry point.

Spaces import a module and look for a top-level `demo`; keeping that here rather
than in app.py leaves app.py runnable locally with `python -m app.app`.
"""

from app.app import build

demo = build()

if __name__ == "__main__":
    demo.launch()
