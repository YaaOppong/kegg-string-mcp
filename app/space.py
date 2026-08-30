"""Hugging Face Space entry point. Spaces import a module and look for `demo`."""

from app.ui import build

demo = build()

if __name__ == "__main__":
    demo.launch()
