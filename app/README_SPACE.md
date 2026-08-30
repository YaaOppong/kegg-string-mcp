---
title: Gene Annotation With Checked Citations
emoji: 🧬
colorFrom: indigo
colorTo: green
sdk: gradio
app_file: app/space.py
pinned: false
license: mit
short_description: Watch an LLM annotate a TB gene, then watch its citations get checked
---

# Deploying this demo to a Hugging Face Space

The Space is a git repo containing this repository's files.

```bash
pip install huggingface_hub
huggingface-cli login
huggingface-cli repo create gene-annotation-checked-citations --type space --space_sdk gradio

git remote add space https://huggingface.co/spaces/<your-username>/gene-annotation-checked-citations
cp app/README_SPACE.md README.md   # the Space reads its config from the root README
git push space main
```

Two files must sit at the **repository root**, because that is where Spaces looks:

* `README.md` carrying the front matter above (copy it from `app/README_SPACE.md`).
* `requirements.txt` — already present at the root for this reason; `app/requirements.txt`
  is a copy kept beside the app for readability and is **not** the one Spaces reads.

Everything else — `app/`, `demo/runs/`, `src/` — is used as-is.

If you would rather not run a server at all, the same demo is published as a static
GitHub Pages build; see `demo/build_pages.py`.

No secrets are required. The demo makes no API calls: it replays runs committed in
`demo/runs/` and re-runs the validator over them.
