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
cp app/README_SPACE.md README_SPACE_HF.md   # the Space needs this front matter as its README.md
git push space main
```

The Space needs `README.md` to carry the front matter above, so on the Space branch
rename `README_SPACE.md` to `README.md`. Everything else — `app/`, `demo/runs/`,
`src/` — is used as-is.

No secrets are required. The demo makes no API calls: it replays runs committed in
`demo/runs/` and re-runs the validator over them.
