# Documentation site

The source for the iblatlas documentation site, built with [Quarto](https://quarto.org) and
[quartodoc](https://machow.github.io/quartodoc/).

## Building the site

You need Quarto and Python 3.10 or newer. From the repository root:

```bash
uv venv --python 3.12 venv
uv pip install --python venv/bin/python -e . -r docs/requirements.txt
source venv/bin/activate

quartodoc build --config docs/_quarto.yml
quarto render docs/
```

Open `docs/_site/index.html`. Run `quartodoc build` before `quarto render`: it generates the
API pages under `docs/reference/api/`. Keep the `griffe<2` pin in `docs/requirements.txt`, as
quartodoc 0.11.1 crashes with griffe 2.

## Adding a page

Put the page in `docs/tutorials/`, `docs/how-to/` or `docs/explanation/` (the Concepts menu),
then add it to the `navbar` in `docs/_quarto.yml`. Citations go in the table on
`docs/citations.qmd`.

## Adding to the API reference

The reference lists only the symbols in `quartodoc.sections` in `docs/_quarto.yml`. Add a new
public symbol there, or record why it is left out in `docs/_exclusions.yml`. Do not edit
`docs/reference/api/`: each build overwrites it.

## Notebook figures

Quarto does not run notebooks when it renders the site; it shows the outputs saved in them.
To add or refresh figures, run the notebook and commit it with its outputs:

```bash
PYTHONWARNINGS=ignore venv/bin/jupyter nbconvert --to notebook --execute --inplace docs/how-to/<page>.ipynb
```

`PYTHONWARNINGS=ignore` keeps warnings, and the local paths they print, out of the page. The
first run downloads atlas data, so run it again once the data is cached to drop the download
progress bars from the outputs.

The figures on the Concepts pages are redrawn with
`PYTHONWARNINGS=ignore venv/bin/python docs/figures/make_figures.py`.
