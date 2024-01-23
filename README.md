# iblatlas
Tools to manipulate hierarchical 3D representations of the mouse brain anatomy for electrophysiolgy experiments.
The tools are mainly using the the Allen CCF although some other atlases can be used.

**This repository uses minimal requirements, based of standard `matplotlib`, `numpy` and `scipy` libraries, to the exclusion
of more complex visualization tools such as `pyqt`.**

## Documentation
[Full package documentation](https://int-brain-lab.github.io/iblenv/_autosummary/iblatlas.html#module-iblatlas)

[Notebooks and examples](https://int-brain-lab.github.io/iblenv/notebooks_external/atlas_working_with_ibllib_atlas.html)

## Installation
`pip install iblatlas`

## Contributing
Changes are merged by pull requests.
Release checklist:
- [x] Update version in `iblatlas/__init__.py`
- [x] Update `CHANGELOG.md`
- [x] Create a pull request to the `main` branch on GitHub
- [x] Once the PR is merged, create a new tag and push the tag

Once a tag is pushed on main the package is uploaded to PyPI using GitHub Actions.
