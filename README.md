# iblatlas
Tools to manipulate hierarchical 3D representations of the mouse brain anatomy for electrophysiolgy experiments.
The tools are mainly using the Allen CCF although some other atlases can be used.

**This repository uses minimal requirements, based of standard `matplotlib`, `numpy` and `scipy` libraries, to the exclusion
of more complex visualization tools such as `pyqt`.**

## Documentation
[Full package documentation](https://int-brain-lab.github.io/iblenv/_autosummary/iblatlas.html#module-iblatlas)

[Notebooks and examples](https://int-brain-lab.github.io/iblenv/notebooks_external/atlas_working_with_ibllib_atlas.html)

## Installation
`pip install iblatlas`

For GUI features, you'll need to install with optional dependencies:
`pip install iblatlas[gui]`

## Contributing
Changes are merged by pull requests.
Release checklist:
- [ ] Update version in `iblatlas/__init__.py`
- [ ] Update `CHANGELOG.md`
- [ ] Merge changes to the `main` branch on GitHub
- [ ] Create a new tag on `main` and push it: `git tag <version> && git push origin <version>`
- [ ] Create a GitHub release from that tag (via the GitHub UI or `gh release create`)

Once a GitHub release is published the package is uploaded to PyPI automatically via GitHub Actions.
