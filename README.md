# iblatlas
This repository contains tools to manipulate 3D representations of the mouse brain anatomy for electrophysiolgy experiments.
The tools are mainly using the the Allen Brain Atlas although some other atlases can be used and are written in Python.

One of the requirment of this repository is to use minimal requirements, based of standard matplotlib, numpy and scipy libraries, and excludes visualisation tool.

The documentation can be found here: https://int-brain-lab.github.io/iblenv/notebooks_external/atlas_working_with_ibllib_atlas.html

## Installation
`pip install iblatlas`


## Contributing
Changes are merged by pull requests.
Release checklist:
- [ ] Update version in `iblatlas/__init__.py`
- [ ] Update `CHANGELOG.md`
- [ ] Create a pull request to the `main` branch on GitHub

Once merged the package is uploaded to PyPI using GitHub Actions.
