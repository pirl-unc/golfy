# Releasing MHCnames

This document explains what do once your [Pull Request](https://www.atlassian.com/git/tutorials/making-a-pull-request/) has been reviewed and all final changes applied. Now you're ready merge your branch into master and release it to the world:

1. Bump the [version](http://semver.org/) on pyproject.toml as part of the PR you want to release.
2. Merge your branch into master.
3. Be sure to have the following tools installed: `build`, `twine`, and `conda-build`
4. Run `deploy.sh` 

