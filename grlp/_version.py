# Single source of truth for the package version, read by pyproject.toml via
# [tool.setuptools.dynamic] (version = {attr = "grlp._version.__version__"})
# and re-exported from grlp/__init__.py.
#
# When bumping, update the version in the other (static) places that external
# tools read and cannot derive from here:
#   - CITATION.cff  (version: and date-released:)  -- GitHub / Zenodo
#   - README.md     (the "(Version x.y.z)" citation line)
#   - CHANGELOG.md  (roll [Unreleased] -> [x.y.z] with date and compare link)
__version__ = "2.1.0"
