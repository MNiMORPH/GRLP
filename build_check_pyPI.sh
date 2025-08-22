#! /bin/sh

# Clean old files in order to not get "400 Client Error: File already exists."
rm dist/*.tar.gz
rm dist/*.whl

# Build the new distribution
python3 -m build
# python3 setup.py sdist bdist_wheel # old way

# Check the new distribution for a successful build
twine check dist/*
