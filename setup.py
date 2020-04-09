import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="GRLP",
    version="1.2.0",
    author="Andrew D. Wickert",
    author_email="awickert@umn.edu",
    description="Evolves gravel-bed river long profiles",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/awickert/GRLP",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Hydrology",
    ],
)
