import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="GRLP",
    version="2.0.0-alpha",
    author="Andrew D. Wickert",
    author_email="awickert@umn.edu",
    description="Evolves gravel-bed river long profiles",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MNiMORPH/GRLP",
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        ],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Hydrology",
        "Intended Audience :: Science/Research",
    ],
    keywords='fluvial geomorphology sediment transport landscape evolution',
    project_urls={
        'Model page': 'https://csdms.colorado.edu/wiki/Model:GRLP',
    },
)
