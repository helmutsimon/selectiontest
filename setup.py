import setuptools
from selectiontest.selectiontest import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="selectiontest-helmutsimon", 
    version=__version__,
    author="Helmut Simon",
    author_email="helmut.simon@anu.edu.au",
    description="A test for selective neutrality using relative likelihood",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/helmutsimon/SelectionTest",
    packages=setuptools.find_packages(),
    scripts=['bin/calc_TajD', 'bin/generate_uniform_variates', 'bin/generate_wf_variates', 'bin/test_neutrality', compute_threshold],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
