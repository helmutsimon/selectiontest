import setuptools
from selectiontest.__init__ import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="selectiontest",
    version=__version__,
    #version="0.1.16",
    author="Helmut Simon",
    author_email="helmut.simon@anu.edu.au",
    description="A test for selective neutrality using relative likelihood",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/helmutsimon/SelectionTest",
    packages=["selectiontest"],
    package_dir={"selectiontest": "selectiontest"},
    install_requires=['click', 'numpy', 'scipy'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    entry_points='''
    [console_scripts]
        cli = stcli:selectiontestcli
        '''
)

