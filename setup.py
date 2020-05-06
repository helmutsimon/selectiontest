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
    py_modules=['cli'],
    install_requires=['Click', 'numpy', 'scipy'],
    #scripts=['bin/calc_TajD', 'bin/sample_wf_distribution', 'bin/sample_uniform_distribution', 'bin/test_neutrality', 'bin/compute_threshold'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    entry_points='''
    [console_scripts]
        cli = cli:selectiontestcli
        '''
)
