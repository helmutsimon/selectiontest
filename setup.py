import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="selectiontest-helmutsimon", 
    version="0.0.1",
    author="Helmut Simon",
    author_email="helmut.simon@anu.edu.au",
    description="A test for selective neutrality using relative likelihood",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/helmutsimon/SelectionTest",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
