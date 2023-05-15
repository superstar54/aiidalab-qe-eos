import pathlib
from setuptools import setup, find_packages

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()


setup(
    name="aiidalab-qe-eos",
    version="0.0.1",
    description="A aiidalab-qe plugin to calculate the equation of state of a crystal.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/superstar54/aiidalab-qe-eos",
    author="Xing Wang",
    author_email="xingwang1991@gmail.com",
    classifiers=[],
    packages=find_packages(),
    entry_points={
        "aiidalab_qe.structure.importer": [
            "eos = aiidalab_qe_eos.structure:EOSStructure",
        ],
        "aiidalab_qe.property": [
            "eos = aiidalab_qe_eos:property",
        ],
    },
    install_requires=[

    ],
    package_data={},
    python_requires=">=3.6",
)
