from setuptools import find_packages, setup

setup(
    name="multimodal_dimelo",
    version="0.0.1alpha",
    packages=find_packages(),
    install_requires=[
        'dimelo @ git+https://github.com/streetslab/dimelo-toolkit.git',
    ],
)
