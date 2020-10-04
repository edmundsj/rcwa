import setuptools
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rcwa",
    version="0.1." + str(os.environ['GITHUB_RUN_NUMBER']),
    author="Jordan Edmunds",
    author_email="jordan.e@berkeley.edu",
    description="Python Implementation of Rigorous Coupled Wave Analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/edmundsj/RCWA",
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={
        "rcwa": ['source/*.py', 'netlist/*', 'nkData/*.csv', 'docs/*', 'examples/*',
            'test/matrixDataOblique/freeSpace/*',
            'test/matrixDataOblique/layer1/*',
            'test/matrixDataOblique/layer2/*',
            'test/matrixDataOblique/reflectionRegion/*',
            'test/matrixDataOblique/transmissionRegion/*',
            'test/matrixDataNormal/freeSpace/*',
            'test/matrixDataNormal/layer1/*',
            'test/matrixDataNormal/layer2/*',
            'test/matrixDataNormal/reflectionRegion/*',
            'test/matrixDataNormal/transmissionRegion/*',
            'test/matrixDataNormal/*',
            'test/matrixDataOblique/*',
            'test/netlists/*',
            'test/*.csv'],
        },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
            'numpy>=1.14.5',
            'matplotlib>=2.0.0',
            'pandas>=0.24.0',
            'scipy>=1.2.2',
        ],
    license="MIT",
)
