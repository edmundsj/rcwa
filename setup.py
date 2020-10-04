import setuptools
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("LICENSE.txt", "r") as fh:
    license = fh.read()

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
        "rcwa": ['source/*.py', 'netlist/*', 'nkData/*.csv', 'docs/*'],
        },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
            'numpy>=1.14.5',
            'matplotlib>=2.2.0',
            'pandas>=0.24.0',
            'scipy>=1.2.2',
        ]
    license=license,
)
