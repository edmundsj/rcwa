import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="RCWATest", # Replace with your own username
    version="0.0.2",
    author="Jordan Edmunds",
    author_email="jordan.e@berkeley.edu",
    description="Python Implementation of Rigorous Coupled Wave Analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/edmundsj/RCWA",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
