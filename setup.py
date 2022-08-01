import setuptools
import os
import glob

nk_data = glob.glob('rcwa/nkData/**', recursive=True)
test_data = glob.glob('rcwa/test/**', recursive=True)
example_data = glob.glob('rcwa/examples/*')
source_data = glob.glob('rcwa/source/*.py')
total_data = nk_data + test_data + example_data + source_data
package_data = [x.strip('rcwa/') for x in total_data]


with open("README.md", "r") as fh:
	long_description = fh.read()

	setuptools.setup(
		name="rcwa",
		version="1.0." + str(5 +  int(os.environ['GITHUB_RUN_NUMBER'])),
		author="Jordan Edmunds",
		author_email="jordan.e@berkeley.edu",
		description="Python Implementation of Rigorous Coupled Wave Analysis",
		long_description=long_description,
		long_description_content_type="text/markdown",
		url="https://github.com/edmundsj/RCWA",
		packages=setuptools.find_packages(),
		include_package_data=True,
		package_data={
			"rcwa": package_data
			},
		classifiers=[
			"Programming Language :: Python :: 3",
			"License :: OSI Approved :: MIT License",
			"Operating System :: OS Independent",
		],
		python_requires='>=3.6',
		install_requires=[
			'numpy>=1.20.0',
			'matplotlib>=2.0.0',
			'pandas>=0.24.0',
			'scipy>=1.2.2',
			'pyyaml>=5.0.0',
			'pytest>6.2.2',
			'progressbar2',
			'autograd',
		],
	license="MIT",
	)
