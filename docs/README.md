Arduino AD7766 Communication Library
====================

Documentation for a library for communicating with the AD7766.

Building Documentation Locally
--------------------------------

**1. Install the toolchain**
Install [Anaconda (miniconda)](https://docs.conda.io/en/latest/miniconda.html). This fantastic little package manager was initially written for python, but has expanded to a huge conglomeration of open-source software, and is by far the easiest way to install it.

Install [Sphinx](https://www.sphinx-doc.org/en/master/). This is the core tool used to generate all the documentation for this project, it's job is to take the source code and a few other files and compile them into HTML. You must first have Sphinx installed on your local system, The install instructions can be found [here](https://www.sphinx-doc.org/en/master/usage/installation.html) if you do not already have Sphinx installed. As previously mentioned, I recommend installing Sphinx with [Anaconda](https://www.anaconda.com/), as it makes the whole process pretty trivial. In this case, you only need to run:
```
conda install sphinx
```

Install [Sphinx_rtd_theme](https://pypi.org/project/sphinx-rtd-theme/) is the theme used for this documentation (contains all the styling and everything). I prefer
```
conda install -c anaconda sphinx_rtd_theme
```

Install [Doxygen](https://www.doxygen.nl/index.html). Old as time itself, this is the intermediary used to generate the documentation from the source code. Unfortanutely, while Sphinx is truly wonderful, it doesn't (as of this writing) directly support C++. Doxygen is a 
```
conda install -c conda-forge doxygen
```

Install [Breathe](https://pypi.org/project/breathe/). This is the link between doxygen and Sphinx. It takes the XML output from Doxygen and pipes it into Sphinx to actually generate the documentation. Strictly, it is an 'extension' for Sphinx.
```
conda install -c conda-forge breathe
```

Install [Exhale](https://pypi.org/project/exhale/). A neat little package that makes Breathe much more usable. In the future I hope this is included in the breathe package itself. This is the only package unavailable via anaconda (currently). Fortunately, this only means one extra command is needed

```
conda install pip
pip install exhale

```

**2. Navigate to the 'docs' directory and run make command**
Make sure you are currently in the 'docs' subdirectory of the project (where this README.md file is located) and run:
```
make html
```

**3. Open the generated index.html file**
The HTML files will be located in a new directory in the 'docs' directory called '\_build'. Navigate to this directory, open the index.html file, and there is the documentation!

Contribute
-------------
To contribute to this documentation or this project, fork [this project](https://github.com/edmundsj/pythonAD7766) on github, make your changes, and initiate a pull request. There's a great tutorial on how to actually do that [here](https://blog.scottlowe.org/2015/01/27/using-fork-branch-git-workflow/).

(Developers) Setting up this Documentation
---------------------------------------------
I created this documentation from [this tutorial](https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/), just without "cmake". Here I describe briefly how to do that, for my reference and for others that might find it useful:

**For projects requiring C++** (if this is not the case, skip to the projects requiring python section):
Once all the dependencies are installed, if you want to set up the documentation from scratch, you will need to set up the appropriate source files. First, build a Doxyfile by running the following command in this directory:
```
doxygen -g
```
This generates a new doxygen file. Since I like commenting stuff with javadoc, I will modify the JAVADOC\_BANNER and JAVADOC\_AUTOBRIEF lines from "NO" to "YES". You also *must* change the GENERATE\_XML line from "NO" to "YES". This should be all we need to set up doxygen. Next, we move on to sphinx.

First, run Sphinx's quickstart utility:
```
sphinx-quickstart
```
Enter the project name (this one was "Arduino AD7766 SCPI"). This generates a couple default files, including conf.py and index.rst, as well as a Makefile and some additional directories. To get this documentation to work with [readthedocs](https://readthedocs.org/) (RTD), there are a couple of tweaks we have to make to conf.py. First, we need to add a line (I put mine right under the "author" line), which will tell RTD where to look for our homepage:
```
master_doc = 'index'
```

Next, we need to add the extensions (including doxygen) that we need to Sphinx. There should be a line with an empty extensions list. Change it to:
```
extensions = [
	'breathe',
	'exhale',
	'sphinx_rtd_theme',
]
```
This will add the rtd theme we installed previously, as well as the "breathe" and "exhale" extensions. We also need to add a few lines to setup the exhale and breathe extensions:

```
# Setup the breathe extension
breathe_projects = {
    project: "./doxyoutput/xml"
}
breathe_default_project = project

# Setup the exhale extension
exhale_args = {
    # These arguments are required
    "containmentFolder":     "./api",
    "rootFileName":          "library_root.rst",
    "rootFileTitle":         "Library API",
    "doxygenStripFromPath":  "..",
    # Suggested optional arguments
    "createTreeView":        True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin":    "INPUT = ../src"
}

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'
```

To actually use the pretty RTD theme we previously installed, we need to change the default html theme at the bottom of the page:
```
html_theme = 'sphinx_rtd_theme'
```

And, to insert all the generated documentation indo our index file, we need only add the following in our index.rst file (the newline between the contents and the api line is important! Indentation is also critical):
```
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api/library_root
```
Finally, if you want to upload documentation to the web via RTD (which I do), create a new project on [RTD](https://readthedocs.org/) with the git repository. Make sure the repository is public. Then, go to "Settings -> Admin -> Integrations", and make sure the github webhook is there.

After that, you need to install the plugins on the RTD server. Fortunately (after you figure out how to do it) this is relatively simple: it involves putting two files in the root directory *of the github repository*: environment.yml and readthedocs.yml. The .readthedocs.yml should contain:
```
version: 2

# Build documentation in the docs/ directory with Sphinx

conda:
  environment: environment.yml
```
This tells our RTD server to use conda as our package manager, and points to yet *another* file, "environment.yml", which is a configuration file for conda. This file should have the following:
```
name: base
channels:
  - conda-forge
  - defaults
dependencies:
  - breathe
  - doxygen
  - pip
  - sphinx
  - sphinx_rtd_theme
  - sphinxcontrib-applehelp=1.0.2=py_0
  - sphinxcontrib-devhelp=1.0.2=py_0
  - sphinxcontrib-htmlhelp=1.0.3=py_0
  - sphinxcontrib-jsmath=1.0.1=py_0
  - sphinxcontrib-qthelp=1.0.3=py_0
  - sphinxcontrib-serializinghtml=1.1.4=py_0
  - pip:
    - exhale
```
Upload this to your repository, and you should be done. Or just clone this repository. Much less painful. If the project is a python project, make sure the environment.yml file contains all the dependencies of the modules documented in that project.
