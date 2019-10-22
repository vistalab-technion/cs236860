---
title: Working locally with the tutorial notebooks
permalink: tutorials/working-locally
toc: true
toc_label: Contents
toc_sticky: true
date: 2019-10-22
---

This document will help you install the required environment for working locally with the tutorial notebooks.

## General

The tutorial notebooks are implemented using a platform called Jupyter notebooks.
[Jupyter](http://jupyter.org/) is a widely-used tool in the machine-learning
ecosystem which allows us to create interactive notebooks containing live code,
equations and text.

To install and manage all the necessary packages and dependencies for the
tutorials, we use [conda](https://conda.io), a popular package-manager for
python.  The tutorials come with an `environment.yml` file which defines
what third-party libraries we depend on. Conda will use this file to create
a virtual environment for you. This virtual environment includes python and
all other packages and tools we specified, separated from any preexisting
python installation you may have. Detailed installation instructions are below.


### Obtaining the tutorial code

The tutorials will be made available on the
[course repository](https://github.com/vistalab-technion/cs236860-tutorials).
You can either download a zip file of the repository or, preferably, use
git to clone it.

In case we need to update the tutorials or make corrections, the repository
will be updated and notice will be given.

### Repository structure

The repository's root directory contains the following files and folders:

- `tutN` where `N` is the tutorial number: a directory containing the tutorial notebook
  ``tutN.ipynb`` and extra needed files.
- `semesters`: a directory containing tutorials from previous semesters.
- `environment.yml`: A file for `conda`, specifying the third-party packages it
  should install into the virtual environment it creates.

## Environment set-up

1. Install the python3 version of [miniconda](https://conda.io/miniconda.html).
   Follow the [installation instructions](https://conda.io/docs/user-guide/install/index.html)
   for your platform.

   For example, on linux you should do:
   ```shell
   curl -fsSLO https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   # Accept EULA
   # Install in default directory
   # Select no for editing .bashrc

   # Update your bashrc like so:
   echo "source $HOME/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc
   ```

   On macOS it's similar but with a different script URL
   ```shell
   curl -fsSLO https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh 
   bash Miniconda3-latest-MacOSX-x86_64.sh
   # Rest is the same
   ```

2. Use conda to create a virtual environment for the assignment.
   From the assignment's root directory, run

   ```shell
   conda env create -f environment.yml
   ```

   This will install all the necessary packages into a new conda virtual
   environment named `cs236860-tutorials`.

3. Activate the new environment by running

   ```shell
   conda activate cs236860-tutorials
   ```

   *Activating* an environment simply means that the path to it's python binaries
   (and packages) is placed at the beginning of your `$PATH` shell variable.
   Therefore, running programs installed into the conda env (e.g. `python`) will
   run the version from the env since it appears in the `$PATH` before any other
   installed version.

   To check what conda environments you have and which is active, run

   ```shell
   conda env list
   ```

   or, you can run `which python` and you should see the python binary is in a
   subfolder of `~/miniconda3/envs/cs236860-tutorials/`.

   You can find more useful info about conda environments
   [here](https://conda.io/docs/user-guide/tasks/manage-environments.html).

Notes: 

- You only need to do steps 1 and 2 above once.
  However, the third-party package dependencies (in the `environment.yml` file)
  might slightly change when we update the repository. To make sure you have
  the correct versions run
  ```shell
  conda env update
  ```
  from the repository root directory every time a new tutorial is published.

- Always make sure the correct environment is active. It will revert to it's
  default each new terminal session. If you want to change the default env you
  can add a `conda activate` in your `~/.bashrc`.

## Running Jupyter

Make sure that the active conda environment is `cs236860-tutorials`, and run

```shell
jupyter lab
```

This will start a [jupyter lab](https://jupyterlab.readthedocs.io/en/stable/)
server and open your browser at the local server's url. You can now start working.
Open the tutorial notebook (`tutN.ipynb`) and run the code.

If you're new to jupyter notebooks, you can get started by reading the
[UI guide](https://jupyter-notebook.readthedocs.io/en/stable/notebook.html#notebook-user-interface)
and also about how to use notebooks in
[JupyterLab](https://jupyterlab.readthedocs.io/en/latest/user/notebook.html).

Note that if you are familiar with and prefer the regular `jupyter notebook` you
can use that instead of `jupyter lab`.
