# Installation

*CristalX* is easy to install; you can get started in minutes.


## Dependencies

*CristalX* is a written in pure Python, although it relies on packages that use other languages (mostly C and C++). However, it is not a problem from the user's perspective as no manual compilation is needed. Those who are interested can see the dependencies in the `environment.yml` file in the root directory.

Whether you install *CristalX* by *conda* or try it online, all but one dependency is installed. In a common workflow, that single missing dependency will not affect you. For details, see ...



## Try it online

If you just want to get a taste of *CristalX*, you can try it online without installing anything.

[Click here](https://mybinder.org/v2/gh/CsatiZoltan/CristalX/master) to jump to the root directory. Existing Jupyter notebooks are in the `notebooks/` directory. You can modify them and create new ones. If you directly want to open the main tutorial, which describes a real-world application, [click here](https://mybinder.org/v2/gh/CsatiZoltan/CristalX/master?filepath=notebooks%2Fexample_application.ipynb).


## Install it locally

If you want more control (debugging, inspecting variables, etc.) over *CristalX* or if you wish to contribute to the project, it is recommended to install it on your machine.



### Obtain the source

If you are a user of *CristalX*, the best is to get the latest release:

-  download it from [GitHub](https://github.com/CsatiZoltan/CristalX/releases)
-  or clone it with Git:

   ```bash
   git clone https://github.com/CsatiZoltan/CristalX.git
   cd CristalX
   git checkout v<version_number>
   ```

   where `<version_number>` is the version you want to use. E.g. if you want to use version 1.0.1, you need to type `git checkout v1.0.1`. See the [available tags](https://github.com/CsatiZoltan/CristalX/tags), corresponding to the published releases, for the possibilities.

   Note that in this case, you will be in a "detached HEAD" state, meaning that the HEAD does not point to a branch but to the specific tag. Any commit you make in this state will not be associated with a branch. Therefore, if you want to develop or contribute to *CristalX*, check out the *master* branch (see the next paragraph).

If you want to develop *CristalX* or simply want to have access to the latest features, you need to fetch the latest state:

-  download it from [GitHub](https://github.com/CsatiZoltan/CristalX/archive/master.zip)
-  or clone it with Git:

   ```bash
   git clone https://github.com/CsatiZoltan/CristalX.git
   ```



### Install with *conda*

All you need to have is the [*conda*](https://docs.conda.io/en/latest/) package manager. Open a terminal (or a command prompt if you are under Windows) in the root directory and type

```bash
conda env create -f environment.yml
```

*CristalX* has been installed to a separate environment, so you can safely work inside it. Activate that environment:

```bash
conda activate CristalX
```

Once you have finished working with *CristalX*, either close the terminal or type

```bash
conda deactivate
```
to return to your default *conda* environment.

If you want to uninstall *CristalX*, make sure that the `CristalX` conda environment is not active and then

```bash
conda env remove -n CristalX
```

When uninstalled, `conda env list` will not show it.



### Install with *pip*

If you do not have *conda* or you prefer *pip*, you can also install *CristalX* by typing

```bash
pip install -r requirements.txt
```

assuming that you have Python installed. It is highly recommended that you first create a new virtual environment to make sure you do not break your Python installation. However, an important component of *CristalX* (functions that rely on *PythonOCC*) will **not** be installed if you choose *pip*. The reason for this is that [*PythonOCC*](https://github.com/tpaviot/pythonocc-core) is not (yet) available on [PyPI](https://pypi.org/).
