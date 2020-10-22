# Development

Intro ...



## General workflow

The development workflow can be followed in the following figure.



![Alt text](images/development_workflow.svg)

1. You write the Python source codes (.py) and the Jupyter notebooks (.ipynb) on your local machine. If you want, you can [create the documentation](documentation) locally using Sphinx.
2. When a logical unit has been finished, you commit the changes with `git push`. This will upload the new version of the modified files to a remote repository (currently GitHub). If you edit files in the remote repository and you want those changes to be present in your local copy, you can use `git pull`. For more details on Git, read the [manual](https://git-scm.com/docs/git-pull).
3. Several commit hooks are attached to the remote repository. When they detect a change, certain actions are activated. One the one hand, the online version of the documentation, hosted on [Read the Docs](https://cristalx.readthedocs.io/) will be updated. On the other hand, static analyzers will reanalyze the new version of the code. Some of them may create a pull request based on their recommendations.



## Profiling


I use [Pyinstrument](https://github.com/joerick/pyinstrument) for profiling the code. It comes with the *CristalX* installation.

Otherwise, you can install it with `pip install pyinstrument` or by *conda* after the *conda-forge* channel has been [activated](https://conda-forge.org/docs/user/introduction.html#how-can-i-install-packages-from-conda-forge).

```bash
conda install -c conda-forge pyinstrument
```

As Pyinstrument has [no dependencies](https://github.com/joerick/pyinstrument/issues/102), you can safely install it to your current environment. If you do not want to take a risk, create a new environment. With *conda*, you can do e.g.

```bash
conda create -n pyinstrument python=3.7 scipy matplotlib
conda activate pyinstrument
conda install -c conda-forge pyinstrument
```



The *profiling* module provides a wrapper around Pyinstrument. Put the code you want to profile in the `profile` context manager, e.g.

```python
>>> import random
>>> from grains.profiling import profile
>>> with profile('html') as p:
...    for _ in range(1000000):
...        rand_num = random.uniform(1, 2.2)
```

For more details see the [API reference](profiling).