# How to use the codes

If you want to get to know the classes and functions, the best is to check the [API documentation](api_doc). Most components are amply documented and contain one or more examples. If you want to see them in action, check out our [workflow](workflow) for a real-world application. If you are curious about the underlying algorithms, consult with the [Algorithms](algorithms) section.



## Tutorials

The tutorials are available in two forms: Python scripts and Jupyter notebooks.



### Python scripts

The scripts are intended to be run in batch mode. They need to be run in a specific order as the subsequent scripts rely on the output of the previous ones. These scripts are useful if you are already a bit familiar with *CristalX* or if you wish to adapt the scripts to your needs. Indeed, if you intend to extend the package, the safest way is to copy the relevant scripts and modify them. The scripts are named as `run_moduleName`, where `moduleName` is the name of the module in the `grains/` directory that the script mainly relies on. In some way, `run_moduleName` acts as a demonstration of what the `moduleName` module is used for. The scripts are meant to be run as modules, i.e. navigate to the root of the project and type

```python
python -m scripts.moduleName
```



### Jupyter notebooks

Almost the same workflow is available as a single Jupyter notebook, located at the `nootebooks/` directory. this allows you to interactively discover *CristalX* through an example. You can easily rerun parts of the code to see the effects of the parameters, and the rich output is embedded into the same document. If you are a novice Python user or just want to have a taste of *CristalX*, notebooks are the recommended way to get started.





## Undo changes

When you experiment with *CristalX*, you will probably change parameter values. As the scripts communicate by reading and writing data, the original data that come with *CristalX* will be overwritten. There are multiple ways to undo the changes.



### In Binder

As written in the [Installation instructions](installation.md#try-it-online), you can try *CristalX* online without the need to install anything on your computer. Then a separate virtual environment is created. Whatever changes you make there, they will not influence your local installation (if you have) or the data on the GitHub repository.



### In a local installation

If you downloaded *CristalX* from GitHub, you can simply replace the new files with the original ones. In case you cloned *CristalX* with Git, you can easily discard the changes you make.