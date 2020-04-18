# Scripts

Demonstrating typical usage



## How to run

The scripts can be run as modules from the root directory, e.g.

```python
python -m scripts.run_segmentation
```

possibly adding the `-d` flag to debug.
You can also run them from within the IPython console:

```python
>>> %run -m scripts.run_segmentation
```

or if you want to enter debug mode, specify the `-d` flag.



## Purpose

Most scripts read and write text files when executed. The files are read from and written to */scripts/data*. The nature of the problem is such that full automation is not possible. Some of the scripts rely on files already created by another script or an external program. Therefore, a sample data set is provided. When the scripts are executed, the original data will not be overwritten. The purpose of the scripts is to show how to use the **grains** package, by performing a typical workflow. It is the role of the tests (residing in the */tests* directory) to thoroughly examine whether the functions work correctly.

| Script name           | Purpose                                                      |
| --------------------- | ------------------------------------------------------------ |
| *run_segmentation.py* | Creates a segmented image from the image of a microstructure.|
| *run_analysis.py*     | Performs measurements on segmented images.                   |
| *run_abaqus.py*       | Mainly material management for Abaqus.                       |
| *run_simulation.py*   | Set up tools for an elastoplastic simulation.                |
