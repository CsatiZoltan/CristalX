# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))
from grains import HAS_PYINSTRUMENT, HAS_MED, HAS_OCCT


# -- Project information -----------------------------------------------------

project = 'CristalX'
copyright = '2020, Zoltan Csati'
author = 'Zoltan Csati'

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.autosummary',
    'sphinx.ext.inheritance_diagram',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinxcontrib.apidoc',
    'recommonmark',
    'sphinx_git',
    'sphinx_copybutton',
    'sphinx_toggleprompt',
    'hoverxref.extension',
    'sphinx_tabs.tabs',
    'matplotlib.sphinxext.plot_directive',
    'nbsphinx',
    'nbsphinx_link'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# Explicitly define the main file (https://stackoverflow.com/a/56448499/4892892)
# so that it works with older Sphinx versions as well
# https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-master_doc
master_doc = 'index'

# Napoleon settings
napoleon_google_docstring = False
napoleon_include_private_with_doc = True
napoleon_use_param = False
napoleon_use_rtype = False

# apidoc settings
apidoc_module_dir = '../'
apidoc_output_dir = '.'
apidoc_excluded_paths = ['../scripts/']
apidoc_separate_modules = True

# autodoc settings
autodoc_mock_imports = []
if not HAS_PYINSTRUMENT:
    autodoc_mock_imports.append('pyinstrument')
if not HAS_MED:
    autodoc_mock_imports.append('MEDLoader')
if not HAS_OCCT:
    autodoc_mock_imports.append('OCC')

# intersphinx settings
intersphinx_mapping = {'python': ('http://docs.python.org/3', None),
                       'numpy': ('https://numpy.org/doc/stable/', None),
                       'scipy': ('http://docs.scipy.org/doc/scipy/reference/', None),
                       'matplotlib': ('http://matplotlib.org', None),
                       'skimage': ('https://scikit-image.org/docs/stable/', None),
                       'sphinx': ('https://www.sphinx-doc.org/en/stable/', None),
                       'skan': ('https://jni.github.io/skan/', None)}

# toggleprompt settings
toggleprompt_offset_right = 25

# todo settings
todo_include_todos = True

# matplotlib plot_directive
plot_include_source = True
plot_html_show_source_link = True

# autosummary settings
autosummary_generate = True

# inheritance_diagram settings
# See the available settings at https://graphviz.org/doc/info/attrs.html
#inheritance_graph_attrs = dict(rankdir="LR", size='"6.0, 8.0"',
#                               fontsize=24, ratio='compress')
#inheritance_node_attrs = dict(shape='ellipse', fontsize=14, height=0.75,
#                              color='dodgerblue1', style='filled')
#inheritance_edge_attrs

# For conversion from markdown to html
import recommonmark.parser
from recommonmark.transform import AutoStructify

# Add type of source files
source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

def setup(app):
    app.add_config_value('recommonmark_config', {
            'url_resolver': lambda url: github_doc_root + url,
            'auto_toc_tree_section': 'Contents',
            'enable_inline_math': True,
            'enable_math': True,
            }, True)
    app.add_transform(AutoStructify)
    app.add_stylesheet('css/custom.css')


