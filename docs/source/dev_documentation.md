# Documentation

Writing documentation is necessary when you contribute to this project, either by writing code or by changing/extending the external documents. In this section, we give tips how you can write them. This is an opinionated topic and it lays down how it is currently done. Feel free to suggest new ideas.

Whenever you write the documentation, first [test it locally](user_documentation.md#local-documentation) before pushing changes: Read the Docs spends about 1000 seconds to build the documentation.

It took me a long time to experiment with the Sphinx settings that provide the output what you can see in the rendered documentation. Some notes concerning these efforts can be found in the [Notes](notes) section. Another way to learn about Sphinx documentation is by reading the source of the existing documentation. On each page of the HTML documentation, the header contains a *View page source* hyperlink.



## Code documentation

As shown in the [figure](development), codes are written either as Python files or as Jupyter notebooks. With the proper extensions (see the `docs/source/conf.py` file), both of them are automatically included in the documentation by Sphinx. In what follows, we concentrate on documenting Python files.

The docstings are written in RST, following the [numpydoc](https://numpydoc.readthedocs.io/en/latest/format.html) style, using the [Napoleon](https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html) Sphinx extension.

Provide [doctests](https://docs.python.org/3/library/doctest.html) to demonstrate the use of the function you write and to provide minimal testing.

The line length is 100 characters, as defined in the `/.editorconfig` file.



## External documentation

The external documentation is written in reStructuredText (RST) and in Markdown. Sphinx, by default, uses RST, extending it with more capabilities. However, [recommonmark](https://github.com/readthedocs/recommonmark) can parse Markdown files and automatically convert them to RST at documentation build time. This is a great help as more people are used to the simple Markdown than to the more complex (and more capable) RST. Actually, the current document you read was also written in Markdown. Of course, it is perfectly fine if you write the documentation exclusively in RST.

Currently, the structure of the documentation is written in RST, and most of the external documentation in Markdown.