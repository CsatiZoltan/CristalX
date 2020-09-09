# -*- coding: utf-8 -*-
"""
This module implements profiling facilities so that user code can be tested for speed.

Functions
---------
.. autosummary::
    :toctree: functions/

    profile

"""
import os
import codecs
import inspect
from contextlib import contextmanager

from grains import HAS_PYINSTRUMENT
if HAS_PYINSTRUMENT:
    from pyinstrument import Profiler
    from pyinstrument.renderers.html import HTMLRenderer
else:
    raise ImportError('This module requires the pyinstrument module.')


@contextmanager
def profile(output_format='html', output_filename=None, html_open=False):
    """Profiles a block of code with Pyinstrument.

    Intended to be used in a `with` statement.

    Parameters
    ----------
    output_format : {'html', 'text'}, optional
        Shows the result either as text or in an HTML file. If :code:`output_format = 'html'`,
        the file is saved according to the :code:`output_filename` parameter.
        The default is 'html'.
    output_filename : str, optional
        Only taken into account if :code:`output_format = 'html'`.
        If not given (default), the html output is saved to the same directory the caller resides.
        The name of the html file is the same as that of the caller.
    html_open : bool, optional
        Only taken into account if :code:`output_format = 'html'`.
        If True, the generated HTML file is opened in the default browser.
        The default is False.

    Yields
    ------

    Notes
    -----
    In the implementation, we move two levels up in the stack frame, one for exiting the context
    manager and one for exiting this generator. This assumes that :meth:`profile` was called as
    a context manager. As a good practice, provide the :code:`output_filename` input argument.

    Examples
    --------
    Measure the time needed to generate 1 million uniformly distributed random numbers.

    >>> import random
    >>> with profile('html') as p:
    ...    for _ in range(1000000):
    ...        rand_num = random.uniform(1, 2.2)

    """
    # If required, guess the path where the results will be saved
    if output_format == 'html' and not output_filename:
        caller = inspect.currentframe().f_back.f_back.f_locals['__file__']
        output_filename = os.path.splitext(caller)[0] + '.html'
    # Start the profiler
    profiler = Profiler()
    profiler.start()
    # Give back the execution to the caller function
    yield
    # Finish profiling and show the results
    profiler.stop()
    if output_format == 'html':
        if html_open:
            HTMLRenderer().open_in_browser(profiler.last_session, output_filename=output_filename)
        else:
            with codecs.open(output_filename, 'w', 'utf-8') as f:
                f.write(HTMLRenderer().render(profiler.last_session))
    elif output_format == 'text':
        print(profiler.output_text(unicode=True, color=True, show_all=True))
