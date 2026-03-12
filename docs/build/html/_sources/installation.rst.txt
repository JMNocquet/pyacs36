Installation
============

It is recommended to install **PYACS** in a dedicated Python environment.

Recommended setup using mamba
-----------------------------

The recommended environment manager is **mamba**, which provides a faster
dependency resolver than conda.

Install Miniforge / mamba::

    curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh
    bash Miniforge3-$(uname)-$(uname -m).sh

Open a new terminal after installation.

Download the environment configuration file from:

https://github.com/JMNocquet/pyacs36/tree/master/environment.yaml

Create the PYACS environment::

    mamba env create -f environment.yaml
    mamba activate pyacs


Installation from an existing Anaconda environment
--------------------------------------------------

If you already use **Anaconda** or **Miniconda**, you can install mamba
in the base environment::

    conda install -n base -c conda-forge mamba
    mamba env create -f environment.yaml

Download the environment configuration file from:

https://github.com/JMNocquet/pyacs36/tree/master/environment.yaml

Note
----

The standard **conda** dependency solver may be significantly slower
than **mamba**.


Getting PYACS
-------------

Download the latest stable version from:

https://github.com/JMNocquet/pyacs36/tree/master/dist

Install using pip::

    pip install pyacs-X.XX.XX.tar.gz


Development installation
------------------------

If you plan to modify the code::

    tar xvfz pyacs-X.XX.XX.tar.gz
    cd pyacs-X.XX.XX
    pip install .


Alternatively, clone the full repository::

    git clone https://github.com/JMNocquet/pyacs36.git
    cd pyacs36
    pip install .


Note
----

The latest development version of the repository may include experimental
changes. Tagged releases are expected to be more stable.


Running tests
-------------

From the directory containing the ``pyacs`` package::

    pytest pyacs/tests


Interactive use
---------------

``ipyacs.py`` is a convenient script that loads the main PYACS libraries and
automatically loads time series located in the current directory (if available).
It allows interactive time-series visualization and analysis using IPython.

You may define a shell alias::

    alias ipyacs='mamba activate pyacs && ipython $(which ipyacs.py) -i'


Working with Jupyter notebooks
------------------------------

Time-series analysis is often conveniently performed using a Jupyter notebook.

Start a notebook, select the kernel corresponding to your PYACS environment,
and run::

    import pyacs
    print(pyacs.__version__)

    import numpy as np
    from pyacs.gts.Sgts import Sgts
    import pyacs.lib.astrotime as at
    from datetime import datetime

    pyacs.verbose("SILENT")

    # data directory containing time series
    ts_dir = "your_time_series_path_directory"

    ts = Sgts(ts_dir, verbose=False)
    print(f"{ts.n()} time series loaded")
    ts.info()


Building a PYACS distribution
-----------------------------

Advanced users can build their own distribution using::

    python -m build

The source and wheel distributions will be generated in the ``dist/`` directory.


Documentation
-------------

An HTML documentation is available online:

https://jmnocquet.github.io/pyacs_docs/pyacs

Alternatively, the documentation can be generated locally::

    ./make_pyacs_doc_html_sphinx.sh