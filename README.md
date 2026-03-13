# PYACS

![Python](https://img.shields.io/badge/python-3.8%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)
[![Documentation](https://img.shields.io/badge/docs-online-blue)](https://jmnocquet.github.io/pyacs_docs/pyacs)

**PYACS** is a collection of Python modules and scripts for analyzing and modeling geodetic data, with a focus on tectonic applications.


Python ≥ 3.8 (current version tested with Python 3.12)

Documentation:  
https://jmnocquet.github.io/pyacs36

Source code:  
https://github.com/JMNocquet/pyacs36

---

# Overview

PYACS includes four main components:

- **Core libraries**  
  A set of Python modules handling coordinates, time representations, estimators, and format conversions.

- **GNSS time series generation**  
  The script `pyacs_make_time_series.py` expresses free GNSS solutions in a single reference frame and produces time series.

- **Geodetic Time Series module (`Gts`)**  
  Provides tools for analyzing individual or collections of time series.

- **Velocity field analysis tools**  
  Tools for GNSS velocity field analysis including Euler pole estimation and basic strain rate analysis.

---

# Related software

PYACS also serves as a core module for several related software tools:

- **pygeca** – processing of large GNSS networks using GAMIT, including in HPC environment.
- **pyeblock** – elastic block modelling
- **pyeq / pyaks** – full time-dependent slip inversion on faults

---

# Typical workflow

A typical PYACS workflow for GNSS time series analysis may involve:

1. Generate station time series from free GNSS solutions in SINEX format
```
pyacs_make_time_series.py
```
2. Load time series and work interactively in Jupyter notebook environment

```ipython
from pyacs.gts.Sgts import Sgts
ts = Sgts("timeseries_directory")
ts.CODE.plot()
ts.CODE.add_offsets_dates([2016.29, 2018.315]).plot()
...
dts = ts.gts('detrend')
```
3. Export velocity field and perform Euler pole calculation or strain rate analysis

---

# Installation

It is recommended to install PYACS in a dedicated Python environment.

## Recommended setup using mamba

The recommended environment manager is mamba

Install Miniforge / mamba:
```
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh
bash Miniforge3-$(uname)-$(uname -m).sh
```
Open a new terminal after installation.

Download the environment configuration file from https://github.com/JMNocquet/pyacs36/tree/master/environment.yaml.

Create the PYACS environment:
```
mamba env create -f environment.yaml
mamba activate pyacs
```

## Installation from an existing Anaconda environment

If you already use Anaconda or Miniconda, you can install mamba in the base environment:

```
Download the environment configuration file from https://github.com/JMNocquet/pyacs36/tree/master/environment.yaml
conda install -n base -c conda-forge mamba
mamba env create -f environment.yaml
```
Note: the standard conda dependency solver may be significantly slower.

## Getting PYACS

- get the latest version from: https://github.com/JMNocquet/pyacs36/tree/master/dist

```
pip install pyacs-X.XX.XX.tar.gz
```
or if you plan to modify the code, select an installation directory and:
```
tar xvfz pyacs-X.XX.XX.tar.gz
cd pyacs-X.XX.XX
pip install .
```

- or use git to clone the whole project:
```
git clone https://github.com/JMNocquet/pyacs36.git
pip install .
```

Note: The latest version might be for development, while commited version are expected to be more stable release.


### Running tests
From the directory containing the pyacs package:
```
pytest pyacs/tests
```
### Interactive use

**ipyacs.py** is a convenient script that load the main PYACS librairies and automatically loads time series in the current directory if available. It then allows interactive time series visualization and analysis under ipython. 
You may define a shell alias:
```
alias ipyacs='mamba activate pyacs && ipython $(which ipyacs.py) -i'
```

### Working with Jupyter notebooks

Time series analysis is more conveniently performed using a Jupyter notebook.
Open a jupyter notebook, select the kernel corresponding to your pyacs installation and paste the following to start:

```
# import
import pyacs
print(pyacs.__version__)
import numpy as np
from pyacs.gts.Sgts import Sgts
import pyacs.lib.astrotime as at
import pyacs
from datetime import datetime

# define your backend
#%matplotlib inline
#%matplotlib qt

# verbose level
pyacs.verbose('SILENT')
# data to analyze
ts_dir = 'your_time_series_path_directory'
# load data
ts = Sgts( ts_dir, verbose=False)
print("%d time series loaded " % ts.n() )
ts.info()
```

# Building a PYACS distribution
Advanced users can build their own distribution using
```
python -m build
```
The source and wheel distributions will be generated in the dist/ directory.

### Documentation

A html documentation is available [here][pyacs_doc],or in pyacs/docs/build/html/index.html. Alternatively, the local documentation cab be generated locally:
```
./make_pyacs_doc_html_sphinx.sh
```

# Reporting issues

Report bugs or desired enhancement in the [Github issues][github_issues].

## Authors
PYACS has been developed by [Jean-Mathieu Nocquet][web_nocquet] and [Dinh Trong Tran][tran_researchgate] as part of his PhD.  

[latest_pyacs]:https://github.com/JMNocquet/pyacs36/tree/master/dist
[pyacs_doc]:https://jmnocquet.github.io/pyacs36
[web_nocquet]:https://jmnocquet.github.io/
[tran_researchgate]:https://www.researchgate.net/profile/Dinh-Tran-14
[github_issues]:https://github.com/JMNocquet/pyacs36/issues

## A few studies which have used PYACS

- Cruz-Atienza, Víctor M., et al. "Short-term interaction between silent and devastating earthquakes in Mexico." Nature communications 12.1 (2021): 1-14.
- Bontemps, Noélie, et al. "Rain and small earthquakes maintain a slow-moving landslide in a persistent critical state." Nature communications 11.1 (2020): 1-10.
- Klein, Emilie, et al. "Interplay of seismic and a-seismic deformation during the 2020 sequence of Atacama, Chile." Earth and Planetary Science Letters 570 (2021): 117081.
- Bletery, Quentin, and Jean-Mathieu Nocquet. "Slip bursts during coalescence of slow slip events in Cascadia." Nature communications 11.1 (2020): 1-6.
- Rolandone, Frédérique, et al. "Areas prone to slow slip events impede earthquake rupture propagation and promote afterslip." Science advances 4.1 (2018): eaao6596.
- Socquet, Anne, et al. "An 8 month slow slip event triggers progressive nucleation of the 2014 Chile megathrust." Geophysical Research Letters 44.9 (2017): 4046-4053.
