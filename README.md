# PYACS
PYACS is a set of scripts and python modules for analyzing and modeling geodetic data.
PYACS includes 4 components:
- a set a core packages handling coordinates, time, estimators, format conversions
- pyacs_make_time_series.py is a script to express free GNSS solutions 
  in a single reference frame and produce time series
  
- the Gts module allows versatile time series analysis either for individual or for a group of time series.  
- a package dedicated to GNSS velocity fields analysis, including Euler pole calculation and basic strain rate analysis

## Installation

PYACS works well under Anaconda with python >= 3.8. 
Get the latest version of PYACS from the [release directory][latest_pyacs] and install it with pip.

```
pip install pyacs-XXX.tar.gz
```
You might prefer to install it under a dedicated conda virtual env
```
conda create -n env_pyacs python=3.9 ipython
conda activate env_pyacs
pip install pyacs-XXX.tar.gz
```
**ipyacs.py** is a useful script for interactive time series visualization and 
analysis under ipython. It is convenient to add an alias in your shell 
configuration file

```
alias ipyacs='ipython `which ipyacs.py` -i'
```
## Documentation

A documentation is available [here][pyacs_doc].
Report bug or desired enhancement in the [Github issues][github_issues].

## Authors
PYACS has been developed by [Jean-Mathieu Nocquet][web_nocquet] and [Dinh Trong Tran][tran_researchgate] as part of his PhD.  

[latest_pyacs]:https://github.com/JMNocquet/pyacs36/tree/master/dist
[pyacs_doc]:https://jmnocquet.github.io/pyacs_docs/pyacs
[web_nocquet]:https://jmnocquet.github.io/
[tran_researchgate]:https://www.researchgate.net/profile/Dinh-Tran-14
[github_issues]:https://github.com/JMNocquet/pyacs36/issues

## A few studies which have used PYACS

- Cruz-Atienza, Víctor M., et al. "Short-term interaction between silent and devastating earthquakes in Mexico." Nature communications 12.1 (2021): 1-14.
- Bontemps, Noélie, et al. "Rain and small earthquakes maintain a slow-moving landslide in a persistent critical state." Nature communications 11.1 (2020): 1-10.
- Klein, Emilied, et al. "Interplay of seismic and a-seismic deformation during the 2020 sequence of Atacama, Chile." Earth and Planetary Science Letters 570 (2021): 117081.
- Bletery, Quentin, and Jean-Mathieu Nocquet. "Slip bursts during coalescence of slow slip events in Cascadia." Nature communications 11.1 (2020): 1-6.
- Rolandone, Frédérique, et al. "Areas prone to slow slip events impede earthquake rupture propagation and promote afterslip." Science advances 4.1 (2018): eaao6596.
- Socquet, Anne, et al. "An 8 month slow slip event triggers progressive nucleation of the 2014 Chile megathrust." Geophysical Research Letters 44.9 (2017): 4046-4053.
