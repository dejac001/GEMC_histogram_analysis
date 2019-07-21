.. GEMC_histogram_analysis documentation master file, created by
   sphinx-quickstart on Sun Jul  7 05:58:17 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to GEMC_histogram_analysis's documentation!
===================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

GEMC_histogram_analysis
=======================

An analysis package utilizing histograms to calculate accurate
vapor and liquid coexistence densities, saturated vapor pressure,
and compressibility factor from a Gibbs ensemble Monte Carlo trajectory
for unary vapor--liquid equilibria near the critical point.


Citation
========
Users of this package should cite:

B.L. Eggimann, Y.Z.S. Sun, R.F. DeJaco, R. Singh, M. Ahsan, T.R. Josephson, and J.I. Siepmann. Assessing the quality of molecular simulations for vapor–liquid equilibria: An analysis of the TraPPE database. *J. Chem. Eng. Data*, submitted for publication on 7/22/2019.


Installation
============

This program is intended to be used on python3.X. It won't work on python2.X
because it uses new types of importing that are only available in python3.X.
The python packages can be installed by the usual route like

.. code-block:: bash

   pip install requirements.txt

but use of a virtualenvironment
or a conda environment is recommended. If the user is unfamiliar with
these approaches, the necessary packages will come with almost
all anaconda python installations which are available on supercomputers.

To test your installation of the program, look through the
functions in `examples.py` and run one of the examples.


GEMC_histogram_analysis runs on Python 3.x only (not Python 2.x).


Usage
============

The program can generate two types of histograms.
A histogram of densities for determining
bulk properties is generated with the file
`run_histogram.py`. A histogram of the sampling
trajectory can be generated with the file
`plot_trajectory.py`.
Either program can be run by

.. code-block:: bash

   python eitherprogram -fi file1 file2 file3 ... -T T1

where `eitherprogram` is either
`run_histogram` or `plot_trajectory.py`,  `file1`, `file2`, `file3`,  and `...`
are the relative paths to all the trajectory
files from one independent simulation, and `T1`
is the temperature at which the simulation was
conducted in Kelvin. All
arguments (including those which are optional)
can be seen by appending `-h` or `--help`
to the line above.

It is not necessary to use this command line approach to run
these programs. To see how
to do this, take a look at the file `examples.py`
or the functions `example_butane()`
and `example_pentanal()` within the file `plot_trajectory.py`.


Useful Scripts
==============
run_histogram.py
----------------
The `run_histogram.py` program will calculate
mean and errors for pressures, densities, and compressibilities,
and compare to the block averages approach.

The mean and error of the density values with
the histogram approach are obtained from the
density probability distributions.

The pressure and compressibilities obtained from
the histogram approach are calculated as the average of all
pressures and compressibilities at MCCs where the pressure was calculated
and the density was in the high density or low density regime
(as defined by 75% of the respective peak height).
The errors are the standard deviations of those set of values.
These values for one independent simulation correspond to the
fluctuations.

plot_trajectory.py
------------------

Methodology
===========
The statistics from the production trajectory of a
:math:`NVT` Gibbs ensemble Monte Carlo
simulation are collected.
These statistics entail :math:`\mathcal{D}`, the distribution
of total box densities :math:`\rho_{i,j}`
calculated in each simulation
box :math:`i` at the end of the :math:`j`-th Monte Carlo Cycle (MCC),
and :math:`\mathcal{P}`, the distribution of box pressures
:math:`P_{i,k}`
calculated in each simulation box
:math:`i` at the end of the :math:`k`-th
MCC where the pressure was calculated.
Usually, the box pressures are calculated once every five to
ten MCCs, so that the number of entries in
:math:`\mathcal{P}` is less than those in :math:`\mathcal{D}`.

The distribution of all observed box densities :math:`\mathcal{D}`
is partitioned into subsets of low
(:math:`\mathcal{L}`) and high (:math:`\mathcal{H}`) densities,
such that
:math:`\rho_{i,l}<\rho_m  \forall i \forall l \in\mathcal{L}`
and
:math:`\rho_{i,h}\ge\rho_m  \forall i \forall h \in\mathcal{H}`
where :math:`\rho_m` is the mean of all density values in :math:`mathcal{D}`.

For each distribution, :math:`\mathcal{L}` and :math:`\mathcal{H}`,
a histogram of the densities is performed with :math:`N_{\mathrm{B}}`
equally-spaced bins (edges) from the minimum to maximum value of density using [NUMPY]_.
The probabilities :math:`Y_i\forall i \in \mathcal{Y}` of observing each density within
each :math:`i` nodes (adjacent bins)
are normalized so that the sum is equal to unity.
The :math:`N_{\mathrm{B}}-1` density values :math:`\rho_i \forall i \in \mathcal{Y}`
for each node (within adjacent edges) are assumed
to be the midpoint between adjacent edges.
A subset of representative probabilities,
:math:`\mathcal{Y}^0`, is chosen as follows

.. math::

    Y_j\ge f\max_{i\in\mathcal{Y}}Y_i\;\;\;\forall j \in \mathcal{Y}^*

where :math:`f` is the fraction of peak maximum.
The fraction of peak maximum is typically 0.75
[DINPAJOOH2015]_.
These probabilities and associated densities,
:math:`\rho_i\forall i \in \mathcal{Y}^*`
are fit to a gaussian function.
The initial guess for the peak height, mean,
and standard deviation of the gaussian are obtained from
:math:`\max_{\forall i \in \mathcal{Y}^*}Y_i`,
:math:`\mu_{\forall i \in \mathcal{Y}^*}\rho_i`,
and the standard deviation of of the distribution :math:`\sigma_{\mathcal{Y}^*}`,
respectively.
The lower bounds allowed for the peak height, mean,
and standard deviation are zero, the minimum density,
and zero, respectively.
The upper bounds allowed for the peak height and mean
are :math:`2\max_{\forall i \in \mathcal{Y}^*}Y_i` and
:math:`\max_{\forall i \in \mathcal{Y}^*}\rho_i`,
respectively.
The maximum in the standard deviation parameter
obtained from the fitting of the gaussian
is not constrained to a maximum finite value.
The constrained parameter fitting is performed using
the Trust Region Reflective algorithm as
implemented in [SCIPY]_.

Finally, the minimum density, :math:`\rho_{\chi,\mathrm{min}}`, and maximum density
:math:`\rho_{\chi,\mathrm{max}}`,
for each distribution :math:`\chi` (:math:`\mathcal{H}` or :math:`\mathcal{L}`)
corresponding to the chosen fraction of peak maximum
:math:`f` are calculated from the gaussian function and the optimized parameters for each peak.

The final high (:math:`\mathcal{H}^*`) and low (:math:`\mathcal{L}^*`)
distributions of densities are determined as follows

.. math::

    \rho_{\mathcal{H},\mathrm{min}}\le \rho_i\le \rho_{\mathcal{H},\mathrm{max}} \;\;\;\;\forall i \in \mathcal{H}^*, \ni i \in \mathcal{H} \\
    \rho_{\mathcal{L},\mathrm{min}}\le \rho_i\le \rho_{\mathcal{L},\mathrm{max}} \;\;\;\;\forall i \in \mathcal{L}^*, \ni i \in \mathcal{L}


Similarly, the pressures corresponding to the high (:math:`\mathcal{P}_{\mathcal{H}^*}`)
and low (:math:`\mathcal{P}_{\mathcal{L}^*}`)distributions of densities are determined as

.. math::

    \rho_{\mathcal{H},\mathrm{min}}\le \rho_i\le \rho_{\mathcal{H},\mathrm{max}} \;\;\;\;\forall i \in \mathcal{P}_{\mathcal{H}^*}, \ni i \in \mathcal{P} \\
    \rho_{\mathcal{L},\mathrm{min}}\le \rho_i\le \rho_{\mathcal{L},\mathrm{max}} \;\;\;\;\forall i \in \mathcal{P}_{\mathcal{L}^*}, \ni i \in \mathcal{P} \\

The compressibility distributions at the high (:math:`\mathcal{Z}_{\mathcal{H}^*}`)
and low (:math:`\mathcal{Z}_{\mathcal{L}^*}`) densities
are calculated at the densities and pressures within :math:`\mathcal{P}_{\mathcal{H}^*}`
and :math:`\mathcal{P}_{\mathcal{L}^*}`, respectively.
The mean and standard deviation of each density regime (:math:`\chi\in\{\mathcal{H}^*,\mathcal{L}^*\}`)
for pressure, compressibility, and density
are determined from the :math:`\chi`, :math:`\mathcal{P}_{\chi}`, and :math:`\mathcal{Z}_{\chi}`
distributions, respectively.

.. [NUMPY] Oliphant, T. NumPy: A guide to NumPy. 2006--2019. USA: Trelgol Publishing. http://www.numpy.org/
.. [DINPAJOOH2015] Dinpajooh M, Bai P, Allan DA, and Siepmann JI. Accurate and precise determination
    of critical properties from Gibbs ensemble Monte Carlo simulations. J. Chem. Phys. 143, 114113 (2015),
    https://doi.org/10.1063/1.4930848
.. [SCIPY] Jones E, Oliphant E, Peterson P, et al. SciPy: Open Source Scientific Tools for Python. 2001--2019. http://www.scipy.org/

histogram
=========
.. automodule:: histogram.chem_constants
   :members:
.. automodule:: histogram.density_histogram
   :members:
.. automodule:: histogram.reader
   :members:
.. automodule:: histogram.trajectory_histogram
   :members:
.. automodule:: histogram.writer
   :members:

scripting files
===============
.. automodule:: run_histogram
   :members:
.. automodule:: plotting_util
   :members:
.. automodule:: plot_trajectory
   :members:
.. automodule:: examples
   :members:
.. automodule:: analysis_factory
   :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
