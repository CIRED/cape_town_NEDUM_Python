====================
Introducing NEDUM-2D
====================

--------
Overview
--------

**NEDUM-2D** (Non-Equilibrium Dynamic Urban Model) is a tool for simulating land-use patterns across a city in two dimensions. The initial rationale for its development was to make urban simulations more realistic by modelling deviations from theoretical static equilibria caused by inertia in the adaptive behaviour of households and developers over time. The current version developed for the *City of Cape Town* (CoCT) incorporates several transportation modes, employment centres, income groups, and housing types, to provide even more realistic prospective scenarios, while remaining tractable in terms of causal mechanisms at play. It allows to model the spatial impact of policies with a special interest for cities in developing countries, such as informal housing regulation and localized flood protection investments.

---------
Resources
---------

Documentation is freely available `here <https://cired.github.io/cape_town_NEDUM_Python/html/index.html>`__.

A simple user interface is available `here <https://kristoffpotgieter-nedumapp-app-f2rto5.streamlitapp.com/>`__.

The reference working paper used along the documentation is available `here <https://openknowledge.worldbank.org/handle/10986/31987?locale-attribute=fr>`__.

------------
Installation
------------

**Step 1**: Git clone **NEDUM-2D** repository in your computer

* Use your terminal and go to a location where you want to store the **NEDUM-2D** model
* Type: ``git clone https://github.com/CIRED/cape_town_NEDUM_Python.git``
* The tree structure of the repository should correspond to (comments in parentheses)::

	code capetown python (name of the repo)
	├── _doc_source (source code for documentation)
	├── _flood_processing (notebook used to pre-process Deltares coastal flood maps)
	├── _research (research articles referenced in the code and the documentation)
	├── calibration (package used for calibration)
	├── docs (html output for documentation)
	├── equilibrium (package used to compute equilibrium outcomes)
	├── inputs (package used to import and process inputs)
	├── outputs (package used for output plots and tables)
	├── .gitignore (defines files to ignore when pushing commits online)
	├── LICENSE (open source license file)
	├── README.rst (introduction file)
	├── calib_nb.ipynb (Jupyter notebook that runs calibration)
	├── calib_nb.py (paired Python script that runs calibration)
	├── main_nb.ipynb (Jupyter notebook that runs the model)
	├── main_nb.py (paired Python script that runs the model)
	├── plots_equil.py (exports plots and tables for initial state static equilibrium)
	├── plots_inputs.py (exports plots and tables for input data)
	└── plots_simul.py (exports plots and tables for subsequent dynamic simulations)

.. **Step 2**: Create a conda environment from the *nedum-2d-env.yml* file

.. ..
.. 	Create the environment file

.. * The *nedum-2d-env.yml* file is in the **NEDUM-2D** repository
.. * Use the terminal and go to the **NEDUM-2D** repository stored on your computer
.. * Type: ``conda env create -f nedum-2d-env.yml``

.. **Step 3**: Activate the new environment

.. * The first line of the *.yml* file sets the new environment’s name
.. * Type: ``conda activate NEDUM-2D``

**Step 2**: Set project directory

* To run properly, the **NEDUM-2D** repository (here, ``code capetown python``) should be included in a project folder that also contains input data (and an empty output folder), according to the following tree structure (comments in parentheses)::

	.
	├── Data
	│   ├── Aux data (auxiliary data used to reference raw inputs in the code)
	│   ├── data_Cape_Town (input data provided by the CoCT)
	│   ├── flood_maps (pre-processed flood maps)
	│   ├── occupation_maps (maps for informal settlement expansion scenarios)
	│   ├── precalculated_inputs (calibrated parameters and data)
	│   ├── precalculated_transport (intermediate outputs of commuting choice model)
	│   ├── housing_types_sal_analysis.xlsx (Small-Area-Level data on housing types)
	│   └── housing_types_grid_sal.xlsx (SAL data transposed to CoCT's grid level)
	├── Output
	└── code capetown python

..
	Do we need to set the repo as a project in Spyder?

**Step 3**: Launch **NEDUM-2D**

* From ``code capetown python`` root, execute the ``main_nb`` notebook (either in .py or .ipynb format) to run the simulations and obtain a preview of results. A non-interactive copy is shown in the documentation for illustrative purposes
* Run one of the ``plots`` scripts to export tables and figures in dedicated subfolders (under the ``Output`` directory)
* If needed, run the ``calib_nb`` notebook (either in .py or .ipynb format) to calibrate parameters again if underlying data has changed
* See :doc:`../technical_doc` for more details on running custom simulations / calibration. Note that when modifying ``main_nb`` or ``calib_nb`` notebooks (either in .py or .ipynb format), the associated paired file needs to be updated using Jupytext.

----------
Versioning
----------

..
	Set as default branch

* The ``gh_pages`` branch contains the latest update of the code and is set as default. If you want to modify the code, please fork the repository and start from this branch, as this is the one used in this documentation.
* The ``main`` branch contains the original code (with some extra features) from our last paper (ported from Matlab to Python)
* The ``TLM-edits`` branch contains some code and folder reorganization without rewriting anything
* The ``TLM-write`` branch contains some rewriting and commenting
* The branches ending in ``_specif`` are tests for several specifications of the code
* All other branches are deprecated

-----------------
About the authors
-----------------

The development of the **NEDUM-2D** model was initiated at *CIRED* in 2014. Coordinated by Vincent Viguié, it involved over the years, in alphabetic order, Paolo Avner, Stéphane Hallegattte, Charlotte Liotta, Thomas Monnier, Basile Pfeiffer, Claus Rabe, Julie Rozenberg, and Harris Selod.

.. _meta_link:

----
Meta
----

If you find **NEDUM-2D** useful, please kindly cite our last paper:

.. code-block:: latex

	@techreport{
	  author      = {Pfeiffer, Basile and Rabe, Claus and Selod, Harris and Viguié, Vincent},
	  title       = {Assessing Urban Policies Using a Simulation Model with Formal and Informal Housing:
	  Application to Cape Town, South Africa},
	  year        = {2019},
	  institution = {World Bank},
	  address     = {Washington, DC},
	  series      = {Policy Research Working Paper},
	  type        = {Working Paper},
	  number      = {8921},
	  url         = {https://openknowledge.worldbank.org/handle/10986/31987}
	}

Thomas Monnier - `tlmonnier.github.io <https://tlmonnier.github.io>`_ - `Github <https://github.com/TLMonnier>`_ - `@TLMonnier <https://twitter.com/TLMonnier>`_ - thomas.monnier@ensae.fr

Distributed under the GNU GENERAL PUBLIC LICENSE.

https://github.com/CIRED/cape_town_NEDUM_Python