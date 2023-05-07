..
	.. image:: images/GFDRR-logo.png
		:width: 45.75 %
		:alt: Global Facility for Disaster Reduction and Recovery logo
	.. image:: images/NT-Logo.jpg
		:width: 24.75 %
		:alt: National Treasury logo
	.. image:: images/CSP-Logo-CMYK.jpg
		:width: 24.95 %
		:alt: Cities Support Programme logo

	.. image:: images/SECO-logo.jpg
		:width: 48 %
		:alt: State Secretariat for Economic Affairs logo
	.. image:: images/COCT-logo.jpeg
		:width: 38 %
		:alt: City of Cape Town logo
	.. image:: images/Logo_CIRED.png
		:width: 12 %
		:alt: Centre International de Recherche sur l’Environnement et le Développement logo

.. image:: images/banner_crop.png
	:width: 100 %
	:alt: Sponsors' logo banner

|

=====
Setup
=====

--------
Overview
--------

**NEDUM-2D** (Non-Equilibrium Dynamic Urban Model) is a tool for simulating land-use patterns across a city in two dimensions. The initial rationale for its development was to have a simple and tractable urban simulation model, based on core urban economic principles, able to reproduce some of the most striking stylized facts observed across cities worldwide, with the goal of studying the impacts of land-use and transport policies on environmental and welfare outcomes. One core addition of **NEDUM-2D** was the introduction of inertia in adaptative responses from households and developers to policy shocks, to study their consequences on transition costs. The current version developed for the *City of Cape Town* (CoCT) incorporates several transportation modes, employment centres, income groups, and housing types, to provide even more realistic prospective scenarios, while remaining tractable in terms of causal mechanisms at play. It allows to model the spatial impact of policies with a special interest for cities in developing countries, such as informal housing regulation and localized flood protection investments.

---------
Resources
---------

Documentation is freely available at `https://cired.github.io/cape_town_NEDUM_Python/ <https://cired.github.io/cape_town_NEDUM_Python/>`_.

The reference working paper used along the documentation is available at `https://openknowledge.worldbank.org/entities/publication/f599340d-9916-5611-8f54-e9dda7eef282 <https://openknowledge.worldbank.org/entities/publication/f599340d-9916-5611-8f54-e9dda7eef282>`_.

------------
Installation
------------

**Step 1**: Git clone **NEDUM-2D** repository in your computer

* Go to `https://github.com/CIRED/cape_town_NEDUM_Python <https://github.com/CIRED/cape_town_NEDUM_Python>`_ and set up the repository on your local machine (choose your preferred option clicking on the ``code`` button).
* The tree structure of the repository should correspond to (comments in parentheses)::

	code capetown python (name of the repo)
	├── _doc_source (source code for documentation)
	├── _ex (exercise notebooks for training sessions)
	├── _flood_processing (notebook + data used to pre-process Deltares coastal flood maps)
	├── _old (old Python scripts used for development)
	├── _research (research articles referenced in the code)
	├── docs (html output for documentation)
	├── calibration (package used for calibration)
	├── equilibrium (package used to compute equilibrium outcomes)
	├── inputs (package used to import and process inputs)
	├── outputs (package used for output plots and tables)
	├── .gitignore (defines files to ignore when pushing commits online)
	├── LICENSE (open source license file)
	├── nedum.yml (packages for conda)
	├── README.rst (introduction file)
	├── requirements.txt (packages for pip)
	├── 0_calib_nb.ipynb (Jupyter notebook that runs calibration)
	├── 1_main_nb.ipynb (Jupyter notebook that runs the model)
	├── 2_plots_equil.py (exports basic plots and tables for initial state static equilibrium)
	├── 2_plots_inputs.py (exports basic plots and tables for input data)
	├── 2_plots_simul.py (exports basic plots and tables for subsequent dynamic simulations)
	├── 3_plots_use_case_cchange.py (exports interactive plots and tables for climate change use case)
	├── 3_plots_use_case_anticip.py (exports interactive plots and tables for anticipation use case)
	├── 4_use_case_nb_empty.ipynb (Jupyter notebook that loads and comments on use case outputs)
	└── 4_use_case_nb_full.ipynb (Jupyter notebook that loads and comments on use case outputs)

**Step 2**: Set up project directory

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

Note that we do not host the data on GitHub, but that it is available upon request.

**Step 3**: Set up Python environment

* If you are not familiar with Python, we recommend using the `Anaconda <https://www.anaconda.com/download/>`_ distribution and setting up a dedicated `virtual environment <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_ for **NEDUM-2D**.
* Open your terminal (Anaconda prompt) and go where you located the ``code capetown python`` folder (by pasting its local path after the command ``cd`` in Anaconda prompt, then press enter)
* Make sure that all the packages indicated in the ``requirements.txt`` file are installed. Alternatively, if you are using Anaconda, you can directly create the adequate environment by typing ``conda env create -f nedum.yml`` from your base environment in the Anaconda terminal (this may take some time to run).

**Step 4**: Launch **NEDUM-2D**

* If you are using Anaconda, activate the adequate environment by typing ``conda activate nedum`` in the Anaconda terminal. Open the integrated development environment (IDE) of your choice by typing its name in the terminal: we recommend using ``jupyter-notebook`` for ``.ipynb`` scripts and ``spyder`` for ``.py`` scripts, as they are the two IDEs set up with the ``nedum`` environment.
* Note that a set of options is left for the user to change in the preamble of most scripts. Each combination corresponds to a specific simulation run, that is saved according to a given naming convention. Before running the scripts, please make sure that the options correspond to the simulation run you are interested in.
* If needed, run the ``0_calib_nb`` notebook to calibrate parameters again (under ``precalculated_inputs``) if underlying data (in ``data_Cape_Town``) has changed. A static copy is shown in the documentation for illustrative purposes.
* Execute the ``1_main_nb`` notebook to run the simulations and obtain a preview of results. Outcomes will be automatically saved in a dedicated subfolder (according to the naming convention defined in the preamble of the script) under the ``Output`` directory. A static copy is shown in the documentation for illustrative purposes.
* Run the ``2_plots`` scripts to export static tables and figures in dedicated subfolders. ``2_plots_inputs.py`` plots input data that does not change across scenarios (output is saved in dedicated subfolders under the ``Output`` directory). ``2_plots_equil.py`` plots outcomes for the initial state static equilibrium. Output changes across scenarios and is saved under the specific subfolder created at previous step with ``1_main_nb``. ``2_plots_simul.py`` plots outcomes for dynamic simulations over subsequent periods. Output is saved in same directory as for ``2_plots_equil.py``: we have created a separate script in case the end user is not interested in the dynamics (which are long to loop over). In case of discrepancies, we recommend using ``2_plots_equil.py`` for each period individually.
* Run the ``3_plots_use_case_anticip.py`` and ``3_plots_use_case_cchange.py`` to export (in a specific subfolder) the interactive plots and tables associated with respectively the anticipation and climate change use cases that we developed for this version of the model. This requires to run both ``1_main_nb`` and ``2_plots_equil.py`` for the scenarios used in each script (with and without anticipation, with and without climate change). The options that define which scenario to run can be modified in the preamble of the scripts (see :doc:`../technical_doc` for more details).
* Run the ``4_use_case_nb_empty.ipynb`` notebook to recover key interactive plots from ``3_plots_use_case_anticip.py`` and ``3_plots_use_case_cchange.py`` with associated comments and interpretation. As the interactive plots are too heavy to save and load as a ``.html`` page, we saved the notebook without the associated cached output. ``4_use_case_nb_full.ipynb`` provides a static version with cached output (without the possibility to zoom in and out or display information by hovering over the plots), that is shown in the documentation for illustrative purposes.
* See :doc:`../technical_doc` for more details on running custom simulations.

----------
Versioning
----------

* The ``gh_pages`` branch of the GitHub repository contains the latest update of the code and is set as default. If you want to modify the code, please fork the repository and start from this branch, as this is the one used in this documentation.
* All other branches are deprecated.

-----------------
About the authors
-----------------

The development of the **NEDUM-2D** model was initiated at *CIRED* in 2009. Coordinated by Vincent Viguié, it involved over the years, in alphabetic order: Paolo Avner, Stéphane Hallegattte, Charlotte Liotta, Thomas Monnier, Basile Pfeiffer, Claus Rabe, Julie Rozenberg, and Harris Selod.

.. _meta_link:

----
Meta
----

If you find **NEDUM-2D** useful, please kindly cite our paper:

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

For internal reference within the CoCT, please contact kristoff.potgieter@capetown.gov.za.

|

User guide conception: Thomas Monnier - `Website <https://tlmonnier.github.io>`_ - `Twitter <https://twitter.com/TLMonnier>`_ - `Github <https://github.com/TLMonnier>`_ - thomas.monnier@ensae.fr

Distributed under the GNU GENERAL PUBLIC LICENSE.

https://github.com/CIRED/cape_town_NEDUM_Python