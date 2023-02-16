=====
Setup
=====

--------
Overview
--------

**NEDUM-2D** (Non-Equilibrium Dynamic Urban Model) is a tool for simulating land-use patterns across a city in two dimensions. The initial rationale for its development was to make urban simulations more realistic by modelling deviations from theoretical static equilibria caused by inertia in the adaptive behaviour of households and developers over time. The current version developed for the *City of Cape Town* (CoCT) incorporates several transportation modes, employment centres, income groups, and housing types, to provide even more realistic prospective scenarios, while remaining tractable in terms of causal mechanisms at play. It allows to model the spatial impact of policies with a special interest for cities in developing countries, such as informal housing regulation and localized flood protection investments.

---------
Resources
---------

Documentation is freely available `here <https://cired.github.io/cape_town_NEDUM_Python/>`__.

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
	├── _flood_processing (notebook + data used to pre-process Deltares coastal flood maps)
	├── _research (research articles referenced in the code)
	├── docs (html output for documentation)
	├── calibration (package used for calibration)
	├── equilibrium (package used to compute equilibrium outcomes)
	├── inputs (package used to import and process inputs)
	├── outputs (package used for output plots and tables)
	├── .gitignore (defines files to ignore when pushing commits online)
	├── LICENSE (open source license file)
	├── README.rst (introduction file)
	├── 0_calib_nb.ipynb (Jupyter notebook that runs calibration)
	├── 0_calib_nb.py (Jupytext-paired Python script that runs calibration)
	├── 1_main_nb.ipynb (Jupyter notebook that runs the model)
	├── 1_main_nb.py (Jupytext-paired Python script that runs the model)
	├── 2_plots_equil.py (exports basic plots and tables for initial state static equilibrium)
	├── 2_plots_inputs.py (exports basic plots and tables for input data)
	├── 2_plots_simul.py (exports basic plots and tables for subsequent dynamic simulations)
	├── 3_plots_use_case_cchange.py (exports interactive plots and tables for c.change use case)
	├── 3_plots_use_case_insur.py (exports interactive plots and tables for insurance use case)
	├── 4_use_case_nb_empty.ipynb (Jupyter notebook that loads and comments on use case outputs)
	└── 4_use_case_nb_full.ipynb (Jupyter notebook that loads and comments on use case outputs)

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

**Step 3**: Launch **NEDUM-2D**

* If needed, run the ``0_calib_nb`` notebook (either in .py or .ipynb format) to calibrate parameters again (under ``precalculated_inputs``) if underlying data (in ``data_Cape_Town``) has changed. A static copy is shown in the documentation for illustrative purposes.
* Execute the ``1_main_nb`` notebook (either in .py or .ipynb format) to run the simulations and obtain a preview of results. Outcomes will be automatically saved in a dedicated subfolder (according to a naming convention defined in the preamble of the script) under the ``Output`` directory. A static copy is shown in the documentation for illustrative purposes.
* Run the ``2_plots`` scripts to export static tables and figures in dedicated subfolders. ``2_plots_inputs.py`` plots input data that does not change across scenarios (output is saved in dedicated subfolders under the ``Output`` directory). ``2_plots_equil.py`` plots outcomes for the initial state static equilibrium. Output changes across scenarios and is saved under the specific subfolder created at previous step with ``1_main_nb``. ``2_plots_simul.py`` plots outcomes for dynamic simulations over subsequent periods. Output is saved in same directory as for ``2_plots_equil.py``: we have created a separate script in case the end user is not interested in the dynamics (which are long to loop over). In case of discrepancies, we recommend using ``2_plots_equil.py`` for each period individually.
* Run the ``3_plots_use_case_anticip.py`` and ``3_plots_use_case_cchange.py`` to export the interactive plots and tables associated with respectively the insurance and climate change use cases that we developed for this iteration of the model. This requires to run both ``1_main_nb`` and ``2_plots_equil.py`` for the scenarios used in each script (with and without insurance, with and without climate change). Those options can be modified in the preamble of the scripts (see :doc:`../technical_doc` for more details).
* Run the ``4_use_case_nb_empty.ipynb`` notebook to recover key interactive plots from ``3_plots_use_case_anticip.py`` and ``3_plots_use_case_cchange.py`` with associated comments and interpretation. As the interactive plots are too heavy to save or load as a ``.html`` page, we save the file without the associated output. ``4_use_case_nb_full.ipynb`` provides a static version with cached output (without the possibility to zoom or hover over the plots), that is shown in the documentation for illustrative purposes.
* See :doc:`../technical_doc` for more details on running custom simulations. Note that to keep ``.py`` and ``.ipynb`` versions of the same script in sync, one needs to pair them by setting up Jupytext locally.

----------
Versioning
----------

* The ``gh_pages`` branch contains the latest update of the code and is set as default. If you want to modify the code, please fork the repository and start from this branch, as this is the one used in this documentation.
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

For internal reference within the CoCT, please contact kristoff.potgieter@capetown.gov.za.

|

Thomas Monnier - `Website <https://tlmonnier.github.io>`_ - `Github <https://github.com/TLMonnier>`_ - `Twitter <https://twitter.com/TLMonnier>`_ - thomas.monnier@ensae.fr

Distributed under the GNU GENERAL PUBLIC LICENSE.

https://github.com/CIRED/cape_town_NEDUM_Python