# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 12:22:55 2020.

@author: Charlotte Liotta
"""

import numpy as np
import pandas as pd


def import_employment_data(households_per_income_class, param, path_data):
    """
    Import number of jobs per income group in each selected employment center.

    Parameters
    ----------
    households_per_income_class : ndarray(float64)
        Exogenous total number of households per income group (excluding people
        out of employment, for 4 groups)
    param : dict
        Dictionary of default parameters
    path_data : str
        Path towards data used in the model

    Returns
    -------
    jobsCentersNGroupRescaled : ndarray(float64, ndim=2)
        Number of jobs in each selected job center (185) per income group (4).
        Remember that we rescale the number of individual jobs to reflect total
        household employment, as our income and population data are for
        households only: one job basically provides employment for two people.
        This simplification allows to model households as a single
        representative agent and to abstract from a two-body problem.
        Empirically, this holds on aggregate as households' position on the
        labor market is often determined by one household head.

    """
    # We import number of jobs per Transport Zone (1,787)
    TAZ = pd.read_csv(path_data + 'TAZ_amp_2013_proj_centro2.csv')

    # Number of employees in each TAZ for the 12 income classes from data
    # NB: we break income classes from TAZ data to stick better to equivalent
    # from census data. Then, we will apply the proper allocation to go to the
    # 4 income classes used in the model (do not correspond to the 4 income
    # classes used in TAZ data)
    jobsCenters12Class = np.array(
        [np.zeros(len(TAZ.Ink1)), TAZ.Ink1/3, TAZ.Ink1/3, TAZ.Ink1/3,
         TAZ.Ink2/2, TAZ.Ink2/2, TAZ.Ink3/3, TAZ.Ink3/3, TAZ.Ink3/3,
         TAZ.Ink4/3, TAZ.Ink4/3, TAZ.Ink4/3]
        )

    # We get geographic coordinates corresponding to CoCT's grid
    codeCentersInitial = TAZ.TZ2013
    xCoord = TAZ.X / 1000
    yCoord = TAZ.Y / 1000

    # We select a restricted number of employment centers with a minimum number
    # of jobs, for the sake of numerical simplicity
    selectedCenters = (
        sum(jobsCenters12Class, 0) > param["job_center_threshold"])

    # We select out job centers that may fall out of Cape Town's exclusive
    # commuting zone (remember that TAZ data has a wider scope) + semi-rural
    # Durbanville
    selectedCenters[xCoord > -10] = np.zeros(1, 'bool')
    selectedCenters[yCoord > -3719] = np.zeros(1, 'bool')
    selectedCenters[(xCoord > -20) & (yCoord > -3765)] = np.zeros(1, 'bool')
    selectedCenters[codeCentersInitial == 1010] = np.zeros(1, 'bool')
    selectedCenters[codeCentersInitial == 1012] = np.zeros(1, 'bool')
    selectedCenters[codeCentersInitial == 1394] = np.zeros(1, 'bool')
    selectedCenters[codeCentersInitial == 1499] = np.zeros(1, 'bool')
    selectedCenters[codeCentersInitial == 4703] = np.zeros(1, 'bool')

    # Number of workers per income group in the model for the selected centers
    jobsCentersNgroup = np.zeros((len(xCoord), param["nb_of_income_classes"]))
    for j in range(0, param["nb_of_income_classes"]):
        jobsCentersNgroup[:, j] = np.sum(
            jobsCenters12Class[param["income_distribution"] == j + 1, :], 0)
    jobsCentersNgroup = jobsCentersNgroup[selectedCenters, :]

    # We rescale result to keep the correct global income distribution:
    # more specifically, this allows to go from individual to household level
    jobsCentersNGroupRescaled = (
        jobsCentersNgroup * households_per_income_class[None, :]
        / np.nansum(jobsCentersNgroup, 0)
        )

    return jobsCentersNGroupRescaled
