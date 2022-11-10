# -*- coding: utf-8 -*-
"""



"""

import numpy as np
import pandas as pd
import copy
import scipy.io


def import_exog_amenities(path_data, path_precalc_inp, dim):
    """
    Import relevant amenity data at SP level.

    Parameters
    ----------
    path_data : str
        Path towards data used in the model
    path_precalc_inp : str
        Path for precalcuted input data (calibrated parameters)
    dim : str
        Geographic level of analysis at which we want to run the commuting
        choice model: should be set to "grid" or "SP"

    Returns
    -------
    table_amenities : DataFrame
        Table yielding, for each selected geographic unit, a set of dummy
        variables corresponding to available exogenous amenites

    """
    # Load amenity file
    amenity_data = pd.read_csv(path_data + dim + '_amenities.csv', sep=',')

    # We replace values for airport cone amenities with a dummy for being
    # located in the airport cone
    airport_cone = copy.deepcopy(amenity_data.airport_cone)
    airport_cone[airport_cone == 55] = 1
    airport_cone[airport_cone == 60] = 1
    airport_cone[airport_cone == 65] = 1
    airport_cone[airport_cone == 70] = 1
    airport_cone[airport_cone == 75] = 1

    # Load distance to RDP house dummies
    distance_RDP = scipy.io.loadmat(
        path_precalc_inp + dim + 'DistanceRDP.mat'
        )
    distance_RDP = distance_RDP[list(distance_RDP)[-1]].squeeze()

    # We store relevant data in an output table
    # NB: we only consider dummies for amenity data crossing some thresholds.
    # This is done by trial and error to simplify calibration process
    table_amenities = pd.DataFrame(
        data=np.transpose(np.array(
            [amenity_data.distance_distr_parks < 2,
             amenity_data.distance_ocean < 2,
             ((amenity_data.distance_ocean > 2)
              & (amenity_data.distance_ocean < 4)),
             amenity_data.distance_world_herit < 2,
             ((amenity_data.distance_world_herit > 2)
              & (amenity_data.distance_world_herit < 4)),
             amenity_data.distance_urban_herit < 2,
             amenity_data.distance_UCT < 2,
             airport_cone,
             ((amenity_data.slope > 1) & (amenity_data.slope < 5)),
             amenity_data.slope > 5,
             amenity_data.distance_train < 2,
             amenity_data.distance_protected_envir < 2,
             ((amenity_data.distance_protected_envir > 2)
              & (amenity_data.distance_protected_envir < 4)),
             distance_RDP,
             amenity_data.distance_power_station < 2,
             amenity_data.distance_biosphere_reserve < 2])
            ),
        columns=['distance_distr_parks', 'distance_ocean',
                 'distance_ocean_2_4', 'distance_world_herit',
                 'distance_world_herit_2_4', 'distance_urban_herit',
                 'distance_UCT', 'airport_cone2', 'slope_1_5', 'slope_5',
                 'distance_train', 'distance_protected_envir',
                 'distance_protected_envir_2_4', 'RDP_proximity',
                 'distance_power_station', 'distance_biosphere_reserve']
        )
    # these new columns have been added to the amenities for testing
    newAmmenities=['publicHealthCare', 'privateHealthcare', 'sportsgrounds',
       'publicHighSchool', 'independentHighSchool', 'HighmastLight']
    table_amenities[newAmmenities]=amenity_data[newAmmenities].astype(int)

    return table_amenities
