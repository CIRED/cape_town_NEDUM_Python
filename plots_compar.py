# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 16:49:54 2022

@author: monni
"""

# %% Preamble

# IMPORT PACKAGES

import os
import itertools
import numpy as np
import pandas as pd
import geopandas as gpd

import plotly.io as pio
import plotly.express as px
import plotly.graph_objects as go

import inputs.parameters_and_options as inpprm


# SET RELATED OPTIONS

pio.renderers.default = 'svg'
# pio.renderers.default = 'browser'


# DEFINE FILE PATHS

path_code = '..'
path_folder = path_code + '/Data/'
path_precalc_inp = path_folder + 'precalculated_inputs/'
path_data = path_folder + 'data_Cape_Town/'
path_precalc_transp = path_folder + 'precalculated_transport/'
path_scenarios = path_data + 'Scenarios/'
path_outputs = path_code + '/Output/'
path_floods = path_folder + "flood_maps/"


# WE CREATE DIRECTORIES TO STORE OUTPUTS (IF NEEDED)

path_use_case = path_outputs + 'use_case_insur/'
path_maps = path_use_case + 'maps/'
path_charts = path_use_case + 'charts/'
path_damage_distrib = path_charts + 'damage_distrib/'
path_maps_abs_damages = path_maps + 'abs_damages/'
path_maps_rel_damages = path_maps + 'rel_damages/'
path_maps_pop_distrib = path_maps + 'pop_distrib/'
path_maps_rent_distrib = path_maps + 'rent_distrib/'

path_maps_rel_damages_formal = path_maps_rel_damages + 'formal/'
path_maps_rel_damages_subsidized = path_maps_rel_damages + 'subsidized/'
path_maps_rel_damages_informal = path_maps_rel_damages + 'informal/'
path_maps_rel_damages_backyard = path_maps_rel_damages + 'backyard/'

path_maps_rent_distrib_formal = path_maps_rent_distrib + 'formal/'
path_maps_rent_distrib_subsidized = path_maps_rent_distrib + 'subsidized/'
path_maps_rent_distrib_informal = path_maps_rent_distrib + 'informal/'
path_maps_rent_distrib_backyard = path_maps_rent_distrib + 'backyard/'

list_output_paths = [
    path_use_case, path_maps, path_maps_abs_damages, path_maps_rel_damages,
    path_maps_pop_distrib, path_maps_rent_distrib,
    path_maps_rel_damages_formal, path_maps_rel_damages_subsidized,
    path_maps_rel_damages_informal, path_maps_rel_damages_backyard,
    path_maps_rent_distrib_formal, path_maps_rent_distrib_subsidized,
    path_maps_rent_distrib_informal, path_maps_rent_distrib_backyard,
    path_charts, path_damage_distrib]

for path in list_output_paths:
    try:
        os.mkdir(path)
    except OSError as error:
        print(error)


# IMPORT PARAMETER AND OPTIONS

options = inpprm.import_options()
param = inpprm.import_param(
    path_precalc_inp, options)


# DEFINE SCENARIOS TO IMPORT

list_scenarios = ['floods00_F0_P11_C10_scenario232',
                  'floods10_F0_P11_C10_scenario232']


# %% Imports

print("Import data")

# IMPORT DAMAGE TABLES

# #For aggregate damages

# NB: We do not import "data" tables since it is not validation data per se,
# only damages computed from the spatial distribution of households we imputed
# from SAL data. We therefore rely on non-flood related validation exercises,
# not only of population distribution but also other relevant endogenous
# variables, to proxy for flood damage validation

# ##Without insurance
damages_fluvialu_sim_noinsur = pd.read_excel(
    path_outputs + 'floods00_F0_P11_C10_scenario232' + '/damage_data.xlsx',
    'fluvialu_sim', index_col=0)
damages_pluvial_sim_noinsur = pd.read_excel(
    path_outputs + 'floods00_F0_P11_C10_scenario232' + '/damage_data.xlsx',
    'pluvial_sim', index_col=0)
damages_coastal_sim_noinsur = pd.read_excel(
    path_outputs + 'floods00_F0_P11_C10_scenario232' + '/damage_data.xlsx',
    'coastal_sim', index_col=0)

# ##With insurance
damages_fluvialu_sim_insur = pd.read_excel(
    path_outputs + 'floods10_F0_P11_C10_scenario232' + '/damage_data.xlsx',
    'fluvialu_sim', index_col=0)
damages_pluvial_sim_insur = pd.read_excel(
    path_outputs + 'floods10_F0_P11_C10_scenario232' + '/damage_data.xlsx',
    'pluvial_sim', index_col=0)
damages_coastal_sim_insur = pd.read_excel(
    path_outputs + 'floods10_F0_P11_C10_scenario232' + '/damage_data.xlsx',
    'coastal_sim', index_col=0)


# #For spatial distribution

# We retrieve all available dimension combinations in the data
flood_types = ['fluvialu', 'pluvial', 'coastal']
housing_types = ['formal', 'subsidized', 'informal', 'backyard']
damage_types = ['structure', 'content']
list_dim = [flood_types, housing_types, damage_types]
all_dim = list(itertools.product(*list_dim))

# We store them in a dictionary

dict_damage_map_noinsur = {}
for dim in all_dim:
    table = pd.read_csv(
        path_outputs + 'floods00_F0_P11_C10_scenario232/tables/floods/'
        + dim[0] + '_' + dim[1] + '_' + dim[2] + '_2d_sim.csv',
        names=['damage'], header=None)
    dict_damage_map_noinsur[dim[0] + '_' + dim[1] + '_' + dim[2]] = table

dict_damage_map_insur = {}
for dim in all_dim:
    table = pd.read_csv(
        path_outputs + 'floods10_F0_P11_C10_scenario232/tables/floods/'
        + dim[0] + '_' + dim[1] + '_' + dim[2] + '_2d_sim.csv',
        names=['damage'], header=None)
    dict_damage_map_insur[dim[0] + '_' + dim[1] + '_' + dim[2]] = table

dict_damage_map_noinsur_shareinc = {}
for dim in all_dim:
    table = pd.read_csv(
        path_outputs + 'floods00_F0_P11_C10_scenario232/tables/floods/'
        + dim[0] + '_' + dim[1] + '_' + dim[2] + '_2d_sim_shareinc.csv',
        names=['damage'], header=None)
    dict_damage_map_noinsur_shareinc[dim[0] +
                                     '_' + dim[1] + '_' + dim[2]] = table

dict_damage_map_insur_shareinc = {}
for dim in all_dim:
    table = pd.read_csv(
        path_outputs + 'floods10_F0_P11_C10_scenario232/tables/floods/'
        + dim[0] + '_' + dim[1] + '_' + dim[2] + '_2d_sim_shareinc.csv',
        names=['damage'], header=None)
    dict_damage_map_insur_shareinc[dim[0] +
                                   '_' + dim[1] + '_' + dim[2]] = table

dict_damage_map = {}
dict_damage_map["noinsur"] = dict_damage_map_noinsur
dict_damage_map["insur"] = dict_damage_map_insur
dict_damage_map["noinsur_shareinc"] = dict_damage_map_noinsur_shareinc
dict_damage_map["insur_shareinc"] = dict_damage_map_insur_shareinc


# Finally, we import key equilibrium outcomes

# We import calibrated income net of commuting costs
income_net_of_commuting_costs = np.load(
    path_precalc_transp + 'GRID_incomeNetOfCommuting_0.npy')

# We retrieve all available dimension combinations for calibrated fraction
# of capital destroyed by floods
housing_types_2 = [
    'formal', 'backyard', 'informal', 'subsidized',
    'formal_1', 'formal_2', 'subsidized_1',
    'informal_settlements', 'informal_backyards']
damage_types_2 = ['structure', 'contents']
list_dim = [damage_types_2, housing_types_2]
all_dim = list(itertools.product(*list_dim))

# We store them in a dictionary
dict_fract_K_destroyed = {}
for dim in all_dim:
    try:
        table = pd.read_csv(
            path_outputs + 'input_tables/'
            + dim[0] + '_' + dim[1] + '_fract_K_destroyed.csv',
            usecols=[1])
        dict_fract_K_destroyed[dim[0] + '_' + dim[1]] = table
    except FileNotFoundError:
        pass

# Then we get endogenous outcomes that change across scenarios
dict_scenario_outcomes = {}

for scenario in list_scenarios:
    initial_state_utility = np.load(
        path_outputs + scenario + '/initial_state_utility.npy')
    initial_state_households_housing_types = np.load(
        path_outputs + scenario
        + '/initial_state_households_housing_types.npy')
    initial_state_household_centers = np.load(
        path_outputs + scenario + '/initial_state_household_centers.npy')
    initial_state_households = np.load(
        path_outputs + scenario + '/initial_state_households.npy')
    initial_state_dwelling_size = np.load(
        path_outputs + scenario + '/initial_state_dwelling_size.npy')
    initial_state_housing_supply = np.load(
        path_outputs + scenario + '/initial_state_housing_supply.npy')
    initial_state_rent = np.load(
        path_outputs + scenario + '/initial_state_rent.npy')

    dict_scenario_outcomes[scenario] = {}
    dict_scenario_outcomes[scenario]['utility'] = initial_state_utility
    dict_scenario_outcomes[scenario]['hh_per_htype'
                                     ] = initial_state_households_housing_types
    dict_scenario_outcomes[scenario]['hh_per_incgroup'
                                     ] = initial_state_household_centers
    dict_scenario_outcomes[scenario]['hh_full_dist'] = initial_state_households
    dict_scenario_outcomes[scenario]['dwelling_size'
                                     ] = initial_state_dwelling_size
    dict_scenario_outcomes[scenario]['hsupply'] = initial_state_housing_supply
    dict_scenario_outcomes[scenario]['rent'] = initial_state_rent

    households_formal_inc = pd.DataFrame(initial_state_households[0, :, :].T)
    households_backyard_inc = pd.DataFrame(initial_state_households[1, :, :].T)
    households_informal_inc = pd.DataFrame(initial_state_households[2, :, :].T)
    households_subsidized_inc = pd.DataFrame(
        initial_state_households[3, :, :].T)
    list_var_temp = [households_formal_inc, households_backyard_inc,
                     households_informal_inc, households_subsidized_inc]

    for var in list_var_temp:
        var.loc[var[0] != 0, 0] = 1
        var.loc[var[1] != 0, 1] = 2
        var.loc[var[2] != 0, 2] = 3
        var.loc[var[3] != 0, 3] = 4

    households_formal_inc_flat = households_formal_inc.sum(axis=1)
    households_backyard_inc_flat = households_backyard_inc.sum(axis=1)
    households_informal_inc_flat = households_informal_inc.sum(axis=1)
    households_subsidized_inc_flat = households_subsidized_inc.sum(axis=1)

    dict_scenario_outcomes[scenario]['formal_incgroup'
                                     ] = households_formal_inc_flat
    dict_scenario_outcomes[scenario]['backyard_incgroup'
                                     ] = households_backyard_inc_flat
    dict_scenario_outcomes[scenario]['informal_incgroup'
                                     ] = households_informal_inc_flat
    dict_scenario_outcomes[scenario]['subsidized_incgroup'
                                     ] = households_subsidized_inc_flat

    for var in list_var_temp:
        var.loc[var[0] != 0, 0] = 1
        var.loc[var[1] != 0, 1] = 1
        var.loc[var[2] != 0, 2] = 1
        var.loc[var[3] != 0, 3] = 1

    net_income_formal = np.nansum(
        income_net_of_commuting_costs * households_formal_inc.T, 0)
    net_income_subsidized = np.nansum(
        income_net_of_commuting_costs * households_subsidized_inc.T, 0)
    net_income_informal = np.nansum(
        income_net_of_commuting_costs * households_informal_inc.T, 0)
    net_income_backyard = np.nansum(
        income_net_of_commuting_costs * households_backyard_inc.T, 0)

    dict_scenario_outcomes[scenario]['net_income_formal'
                                     ] = net_income_formal
    dict_scenario_outcomes[scenario]['net_income_subsidized'
                                     ] = net_income_subsidized
    dict_scenario_outcomes[scenario]['net_income_informal'
                                     ] = net_income_informal
    dict_scenario_outcomes[scenario]['net_income_backyard'
                                     ] = net_income_backyard


# We also import the geography reference grid for plots
geo_grid = gpd.read_file(path_data + "grid_reference_500.shp")
geo_grid.to_crs(4326, inplace=True)
geo_grid['lon'] = geo_grid.centroid.x
geo_grid['lat'] = geo_grid.centroid.y

# NB: Should we show the distribution of exposure to flood risks by income
# group, or is it redundant with the distribution of damages? Not necessarily,
# since equivalent damages can be obtained in high and low-risk areas
# depending on the value of the underlying capital


# %% Data processing

print("Data processing")

# REARRANGE TABLES

# #For aggregate damages

# We store all tables in a list for looping over
table_list = [damages_fluvialu_sim_noinsur, damages_pluvial_sim_noinsur,
              damages_coastal_sim_noinsur, damages_fluvialu_sim_insur,
              damages_pluvial_sim_insur, damages_coastal_sim_insur]

# We give labels to index and variables corresponding to subdimensions in
# combined graphs that are not easily modifiable ex post
for table in table_list:
    table.rename(
        index={'formal': 'Formal private',
               'subsidized': 'Formal subsidized',
               'informal': 'Informal settlements',
               'backyard': 'Informal backyards'},
        columns={'struct_damage': 'Structures', 'content_damage': 'Contents'},
        inplace=True)

# We also adapt naming conventions to subsequent combined graphs
newnames_noinsur = {'Structures': 'Structures (w/o/ insurance)',
                    'Contents': 'Contents (w/o/ insurance)'}
newnames_insur = {'Structures': 'Structures (w/ insurance)',
                  'Contents': 'Contents (w/ insurance)'}


# #For spatial distribution

# We build an aggregate damage table to show on the map, and will display all
# the additional information as hover data

list_maps = ['noinsur', 'insur']
damage_maps = {}
list_maps_shareinc = ['noinsur_shareinc', 'insur_shareinc']
damage_maps_shareinc = {}

for item in list_maps:
    # We first concatenate all damage tables into one data frame
    concat_damage_map = (
        pd.concat(dict_damage_map[item], axis=1).groupby(axis=1, level=0).sum()
    )

    # We also define filters for column names to be used as part of processing
    filter_names = ['fluvialu', 'pluvial', 'coastal', 'fluvialu_formal',
                    'fluvialu_subsidized', 'fluvialu_informal',
                    'fluvialu_backyard', 'pluvial_formal',
                    'pluvial_subsidized', 'pluvial_informal',
                    'pluvial_backyard', 'coastal_formal', 'coastal_subsidized',
                    'coastal_informal', 'coastal_backyard']
    subfilter_names = ['fluvialu_formal', 'fluvialu_subsidized',
                       'fluvialu_informal', 'fluvialu_backyard',
                       'pluvial_formal', 'pluvial_subsidized',
                       'pluvial_informal', 'pluvial_backyard',
                       'coastal_formal', 'coastal_subsidized',
                       'coastal_informal', 'coastal_backyard']

    # We get the columns associated with each column name and store them in a
    # dictionary
    dict_filters = {}
    for filtre in filter_names:
        dict_filters[filtre] = [
            col for col in concat_damage_map if col.startswith(filtre)]

    # We obtain the sum of all damage values corresponding to each data categ
    for key, value in dict_filters.items():
        concat_damage_map['sum_' + key] = (
            concat_damage_map[dict_filters[key]].sum(axis=1))

    # We take the maximum across aggregate flood type damages to display as our
    # main outcome of interest (cf. bath-tub model)
    concat_damage_map['max_val'] = concat_damage_map[
        ['sum_fluvialu', 'sum_pluvial', 'sum_coastal']].max(axis=1)

    # We also define the dominant flood type to display as hover data
    concat_damage_map['flood_type'] = 'None'
    concat_damage_map.loc[
        (concat_damage_map['sum_fluvialu'] == concat_damage_map['max_val'])
        & (concat_damage_map['max_val'] != 0), 'flood_type'] = 'fluvial'
    concat_damage_map.loc[
        (concat_damage_map['sum_pluvial'] == concat_damage_map['max_val'])
        & (concat_damage_map['max_val'] != 0), 'flood_type'] = 'pluvial'
    concat_damage_map.loc[
        (concat_damage_map['sum_coastal'] == concat_damage_map['max_val'])
        & (concat_damage_map['max_val'] != 0), 'flood_type'] = 'coastal'

    # Then, we recover the breakdown of damages across housing types, to show
    # as hover data as well

    # We initialize values for each housing type
    concat_damage_map['damage_formal'] = 0
    concat_damage_map['damage_subsidized'] = 0
    concat_damage_map['damage_informal'] = 0
    concat_damage_map['damage_backyard'] = 0

    # We associate each value to the corresponding sum of damages for the main
    # flood type only
    for housing_type in housing_types:
        concat_damage_map.loc[
            concat_damage_map['flood_type']
            == 'fluvial', 'damage_' + housing_type
        ] = concat_damage_map['sum_fluvialu_' + housing_type]
        concat_damage_map.loc[
            concat_damage_map['flood_type']
            == 'pluvial', 'damage_' + housing_type
        ] = concat_damage_map['sum_pluvial_' + housing_type]
        concat_damage_map.loc[
            concat_damage_map['flood_type']
            == 'coastal', 'damage_' + housing_type
        ] = concat_damage_map['sum_coastal_' + housing_type]

    # Finally, we also define the share of content (vs. structures) damages in
    # overall housing-type-specific damages, to display as hover data

    # We start by defining the ratio for all flood + housing type combinations
    for subfilter in subfilter_names:
        concat_damage_map[subfilter + '_content_share'] = (
            concat_damage_map[subfilter + '_content']
            / concat_damage_map['sum_' + subfilter])

    # We initialize values for housing types only
    concat_damage_map['content_share_formal'] = 0
    concat_damage_map['content_share_subsidized'] = 0
    concat_damage_map['content_share_informal'] = 0
    concat_damage_map['content_share_backyard'] = 0

    # Then, we associate each value to the associated content ratio for the
    # dominant flood type only
    for housing_type in housing_types:
        concat_damage_map.loc[
            concat_damage_map['flood_type'] == 'fluvial',
            'content_share_' + housing_type
        ] = concat_damage_map['fluvialu_' + housing_type +
                              '_content_share']
        concat_damage_map.loc[
            concat_damage_map['flood_type'] == 'pluvial',
            'content_share_' + housing_type
        ] = concat_damage_map['pluvial_' + housing_type + '_content_share']
        concat_damage_map.loc[
            concat_damage_map['flood_type'] == 'coastal',
            'content_share_' + housing_type
        ] = concat_damage_map['coastal_' + housing_type + '_content_share']

    # We convert nans to zeros for proper rendering in output map
    concat_damage_map = concat_damage_map.fillna(0)
    concat_damage_map = concat_damage_map.replace([np.inf, -np.inf], 0)

    # We store the output in the corresponding dictionary entry
    damage_maps[item] = concat_damage_map


# We do the same expressing damages as a share of net income

for item in list_maps_shareinc:
    # We first concatenate all damage tables into one data frame
    concat_damage_map = (
        pd.concat(dict_damage_map[item], axis=1).groupby(axis=1, level=0).sum()
    )

    # We also define filters for column names to be used as part of processing
    subfilter_names = ['fluvialu_formal', 'fluvialu_subsidized',
                       'fluvialu_informal', 'fluvialu_backyard',
                       'pluvial_formal', 'pluvial_subsidized',
                       'pluvial_informal', 'pluvial_backyard',
                       'coastal_formal', 'coastal_subsidized',
                       'coastal_informal', 'coastal_backyard']

    # We get the columns associated with each column name and store them in a
    # dictionary
    dict_filters = {}
    for filtre in subfilter_names:
        dict_filters[filtre] = [
            col for col in concat_damage_map if col.startswith(filtre)]

    # We obtain the sum of all damage values corresponding to each data categ
    for key, value in dict_filters.items():
        concat_damage_map['sum_' + key] = (
            concat_damage_map[dict_filters[key]].sum(axis=1))

    # We take the maximum across aggregate flood type damages to display as our
    # main outcome of interest (cf. bath-tub model)
    concat_damage_map['max_val_formal'] = concat_damage_map[
        ['sum_fluvialu_formal', 'sum_pluvial_formal', 'sum_coastal_formal']
    ].max(axis=1)
    concat_damage_map['max_val_subsidized'] = concat_damage_map[
        ['sum_fluvialu_subsidized', 'sum_pluvial_subsidized',
         'sum_coastal_subsidized']
    ].max(axis=1)
    concat_damage_map['max_val_informal'] = concat_damage_map[
        ['sum_fluvialu_informal', 'sum_pluvial_informal',
         'sum_coastal_informal']
    ].max(axis=1)
    concat_damage_map['max_val_backyard'] = concat_damage_map[
        ['sum_fluvialu_backyard', 'sum_pluvial_backyard',
         'sum_coastal_backyard']
    ].max(axis=1)

    # We do the same for content share only (to display as hover data)
    concat_damage_map['max_content_formal'] = concat_damage_map[
        ['fluvialu_formal_content', 'pluvial_formal_content',
         'coastal_formal_content']
    ].max(axis=1)
    concat_damage_map['max_content_subsidized'] = concat_damage_map[
        ['fluvialu_subsidized_content', 'pluvial_subsidized_content',
            'coastal_subsidized_content']
    ].max(axis=1)
    concat_damage_map['max_content_informal'] = concat_damage_map[
        ['fluvialu_informal_content', 'pluvial_informal_content',
            'coastal_informal_content']
    ].max(axis=1)
    concat_damage_map['max_content_backyard'] = concat_damage_map[
        ['fluvialu_backyard_content', 'pluvial_backyard_content',
            'coastal_backyard_content']
    ].max(axis=1)

    # We also define the dominant flood type to display as hover data
    for housing_type in housing_types:
        concat_damage_map['flood_type_' + housing_type] = 'None'
        concat_damage_map.loc[
            (concat_damage_map['sum_fluvialu_' + housing_type]
             == concat_damage_map['max_val_' + housing_type])
            & (concat_damage_map['max_val_' + housing_type] != 0),
            'flood_type_' + housing_type] = 'fluvial'
        concat_damage_map.loc[
            (concat_damage_map['sum_pluvial_' + housing_type]
             == concat_damage_map['max_val_' + housing_type])
            & (concat_damage_map['max_val_' + housing_type] != 0),
            'flood_type_' + housing_type] = 'pluvial'
        concat_damage_map.loc[
            (concat_damage_map['sum_coastal_' + housing_type]
             == concat_damage_map['max_val_' + housing_type])
            & (concat_damage_map['max_val_' + housing_type] != 0),
            'flood_type_' + housing_type] = 'coastal'

    # We also add some population info to display as hover data
    if item == 'noinsur_shareinc':
        concat_damage_map['nb_households_formal'] = (
            dict_scenario_outcomes[list_scenarios[0]]['hh_per_htype'][0, :])
        concat_damage_map['nb_households_subsidized'] = (
            dict_scenario_outcomes[list_scenarios[0]]['hh_per_htype'][3, :])
        concat_damage_map['nb_households_informal'] = (
            dict_scenario_outcomes[list_scenarios[0]]['hh_per_htype'][2, :])
        concat_damage_map['nb_households_backyard'] = (
            dict_scenario_outcomes[list_scenarios[0]]['hh_per_htype'][1, :])
        concat_damage_map['rent_formal'] = (
            dict_scenario_outcomes[list_scenarios[0]]['rent'][0, :])
        concat_damage_map['rent_formal'][
            concat_damage_map['nb_households_formal'] == 0] = 0
        concat_damage_map['rent_subsidized'] = (
            dict_scenario_outcomes[list_scenarios[0]]['rent'][3, :])
        concat_damage_map['rent_subsidized'][
            concat_damage_map['nb_households_subsidized'] == 0] = 0
        concat_damage_map['rent_informal'] = (
            dict_scenario_outcomes[list_scenarios[0]]['rent'][2, :])
        concat_damage_map['rent_informal'][
            concat_damage_map['nb_households_informal'] == 0] = 0
        concat_damage_map['rent_backyard'] = (
            dict_scenario_outcomes[list_scenarios[0]]['rent'][1, :])
        concat_damage_map['rent_subsidized'][
            concat_damage_map['nb_households_subsidized'] == 0] = 0
        for housing_type in housing_types:
            concat_damage_map['incgroup_' + housing_type] = (
                dict_scenario_outcomes[list_scenarios[0]
                                       ][housing_type + '_incgroup'])
            concat_damage_map['net_income_' + housing_type] = (
                dict_scenario_outcomes[list_scenarios[0]
                                       ]['net_income_' + housing_type])

    elif item == 'insur_shareinc':
        concat_damage_map['nb_households_formal'] = (
            dict_scenario_outcomes[list_scenarios[1]]['hh_per_htype'][0, :])
        concat_damage_map['nb_households_subsidized'] = (
            dict_scenario_outcomes[list_scenarios[1]]['hh_per_htype'][3, :])
        concat_damage_map['nb_households_informal'] = (
            dict_scenario_outcomes[list_scenarios[1]]['hh_per_htype'][2, :])
        concat_damage_map['nb_households_backyard'] = (
            dict_scenario_outcomes[list_scenarios[1]]['hh_per_htype'][1, :])
        concat_damage_map['rent_formal'] = (
            dict_scenario_outcomes[list_scenarios[1]]['rent'][0, :])
        concat_damage_map['rent_formal'][
            concat_damage_map['nb_households_formal'] == 0] = 0
        concat_damage_map['rent_subsidized'] = (
            dict_scenario_outcomes[list_scenarios[1]]['rent'][3, :])
        concat_damage_map['rent_subsidized'][
            concat_damage_map['nb_households_subsidized'] == 0] = 0
        concat_damage_map['rent_informal'] = (
            dict_scenario_outcomes[list_scenarios[1]]['rent'][2, :])
        concat_damage_map['rent_informal'][
            concat_damage_map['nb_households_informal'] == 0] = 0
        concat_damage_map['rent_backyard'] = (
            dict_scenario_outcomes[list_scenarios[1]]['rent'][1, :])
        concat_damage_map['rent_backyard'][
            concat_damage_map['nb_households_backyard'] == 0] = 0
        for housing_type in housing_types:
            concat_damage_map['incgroup_' + housing_type] = (
                dict_scenario_outcomes[list_scenarios[1]
                                       ][housing_type + '_incgroup'])
            concat_damage_map['net_income_' + housing_type] = (
                dict_scenario_outcomes[list_scenarios[1]
                                       ]['net_income_' + housing_type])

    # We convert nans to zeros for proper rendering in output map
    concat_damage_map = concat_damage_map.fillna(0)
    concat_damage_map = concat_damage_map.replace([np.inf, -np.inf], 0)

    # We store the output in the corresponding dictionary entry
    damage_maps_shareinc[item] = concat_damage_map


# We also create a data table for other equilibrium outcomes
equil_maps = {}

for item in list_scenarios:

    # For overall distribution

    hh_tot = np.nansum(dict_scenario_outcomes[item]['hh_per_htype'], 0)
    hh_formal = dict_scenario_outcomes[item]['hh_per_htype'][0, :]
    hh_backyard = dict_scenario_outcomes[item]['hh_per_htype'][1, :]
    hh_informal = dict_scenario_outcomes[item]['hh_per_htype'][2, :]
    hh_subsidized = dict_scenario_outcomes[item]['hh_per_htype'][3, :]

    formal_incgroup = dict_scenario_outcomes[item]['formal_incgroup']
    backyard_incgroup = dict_scenario_outcomes[item]['backyard_incgroup']
    informal_incgroup = dict_scenario_outcomes[item]['informal_incgroup']
    subsidized_incgroup = dict_scenario_outcomes[item]['subsidized_incgroup']

    # + for rent distribution (per housing type)

    hsupply_formal = dict_scenario_outcomes[item]['hsupply'][0, :]
    hsupply_backyard = dict_scenario_outcomes[item]['hsupply'][1, :]
    hsupply_informal = dict_scenario_outcomes[item]['hsupply'][2, :]
    hsupply_subsidized = dict_scenario_outcomes[item]['hsupply'][3, :]

    rent_formal = dict_scenario_outcomes[item]['rent'][0, :]
    rent_formal[hh_formal == 0] = 0
    rent_backyard = dict_scenario_outcomes[item]['rent'][1, :]
    rent_backyard[hh_backyard == 0] = 0
    rent_informal = dict_scenario_outcomes[item]['rent'][2, :]
    rent_informal[hh_informal == 0] = 0
    rent_subsidized = dict_scenario_outcomes[item]['rent'][3, :]
    rent_subsidized[hh_subsidized == 0] = 0

    dwelling_size_formal = dict_scenario_outcomes[item]['dwelling_size'][0, :]
    dwelling_size_formal[hh_formal == 0] = 0
    dwelling_size_backyard = (
        dict_scenario_outcomes[item]['dwelling_size'][1, :])
    dwelling_size_backyard[hh_backyard == 0] = 0
    dwelling_size_informal = (
        dict_scenario_outcomes[item]['dwelling_size'][2, :])
    dwelling_size_informal[hh_informal == 0] = 0
    dwelling_size_subsidized = (
        dict_scenario_outcomes[item]['dwelling_size'][3, :])
    dwelling_size_subsidized[hh_subsidized == 0] = 0

    # Note that content depreciation is in fact identical across htypes
    flood_deprec_content_formal = dict_fract_K_destroyed['contents_formal']
    flood_deprec_content_backyard = dict_fract_K_destroyed['contents_backyard']
    flood_deprec_content_informal = dict_fract_K_destroyed['contents_informal']
    flood_deprec_content_subsidized = (
        dict_fract_K_destroyed['contents_subsidized'])

    flood_deprec_struct_formal = dict_fract_K_destroyed['structure_formal_1']
    flood_deprec_struct_formal[dwelling_size_formal > param["threshold"]] = (
        dict_fract_K_destroyed['structure_formal_2'][
            dwelling_size_formal > param["threshold"]])
    flood_deprec_struct_backyard = dict_fract_K_destroyed[
        'structure_informal_backyards']
    flood_deprec_struct_informal = dict_fract_K_destroyed[
        'structure_informal_settlements']
    flood_deprec_struct_subsidized = dict_fract_K_destroyed[
        'structure_subsidized_1']

    # We store everything into one data frame

    concat_equil_map = pd.DataFrame({
        'hh_tot': hh_tot,
        'hh_formal': hh_formal,
        'hh_backyard': hh_backyard,
        'hh_informal': hh_informal,
        'hh_subsidized': hh_subsidized,
        'formal_incgroup': formal_incgroup,
        'backyard_incgroup': backyard_incgroup,
        'informal_incgroup': informal_incgroup,
        'subsidized_incgroup': subsidized_incgroup,
        'hsupply_formal': hsupply_formal,
        'hsupply_backyard': hsupply_backyard,
        'hsupply_informal': hsupply_informal,
        'hsupply_subsidized': hsupply_subsidized,
        'rent_formal': rent_formal,
        'rent_backyard': rent_backyard,
        'rent_informal': rent_informal,
        'rent_subsidized': rent_subsidized,
        'dwelling_size_formal': dwelling_size_formal,
        'dwelling_size_backyard': dwelling_size_backyard,
        'dwelling_size_informal': dwelling_size_informal,
        'dwelling_size_subsidized': dwelling_size_subsidized,
        'flood_deprec_content_formal':
            flood_deprec_content_formal.to_numpy().flatten(),
        'flood_deprec_content_backyard':
            flood_deprec_content_backyard.to_numpy().flatten(),
        'flood_deprec_content_informal':
            flood_deprec_content_informal.to_numpy().flatten(),
        'flood_deprec_content_subsidized':
            flood_deprec_content_subsidized.to_numpy().flatten(),
        'flood_deprec_struct_formal':
            flood_deprec_struct_formal.to_numpy().flatten(),
        'flood_deprec_struct_backyard':
            flood_deprec_struct_backyard.to_numpy().flatten(),
        'flood_deprec_struct_informal':
            flood_deprec_struct_informal.to_numpy().flatten(),
        'flood_deprec_struct_subsidized':
            flood_deprec_struct_subsidized.to_numpy().flatten()
    })

    # We convert nans to zeros for proper rendering in output map
    concat_equil_map = concat_equil_map.fillna(0)
    concat_equil_map = concat_equil_map.replace([np.inf, -np.inf], 0)

    # We store the output in the corresponding dictionary entry
    equil_maps[item] = concat_equil_map


# #Then, we create comparative maps across scenarios to better see changes

# We store name of variables to be created
list_num_var = ['max_val', 'damage_formal', 'damage_subsidized',
                'damage_informal', 'damage_backyard']
list_num_pct = ['max_val_pct', 'damage_formal_pct', 'damage_subsidized_pct',
                'damage_informal_pct', 'damage_backyard_pct']

# We create the table for changes
damage_map_compar = pd.DataFrame()
damage_map_compar[list_num_var] = (
    damage_maps['noinsur'][list_num_var] - damage_maps['insur'][list_num_var])
damage_map_compar[list_num_pct] = (
    damage_maps['noinsur'][list_num_var] / damage_maps['insur'][list_num_var]
    - 1)
damage_map_compar['lon'] = geo_grid.lon
damage_map_compar['lat'] = geo_grid.lat
damage_map_compar['flood_type'] = damage_maps['noinsur']['flood_type']
damage_map_compar.loc[
    damage_map_compar['flood_type'] == 'None', 'flood_type'
] = damage_maps['insur']['flood_type']

damage_map_compar = damage_map_compar.fillna(0)
damage_map_compar = damage_map_compar.replace([np.inf, -np.inf], 0)


# We do the same for damages expressed as a share of income

# We store name of variables to be created
list_num_var = [
    'max_val_formal', 'max_val_subsidized', 'max_val_informal',
    'max_val_backyard', 'max_content_formal', 'max_content_subsidized',
    'max_content_informal', 'max_content_backyard', 'rent_formal',
    'rent_subsidized', 'rent_informal', 'rent_backyard']

list_num_pct = [
    'max_val_formal_pct', 'max_val_subsidized_pct', 'max_val_informal_pct',
    'max_val_backyard_pct', 'max_content_formal_pct',
    'max_content_subsidized_pct', 'max_content_informal_pct',
    'max_content_backyard_pct', 'rent_formal_pct', 'rent_subsidized_pct',
    'rent_informal_pct', 'rent_backyard_pct']

# We create the table for changes

damage_map_compar_shareinc = pd.DataFrame()
damage_map_compar_shareinc[list_num_var] = (
    damage_maps_shareinc['noinsur_shareinc'][list_num_var]
    - damage_maps_shareinc['insur_shareinc'][list_num_var])
damage_map_compar_shareinc[list_num_pct] = (
    damage_maps_shareinc['noinsur_shareinc'][list_num_var] /
    damage_maps_shareinc['insur_shareinc'][list_num_var]
    - 1)
damage_map_compar_shareinc['lon'] = geo_grid.lon
damage_map_compar_shareinc['lat'] = geo_grid.lat

for housing_type in housing_types:
    damage_map_compar_shareinc['flood_type_' + housing_type] = (
        damage_maps_shareinc['noinsur_shareinc']['flood_type_' + housing_type])
    damage_map_compar_shareinc.loc[
        damage_map_compar_shareinc['flood_type_' + housing_type] == 'None',
        'flood_type_' + housing_type
    ] = damage_maps_shareinc['insur_shareinc']['flood_type_' + housing_type]
    # NB: we take the no-insurance case as a benchmark, since we want to show
    # surplus damages from no insurance. We'll rely on other maps to show
    # population moves and composition effects, and the extent to which they
    # explain what we observe
    damage_map_compar_shareinc['nb_households_' + housing_type] = (
        damage_maps_shareinc['noinsur_shareinc'][
            'nb_households_' + housing_type])
    damage_map_compar_shareinc['incgroup_' + housing_type] = (
        damage_maps_shareinc['noinsur_shareinc']['incgroup_' + housing_type])
    damage_map_compar_shareinc['net_income_' + housing_type] = (
        damage_maps_shareinc['noinsur_shareinc']['net_income_' + housing_type])

damage_map_compar_shareinc = damage_map_compar_shareinc.fillna(0)
damage_map_compar_shareinc = damage_map_compar_shareinc.replace(
    [np.inf, -np.inf], 0)


# Then for other equilibrium outcomes

# We store name of variables to be created
list_num_var = [
    'hh_tot', 'hh_formal', 'hh_backyard', 'hh_informal', 'hh_subsidized',
    'hsupply_formal', 'hsupply_backyard', 'hsupply_informal',
    'hsupply_subsidized',
    'rent_formal', 'rent_backyard', 'rent_informal', 'rent_subsidized',
    'dwelling_size_formal', 'dwelling_size_backyard', 'dwelling_size_informal',
    'dwelling_size_subsidized']
list_num_pct = [
    'hh_tot_pct', 'hh_formal_pct', 'hh_backyard_pct', 'hh_informal_pct',
    'hh_subsidized_pct',
    'hsupply_formal_pct', 'hsupply_backyard_pct', 'hsupply_informal_pct',
    'hsupply_subsidized_pct',
    'rent_formal_pct', 'rent_backyard_pct', 'rent_informal_pct',
    'rent_subsidized_pct',
    'dwelling_size_formal_pct', 'dwelling_size_backyard_pct',
    'dwelling_size_informal_pct', 'dwelling_size_subsidized_pct']

# We create the table for changes
equil_map_compar = pd.DataFrame()
equil_map_compar[list_num_var] = (
    equil_maps[list_scenarios[0]][list_num_var]
    - equil_maps[list_scenarios[1]][list_num_var])
equil_map_compar[list_num_pct] = (
    equil_maps[list_scenarios[0]][list_num_var]
    / equil_maps[list_scenarios[1]][list_num_var]
    - 1)
equil_map_compar['lon'] = geo_grid.lon
equil_map_compar['lat'] = geo_grid.lat

# NB: Again, we take the no-insurance case as a benchmark, while considering
# the income group of households who left when the area is left unoccupied
# in the benchmark case
for housing_type in housing_types:
    equil_map_compar['flood_deprec_content_' + housing_type] = (
        equil_maps[list_scenarios[0]]['flood_deprec_content_' + housing_type])
    equil_map_compar['flood_deprec_struct_' + housing_type] = (
        equil_maps[list_scenarios[0]]['flood_deprec_struct_' + housing_type])
    equil_map_compar[housing_type + '_incgroup'] = (
        equil_maps[list_scenarios[0]][housing_type + '_incgroup'])
    equil_map_compar.loc[
        equil_map_compar[housing_type + '_incgroup'] == 0,
        housing_type + '_incgroup'
        ] = equil_maps[list_scenarios[1]][housing_type + '_incgroup']

equil_map_compar = equil_map_compar.fillna(0)
equil_map_compar = equil_map_compar.replace([np.inf, -np.inf], 0)


# %% Output graphs

print("Output graphs")

# PLOT AGGREGATE DAMAGES

# #Fluvial floods

# We create the first plot with insurance
fig_fluvialu_insur = px.bar(
    damages_fluvialu_sim_insur,
    x=damages_fluvialu_sim_insur.index,
    y=['Structures', 'Contents'],
    barmode='group',
    color_discrete_sequence=px.colors.qualitative.Set1,
    title='Estimated annual damages from fluvial floods (in M rands, 2011)',
    labels={'index': 'Housing type', 'value': 'Annual damages',
            'variable': 'Damage type'},
    template='plotly_white')

# We adapt labels accordingly
fig_fluvialu_insur.for_each_trace(
    lambda t: t.update(
        name=newnames_insur[t.name],
        legendgroup=newnames_insur[t.name],
        hovertemplate=t.hovertemplate.replace(t.name, newnames_insur[t.name])
    ))

# fig_fluvialu_insur.show()

# We create the second plot without insurance (with a less opaque color for
# later superimposition, as this takes bigger values)
fig_fluvialu_noinsur = px.bar(
    damages_fluvialu_sim_noinsur,
    x=damages_fluvialu_sim_noinsur.index,
    y=['Structures', 'Contents'],
    barmode='group',
    color_discrete_sequence=px.colors.qualitative.Pastel1,
    title='Estimated annual damages from fluvial floods (in M rands, 2011)',
    labels={'index': 'Housing type', 'value': 'Annual damages',
            'variable': 'Damage type'},
    template='plotly_white')

# We adapt labels accordingly
fig_fluvialu_noinsur.for_each_trace(
    lambda t: t.update(
        name=newnames_noinsur[t.name],
        legendgroup=newnames_noinsur[t.name],
        hovertemplate=t.hovertemplate.replace(
            t.name, newnames_noinsur[t.name])
    ))

# fig_fluvialu_noinsur.show()

# We combine the two graphs and store the output both as an interactive and a
# static image
fig_fluvialu = go.Figure(fig_fluvialu_noinsur)
fig_fluvialu = fig_fluvialu.add_traces(fig_fluvialu_insur.data)
# fig_fluvialu.show()
fig_fluvialu.write_html(path_charts + "fluvialu_sim_damage_sum_insur.html")
fig_fluvialu.write_image(path_charts + "fluvialu_sim_damage_sum_insur.png")

print("Aggregate fluvial damages done")


# #Pluvial floods

# We create the first plot with insurance
fig_pluvial_insur = px.bar(
    damages_pluvial_sim_insur,
    x=damages_pluvial_sim_insur.index,
    y=['Structures', 'Contents'],
    barmode='group',
    color_discrete_sequence=px.colors.qualitative.Set1,
    title='Estimated annual damages from pluvial floods (in M rands, 2011)',
    labels={'index': 'Housing type', 'value': 'Annual damages',
            'variable': 'Damage type'},
    template='plotly_white')

# We adapt labels accordingly
fig_pluvial_insur.for_each_trace(
    lambda t: t.update(
        name=newnames_insur[t.name],
        legendgroup=newnames_insur[t.name],
        hovertemplate=t.hovertemplate.replace(t.name, newnames_insur[t.name])
    ))

# fig_pluvial_insur.show()

# We create the second plot without insurance (with a less opaque color for
# later superimposition, as this takes bigger values)
fig_pluvial_noinsur = px.bar(
    damages_pluvial_sim_noinsur,
    x=damages_pluvial_sim_noinsur.index,
    y=['Structures', 'Contents'],
    barmode='group',
    color_discrete_sequence=px.colors.qualitative.Pastel1,
    title='Estimated annual damages from pluvial floods (in M rands, 2011)',
    labels={'index': 'Housing type', 'value': 'Annual damages',
            'variable': 'Damage type'},
    template='plotly_white')

# We adapt labels accordingly
fig_pluvial_noinsur.for_each_trace(
    lambda t: t.update(
        name=newnames_noinsur[t.name],
        legendgroup=newnames_noinsur[t.name],
        hovertemplate=t.hovertemplate.replace(
            t.name, newnames_noinsur[t.name])
    ))

# fig_pluvial_noinsur.show()

# We combine the two graphs and store the output both as an interactive and a
# static image
fig_pluvial = go.Figure(fig_pluvial_noinsur)
fig_pluvial = fig_pluvial.add_traces(fig_pluvial_insur.data)
# fig_pluvial.show()
fig_pluvial.write_html(path_charts + "pluvial_sim_damage_sum_insur.html")
fig_pluvial.write_image(path_charts + "pluvial_sim_damage_sum_insur.png")

print("Aggregate pluvial damages done")


# #Coastal floods

# We create the first plot with insurance
fig_coastal_insur = px.bar(
    damages_coastal_sim_insur,
    x=damages_coastal_sim_insur.index,
    y=['Structures', 'Contents'],
    barmode='group',
    color_discrete_sequence=px.colors.qualitative.Set1,
    title='Estimated annual damages from coastal floods (in M rands, 2011)',
    labels={'index': 'Housing type', 'value': 'Annual damages',
            'variable': 'Damage type'},
    template='plotly_white')

# We adapt labels accordingly
fig_coastal_insur.for_each_trace(
    lambda t: t.update(
        name=newnames_insur[t.name],
        legendgroup=newnames_insur[t.name],
        hovertemplate=t.hovertemplate.replace(t.name, newnames_insur[t.name])
    ))

# fig_coastal_insur.show()

# We create the second plot without insurance (with a less opaque color for
# later superimposition, as this takes bigger values)
fig_coastal_noinsur = px.bar(
    damages_coastal_sim_noinsur,
    x=damages_coastal_sim_noinsur.index,
    y=['Structures', 'Contents'],
    barmode='group',
    color_discrete_sequence=px.colors.qualitative.Pastel1,
    title='Estimated annual damages from coastal floods (in M rands, 2011)',
    labels={'index': 'Housing type', 'value': 'Annual damages',
            'variable': 'Damage type'},
    template='plotly_white')

# We adapt labels accordingly
fig_coastal_noinsur.for_each_trace(
    lambda t: t.update(
        name=newnames_noinsur[t.name],
        legendgroup=newnames_noinsur[t.name],
        hovertemplate=t.hovertemplate.replace(
            t.name, newnames_noinsur[t.name])
    ))

# fig_coastal_noinsur.show()

# We combine the two graphs and store the output both as an interactive and a
# static image
fig_coastal = go.Figure(fig_coastal_noinsur)
fig_coastal = fig_coastal.add_traces(fig_coastal_insur.data)
# fig_coastal.show()
fig_coastal.write_html(path_charts + "coastal_sim_damage_sum_insur.html")
fig_coastal.write_image(path_charts + "coastal_sim_damage_sum_insur.png")

print("Aggregate coastal damages done")


# PLOT SPATIAL DAMAGE DISTRIBUTION

# We plot all the relevant information for a given scenario as a chloropleth
# map

# #First for damages in absolute values

for item in list_maps:

    damage_maps[item]['lon'] = geo_grid.lon
    damage_maps[item]['lat'] = geo_grid.lat

    fig = px.choropleth_mapbox(
        damage_maps[item],
        geojson=geo_grid.geometry,
        locations=geo_grid.index,
        color='max_val',
        center={"lat": -33.92345542582841, "lon": 18.434424141913478},
        zoom=9.25,
        mapbox_style='stamen-toner',
        opacity=0.75,
        labels={'lon': 'Lon.', 'lat': 'Lat.',
                'max_val': 'Total', 'locations': 'Pixel ID',
                'flood_type': 'Flood type',
                'damage_formal': 'Formal private',
                'damage_subsidized': 'Formal subsidized',
                'damage_informal': 'Informal settlements',
                'damage_backyard': 'Informal backyards',
                'content_share_formal': '% of content (formal)',
                'content_share_subsidized': '% of content (subsidized)',
                'content_share_informal': '% of content (informal)',
                'content_share_backyard': '% of content (backyard)'},
        title='Estimated annual flood damages (in rands, 2011)',
        color_continuous_scale="Reds",
        template='plotly_white',
        hover_data={'lon': ':.2f', 'lat': ':.2f', 'flood_type': True,
                    'damage_formal': ':,.0f',
                    'content_share_formal': ':.0%',
                    'damage_subsidized': ':,.0f',
                    'content_share_subsidized': ':.0%',
                    'damage_informal': ':,.0f',
                    'content_share_informal': ':.0%',
                    'damage_backyard': ':,.0f',
                    'content_share_backyard': ':.0%'})
    fig.update_layout(margin={"r": 0, "t": 30, "l": 0, "b": 0})
    fig.update_traces(marker_line_width=0)
    # fig.show()
    fig.write_html(path_maps_abs_damages + "map_damages_" + item + ".html")
    fig.write_image(path_maps_abs_damages + "map_damages_" + item + ".png")

    print("map_damages_" + item + " done")

fig = px.choropleth_mapbox(
    damage_map_compar,
    geojson=geo_grid.geometry,
    locations=geo_grid.index,
    color='max_val',
    center={"lat": -33.92345542582841, "lon": 18.434424141913478},
    zoom=9.25,
    mapbox_style='stamen-toner',
    opacity=0.75,
    labels={'lon': 'Lon.', 'lat': 'Lat.',
            'max_val': 'Total change', 'max_val_pct': '% change',
            'locations': 'Pixel ID',
            'flood_type': 'Flood type',
            'damage_formal': 'Formal private',
            'damage_formal_pct': '% change (formal)',
            'damage_subsidized': 'Formal subsidized',
            'damage_subsidized_pct': '% change (subsidized)',
            'damage_informal': 'Informal settlements',
            'damage_informal_pct': '% change (informal)',
            'damage_backyard': 'Informal backyards',
            'damage_backyard_pct': '% change (backyards)'},
    title='Flood damage surplus from no insurance (in rands, 2011)',
    color_continuous_scale="Picnic",
    color_continuous_midpoint=0,
    template='plotly_white',
    hover_data={'lon': ':.2f', 'lat': ':.2f', 'flood_type': True,
                'damage_formal': ':,.0f', 'damage_formal_pct': ':+.0%',
                'damage_subsidized': ':,.0f', 'damage_subsidized_pct': ':+.0%',
                'damage_informal': ':,.0f', 'damage_informal_pct': ':+.0%',
                'damage_backyard': ':,.0f', 'damage_backyard_pct': ':+.0%',
                'max_val': True, 'max_val_pct': ':+.0%'})
fig.update_layout(margin={"r": 0, "t": 30, "l": 0, "b": 0})
fig.update_traces(marker_line_width=0)
# fig.show()
fig.write_html(path_maps_abs_damages + "map_damages_compar_insur.html")
fig.write_image(path_maps_abs_damages + "map_damages_compar_insur.png")

print("map_damages_compar_insur done")


# #Then for damages in relative terms

for housing_type in housing_types:

    for item in list_maps_shareinc:

        damage_maps_shareinc[item]['lon'] = geo_grid.lon
        damage_maps_shareinc[item]['lat'] = geo_grid.lat

        fig = px.choropleth_mapbox(
            damage_maps_shareinc[item],
            geojson=geo_grid.geometry,
            locations=geo_grid.index,
            color='max_val_' + housing_type,
            center={"lat": -33.92345542582841, "lon": 18.434424141913478},
            zoom=9.25,
            mapbox_style='stamen-toner',
            opacity=0.75,
            labels={'lon': 'Lon.', 'lat': 'Lat.',
                    'max_val_' + housing_type: '% net income',
                    'max_content_' + housing_type: 'o.w. content',
                    'locations': 'Pixel ID',
                    'flood_type_' + housing_type: 'Flood type',
                    'nb_households_' + housing_type: 'Nb of households',
                    'incgroup_' + housing_type: 'Dominant inc. group',
                    'net_income_' + housing_type: 'Net income',
                    'rent_' + housing_type: 'Annual rent / m²'},
            title='Estimated annual flood damages in ' + housing_type
            + ' housing (as % of net income)',
            color_continuous_scale="Reds",
            template='plotly_white',
            hover_data={'lon': ':.2f', 'lat': ':.2f',
                        'flood_type_' + housing_type: True,
                        'nb_households_' + housing_type: ':,.0f',
                        'incgroup_' + housing_type: ':,.0f',
                        'rent_' + housing_type: ':,.0f',
                        'net_income_' + housing_type: ':,.0f',
                        'max_val_' + housing_type: ':.2%',
                        'max_content_' + housing_type: ':.2%'})
        fig.update_layout(margin={"r": 0, "t": 30, "l": 0, "b": 0})
        fig.update_traces(marker_line_width=0)
        # fig.show()
        fig.write_html(path_maps_rel_damages + housing_type
                       + "/map_damages_" + item + '_' + housing_type + ".html")
        fig.write_image(path_maps_rel_damages + housing_type
                        + "/map_damages_" + item + '_' + housing_type + ".png")

        print("map_damages_" + item + "_" + housing_type + " done")

    fig = px.choropleth_mapbox(
        damage_map_compar_shareinc,
        geojson=geo_grid.geometry,
        locations=geo_grid.index,
        color='max_val_' + housing_type,
        center={"lat": -33.92345542582841, "lon": 18.434424141913478},
        zoom=9.25,
        mapbox_style='stamen-toner',
        opacity=0.75,
        labels={'lon': 'Lon.', 'lat': 'Lat.',
                'max_val_' + housing_type: '% surplus',
                'max_content_' + housing_type: 'o.w. content',
                'locations': 'Pixel ID',
                'flood_type_' + housing_type: 'Flood type',
                'nb_households_' + housing_type: 'Nb of households',
                'incgroup_' + housing_type: 'Dominant inc. group',
                'net_income_' + housing_type: 'Net income',
                'rent_' + housing_type + '_pct': '% change in rent'},
        title='Flood damage surplus from no insurance in ' + housing_type
        + ' housing (as % of net income)',
        color_continuous_scale="Picnic",
        color_continuous_midpoint=0,
        template='plotly_white',
        hover_data={'lon': ':.2f', 'lat': ':.2f',
                    'flood_type_' + housing_type: True,
                    'nb_households_' + housing_type: ':,.0f',
                    'incgroup_' + housing_type: ':,.0f',
                    'rent_' + housing_type + '_pct': ':+.2%',
                    'net_income_' + housing_type: ':,.0f',
                    'max_val_' + housing_type: ':+.2%',
                    'max_content_' + housing_type: ':+.2%'})
    fig.update_layout(margin={"r": 0, "t": 30, "l": 0, "b": 0})
    fig.update_traces(marker_line_width=0)
    # fig.show()
    fig.write_html(path_maps_rel_damages + housing_type
                   + "/map_damages_" + item + '_' + housing_type + ".html")
    fig.write_image(path_maps_rel_damages + housing_type
                    + "/map_damages_" + item + '_' + housing_type + ".png")

    print("map_damages_compar_shareinc_insur_" + housing_type + " done")


# PLOT OTHER EQUILIBRIUM OUTCOMES

# We start with population distribution

for item in list_scenarios:

    equil_maps[item]['lon'] = geo_grid.lon
    equil_maps[item]['lat'] = geo_grid.lat

    fig = px.choropleth_mapbox(
        equil_maps[item],
        geojson=geo_grid.geometry,
        locations=geo_grid.index,
        color='hh_tot',
        center={"lat": -33.92345542582841, "lon": 18.434424141913478},
        zoom=9.25,
        mapbox_style='stamen-toner',
        opacity=0.75,
        labels={'lon': 'Lon.', 'lat': 'Lat.',
                'hh_formal': 'Formal private',
                'formal_incgroup': 'Inc. group (formal)',
                'hh_subsidized': 'Formal subsidized',
                'subsidized_incgroup': 'Inc. group (subsidized)',
                'hh_informal': 'Informal settlements',
                'informal_incgroup': 'Inc. group (informal)',
                'hh_backyard': 'Informal backyards',
                'backyard_incgroup': 'Inc. group (backyard)',
                'hh_tot': 'Total', 'locations': 'Pixel ID'},
        title='Estimated number of households',
        color_continuous_scale="Reds",
        template='plotly_white',
        hover_data={'lon': ':.2f', 'lat': ':.2f',
                    'hh_formal': ':,.0f',
                    'formal_incgroup': ':,.0f',
                    'hh_subsidized': ':,.0f',
                    'subsidized_incgroup': ':,.0f',
                    'hh_informal': ':,.0f',
                    'informal_incgroup': ':,.0f',
                    'hh_backyard': ':,.0f',
                    'backyard_incgroup': ':,.0f',
                    'hh_tot': ':,.0f'})
    fig.update_layout(margin={"r": 0, "t": 30, "l": 0, "b": 0})
    fig.update_traces(marker_line_width=0)
    # fig.show()
    type_map = ""
    if item == list_scenarios[0]:
        type_map = "noinsur"
    elif item == list_scenarios[1]:
        type_map = "insur"
    fig.write_html(path_maps_pop_distrib + "map_pop_distrib_" + type_map
                   + ".html")
    fig.write_image(path_maps_pop_distrib + "map_pop_distrib_" + type_map
                    + ".png")

    print("map_pop_distrib_" + type_map + " done")

fig = px.choropleth_mapbox(
    equil_map_compar,
    geojson=geo_grid.geometry,
    locations=geo_grid.index,
    color='hh_tot',
    center={"lat": -33.92345542582841, "lon": 18.434424141913478},
    zoom=9.25,
    mapbox_style='stamen-toner',
    opacity=0.75,
    labels={'lon': 'Lon.', 'lat': 'Lat.',
            'hh_formal': 'Formal private',
            'formal_incgroup': 'Inc. group (formal)',
            'hh_subsidized': 'Formal subsidized',
            'subsidized_incgroup': 'Inc. group (subsidized)',
            'hh_informal': 'Informal settlements',
            'informal_incgroup': 'Inc. group (informal)',
            'hh_backyard': 'Informal backyards',
            'backyard_incgroup': 'Inc. group (backyard)',
            'hh_tot': 'Total change', 'hh_tot_pct': '% change',
            'locations': 'Pixel ID'},
    title='Evolution of number of households under no insurance',
    color_continuous_scale="Picnic",
    color_continuous_midpoint=0,
    template='plotly_white',
    hover_data={'lon': ':.2f', 'lat': ':.2f',
                'hh_formal': ':,.0f',
                'formal_incgroup': ':,.0f',
                'hh_subsidized': ':,.0f',
                'subsidized_incgroup': ':,.0f',
                'hh_informal': ':,.0f',
                'informal_incgroup': ':,.0f',
                'hh_backyard': ':,.0f',
                'backyard_incgroup': ':,.0f',
                'hh_tot': ':,.0f',
                'hh_tot_pct': ':+.2%'})
fig.update_layout(margin={"r": 0, "t": 30, "l": 0, "b": 0})
fig.update_traces(marker_line_width=0)
# fig.show()
fig.write_html(path_maps_pop_distrib + "map_pop_distrib_compar_insur.html")
fig.write_image(path_maps_pop_distrib + "map_pop_distrib_compar_insur.png")

print("map_pop_distrib_compar_insur done")


# Then we plot rent distribution (housing-type-specific)

for housing_type in housing_types:

    for item in list_scenarios:

        equil_maps[item]['lon'] = geo_grid.lon
        equil_maps[item]['lat'] = geo_grid.lat
        # equil_maps[item]['hsupply_' + housing_type] = (
        #     equil_maps[item]['hsupply_' + housing_type] / 1000000)

        fig = px.choropleth_mapbox(
            equil_maps[item],
            geojson=geo_grid.geometry,
            locations=geo_grid.index,
            color='rent_' + housing_type,
            center={"lat": -33.92345542582841, "lon": 18.434424141913478},
            zoom=9.25,
            mapbox_style='stamen-toner',
            opacity=0.75,
            labels={'locations': 'Pixel ID', 'lon': 'Lon.', 'lat': 'Lat.',
                    'hh_' + housing_type: 'Nb of households',
                    housing_type + '_incgroup': 'Inc. group',
                    # 'hsupply_' + housing_type: 'FAR',
                    'dwelling_size_' + housing_type: 'Avg dwelling size',
                    'flood_deprec_struct_' + housing_type:
                        'Struct. deprec. from floods',
                    'flood_deprec_content_' + housing_type:
                        'Content deprec. from floods',
                    'rent_' + housing_type: 'Annual rent / m²'},
            title='Estimated average annual rent / m² (rands, 2011)',
            color_continuous_scale="Reds",
            template='plotly_white',
            hover_data={'hh_' + housing_type: ':,.0f',
                        housing_type + '_incgroup': ':,.0f',
                        # 'hsupply_' + housing_type: ':.2f',
                        'dwelling_size_' + housing_type: ':,.0f',
                        'flood_deprec_struct_' + housing_type: ':.2%',
                        'flood_deprec_content_' + housing_type: ':.2%',
                        'rent_' + housing_type: ':,.0f'})
        fig.update_layout(margin={"r": 0, "t": 30, "l": 0, "b": 0})
        fig.update_traces(marker_line_width=0)
        # fig.show()
        type_map = ""
        if item == list_scenarios[0]:
            type_map = "noinsur"
        elif item == list_scenarios[1]:
            type_map = "insur"
        fig.write_html(path_maps_rent_distrib + housing_type
                       + "/map_rent_distrib_" + type_map + "_" + housing_type
                       + ".html")
        fig.write_image(path_maps_rent_distrib + housing_type
                        + "/map_rent_distrib_" + type_map + "_" + housing_type
                        + ".png")

        print("map_rent_distrib_" + type_map + "_" + housing_type + " done")

    fig = px.choropleth_mapbox(
        equil_map_compar,
        geojson=geo_grid.geometry,
        locations=geo_grid.index,
        color='rent_' + housing_type,
        center={"lat": -33.92345542582841, "lon": 18.434424141913478},
        zoom=9.25,
        mapbox_style='stamen-toner',
        opacity=0.75,
        labels={'locations': 'Pixel ID', 'lon': 'Lon.', 'lat': 'Lat.',
                'hh_' + housing_type + '_pct': 'Nb of HHs (% change)',
                housing_type + '_incgroup': 'Inc. group',
                # 'hsupply_' + housing_type + '_pct': 'FAR (% change)',
                'dwelling_size_' + housing_type + '_pct':
                    'Dwelling size (% change)',
                'flood_deprec_struct_' + housing_type:
                    'Struct. deprec. from floods',
                'flood_deprec_content_' + housing_type:
                    'Content deprec. from floods',
                'rent_' + housing_type: 'Total change',
                'rent_' + housing_type + '_pct': '% change'},
        title='Evolution of annual rent / m² under no insurance (rands, 2011)',
        color_continuous_scale="Picnic",
        color_continuous_midpoint=0,
        template='plotly_white',
        hover_data={'hh_' + housing_type + '_pct': ':+.2%',
                    housing_type + '_incgroup': ':,.0f',
                    # 'hsupply_' + housing_type + '_pct': ':+.2%',
                    'dwelling_size_' + housing_type + '_pct': ':+.2%',
                    'flood_deprec_struct_' + housing_type: ':.2%',
                    'flood_deprec_content_' + housing_type: ':.2%',
                    'rent_' + housing_type: ':,.0f',
                    'rent_' + housing_type + '_pct': ':+.2%'})
    fig.update_layout(margin={"r": 0, "t": 30, "l": 0, "b": 0})
    fig.update_traces(marker_line_width=0)
    # fig.show()
    fig.write_html(path_maps_rent_distrib + housing_type
                   + "/map_rent_distrib_compar_insur_" + housing_type + ".html"
                   )
    fig.write_image(path_maps_rent_distrib + housing_type
                    + "/map_rent_distrib_compar_insur_" + housing_type + ".png"
                    )

    print("map_rent_distrib_compar_insur_" + housing_type + " done")


# NOW LET US PLOT DAMAGE DISTRIBUTION ACROSS INCOME GROUPS IN ONE DIMENSION

# First, we define the budget share spent on flood damages by income group
list_incgroup = ['1', '2', '3', '4']

# TODO: correction
for item in list_maps_shareinc:
    for incgroup in list_incgroup:
        for flood in flood_types:
            damage_maps_shareinc[item]['sum_' + flood + '_' + incgroup] = 0
            for housing in housing_types:
                temp1 = damage_maps_shareinc[item]['incgroup_' + housing_type]
                temp2 = damage_maps_shareinc[item][
                    'sum_' + flood + '_' + housing_type]
                damage_maps_shareinc[item].loc[
                    temp1 == int(incgroup),
                    'sum_' + flood + '_' + incgroup] = (
                    damage_maps_shareinc[item]['sum_' + flood + '_' + incgroup]
                    + temp2)
        if item == 'noinsur_shareinc':
            damage_maps_shareinc[item]['nb_households_' + incgroup] = (
                dict_scenario_outcomes[list_scenarios[0]]['hh_per_incgroup'][
                    int(incgroup) - 1, :])
        elif item == 'insur_shareinc':
            damage_maps_shareinc[item]['nb_households_' + incgroup] = (
                dict_scenario_outcomes[list_scenarios[1]]['hh_per_incgroup'][
                    int(incgroup) - 1, :])

# TODO: Repeat instead of doing subplots?
# NB: add pattern_shape for content?
# Write somewhere where the rest of households is!

damage_maps_shareinc['noinsur_shareinc']['insur'] = 'No'
damage_maps_shareinc['insur_shareinc']['insur'] = 'Yes'
damage_maps_shareinc_concat = pd.concat(
    [damage_maps_shareinc['insur_shareinc'],
     damage_maps_shareinc['noinsur_shareinc']])

for flood in flood_types:
    for incgroup in list_incgroup:
        dist = px.histogram(
            damage_maps_shareinc_concat.loc[
                damage_maps_shareinc_concat[
                    'sum_' + flood + '_' + incgroup] > 0],
            x='sum_' + flood + '_' + incgroup,
            y='nb_households_' + incgroup,
            color='insur',
            labels={'sum_' + flood + '_' + incgroup: '% of net income',
                    'nb_households_' + incgroup: 'nb of households',
                    'insur': 'w/ insurance'},
            barmode='group',
            hover_data={'sum_' + flood + '_' + incgroup: False},
            title='Distribution of flood damages among income group '
            + incgroup + ' (as % of net income)',
            template='plotly_white')

        dist.update_layout(
            yaxis_title='Total nb of households',
            xaxis=dict(tickmode='linear', tick0=0, dtick=0.01)
            )

        dist.write_html(path_damage_distrib + flood + "_damage_dist_"
                        + incgroup + "_insur.html")
        dist.write_image(path_damage_distrib + flood + "_damage_dist_"
                        + incgroup + "_insur.png")

        print(flood + "_damage_dist_" + incgroup + "_insur done")

# TODO: show utility changes in income equiv. What about land value creation?
# Read Paolo's paper
# Also save data
