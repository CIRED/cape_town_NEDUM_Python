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
import inputs.data as inpdt


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

path_use_case = path_outputs + 'use_case_anticip/'
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

path_damage_distrib_total = path_damage_distrib + 'total/'
path_damage_distrib_content = path_damage_distrib + 'content/'

list_output_paths = [
    path_use_case, path_maps, path_maps_abs_damages, path_maps_rel_damages,
    path_maps_pop_distrib, path_maps_rent_distrib,
    path_maps_rel_damages_formal, path_maps_rel_damages_subsidized,
    path_maps_rel_damages_informal, path_maps_rel_damages_backyard,
    path_maps_rent_distrib_formal, path_maps_rent_distrib_subsidized,
    path_maps_rent_distrib_informal, path_maps_rent_distrib_backyard,
    path_charts, path_damage_distrib,
    path_damage_distrib_total, path_damage_distrib_content]

for path in list_output_paths:
    try:
        os.mkdir(path)
    except OSError as error:
        print(error)


# IMPORT PARAMETER AND OPTIONS

options = inpprm.import_options()
param = inpprm.import_param(
    path_precalc_inp, options)


# LAND USE PROJECTIONS + AMENITIES (for input maps)

(interest_rate, population, housing_type_data, total_RDP
 ) = inpdt.import_macro_data(param, path_scenarios, path_folder)
grid, center = inpdt.import_grid(path_data)
(data_rdp, housing_types_sp, data_sp, mitchells_plain_grid_baseline,
 grid_formal_density_HFA, threshold_income_distribution, income_distribution,
 cape_town_limits) = inpdt.import_households_data(path_precalc_inp)
housing_types_data = pd.read_excel(path_folder + 'housing_types_grid_sal.xlsx')
housing_types_data[np.isnan(housing_types_data)] = 0

(spline_RDP, spline_estimate_RDP, spline_land_RDP,
 spline_land_backyard, spline_land_informal, spline_land_constraints,
 number_properties_RDP) = (
     inpdt.import_land_use(grid, options, param, data_rdp, housing_types_data,
                           housing_type_data, path_data, path_folder)
)

coeff_land = inpdt.import_coeff_land(
    spline_land_constraints, spline_land_backyard, spline_land_informal,
    spline_land_RDP, param, 0)

amenities = inpdt.import_amenities(path_precalc_inp, options)


# DEFINE SCENARIOS TO IMPORT

list_scenarios = ['floods00_F0_P11_C10_scenario232',
                  'floods10_F0_P11_C10_scenario232']

# We store other dimension names
housing_types = ['formal', 'subsidized', 'informal', 'backyard']
flood_types = ['fluvialu', 'pluvial', 'coastal']
damage_types = ['structure', 'content']


# %% Imports

print("Import data")

# DAMAGE DATA

# We retrieve all available dimension combinations in the data
list_dim = [flood_types, housing_types, damage_types]
all_dim = list(itertools.product(*list_dim))

# We store them in a dictionary

dict_damage_map_noanticip = {}
for dim in all_dim:
    table = pd.read_csv(
        path_outputs + list_scenarios[0] + '/tables/floods/'
        + dim[0] + '_' + dim[1] + '_' + dim[2] + '_2d_sim.csv',
        names=['damage'], header=None)
    dict_damage_map_noanticip[dim[0] + '_' + dim[1] + '_' + dim[2]] = table

dict_damage_map_anticip = {}
for dim in all_dim:
    table = pd.read_csv(
        path_outputs + list_scenarios[1] + '/tables/floods/'
        + dim[0] + '_' + dim[1] + '_' + dim[2] + '_2d_sim.csv',
        names=['damage'], header=None)
    dict_damage_map_anticip[dim[0] + '_' + dim[1] + '_' + dim[2]] = table

dict_damage_map_noanticip_shareinc = {}
for dim in all_dim:
    table = pd.read_csv(
        path_outputs + list_scenarios[0] + '/tables/floods/'
        + dim[0] + '_' + dim[1] + '_' + dim[2] + '_2d_sim_shareinc.csv',
        names=['damage'], header=None)
    dict_damage_map_noanticip_shareinc[dim[0] +
                                     '_' + dim[1] + '_' + dim[2]] = table

dict_damage_map_anticip_shareinc = {}
for dim in all_dim:
    table = pd.read_csv(
        path_outputs + list_scenarios[1] + '/tables/floods/'
        + dim[0] + '_' + dim[1] + '_' + dim[2] + '_2d_sim_shareinc.csv',
        names=['damage'], header=None)
    dict_damage_map_anticip_shareinc[dim[0] +
                                   '_' + dim[1] + '_' + dim[2]] = table

dict_damage_map = {}
dict_damage_map["noanticip"] = dict_damage_map_noanticip
dict_damage_map["anticip"] = dict_damage_map_anticip
dict_damage_map["noanticip_shareinc"] = dict_damage_map_noanticip_shareinc
dict_damage_map["anticip_shareinc"] = dict_damage_map_anticip_shareinc


# OTHER EQUILIBRIUM OUTCOMES

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
for scenario in list_scenarios:
    for dim in all_dim:
        try:
            table = pd.read_csv(
                path_outputs + scenario + '/tables/floods/'
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

# #For spatial distribution

# We build an aggregate damage table to show on the map, and will display all
# the additional information as hover data

list_maps = ['noanticip', 'anticip']
damage_maps = {}
list_maps_shareinc = ['noanticip_shareinc', 'anticip_shareinc']
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

    # For reference to rainwater floods
    concat_damage_map['max_val_rainwater'] = concat_damage_map[
        ['sum_fluvialu', 'sum_pluvial']].max(axis=1)

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
    if item == 'noanticip_shareinc':
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
        concat_damage_map.loc[
            concat_damage_map['nb_households_formal'] == 0, 'rent_formal'] = 0
        concat_damage_map['rent_subsidized'] = (
            dict_scenario_outcomes[list_scenarios[0]]['rent'][3, :])
        concat_damage_map.loc[
            concat_damage_map['nb_households_subsidized'] == 0,
            'rent_subsidized'] = 0
        concat_damage_map['rent_informal'] = (
            dict_scenario_outcomes[list_scenarios[0]]['rent'][2, :])
        concat_damage_map.loc[
            concat_damage_map['nb_households_informal'] == 0,
            'rent_informal'] = 0
        concat_damage_map['rent_backyard'] = (
            dict_scenario_outcomes[list_scenarios[0]]['rent'][1, :])
        concat_damage_map.loc[
            concat_damage_map['nb_households_subsidized'] == 0,
            'rent_subsidized'] = 0
        for housing_type in housing_types:
            concat_damage_map['incgroup_' + housing_type] = (
                dict_scenario_outcomes[list_scenarios[0]
                                       ][housing_type + '_incgroup'])
            concat_damage_map['net_income_' + housing_type] = (
                dict_scenario_outcomes[list_scenarios[0]
                                       ]['net_income_' + housing_type])

    elif item == 'anticip_shareinc':
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
        concat_damage_map.loc[
            concat_damage_map['nb_households_formal'] == 0, 'rent_formal'] = 0
        concat_damage_map['rent_subsidized'] = (
            dict_scenario_outcomes[list_scenarios[1]]['rent'][3, :])
        concat_damage_map.loc[
            concat_damage_map['nb_households_subsidized'] == 0,
            'rent_subsidized'] = 0
        concat_damage_map['rent_informal'] = (
            dict_scenario_outcomes[list_scenarios[1]]['rent'][2, :])
        concat_damage_map.loc[
            concat_damage_map['nb_households_informal'] == 0,
            'rent_informal'] = 0
        concat_damage_map['rent_backyard'] = (
            dict_scenario_outcomes[list_scenarios[1]]['rent'][1, :])
        concat_damage_map.loc[
            concat_damage_map['nb_households_backyard'] == 0,
            'rent_backyard'] = 0
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
    # (here, value does not change across scenarios as there is no climate
    # change)
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
    damage_maps['noanticip'][list_num_var] - damage_maps['anticip'][list_num_var])
damage_map_compar[list_num_pct] = (
    damage_maps['noanticip'][list_num_var] / damage_maps['anticip'][list_num_var]
    - 1)
damage_map_compar['lon'] = geo_grid.lon
damage_map_compar['lat'] = geo_grid.lat
damage_map_compar['flood_type'] = damage_maps['noanticip']['flood_type']
damage_map_compar.loc[
    damage_map_compar['flood_type'] == 'None', 'flood_type'
] = damage_maps['anticip']['flood_type']

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
    damage_maps_shareinc['noanticip_shareinc'][list_num_var]
    - damage_maps_shareinc['anticip_shareinc'][list_num_var])
damage_map_compar_shareinc[list_num_pct] = (
    damage_maps_shareinc['noanticip_shareinc'][list_num_var] /
    damage_maps_shareinc['anticip_shareinc'][list_num_var]
    - 1)
damage_map_compar_shareinc['lon'] = geo_grid.lon
damage_map_compar_shareinc['lat'] = geo_grid.lat

for housing_type in housing_types:
    damage_map_compar_shareinc['flood_type_' + housing_type] = (
        damage_maps_shareinc['noanticip_shareinc']['flood_type_' + housing_type])
    damage_map_compar_shareinc.loc[
        damage_map_compar_shareinc['flood_type_' + housing_type] == 'None',
        'flood_type_' + housing_type
    ] = damage_maps_shareinc['anticip_shareinc']['flood_type_' + housing_type]
    # NB: we take the no-anticipation case as a benchmark, since we want to show
    # increase damages from no anticipation. We'll rely on other maps to show
    # population moves and composition effects, and the extent to which they
    # explain what we observe
    damage_map_compar_shareinc['nb_households_' + housing_type] = (
        damage_maps_shareinc['noanticip_shareinc'][
            'nb_households_' + housing_type])
    damage_map_compar_shareinc['incgroup_' + housing_type] = (
        damage_maps_shareinc['noanticip_shareinc']['incgroup_' + housing_type])
    damage_map_compar_shareinc['net_income_' + housing_type] = (
        damage_maps_shareinc['noanticip_shareinc']['net_income_' + housing_type])

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

# NB: Again, we take the no-anticipation case as a benchmark, while considering
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

# PLOT SPATIAL DAMAGE DISTRIBUTION

# We plot all the relevant information for a given scenario as a chloropleth
# map

# #First for damages in absolute values

for item in list_maps:

    damage_maps[item]['lon'] = geo_grid.lon
    damage_maps[item]['lat'] = geo_grid.lat
    damage_maps[item].loc[
        damage_maps[item]['max_val'] == 0, 'max_val'] = np.nan

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
    fig.write_image(path_maps_abs_damages + "map_damages_" + item + ".png",
                    height=650, width=1000)

    print("map_damages_" + item + " done")

damage_map_compar.loc[
    damage_map_compar['max_val'] == 0, 'max_val'] = np.nan
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
    title='Flood damage increase from no anticipation (in rands, 2011)',
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
fig.write_html(path_maps_abs_damages + "map_damages_compar_anticip.html")
fig.write_image(path_maps_abs_damages + "map_damages_compar_anticip.png",
                height=650, width=1000)

print("map_damages_compar_anticip done")


# #Then for damages in relative terms

for housing_type in housing_types:

    for item in list_maps_shareinc:

        damage_maps_shareinc[item]['lon'] = geo_grid.lon
        damage_maps_shareinc[item]['lat'] = geo_grid.lat
        damage_maps_shareinc[item].loc[
            damage_maps_shareinc[item]['max_val_' + housing_type] == 0,
            'max_val_' + housing_type] = np.nan

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
            + ' housing (as share of net income)',
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
                        + "/map_damages_" + item + '_' + housing_type + ".png",
                        height=650, width=1000
                        )

        print("map_damages_" + item + "_" + housing_type + " done")

    damage_map_compar_shareinc.loc[
        damage_map_compar_shareinc['max_val_' + housing_type] == 0,
        'max_val_' + housing_type] = np.nan
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
                'max_val_' + housing_type: 'Total change',
                'max_content_' + housing_type: 'o.w. content',
                'locations': 'Pixel ID',
                'flood_type_' + housing_type: 'Flood type',
                'nb_households_' + housing_type: 'Nb of households',
                'incgroup_' + housing_type: 'Dominant inc. group',
                'net_income_' + housing_type: 'Net income',
                'rent_' + housing_type + '_pct': '% change in rent'},
        title='Flood damage increase from no anticipation in ' + housing_type
        + ' housing (as share of net income)',
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
                   + "/map_damages_compar_" + housing_type + ".html")
    fig.write_image(path_maps_rel_damages + housing_type
                    + "/map_damages_compar_" + housing_type + ".png",
                    height=650, width=1000)

    print("map_damages_compar_shareinc_anticip_" + housing_type + " done")


# PLOT OTHER EQUILIBRIUM OUTCOMES

# We start with population distribution

for item in list_scenarios:

    equil_maps[item]['lon'] = geo_grid.lon
    equil_maps[item]['lat'] = geo_grid.lat
    equil_maps[item].loc[
        equil_maps[item]['hh_tot'] == 0, 'hh_tot'] = np.nan

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
        type_map = "noanticip"
    elif item == list_scenarios[1]:
        type_map = "anticip"
    fig.write_html(path_maps_pop_distrib + "map_pop_distrib_" + type_map
                   + ".html")
    fig.write_image(path_maps_pop_distrib + "map_pop_distrib_" + type_map
                    + ".png", height=650, width=1000)

    print("map_pop_distrib_" + type_map + " done")

equil_map_compar.loc[
    equil_map_compar['hh_tot'] == 0, 'hh_tot'] = np.nan
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
    title='Evolution of number of households under no anticipation',
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
fig.write_html(path_maps_pop_distrib + "map_pop_distrib_compar_anticip.html")
fig.write_image(path_maps_pop_distrib + "map_pop_distrib_compar_anticip.png",
                height=650, width=1000)

print("map_pop_distrib_compar_anticip done")


# Then we plot rent distribution (housing-type-specific)

for housing_type in housing_types:

    for item in list_scenarios:

        equil_maps[item]['lon'] = geo_grid.lon
        equil_maps[item]['lat'] = geo_grid.lat
        # equil_maps[item]['hsupply_' + housing_type] = (
        #     equil_maps[item]['hsupply_' + housing_type] / 1000000)
        equil_maps[item].loc[
            equil_maps[item]['rent_' + housing_type] == 0,
            'rent_' + housing_type] = np.nan

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
            title='Estimated average annual rent / m² in ' + housing_type
            + ' housing (in rands, 2011)',
            color_continuous_scale="Reds",
            template='plotly_white',
            hover_data={'lon': ':.2f', 'lat': ':.2f',
                        'hh_' + housing_type: ':,.0f',
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
            type_map = "noanticip"
        elif item == list_scenarios[1]:
            type_map = "anticip"
        fig.write_html(path_maps_rent_distrib + housing_type
                       + "/map_rent_distrib_" + type_map + "_" + housing_type
                       + ".html")
        fig.write_image(path_maps_rent_distrib + housing_type
                        + "/map_rent_distrib_" + type_map + "_" + housing_type
                        + ".png", height=650, width=1000)

        print("map_rent_distrib_" + type_map + "_" + housing_type + " done")

    equil_map_compar.loc[
        equil_map_compar['rent_' + housing_type] == 0,
        'rent_' + housing_type] = np.nan
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
        title='Evolution of annual rent / m² under no anticipation in '
        + housing_type + ' housing (in rands, 2011)',
        color_continuous_scale="Picnic",
        color_continuous_midpoint=0,
        template='plotly_white',
        hover_data={'lon': ':.2f', 'lat': ':.2f',
                    'hh_' + housing_type + '_pct': ':+.2%',
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
                   + "/map_rent_distrib_compar_anticip_" + housing_type + ".html"
                   )
    fig.write_image(path_maps_rent_distrib + housing_type
                    + "/map_rent_distrib_compar_anticip_" + housing_type
                    + ".png", height=650, width=1000
                    )

    print("map_rent_distrib_compar_anticip_" + housing_type + " done")


# NOW LET US PLOT DAMAGE DISTRIBUTION ACROSS INCOME GROUPS IN ONE DIMENSION

# First, we define the budget share spent on flood damages by income group
list_incgroup = ['1', '2', '3', '4']

# Then, we create a variable associating each observed damages (as a share of
# income) to the dominant income group in the location, based on observed
# housing type
for item in list_maps_shareinc:
    for incgroup in list_incgroup:
        if item == 'noanticip_shareinc':
            damage_maps_shareinc[item]['nb_households_' + incgroup] = (
                dict_scenario_outcomes[list_scenarios[0]]['hh_per_incgroup'][
                    int(incgroup) - 1, :])
        elif item == 'anticip_shareinc':
            damage_maps_shareinc[item]['nb_households_' + incgroup] = (
                dict_scenario_outcomes[list_scenarios[1]]['hh_per_incgroup'][
                    int(incgroup) - 1, :])
        for flood in flood_types:
            damage_maps_shareinc[item]['sum_' + flood + '_' + incgroup] = 0
            for housing in housing_types:
                temp1 = damage_maps_shareinc[item]['incgroup_' + housing]
                temp2 = damage_maps_shareinc[item][
                    'sum_' + flood + '_' + housing]
                damage_maps_shareinc[item].loc[
                    temp1 == int(incgroup),
                    'sum_' + flood + '_' + incgroup] = (
                    damage_maps_shareinc[item]['sum_' + flood + '_' + incgroup]
                    + temp2)


# We repeat the observation for content damages only (to display as a share of
# total damages)
for item in list_maps_shareinc:
    for incgroup in list_incgroup:
        for flood in flood_types:
            damage_maps_shareinc[item]['content_' + flood + '_' + incgroup] = 0
            for housing in housing_types:
                temp1 = damage_maps_shareinc[item]['incgroup_' + housing]
                temp2 = damage_maps_shareinc[item][
                    flood + '_' + housing + '_content']
                damage_maps_shareinc[item].loc[
                    temp1 == int(incgroup),
                    'content_' + flood + '_' + incgroup] = (
                    damage_maps_shareinc[item][
                        'content_' + flood + '_' + incgroup]
                    + temp2)

# NB: Write somewhere where the rest of households is (on top of aggregate
# satistics)!

damage_maps_shareinc['noanticip_shareinc']['anticip'] = 'No'
damage_maps_shareinc['anticip_shareinc']['anticip'] = 'Yes'
damage_maps_shareinc_concat = pd.concat(
    [damage_maps_shareinc['anticip_shareinc'],
     damage_maps_shareinc['noanticip_shareinc']])

damage_maps['noanticip']['anticip'] = 'No'
damage_maps['anticip']['anticip'] = 'Yes'
damage_maps_concat = pd.concat(
    [damage_maps['anticip'],
     damage_maps['noanticip']])

equil_maps[list_scenarios[0]]['anticip'] = 'No'
equil_maps[list_scenarios[1]]['anticip'] = 'Yes'
equil_maps_concat = pd.concat(
    [equil_maps[list_scenarios[1]],
     equil_maps[list_scenarios[0]]])


# Then we plot the two variables

agg_stat_rel_damages = pd.DataFrame()

for flood in flood_types:
    for incgroup in list_incgroup:
        dist_sum = px.histogram(
            damage_maps_shareinc_concat.loc[
                damage_maps_shareinc_concat[
                    'sum_' + flood + '_' + incgroup] > 0],
            x='sum_' + flood + '_' + incgroup,
            y='nb_households_' + incgroup,
            color='anticip',
            labels={'sum_' + flood + '_' + incgroup: 'Share of net income',
                    'nb_households_' + incgroup: 'nb of households',
                    'anticip': 'w/ anticipation'},
            barmode='group',
            hover_data={'sum_' + flood + '_' + incgroup: False},
            title='Distribution of damages in flood zones among income group '
            + incgroup + ' (as share of net income)',
            template='plotly_white')

        try:
            noanticip_avg = np.average(
                damage_maps_shareinc['noanticip_shareinc'].loc[
                    damage_maps_shareinc['noanticip_shareinc'
                                         ]['sum_' + flood + '_' + incgroup]
                    > 0, 'sum_' + flood + '_' + incgroup],
                weights=damage_maps_shareinc['noanticip_shareinc'].loc[
                    damage_maps_shareinc['noanticip_shareinc'
                                         ]['sum_' + flood + '_' + incgroup]
                    > 0, 'nb_households_' + incgroup]
            )
        except ZeroDivisionError:
            noanticip_avg = 0

        noanticip_pop = np.nansum(
            damage_maps_shareinc['noanticip_shareinc'].loc[
                damage_maps_shareinc['noanticip_shareinc'
                                     ]['sum_' + flood + '_' + incgroup]
                > 0, 'nb_households_' + incgroup]
        )

        try:
            anticip_avg = np.average(
                damage_maps_shareinc['anticip_shareinc'].loc[
                    damage_maps_shareinc['anticip_shareinc'
                                         ]['sum_' + flood + '_' + incgroup]
                    > 0, 'sum_' + flood + '_' + incgroup],
                weights=damage_maps_shareinc['anticip_shareinc'].loc[
                    damage_maps_shareinc['anticip_shareinc'
                                         ]['sum_' + flood + '_' + incgroup]
                    > 0, 'nb_households_' + incgroup]
            )
        except ZeroDivisionError:
            anticip_avg = 0

        anticip_pop = np.nansum(
            damage_maps_shareinc['anticip_shareinc'].loc[
                damage_maps_shareinc['anticip_shareinc'
                                     ]['sum_' + flood + '_' + incgroup]
                > 0, 'nb_households_' + incgroup]
        )

        agg_stat_rel_damages.at[0, 'anticip_avg_' + flood + '_' + incgroup
                                ] = anticip_avg
        agg_stat_rel_damages.at[0, 'noanticip_avg_' + flood + '_' + incgroup
                                ] = noanticip_avg
        agg_stat_rel_damages.at[0, 'anticip_pop_' + flood + '_' + incgroup
                                ] = anticip_pop
        agg_stat_rel_damages.at[0, 'noanticip_pop_' + flood + '_' + incgroup
                                ] = noanticip_pop

        dist_sum.update_layout(
            # xaxis=dict(tickmode='linear', tick0=0, dtick=0.01),
            yaxis_title='Total nb of households'
        )
        dist_sum.write_html(path_damage_distrib_total + flood + "_damage_dist_"
                            + incgroup + "_anticip_sum.html")
        dist_sum.write_image(path_damage_distrib_total + flood
                             + "_damage_dist_" + incgroup + "_anticip_sum.png",
                             height=650, width=1000)
        print(flood + "_damage_dist_" + incgroup + "_anticip_sum done")

        dist_content = px.histogram(
            damage_maps_shareinc_concat.loc[
                damage_maps_shareinc_concat[
                    'content_' + flood + '_' + incgroup] > 0],
            x='content_' + flood + '_' + incgroup,
            y='nb_households_' + incgroup,
            color='anticip',
            labels={'content_' + flood + '_' + incgroup: '(o.w. content)',
                    'nb_households_' + incgroup: 'nb of households',
                    'anticip': 'w/ anticipation'},
            barmode='group',
            hover_data={'content_' + flood + '_' + incgroup: False},
            title='Distribution of damages in flood zones among income group '
            + incgroup + ' (as share of net income)',
            template='plotly_white')

        dist_content.update_layout(
            # xaxis=dict(tickmode='linear', tick0=0, dtick=0.01),
            yaxis_title='Total nb of households'
        )
        dist_content.write_html(
            path_damage_distrib_content + flood + "_damage_dist_" + incgroup
            + "_anticip_content.html")
        dist_content.write_image(
            path_damage_distrib_content + flood + "_damage_dist_" + incgroup
            + "_anticip_content.png", height=650, width=1000)

        print(flood + "_damage_dist_" + incgroup + "_anticip_content done")


# Since we do not allow for substitution effects, income losses from flood
# damages in the no-anticipation case directly translate into a utility loss
# through an income effect (note that spatial indifference does not hold in
# this case). Realized utility levels cannot be directly compared however since
# they only have an ordinal (as opposed to cardinal) meaning.

# NB: for households living in formal private housing, we can consider that
# they pay all flood damages, even structural ones through increases in rents:
# this is an application of tax equivalence with a perfectly elastic supply

# Land value creation potential is another brick that we will add in another
# use case on mitigation investments


# AGGREGATE DAMAGES

# We first create adequate tables for bar charts

yes_no = ['Yes', 'No']

dict_agg_damages = {}
for indic in yes_no:
    for flood in flood_types:
        for housing in housing_types:
            for damage in damage_types:
                dict_agg_damages[
                    indic + '_' + flood + '_' + housing + '_' + damage
                ] = damage_maps_concat.loc[
                    damage_maps_concat['anticip'] == indic,
                    flood + '_' + housing + '_' + damage].sum()

for flood in flood_types:
    for housing in housing_types:
        for damage in damage_types:
            dict_agg_damages[
                'pct_change_' + flood + '_' + housing + '_' + damage
            ] = (dict_agg_damages[
                'No_' + flood + '_' + housing + '_' + damage]
                / dict_agg_damages[
                'Yes_' + flood + '_' + housing + '_' + damage]
                - 1)

dict_df_agg_damages = {}

for flood in flood_types:

    df_agg_damages_anticip = pd.DataFrame()
    df_agg_damages_anticip.at[0, 'Structures'] = (
        dict_agg_damages['Yes_' + flood + '_formal_structure'])
    df_agg_damages_anticip.at[0, 'Contents'] = (
        dict_agg_damages['Yes_' + flood + '_formal_content'])
    df_agg_damages_anticip.at[1, 'Structures'] = (
        dict_agg_damages['Yes_' + flood + '_subsidized_structure'])
    df_agg_damages_anticip.at[1, 'Contents'] = (
        dict_agg_damages['Yes_' + flood + '_subsidized_content'])
    df_agg_damages_anticip.at[2, 'Structures'] = (
        dict_agg_damages['Yes_' + flood + '_informal_structure'])
    df_agg_damages_anticip.at[2, 'Contents'] = (
        dict_agg_damages['Yes_' + flood + '_informal_content'])
    df_agg_damages_anticip.at[3, 'Structures'] = (
        dict_agg_damages['Yes_' + flood + '_backyard_structure'])
    df_agg_damages_anticip.at[3, 'Contents'] = (
        dict_agg_damages['Yes_' + flood + '_backyard_content'])

    df_agg_damages_anticip.rename(index={
        0: 'Formal private', 1: 'Formal subsidized', 2: 'Informal settlements',
        3: 'Informal backyards'}, inplace=True)
    df_agg_damages_anticip = df_agg_damages_anticip / 1000000

    df_agg_damages_noanticip = pd.DataFrame()
    df_agg_damages_noanticip.at[0, 'Structures'] = (
        dict_agg_damages['No_' + flood + '_formal_structure'])
    df_agg_damages_noanticip.at[0, 'Contents'] = (
        dict_agg_damages['No_' + flood + '_formal_content'])
    df_agg_damages_noanticip.at[1, 'Structures'] = (
        dict_agg_damages['No_' + flood + '_subsidized_structure'])
    df_agg_damages_noanticip.at[1, 'Contents'] = (
        dict_agg_damages['No_' + flood + '_subsidized_content'])
    df_agg_damages_noanticip.at[2, 'Structures'] = (
        dict_agg_damages['No_' + flood + '_informal_structure'])
    df_agg_damages_noanticip.at[2, 'Contents'] = (
        dict_agg_damages['No_' + flood + '_informal_content'])
    df_agg_damages_noanticip.at[3, 'Structures'] = (
        dict_agg_damages['No_' + flood + '_backyard_structure'])
    df_agg_damages_noanticip.at[3, 'Contents'] = (
        dict_agg_damages['No_' + flood + '_backyard_content'])

    df_agg_damages_noanticip.rename(index={
        0: 'Formal private', 1: 'Formal subsidized', 2: 'Informal settlements',
        3: 'Informal backyards'}, inplace=True)
    df_agg_damages_noanticip = df_agg_damages_noanticip / 1000000

    dict_df_agg_damages[flood + '_anticip'] = df_agg_damages_anticip
    dict_df_agg_damages[flood + '_noanticip'] = df_agg_damages_noanticip


# We then do the plots

newnames_noanticip = {'Structures': 'Structures (w/o/ anticipation)',
                    'Contents': 'Contents (w/o/ anticipation)'}
newnames_anticip = {'Structures': 'Structures (w/ anticipation)',
                  'Contents': 'Contents (w/ anticipation)'}

# First for fluvial floods

fig_fluvialu_noanticip = px.bar(
    dict_df_agg_damages['fluvialu_noanticip'],
    x=dict_df_agg_damages['fluvialu_noanticip'].index,
    y=dict_df_agg_damages['fluvialu_noanticip'].columns,
    barmode='group',
    color_discrete_sequence=px.colors.qualitative.Pastel1,
    title='Estimated annual damages from fluvial floods (in M rands, 2011)',
    labels={'index': 'Housing type', 'value': 'Annual damages',
            'variable': 'Damages'},
    opacity=0.8,
    template='plotly_white')

fig_fluvialu_noanticip.for_each_trace(
    lambda t: t.update(
        name=newnames_noanticip[t.name],
        legendgroup=newnames_noanticip[t.name],
        hovertemplate=t.hovertemplate.replace(
         t.name, newnames_noanticip[t.name])
    ))

fig_fluvialu_anticip = px.bar(
    dict_df_agg_damages['fluvialu_anticip'],
    x=dict_df_agg_damages['fluvialu_anticip'].index,
    y=dict_df_agg_damages['fluvialu_anticip'].columns,
    barmode='group',
    color_discrete_sequence=px.colors.qualitative.Set1,
    title='Estimated annual damages from fluvial floods (in M rands, 2011)',
    labels={'index': 'Housing type', 'value': 'Annual damages',
            'variable': 'Damages'},
    opacity=0.8,
    template='plotly_white')

fig_fluvialu_anticip.for_each_trace(
    lambda t: t.update(
        name=newnames_anticip[t.name],
        legendgroup=newnames_anticip[t.name],
        hovertemplate=t.hovertemplate.replace(
         t.name, newnames_anticip[t.name])
    ))

fig_fluvialu = go.Figure(fig_fluvialu_noanticip)
fig_fluvialu = fig_fluvialu.add_traces(fig_fluvialu_anticip.data)
# fig_fluvialu.show()
fig_fluvialu.write_html(path_charts + "fluvialu_sim_damage_sum_anticip.html")
fig_fluvialu.write_image(path_charts + "fluvialu_sim_damage_sum_anticip.png",
                         height=650, width=1000)

# Then for pluvial floods

fig_pluvial_noanticip = px.bar(
    dict_df_agg_damages['pluvial_noanticip'],
    x=dict_df_agg_damages['pluvial_noanticip'].index,
    y=dict_df_agg_damages['pluvial_noanticip'].columns,
    barmode='group',
    color_discrete_sequence=px.colors.qualitative.Pastel1,
    title='Estimated annual damages from pluvial floods (in M rands, 2011)',
    labels={'index': 'Housing type', 'value': 'Annual damages',
            'variable': 'Damages'},
    opacity=0.8,
    template='plotly_white')

fig_pluvial_noanticip.for_each_trace(
    lambda t: t.update(
        name=newnames_noanticip[t.name],
        legendgroup=newnames_noanticip[t.name],
        hovertemplate=t.hovertemplate.replace(
         t.name, newnames_noanticip[t.name])
    ))

fig_pluvial_anticip = px.bar(
    dict_df_agg_damages['pluvial_anticip'],
    x=dict_df_agg_damages['pluvial_anticip'].index,
    y=dict_df_agg_damages['pluvial_anticip'].columns,
    barmode='group',
    color_discrete_sequence=px.colors.qualitative.Set1,
    title='Estimated annual damages from pluvial floods (in M rands, 2011)',
    labels={'index': 'Housing type', 'value': 'Annual damages',
            'variable': 'Damages'},
    opacity=0.8,
    template='plotly_white')

fig_pluvial_anticip.for_each_trace(
    lambda t: t.update(
        name=newnames_anticip[t.name],
        legendgroup=newnames_anticip[t.name],
        hovertemplate=t.hovertemplate.replace(
         t.name, newnames_anticip[t.name])
    ))

fig_pluvial = go.Figure(fig_pluvial_noanticip)
fig_pluvial = fig_pluvial.add_traces(fig_pluvial_anticip.data)
# fig_pluvial.show()
fig_pluvial.write_html(path_charts + "pluvial_sim_damage_sum_anticip.html")
fig_pluvial.write_image(path_charts + "pluvial_sim_damage_sum_anticip.png",
                        height=650, width=1000)

# Finally for coastal floods

fig_coastal_noanticip = px.bar(
    dict_df_agg_damages['coastal_noanticip'],
    x=dict_df_agg_damages['coastal_noanticip'].index,
    y=dict_df_agg_damages['coastal_noanticip'].columns,
    barmode='group',
    color_discrete_sequence=px.colors.qualitative.Pastel1,
    title='Estimated annual damages from coastal floods (in M rands, 2011)',
    labels={'index': 'Housing type', 'value': 'Annual damages',
            'variable': 'Damages'},
    opacity=0.8,
    template='plotly_white')

fig_coastal_noanticip.for_each_trace(
    lambda t: t.update(
        name=newnames_noanticip[t.name],
        legendgroup=newnames_noanticip[t.name],
        hovertemplate=t.hovertemplate.replace(
         t.name, newnames_noanticip[t.name])
    ))

fig_coastal_anticip = px.bar(
    dict_df_agg_damages['coastal_anticip'],
    x=dict_df_agg_damages['coastal_anticip'].index,
    y=dict_df_agg_damages['coastal_anticip'].columns,
    barmode='group',
    color_discrete_sequence=px.colors.qualitative.Set1,
    title='Estimated annual damages from coastal floods (in M rands, 2011)',
    labels={'index': 'Housing type', 'value': 'Annual damages',
            'variable': 'Damages'},
    opacity=0.8,
    template='plotly_white')

fig_coastal_anticip.for_each_trace(
    lambda t: t.update(
        name=newnames_anticip[t.name],
        legendgroup=newnames_anticip[t.name],
        hovertemplate=t.hovertemplate.replace(
         t.name, newnames_anticip[t.name])
    ))

fig_coastal = go.Figure(fig_coastal_noanticip)
fig_coastal = fig_coastal.add_traces(fig_coastal_anticip.data)
# fig_coastal.show()
fig_coastal.write_html(path_charts + "coastal_sim_damage_sum_anticip.html")
fig_coastal.write_image(path_charts + "coastal_sim_damage_sum_anticip.png",
                        height=650, width=1000)


# INPUT MAPS

# We start with flood depreciation

equil_maps[list_scenarios[0]].loc[
    equil_maps[list_scenarios[0]]['flood_deprec_content_formal'] == 0,
    'flood_deprec_content_formal'] = np.nan
fig = px.choropleth_mapbox(
    equil_maps[list_scenarios[0]],
    geojson=geo_grid.geometry,
    locations=geo_grid.index,
    color='flood_deprec_content_formal',
    center={"lat": -33.92345542582841, "lon": 18.434424141913478},
    zoom=9.25,
    mapbox_style='stamen-terrain',
    opacity=0.5,
    labels={'locations': 'Pixel ID', 'lon': 'Lon.', 'lat': 'Lat.',
            'flood_deprec_content_formal': 'Deprec. rate'},
    title='Estimated fraction of capital destroyed due to floods'
    + ' (content damages)',
    color_continuous_scale="Reds",
    template='plotly_white',
    hover_data={'lon': ':.2f', 'lat': ':.2f',
                'flood_deprec_content_formal': ':.2f'})
fig.update_layout(margin={"r": 0, "t": 30, "l": 0, "b": 0})
fig.update_traces(marker_line_width=0)
# fig.show()
fig.write_html(path_maps + "map_fract_K_destroyed.html")
fig.write_image(path_maps + "map_fract_K_destroyed.png",
                height=650, width=1000)

# We then plot estimated amenity index

amenities_df = pd.DataFrame({'amenities': amenities})
amenities_df['lon'] = geo_grid.lon
amenities_df['lat'] = geo_grid.lat
amenities_df.loc[
    (coeff_land[0, :] == 0) & (coeff_land[1, :] == 0)
    & (coeff_land[2, :] == 0) & (coeff_land[3, :] == 0), 'amenities'] = np.nan

fig = px.choropleth_mapbox(
    amenities_df,
    geojson=geo_grid.geometry,
    locations=geo_grid.index,
    color='amenities',
    center={"lat": -33.92345542582841, "lon": 18.434424141913478},
    zoom=9.25,
    mapbox_style='stamen-terrain',
    opacity=0.5,
    labels={'locations': 'Pixel ID', 'lon': 'Lon.', 'lat': 'Lat.',
            'amenities': 'Amenity index'},
    title='Estimated amenity index (in habitable areas)',
    color_continuous_scale="Picnic",
    color_continuous_midpoint=1,
    template='plotly_white',
    hover_data={'lon': ':.2f', 'lat': ':.2f',
                'amenities': ':.2f'})
fig.update_layout(margin={"r": 0, "t": 30, "l": 0, "b": 0})
fig.update_traces(marker_line_width=0)
# fig.show()
fig.write_html(path_maps + "map_amenities.html")
fig.write_image(path_maps + "map_amenities.png", height=650, width=1000)


# Finally, we plot theoretical income net of commuting costs (to picture better
# the choice set faced by households a priori)

income_net_of_commuting_costs_df = pd.DataFrame(
    {'1': income_net_of_commuting_costs[0, :],
     '2': income_net_of_commuting_costs[1, :],
     '3': income_net_of_commuting_costs[2, :],
     '4': income_net_of_commuting_costs[3, :]})
income_net_of_commuting_costs_df.loc[
    np.isnan(amenities_df['amenities']), '1'] = np.nan
income_net_of_commuting_costs_df.loc[
    np.isnan(amenities_df['amenities']), '2'] = np.nan
income_net_of_commuting_costs_df.loc[
    np.isnan(amenities_df['amenities']), '3'] = np.nan
income_net_of_commuting_costs_df.loc[
    np.isnan(amenities_df['amenities']), '4'] = np.nan
income_net_of_commuting_costs_df['lon'] = geo_grid.lon
income_net_of_commuting_costs_df['lat'] = geo_grid.lat

fig = px.choropleth_mapbox(
    income_net_of_commuting_costs_df,
    geojson=geo_grid.geometry,
    locations=geo_grid.index,
    color='1',
    center={"lat": -33.92345542582841, "lon": 18.434424141913478},
    zoom=9.25,
    mapbox_style='stamen-terrain',
    opacity=0.5,
    labels={'locations': 'Pixel ID', 'lon': 'Lon.', 'lat': 'Lat.',
            '1': 'Net income'},
    title='Estimated annual income net of commuting costs for income group 1'
    + ' (in rands, 2011)',
    color_continuous_scale="Reds",
    template='plotly_white',
    hover_data={'lon': ':.2f', 'lat': ':.2f',
                '1': ':,.0f'})
fig.update_layout(margin={"r": 0, "t": 30, "l": 0, "b": 0})
fig.update_traces(marker_line_width=0)
# fig.show()
fig.write_html(path_maps + "map_netincome_1.html")
fig.write_image(path_maps + "map_netincome_1.png", height=650, width=1000)

fig = px.choropleth_mapbox(
    income_net_of_commuting_costs_df,
    geojson=geo_grid.geometry,
    locations=geo_grid.index,
    color='2',
    center={"lat": -33.92345542582841, "lon": 18.434424141913478},
    zoom=9.25,
    mapbox_style='stamen-terrain',
    opacity=0.5,
    labels={'locations': 'Pixel ID', 'lon': 'Lon.', 'lat': 'Lat.',
            '2': 'Net income'},
    title='Estimated annual income net of commuting costs for income group 2'
    + ' (in rands, 2011)',
    color_continuous_scale="Reds",
    template='plotly_white',
    hover_data={'lon': ':.2f', 'lat': ':.2f',
                '2': ':,.0f'})
fig.update_layout(margin={"r": 0, "t": 30, "l": 0, "b": 0})
fig.update_traces(marker_line_width=0)
# fig.show()
fig.write_html(path_maps + "map_netincome_2.html")
fig.write_image(path_maps + "map_netincome_2.png", height=650, width=1000)

fig = px.choropleth_mapbox(
    income_net_of_commuting_costs_df,
    geojson=geo_grid.geometry,
    locations=geo_grid.index,
    color='3',
    center={"lat": -33.92345542582841, "lon": 18.434424141913478},
    zoom=9.25,
    mapbox_style='stamen-terrain',
    opacity=0.5,
    labels={'locations': 'Pixel ID', 'lon': 'Lon.', 'lat': 'Lat.',
            '3': 'Net income'},
    title='Estimated annual income net of commuting costs for income group 3'
    + ' (in rands, 2011)',
    color_continuous_scale="Reds",
    template='plotly_white',
    hover_data={'lon': ':.2f', 'lat': ':.2f',
                '3': ':,.0f'})
fig.update_layout(margin={"r": 0, "t": 30, "l": 0, "b": 0})
fig.update_traces(marker_line_width=0)
# fig.show()
fig.write_html(path_maps + "map_netincome_3.html")
fig.write_image(path_maps + "map_netincome_3.png", height=650, width=1000)

fig = px.choropleth_mapbox(
    income_net_of_commuting_costs_df,
    geojson=geo_grid.geometry,
    locations=geo_grid.index,
    color='4',
    center={"lat": -33.92345542582841, "lon": 18.434424141913478},
    zoom=9.25,
    mapbox_style='stamen-terrain',
    opacity=0.5,
    labels={'locations': 'Pixel ID', 'lon': 'Lon.', 'lat': 'Lat.',
            '4': 'Net income'},
    title='Estimated annual income net of commuting costs for income group 4'
    + ' (in rands, 2011)',
    color_continuous_scale="Reds",
    template='plotly_white',
    hover_data={'lon': ':.2f', 'lat': ':.2f',
                '4': ':,.0f'})
fig.update_layout(margin={"r": 0, "t": 30, "l": 0, "b": 0})
fig.update_traces(marker_line_width=0)
# fig.show()
fig.write_html(path_maps + "map_netincome_4.html")
fig.write_image(path_maps + "map_netincome_4.png", height=650, width=1000)

print("Input maps done")


# %% Save underlying data

# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter(path_use_case + 'agg_stat_abs_damages.xlsx',
                        engine='xlsxwriter')

# Write each dataframe to a different worksheet.
for key, value in dict_df_agg_damages.items():
    value.to_excel(writer, sheet_name=key)

# Close the Pandas Excel writer and output the Excel file.
writer.close()

# Save other data sets
damage_maps_concat.to_csv(path_use_case + 'damage_maps_concat.csv')
damage_map_compar.to_csv(path_use_case + 'damage_map_compar.csv')
damage_maps_shareinc_concat.to_csv(
    path_use_case + 'damage_maps_shareinc_concat.csv')
damage_map_compar_shareinc.to_csv(
    path_use_case + 'damage_map_compar_shareinc.csv')
equil_maps_concat.to_csv(path_use_case + 'equil_maps_concat.csv')
equil_map_compar.to_csv(path_use_case + 'equil_map_compar.csv')
agg_stat_rel_damages.to_csv(path_use_case + 'agg_stat_rel_damages.csv')

utility_noanticip = pd.DataFrame(
    dict_scenario_outcomes[list_scenarios[0]]['utility'])
utility_anticip = pd.DataFrame(
    dict_scenario_outcomes[list_scenarios[1]]['utility'])
utility_noanticip.to_csv(
    path_use_case + 'utility_noanticip.csv')
utility_anticip.to_csv(
    path_use_case + 'utility_anticip.csv')
