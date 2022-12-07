# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 16:49:54 2022

@author: monni
"""

# %% Preamble

# IMPORT PACKAGES

import itertools
import numpy as np
import pandas as pd
import geopandas as gpd

import plotly.io as pio
import plotly.express as px
import plotly.graph_objects as go


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

# TODO: import tables with insurance, but with and without climate change
# (as we need anticipations for spatial sorting to adapt)

# #For spatial distribution

# We retrieve all available dimension combinations in the data
flood_types = ['fluvialu', 'pluvial', 'coastal']
housing_types = ['formal', 'subsidized', 'informal', 'backyard']
# income_groups = ['poor', 'midpoor', 'midrich', 'rich']
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
    dict_damage_map_noinsur_shareinc[dim[0] + '_' + dim[1] + '_' + dim[2]] = table

dict_damage_map_insur_shareinc = {}
for dim in all_dim:
    table = pd.read_csv(
        path_outputs + 'floods10_F0_P11_C10_scenario232/tables/floods/'
        + dim[0] + '_' + dim[1] + '_' + dim[2] + '_2d_sim_shareinc.csv',
        names=['damage'], header=None)
    dict_damage_map_insur_shareinc[dim[0] + '_' + dim[1] + '_' + dim[2]] = table

dict_damage_map = {}
dict_damage_map["noinsur"] = dict_damage_map_noinsur
dict_damage_map["insur"] = dict_damage_map_insur
dict_damage_map["noinsur_shareinc"] = dict_damage_map_noinsur_shareinc
dict_damage_map["insur_shareinc"] = dict_damage_map_insur_shareinc

# Finally, we import key equilibrium outcomes

list_scenarios = ['floods00_F0_P11_C10_scenario232',
                  'floods10_F0_P11_C10_scenario232']
dict_scenario_outcomes = {}

for scenario in list_scenarios:
    initial_state_utility = np.load(
        path_outputs + scenario + '/initial_state_utility.npy')
    initial_state_households_housing_types = np.load(
        path_outputs + scenario + '/initial_state_households_housing_types.npy')
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
    dict_scenario_outcomes[scenario]['hh_per_htype'] = initial_state_households_housing_types
    dict_scenario_outcomes[scenario]['hh_per_incgroup'] = initial_state_household_centers
    dict_scenario_outcomes[scenario]['hh_full_dist'] = initial_state_households
    dict_scenario_outcomes[scenario]['dwelling_size'] = initial_state_dwelling_size
    dict_scenario_outcomes[scenario]['hsupply'] = initial_state_housing_supply
    dict_scenario_outcomes[scenario]['rent'] = initial_state_rent

    households_formal_inc = pd.DataFrame(initial_state_households[0, :, :].T)
    households_backyard_inc = pd.DataFrame(initial_state_households[1, :, :].T)
    households_informal_inc = pd.DataFrame(initial_state_households[2, :, :].T)
    households_subsidized_inc = pd.DataFrame(initial_state_households[3, :, :].T)
    list_var_temp = [households_formal_inc, households_backyard_inc,
                     households_informal_inc, households_subsidized_inc]

    for var in list_var_temp:
        var.loc[var[0] != 0, 0] = 1
        var.loc[var[1] != 0, 1] = 2
        var.loc[var[2] != 0, 2] = 3
        var.loc[var[3] != 0, 3] = 4
        var = var.sum(axis=1)

    dict_scenario_outcomes[scenario]['formal_incgroup'
                                     ] = households_formal_inc
    dict_scenario_outcomes[scenario]['backyard_incgroup'
                                     ] = households_backyard_inc
    dict_scenario_outcomes[scenario]['informal_incgroup'
                                     ] = households_informal_inc
    dict_scenario_outcomes[scenario]['subsidized_incgroup'
                                     ] = households_subsidized_inc


for col in households_formal_inc.columns:
    print(col)

# We also import the geography reference grid for plots
geo_grid = gpd.read_file(path_data + "grid_reference_500.shp")
geo_grid.to_crs(4326, inplace=True)
geo_grid['lon'] = geo_grid.centroid.x
geo_grid['lat'] = geo_grid.centroid.y

# TODO: repeat (more or less) the same exercise with damages as a share of
# income (when relevant). Also show achieved utility levels.

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
        pd.concat(dict_damage_map[item], axis=1).groupby(axis=1, level=0).sum())

    # We also define filters for column names to be used as part of data processing
    filter_names = ['fluvialu', 'pluvial', 'coastal', 'fluvialu_formal',
                    'fluvialu_subsidized', 'fluvialu_informal',
                    'fluvialu_backyard', 'pluvial_formal', 'pluvial_subsidized',
                    'pluvial_informal', 'pluvial_backyard', 'coastal_formal',
                    'coastal_subsidized', 'coastal_informal', 'coastal_backyard']
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

    # We obtain the sum of all damage values corresponding to each data category
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

    # Then, we recover the breakdown of damages across housing types, to display
    # as hover data as well

    # We initialize values for each housing type
    concat_damage_map['damage_formal'] = 0
    concat_damage_map['damage_subsidized'] = 0
    concat_damage_map['damage_informal'] = 0
    concat_damage_map['damage_backyard'] = 0

    # We associate each value to the corresponding sum of damages for the dominant
    # flood type only
    for housing_type in housing_types:
        concat_damage_map.loc[
            concat_damage_map['flood_type'] == 'fluvial', 'damage_' + housing_type
        ] = concat_damage_map['sum_fluvialu_' + housing_type]
        concat_damage_map.loc[
            concat_damage_map['flood_type'] == 'pluvial', 'damage_' + housing_type
        ] = concat_damage_map['sum_pluvial_' + housing_type]
        concat_damage_map.loc[
            concat_damage_map['flood_type'] == 'coastal', 'damage_' + housing_type
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

    # We store the output in the corresponding dictionary entry
    damage_maps[item] = concat_damage_map


for item in list_maps_shareinc:
    # We first concatenate all damage tables into one data frame
    concat_damage_map = (
        pd.concat(dict_damage_map[item], axis=1).groupby(axis=1, level=0).sum())

    # We also define filters for column names to be used as part of data processing
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

    # We obtain the sum of all damage values corresponding to each data category
    for key, value in dict_filters.items():
        concat_damage_map['sum_' + key] = (
            concat_damage_map[dict_filters[key]].sum(axis=1))

    # We take the maximum across aggregate flood type damages to display as our
    # main outcome of interest (cf. bath-tub model)
    concat_damage_map['max_val_formal'] = concat_damage_map[
        ['sum_fluvialu_formal', 'sum_pluvial_formal', 'sum_coastal_formal']
        ].max(axis=1)
    concat_damage_map['max_val_subsidized'] = concat_damage_map[
        ['sum_fluvialu_subsidized', 'sum_pluvial_subsidized', 'sum_coastal_subsidized']
        ].max(axis=1)
    concat_damage_map['max_val_informal'] = concat_damage_map[
        ['sum_fluvialu_informal', 'sum_pluvial_informal', 'sum_coastal_informal']
        ].max(axis=1)
    concat_damage_map['max_val_backyard'] = concat_damage_map[
        ['sum_fluvialu_backyard', 'sum_pluvial_backyard', 'sum_coastal_backyard']
        ].max(axis=1)

    # We do the same for content share only (to display as hover data)
    concat_damage_map['max_content_formal'] = concat_damage_map[
        ['fluvialu_formal_content', 'pluvial_formal_content', 'coastal_formal_content']
        ].max(axis=1)
    concat_damage_map['max_content_subsidized'] = concat_damage_map[
        ['fluvialu_subsidized_content', 'pluvial_subsidized_content', 'coastal_subsidized_content']
        ].max(axis=1)
    concat_damage_map['max_content_informal'] = concat_damage_map[
        ['fluvialu_informal_content', 'pluvial_informal_content', 'coastal_informal_content']
        ].max(axis=1)
    concat_damage_map['max_content_backyard'] = concat_damage_map[
        ['fluvialu_backyard_content', 'pluvial_backyard_content', 'coastal_backyard_content']
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

    # Also add some population info!
    if item == 'noinsur_shareinc':
        concat_damage_map['nb_households_formal'] = (
            dict_scenario_outcomes['floods00_F0_P11_C10_scenario232'
                                   ]['hh_per_htype'][0, :])
        concat_damage_map['nb_households_subsidized'] = (
            dict_scenario_outcomes['floods00_F0_P11_C10_scenario232'
                                   ]['hh_per_htype'][3, :])
        concat_damage_map['nb_households_informal'] = (
            dict_scenario_outcomes['floods00_F0_P11_C10_scenario232'
                                   ]['hh_per_htype'][2, :])
        concat_damage_map['nb_households_backyard'] = (
            dict_scenario_outcomes['floods00_F0_P11_C10_scenario232'
                                   ]['hh_per_htype'][1, :])
        concat_damage_map['incgroup_formal']
        concat_damage_map['incgroup_subsidized']
        concat_damage_map['incgroup_informal']
        concat_damage_map['incgroup_backyard']
    elif item == 'insur_shareinc':
        concat_damage_map['nb_households_formal'] = (
            dict_scenario_outcomes['floods10_F0_P11_C10_scenario232'
                                   ]['hh_per_htype'][0, :])
        concat_damage_map['nb_households_subsidized'] = (
            dict_scenario_outcomes['floods10_F0_P11_C10_scenario232'
                                   ]['hh_per_htype'][3, :])
        concat_damage_map['nb_households_informal'] = (
            dict_scenario_outcomes['floods10_F0_P11_C10_scenario232'
                                   ]['hh_per_htype'][2, :])
        concat_damage_map['nb_households_backyard'] = (
            dict_scenario_outcomes['floods10_F0_P11_C10_scenario232'
                                   ]['hh_per_htype'][1, :])
        concat_damage_map['incgroup_formal']
        concat_damage_map['incgroup_subsidized']
        concat_damage_map['incgroup_informal']
        concat_damage_map['incgroup_backyard']
    

    # We convert nans to zeros for proper rendering in output map
    concat_damage_map = concat_damage_map.fillna(0)

    # We store the output in the corresponding dictionary entry
    damage_maps_shareinc[item] = concat_damage_map

# Then, we create a comparative map across scenarios to better see changes

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


# %% Output graphs

print("Output graphs")

# PLOT AGGREGATE DAMAGES

# TODO: In which unit should we convert outputs?
# NB: We could do a loop for all 3 plots and store them in a dictionary

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
fig_fluvialu.write_html(path_outputs + "fluvialu_sim_damage_sum_insur.html")
fig_fluvialu.write_image(path_outputs + "fluvialu_sim_damage_sum_insur.png")

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
fig_pluvial.write_html(path_outputs + "pluvial_sim_damage_sum_insur.html")
fig_pluvial.write_image(path_outputs + "pluvial_sim_damage_sum_insur.png")

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
fig_coastal.write_html(path_outputs + "coastal_sim_damage_sum_insur.html")
fig_coastal.write_image(path_outputs + "coastal_sim_damage_sum_insur.png")

print("Aggregate coastal damages done")


# PLOT SPATIAL DAMAGE DISTRIBUTION

# We plot all the relevant information for a given scenario as a chloropleth
# map

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
    fig.write_html(path_outputs + "map_damages_" + item + ".html")

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
fig.write_html(path_outputs + "map_damages_compar_insur.html")

print("map_damages_compar_insur done")

###

dist = px.histogram(
    damages_fluvialu_sim_insur,
    x=damages_fluvialu_sim_insur.index,
    y=['Structures', 'Contents'],
    barmode='group',
    color_discrete_sequence=px.colors.qualitative.Set1,
    title='Estimated annual damages from fluvial floods (in M rands, 2011)',
    labels={'index': 'Housing type', 'value': 'Annual damages',
            'variable': 'Damage type'},
    template='plotly_white')
