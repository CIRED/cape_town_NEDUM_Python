# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 15:31:33 2022

@author: monni
"""

# %% Preamble

# IMPORT PACKAGES

import pandas as pd
# import xslxwriter

import inputs.parameters_and_options as inpprm
import outputs.flood_outputs as outfld


# DEFINE FILE PATHS

path_code = '..'
path_folder = path_code + '/Data/'
path_precalc_inp = path_folder + 'precalculated_inputs/'
path_data = path_folder + 'data_Cape_Town/'
path_precalc_transp = path_folder + 'precalculated_transport/'
path_scenarios = path_data + 'Scenarios/'
path_outputs = path_code + '/Output/'
path_floods = path_folder + "flood_maps/"


# IMPORT PARAMETERS AND OPTIONS

options = inpprm.import_options()
param = inpprm.import_param(
    path_precalc_inp, options)

# Set custom options for this simulation
#  Dummy for taking floods into account in the utility function
options["agents_anticipate_floods"] = 1
#  Dummy for preventing new informal settlement development
options["informal_land_constrained"] = 0

# More custom options regarding flood model
#  Dummy for taking pluvial floods into account (on top of fluvial floods)
options["pluvial"] = 1
#  Dummy for reducing pluvial risk for (better protected) formal structures
options["correct_pluvial"] = 1
#  Dummy for taking coastal floods into account (on top of fluvial floods)
options["coastal"] = 1
#  Digital elevation to be used with coastal flood data (MERITDEM or NASADEM)
#  NB: MERITDEM is also the DEM used for fluvial and pluvial flood data
options["dem"] = "MERITDEM"
#  We consider undefended flood maps as our default because they are more
#  reliable
options["defended"] = 0
#  Dummy for taking sea-level rise into account in coastal flood data
#  NB: Projections are up to 2050, based upon IPCC AR5 assessment for the
#  RCP 8.5 scenario
options["climate_change"] = 0

# More custom options regarding scenarios
options["inc_ineq_scenario"] = 2
options["pop_growth_scenario"] = 3
options["fuel_price_scenario"] = 2

# Processing options for this simulation
options["convert_sp_data"] = 0


# GIVE NAME TO SIMULATION TO EXPORT THE RESULTS
# (change according to custom parameters to be included)

name = ('floods' + str(options["agents_anticipate_floods"])
        + str(options["informal_land_constrained"])
        + '_F' + str(options["defended"])
        + '_P' + str(options["pluvial"]) + str(options["correct_pluvial"])
        + '_C' + str(options["coastal"]) + str(options["climate_change"])
        + '_scenario' + str(options["inc_ineq_scenario"])
        + str(options["pop_growth_scenario"])
        + str(options["fuel_price_scenario"]))

path_plots = path_outputs + name + '/plots/'
path_tables = path_outputs + name + '/tables/'
path_plots_floods = path_plots + 'floods/'
path_tables_floods = path_tables + 'floods/'

# CREATE A DICTIONARY OF FLOOD DAMAGE DATA

flood_categ = ['fluvialu_data', 'fluvialu_sim',
               'pluvial_data', 'pluvial_sim',
               'coastal_data', 'coastal_sim']

housing_types = ['formal', 'subsidized', 'informal', 'backyard']

damage_data_dict = {key: pd.read_csv(path_tables_floods + key + '_damages.csv')
                    for key in flood_categ}

# ANNUALIZE DAMAGES

# First for fluvial floods

fluvialu_struct = [
    outfld.annualize_damages(
        damage_data_dict['fluvialu_sim'].formal_structure_damages,
        'fluvialu', 'formal', options) / 1000000,
    outfld.annualize_damages(
        damage_data_dict['fluvialu_sim'].subsidized_structure_damages,
        'fluvialu', 'subsidized', options) / 1000000,
    outfld.annualize_damages(
        damage_data_dict['fluvialu_sim'].informal_structure_damages,
        'fluvialu', 'informal', options) / 1000000,
    outfld.annualize_damages(
        damage_data_dict['fluvialu_sim'].backyard_structure_damages,
        'fluvialu', 'backyard', options) / 1000000]
fluvialu_struct = pd.DataFrame(
    fluvialu_struct, housing_types, ['struct_damage'])

fluvialu_content = [
    outfld.annualize_damages(
        damage_data_dict['fluvialu_sim'].formal_content_damages,
        'fluvialu', 'formal', options) / 1000000,
    outfld.annualize_damages(
         damage_data_dict['fluvialu_sim'].subsidized_content_damages,
         'fluvialu', 'subsidized',
         options) / 1000000,
    outfld.annualize_damages(
         damage_data_dict['fluvialu_sim'].informal_content_damages,
         'fluvialu', 'informal',
         options) / 1000000,
    outfld.annualize_damages(
         damage_data_dict['fluvialu_sim'].backyard_content_damages,
         'fluvialu', 'backyard',
         options) / 1000000
    ]
fluvialu_content = pd.DataFrame(
    fluvialu_content, housing_types, ['content_damage'])

fluvialu_sim = pd.merge(fluvialu_struct, fluvialu_content,
                        left_index=True, right_index=True)

# Then for pluvial floods

pluvial_struct = [
    outfld.annualize_damages(
        damage_data_dict['pluvial_sim'].formal_structure_damages,
        'pluvial', 'formal', options) / 1000000,
    outfld.annualize_damages(
        damage_data_dict['pluvial_sim'].subsidized_structure_damages,
        'pluvial', 'subsidized', options) / 1000000,
    outfld.annualize_damages(
        damage_data_dict['pluvial_sim'].informal_structure_damages,
        'pluvial', 'informal', options) / 1000000,
    outfld.annualize_damages(
        damage_data_dict['pluvial_sim'].backyard_structure_damages,
        'pluvial', 'backyard', options) / 1000000
    ]
pluvial_struct = pd.DataFrame(pluvial_struct, housing_types, ['struct_damage'])

pluvial_content = [
    outfld.annualize_damages(
        damage_data_dict['pluvial_sim'].formal_content_damages,
        'pluvial', 'formal', options) / 1000000,
    outfld.annualize_damages(
        damage_data_dict['pluvial_sim'].subsidized_content_damages,
        'pluvial', 'subsidized', options) / 1000000,
    outfld.annualize_damages(
        damage_data_dict['pluvial_sim'].informal_content_damages,
        'pluvial', 'informal', options) / 1000000,
    outfld.annualize_damages(
        damage_data_dict['pluvial_sim'].backyard_content_damages,
        'pluvial', 'backyard', options) / 1000000
    ]
pluvial_content = pd.DataFrame(
    pluvial_content, housing_types, ['content_damage'])

pluvial_sim = pd.merge(pluvial_struct, pluvial_content,
                       left_index=True, right_index=True)

# Finally for coastal floods

coastal_struct = [
    outfld.annualize_damages(
        damage_data_dict['coastal_sim'].formal_structure_damages,
        'coastal', 'formal', options) / 1000000,
    outfld.annualize_damages(
        damage_data_dict['coastal_sim'].subsidized_structure_damages,
        'coastal', 'subsidized', options) / 1000000,
    outfld.annualize_damages(
        damage_data_dict['coastal_sim'].informal_structure_damages,
        'coastal', 'informal', options) / 1000000,
    outfld.annualize_damages(
        damage_data_dict['coastal_sim'].backyard_structure_damages,
        'coastal', 'backyard', options) / 1000000
    ]
coastal_struct = pd.DataFrame(coastal_struct, housing_types, ['struct_damage'])

coastal_content = [
    outfld.annualize_damages(
        damage_data_dict['coastal_sim'].formal_content_damages,
        'coastal', 'formal', options) / 1000000,
    outfld.annualize_damages(
        damage_data_dict['coastal_sim'].subsidized_content_damages,
        'coastal', 'subsidized', options) / 1000000,
    outfld.annualize_damages(
        damage_data_dict['coastal_sim'].informal_content_damages,
        'coastal', 'informal', options) / 1000000,
    outfld.annualize_damages(
        damage_data_dict['coastal_sim'].backyard_content_damages,
        'coastal', 'backyard', options) / 1000000
    ]
coastal_content = pd.DataFrame(
    coastal_content, housing_types, ['content_damage'])

coastal_sim = pd.merge(coastal_struct, coastal_content,
                       left_index=True, right_index=True)


# SAVE OUTPUTS FOR THE SIMULATION

df_dict = {'fluvialu_sim': fluvialu_sim,
           'pluvial_sim': pluvial_sim,
           'coastal_sim': coastal_sim}

# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter(path_outputs + name + '/damage_data.xlsx',
                        engine='xlsxwriter')

# Write each dataframe to a different worksheet.
for key, value in df_dict.items():
    value.to_excel(writer, sheet_name=key)

# Close the Pandas Excel writer and output the Excel file.
writer.save()
