# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 15:33:37 2020.

@author: Charlotte Liotta
"""


# %% Preamble


# IMPORT PACKAGES

import numpy as np
import pandas as pd
import time
import os

import inputs.data as inpdt
import inputs.parameters_and_options as inpprm

import equilibrium.compute_equilibrium as eqcmp
import equilibrium.run_simulations as eqsim
import equilibrium.functions_dynamic as eqdyn

import calibration.calib_main_func as calmain


# DEFINE FILE PATHS

path_code = '..'
path_folder = path_code + '/2. Data/'
path_precalc_inp = path_folder + '0. Precalculated inputs/'
path_data = path_folder + 'data_Cape_Town/'
path_precalc_transp = path_folder + 'precalculated_transport/'
path_scenarios = path_folder + 'data_Cape_Town/Scenarios/'
path_outputs = path_code + '/4. Sorties/'
path_floods = path_folder + "FATHOM/"


# START TIMER FOR CODE OPTIMIZATION

start = time.process_time()


# %% Import parameters and options


# IMPORT PARAMETERS AND OPTIONS

options = inpprm.import_options()
param = inpprm.import_param(path_precalc_inp, path_outputs)

#  Set custom options for this simulation
options["agents_anticipate_floods"] = 1
options["informal_land_constrained"] = 0

#  More custom options regarding flood model
options["pluvial"] = 1
options["correct_pluvial"] = 1
options["coastal"] = 1
# This is in line with the DEM used in FATHOM data for fluvial and pluvial
options["dem"] = "MERITDEM"
options["slr"] = 1

#  Re-processing options
options["convert_sal_data"] = 0
options["compute_net_income"] = 0

#  Code correction options
options["actual_backyards"] = 0
options["unempl_reweight"] = 1
# TODO: recalibrate incomes net of commuting costs using implicit empl rate?
# implicit_empl_rate = 0.74/0.99/0.98/0.99
options["correct_agri_rent"] = 1

#  Options for calibration code correction
options["run_calib"] = 1
options["correct_dominant_incgrp"] = 0
options["substract_RDP_from_formal"] = 1
options["correct_mitchells_plain"] = 0
options["correct_selected_density"] = 1
options["correct_kappa"] = 1
#  TODO: set ranges for parameter scanning?

#  Set timeline for simulations
t = np.arange(0, 30)

# GIVE NAME TO SIMULATION TO EXPORT THE RESULTS
# (change according to custom parameters to be included)

name = ('allfloods_precal_modif')


# %% Load data


# BASIC GEOGRAPHIC DATA

grid, center = inpdt.import_grid(path_data)
amenities = inpdt.import_amenities(path_precalc_inp)


# MACRO DATA

(interest_rate, population, housing_type_data, total_RDP
 ) = inpdt.import_macro_data(param, path_scenarios)


# HOUSEHOLDS AND INCOME DATA

income_class_by_housing_type = inpdt.import_hypothesis_housing_type()

# See appendix A1
(mean_income, households_per_income_class, average_income, income_mult,
 income_2011, households_per_income_and_housing
 ) = inpdt.import_income_classes_data(param, path_data)

#  We create this parameter to maintain money illusion in simulations
#  (see eqsim.run_simulation)
#  TODO: Set as a variable, not a parameter
param["income_year_reference"] = mean_income

(data_rdp, housing_types_sp, data_sp, mitchells_plain_grid_2011,
 grid_formal_density_HFA, threshold_income_distribution, income_distribution,
 cape_town_limits) = inpdt.import_households_data(path_precalc_inp)

#  Import nb of households per pixel, by housing type.
#  Note that RDP is included in formal, and there are both formal and informal
#  backyards

if options["convert_sal_data"] == 1:
    housing_types = inpdt.import_sal_data(grid, path_folder, path_data,
                                          housing_type_data)

housing_types = pd.read_excel(path_folder + 'housing_types_grid_sal.xlsx')


# LAND USE PROJECTIONS

(spline_RDP, spline_estimate_RDP, spline_land_RDP,
 spline_land_backyard, spline_land_informal, spline_land_constraints,
 number_properties_RDP) = (
     inpdt.import_land_use(grid, options, param, data_rdp, housing_types,
                           housing_type_data, path_data, path_folder)
     )

# Correction needed with Charlotte's calibration
# TODO: check if still needed after recalibration
param["pockets"][
    (spline_land_informal(29) > 0) & (spline_land_informal(0) == 0)
    ] = 0.79

#  We correct areas for each housing type at baseline year for the amount of
#  constructible land in each type
coeff_land = inpdt.import_coeff_land(
    spline_land_constraints, spline_land_backyard, spline_land_informal,
    spline_land_RDP, param, 0)

#  We update land use parameters at baseline (relies on loaded data)
housing_limit = inpdt.import_housing_limit(grid, param)

#  TODO: plug outputs in a new variable (not param) and adapt linked functions
(param, minimum_housing_supply, agricultural_rent
 ) = inpprm.import_construction_parameters(
    param, grid, housing_types_sp, data_sp["dwelling_size"],
    mitchells_plain_grid_2011, grid_formal_density_HFA, coeff_land,
    interest_rate, options
    )

# FLOOD DATA (takes some time when agents anticipate floods)
#  TODO: create a new variable instead of storing in param
#  TODO: check if WBUS2 data is indeed deprecated
#  param = inpdt.infer_WBUS2_depth(housing_types, param, path_floods)
if options["agents_anticipate_floods"] == 1:
    (fraction_capital_destroyed, structural_damages_small_houses,
     structural_damages_medium_houses, structural_damages_large_houses,
     content_damages, structural_damages_type1, structural_damages_type2,
     structural_damages_type3a, structural_damages_type3b,
     structural_damages_type4a, structural_damages_type4b
     ) = inpdt.import_full_floods_data(options, param, path_folder,
                                       housing_type_data)
elif options["agents_anticipate_floods"] == 0:
    fraction_capital_destroyed = pd.DataFrame()
    fraction_capital_destroyed["structure_formal_2"] = np.zeros(24014)
    fraction_capital_destroyed["structure_formal_1"] = np.zeros(24014)
    fraction_capital_destroyed["structure_subsidized_2"] = np.zeros(24014)
    fraction_capital_destroyed["structure_subsidized_1"] = np.zeros(24014)
    fraction_capital_destroyed["contents_formal"] = np.zeros(24014)
    fraction_capital_destroyed["contents_informal"] = np.zeros(24014)
    fraction_capital_destroyed["contents_subsidized"] = np.zeros(24014)
    fraction_capital_destroyed["contents_backyard"] = np.zeros(24014)
    fraction_capital_destroyed["structure_backyards"] = np.zeros(24014)
    fraction_capital_destroyed["structure_informal_settlements"
                               ] = np.zeros(24014)

# SCENARIOS

(spline_agricultural_rent, spline_interest_rate,
 spline_population_income_distribution, spline_inflation,
 spline_income_distribution, spline_population,
 spline_income, spline_minimum_housing_supply, spline_fuel
 ) = eqdyn.import_scenarios(income_2011, param, grid, path_scenarios)

#  Import income net of commuting costs, as calibrated in Pfeiffer et al.
#  (see part 3.1 or appendix C3)

if options["compute_net_income"] == 1:
    for t_temp in t:
        print(t_temp)
        (incomeNetOfCommuting, modalShares, ODflows, averageIncome
         ) = inpdt.import_transport_data(
             grid, param, t_temp, households_per_income_class, average_income,
             spline_inflation, spline_fuel,
             spline_population_income_distribution, spline_income_distribution,
             path_precalc_inp, path_precalc_transp, 'GRID')

income_net_of_commuting_costs = np.load(
    path_precalc_transp + 'GRID_incomeNetOfCommuting_0.npy')


# %% Re-run calibration (takes time, only if needed)

if options["run_calib"] == 1:

    # We estimate the coefficients of construction function
    # Note that scale factor is significantly smaller than in paper
    coeff_b, coeff_a, coeffKappa = calmain.estim_construct_func_param(
        options, param, data_sp, threshold_income_distribution,
        income_distribution, data_rdp, housing_types_sp,
        path_data, path_precalc_inp)

    # We scan values for the gravity parameter to estimate incomes as a
    # function of it.
    # The value range is set by trial and error: the wider the range you want
    # to test, the longer.
    list_lambda = 10 ** np.arange(0.5, 0.7, 0.05)
    # list_lambda = 10 ** np.arange(0.6, 0.65, 0.01)
    # list_lambda = 10 ** np.arange(0.6, 0.61, 0.005)

    # TODO: update param vector


# %% Compute initial state

# TODO: Note that we use a Cobb-Douglas production function (with rho+delta)
# all along!
# Also note that we simulate households as two representative agents
# (not as in the paper)

# TODO: create option to run on old or new calibrated parameters

# population = sum(income_2011.Households_nb)
# param["max_iter"] = 10000

(initial_state_utility,
 initial_state_error,
 initial_state_simulated_jobs,
 initial_state_households_housing_types,
 initial_state_household_centers,
 initial_state_households,
 initial_state_dwelling_size,
 initial_state_housing_supply,
 initial_state_rent,
 initial_state_rent_matrix,
 initial_state_capital_land,
 initial_state_average_income,
 initial_state_limit_city) = eqcmp.compute_equilibrium(
     fraction_capital_destroyed,
     amenities,
     param,
     housing_limit,
     population,
     households_per_income_class,
     total_RDP,
     coeff_land,
     income_net_of_commuting_costs,
     grid,
     options,
     agricultural_rent,
     interest_rate,
     number_properties_RDP,
     average_income,
     mean_income,
     income_class_by_housing_type,
     minimum_housing_supply,
     param["coeff_A"])


# Reminder: income groups are ranked from poorer to richer, and housing types
# follow the following order: formal-backyard-informal-RDP

# Note on outputs (with dimensions in same order as axes):
# initial_state_utility = utility for each income group (no RDP)
#   after optimization
# initial_state_error = value of error term for each group after optimization
# initial_state_simulated_jobs = total number of households per housing type
#   (no RDP) and income group
# initial_state_households_housing_types = number of households
#   per housing type (with RDP) per pixel
# initial_state_household_centers = number of households per income group
#   per pixel
# initial_state_households = number of households in each housing type
#   and income group per pixel
# initial_state_dwelling_size = dwelling size (in m²) for each housing type
#   per pixel
# initial_state_housing_supply = housing surface built (in m²) per unit of
#   available land (in km²) for each housing type in each pixel
# initial_state_rent = average rent (in rands/m²) for each housing type
#   in each pixel
# initial_state_rent_matrix = average willingness to pay (in rands)
#   for each housing type (no RDP) and each income group in each pixel
# initial_state_capital_land = value of the (housing construction sector)
#   capital stock (in available-land unit equivalent) per unit of available
#   land (in km²) in each housing type (no RDP) and each selected pixel
# initial_state_average_income = average income per income group
#   (not an output of the model)
# initial_state_limit_city = indicator dummy for having strictly more
#   than one household per housing type and income group in each pixel

# Save outputs

try:
    os.mkdir(path_outputs + name)
except OSError as error:
    print(error)


np.save(path_outputs + name + '/initial_state_utility.npy',
        initial_state_utility)
np.save(path_outputs + name + '/initial_state_error.npy',
        initial_state_error)
np.save(path_outputs + name + '/initial_state_simulated_jobs.npy',
        initial_state_simulated_jobs)
np.save(path_outputs + name + '/initial_state_households_housing_types.npy',
        initial_state_households_housing_types)
np.save(path_outputs + name + '/initial_state_household_centers.npy',
        initial_state_household_centers)
np.save(path_outputs + name + '/initial_state_households.npy',
        initial_state_households)
np.save(path_outputs + name + '/initial_state_dwelling_size.npy',
        initial_state_dwelling_size)
np.save(path_outputs + name + '/initial_state_housing_supply.npy',
        initial_state_housing_supply)
np.save(path_outputs + name + '/initial_state_rent.npy',
        initial_state_rent)
np.save(path_outputs + name + '/initial_state_rent_matrix.npy',
        initial_state_rent_matrix)
np.save(path_outputs + name + '/initial_state_capital_land.npy',
        initial_state_capital_land)
np.save(path_outputs + name + '/initial_state_average_income.npy',
        initial_state_average_income)
np.save(path_outputs + name + '/initial_state_limit_city.npy',
        initial_state_limit_city)

# %% Scenarios

# NB: From simulation 22 onwards (with constraint), algorithm does not converge
# Note that this does not depend on calibration used!

# TODO: choose between right and original specification
# from scipy.interpolate import interp1d
# RDP_2011 = 2.2666e+05
# RDP_2001 = 1.1718e+05
# spline_RDP = interp1d(
#     [2001 - param["baseline_year"], 2011 - param["baseline_year"],
#      2018 - param["baseline_year"], 2041 - param["baseline_year"]],
#     [RDP_2001, RDP_2011, RDP_2011 + 7*5000,
#      RDP_2011 + 7*5000 + 23 * param["future_rate_public_housing"]], 'linear'
#     )

# RUN SIMULATION: time depends on the timeline (takes hours with 30 years)
(simulation_households_center,
 simulation_households_housing_type,
 simulation_dwelling_size,
 simulation_rent,
 simulation_households,
 simulation_error,
 simulation_housing_supply,
 simulation_utility,
 simulation_deriv_housing,
 simulation_T) = eqsim.run_simulation(
     t,
     options,
     param,
     grid,
     initial_state_utility,
     initial_state_error,
     initial_state_households,
     initial_state_households_housing_types,
     initial_state_housing_supply,
     initial_state_household_centers,
     initial_state_average_income,
     initial_state_rent,
     initial_state_dwelling_size,
     fraction_capital_destroyed,
     amenities,
     housing_limit,
     spline_estimate_RDP,
     spline_land_constraints,
     spline_land_backyard,
     spline_land_RDP,
     spline_land_informal,
     income_class_by_housing_type,
     path_precalc_transp,
     spline_RDP,
     spline_agricultural_rent,
     spline_interest_rate,
     spline_population_income_distribution,
     spline_inflation,
     spline_income_distribution,
     spline_population,
     spline_income,
     spline_minimum_housing_supply,
     spline_fuel
     )

# Save outputs

try:
    os.mkdir(path_outputs + name)
except OSError as error:
    print(error)

np.save(path_outputs + name + '/simulation_households_center.npy',
        simulation_households_center)
np.save(path_outputs + name + '/simulation_households_housing_type.npy',
        simulation_households_housing_type)
np.save(path_outputs + name + '/simulation_dwelling_size.npy',
        simulation_dwelling_size)
np.save(path_outputs + name + '/simulation_rent.npy',
        simulation_rent)
np.save(path_outputs + name + '/simulation_households.npy',
        simulation_households)
np.save(path_outputs + name + '/simulation_error.npy',
        simulation_error)
np.save(path_outputs + name + '/simulation_housing_supply.npy',
        simulation_housing_supply)
np.save(path_outputs + name + '/simulation_utility.npy',
        simulation_utility)
np.save(path_outputs + name + '/simulation_deriv_housing.npy',
        simulation_deriv_housing)
np.save(path_outputs + name + '/simulation_T.npy',
        simulation_T)

# NB: how do we model income and amenity changes with floods?
