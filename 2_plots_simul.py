# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 16:21:32 2022.

@author: monni
"""

# %% Preamble

# IMPORT PACKAGES

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl

import inputs.parameters_and_options as inpprm
import inputs.data as inpdt
import equilibrium.functions_dynamic as eqdyn
import outputs.export_outputs as outexp
import outputs.flood_outputs as outfld
import outputs.export_outputs_floods as outval

print("Import information to be used in the simulation")


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
options["defended"] = 1
#  Dummy for taking sea-level rise into account in coastal flood data
#  NB: Projections are up to 2050, based upon IPCC AR5 assessment for the
#  RCP 8.5 scenario
options["climate_change"] = 1

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


# %% Load data

print("Load data and results to be plotted as outputs")


# BASIC GEOGRAPHIC DATA

grid, center = inpdt.import_grid(path_data)
geo_grid = gpd.read_file(path_data + "grid_reference_500.shp")
geo_TAZ = gpd.read_file(path_data + "TAZ_ampp_prod_attr_2013_2032.shp")
amenities = inpdt.import_amenities(path_precalc_inp, options)


# MACRO DATA

(interest_rate, population, housing_type_data, total_RDP
 ) = inpdt.import_macro_data(param, path_scenarios, path_folder)


# HOUSEHOLDS AND INCOME DATA

(mean_income, households_per_income_class, average_income, income_mult,
 income_baseline, households_per_income_and_housing
 ) = inpdt.import_income_classes_data(param, path_data)

(data_rdp, housing_types_sp, data_sp, mitchells_plain_grid_baseline,
 grid_formal_density_HFA, threshold_income_distribution, income_distribution,
 cape_town_limits) = inpdt.import_households_data(path_precalc_inp)

housing_types = pd.read_excel(path_folder + 'housing_types_grid_sal.xlsx')
housing_types[np.isnan(housing_types)] = 0

# We convert income distribution data (at SP level) to grid dimensions for use
# in income calibration: long to run, uncomment only if needed

if options["convert_sp_data"] == 1:
    print("Convert SP data to grid dimensions - start")
    income_distribution_grid = inpdt.convert_income_distribution(
        income_distribution, grid, path_data, data_sp)
    print("Convert SP data to grid dimensions - end")

income_distribution_grid = np.load(path_data + "income_distrib_grid.npy")


# LAND USE PROJECTIONS

(spline_RDP, spline_estimate_RDP, spline_land_RDP,
 spline_land_backyard, spline_land_informal, spline_land_constraints,
 number_properties_RDP) = (
     inpdt.import_land_use(grid, options, param, data_rdp, housing_types,
                           housing_type_data, path_data, path_folder)
)

#  We correct areas for each housing type at baseline year for the amount of
#  constructible land in each type
coeff_land = inpdt.import_coeff_land(
    spline_land_constraints, spline_land_backyard, spline_land_informal,
    spline_land_RDP, param, 0)

#  We update parameter vector with construction parameters
(param, minimum_housing_supply, agricultural_rent
 ) = inpprm.import_construction_parameters(
    param, grid, housing_types_sp, data_sp["dwelling_size"],
    mitchells_plain_grid_baseline, grid_formal_density_HFA, coeff_land,
    interest_rate, options
    )

# LOAD EQUILIBRIUM DATA

initial_state_utility = np.load(
    path_outputs + name + '/initial_state_utility.npy')
initial_state_error = np.load(
    path_outputs + name + '/initial_state_error.npy')
initial_state_simulated_jobs = np.load(
    path_outputs + name + '/initial_state_simulated_jobs.npy')
initial_state_households_housing_types = np.load(
    path_outputs + name + '/initial_state_households_housing_types.npy')
initial_state_household_centers = np.load(
    path_outputs + name + '/initial_state_household_centers.npy')
initial_state_households = np.load(
    path_outputs + name + '/initial_state_households.npy')
initial_state_dwelling_size = np.load(
    path_outputs + name + '/initial_state_dwelling_size.npy')
initial_state_housing_supply = np.load(
    path_outputs + name + '/initial_state_housing_supply.npy')
initial_state_rent = np.load(
    path_outputs + name + '/initial_state_rent.npy')
initial_state_rent_matrix = np.load(
    path_outputs + name + '/initial_state_rent_matrix.npy')
initial_state_capital_land = np.load(
    path_outputs + name + '/initial_state_capital_land.npy')
initial_state_average_income = np.load(
    path_outputs + name + '/initial_state_average_income.npy')
initial_state_limit_city = np.load(
    path_outputs + name + '/initial_state_limit_city.npy')


# LOAD SIMULATION DATA (from main.py)

simulation_households_center = np.load(
    path_outputs + name + '/simulation_households_center.npy')
simulation_households_housing_type = np.load(
    path_outputs + name + '/simulation_households_housing_type.npy')
simulation_dwelling_size = np.load(
    path_outputs + name + '/simulation_dwelling_size.npy')
simulation_rent = np.load(
    path_outputs + name + '/simulation_rent.npy')
simulation_households_housing_type = np.load(
    path_outputs + name + '/simulation_households_housing_type.npy')
simulation_households = np.load(
    path_outputs + name + '/simulation_households.npy')
simulation_error = np.load(
    path_outputs + name + '/simulation_error.npy')
simulation_utility = np.load(
    path_outputs + name + '/simulation_utility.npy')
simulation_housing_supply = np.load(
    path_outputs + name + '/simulation_housing_supply.npy')
simulation_deriv_housing = np.load(
    path_outputs + name + '/simulation_deriv_housing.npy')
simulation_T = np.load(
    path_outputs + name + '/simulation_T.npy')


# LOAD FLOOD DATA

if options["agents_anticipate_floods"] == 1:
    (fraction_capital_destroyed, structural_damages_small_houses,
     structural_damages_medium_houses, structural_damages_large_houses,
     content_damages, structural_damages_type1, structural_damages_type2,
     structural_damages_type3a, structural_damages_type3b,
     structural_damages_type4a, structural_damages_type4b
     ) = inpdt.import_full_floods_data(
         options, param, path_folder)

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
    fraction_capital_destroyed["structure_formal_backyards"] = np.zeros(24014)
    fraction_capital_destroyed["structure_informal_backyards"
                               ] = np.zeros(24014)
    fraction_capital_destroyed["structure_informal_settlements"
                               ] = np.zeros(24014)


# SCENARIOS

#  We create this parameter to maintain money illusion in simulations
#  (see eqsim.run_simulation)
param["income_year_reference"] = mean_income

(spline_agricultural_price, spline_interest_rate,
 spline_population_income_distribution, spline_inflation,
 spline_income_distribution, spline_population,
 spline_income, spline_minimum_housing_supply, spline_fuel
 ) = eqdyn.import_scenarios(income_baseline, param, grid, path_scenarios,
                            options)

fluvialu_damages_2d_dyn = []
pluvial_damages_2d_dyn = []
coastal_damages_2d_dyn = []


# %% DYNAMICS

years_simul = np.arange(2011, 2011 + 30)

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(years_simul, simulation_utility[:, 0],
        color="maroon", label="Poor")
ax.plot(years_simul, simulation_utility[:, 1],
        color="red", label="Mid-poor")
ax.plot(years_simul, simulation_utility[:, 2],
        color="darkorange", label="Mid-rich")
ax.plot(years_simul, simulation_utility[:, 3],
        color="gold", label="Rich")
ax.set_ylim(0)
ax.yaxis.set_major_formatter(
    mpl.ticker.StrMethodFormatter('{x:,.0f}'))
plt.legend()
plt.tick_params(labelbottom=True)
plt.ylabel("Utility levels", labelpad=15)
plt.savefig(path_plots + 'evol_util_levels.png')
plt.close()

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(years_simul, np.nansum(simulation_households_center, 2)[:, 0],
        color="maroon", label="Poor")
ax.plot(years_simul, np.nansum(simulation_households_center, 2)[:, 1],
        color="red", label="Mid-poor")
ax.plot(years_simul, np.nansum(simulation_households_center, 2)[:, 2],
        color="darkorange", label="Mid-rich")
ax.plot(years_simul, np.nansum(simulation_households_center, 2)[:, 3],
        color="gold", label="Rich")
ax.set_ylim(0)
ax.yaxis.set_major_formatter(
    mpl.ticker.StrMethodFormatter('{x:,.0f}'))
plt.legend()
plt.tick_params(labelbottom=True)
plt.ylabel("Total number of households per income group", labelpad=15)
plt.savefig(path_plots + 'evol_nb_households_incgroup.png')
plt.close()

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(years_simul, np.nansum(simulation_households_housing_type, 2)[:, 0],
        color="gold", label="Formal")
ax.plot(years_simul, np.nansum(simulation_households_housing_type, 2)[:, 1],
        color="darkorange", label="Backyard")
ax.plot(years_simul, np.nansum(simulation_households_housing_type, 2)[:, 2],
        color="red", label="Informal")
ax.plot(years_simul, np.nansum(simulation_households_housing_type, 2)[:, 3],
        color="maroon", label="Subsidized")
ax.set_ylim(0)
ax.yaxis.set_major_formatter(
    mpl.ticker.StrMethodFormatter('{x:,.0f}'))
plt.legend()
plt.tick_params(labelbottom=True)
plt.ylabel("Total number of households per housing type", labelpad=15)
plt.savefig(path_plots + 'evol_nb_households_htype.png')
plt.close()


# %% BEGIN THE LOOP AFTER CREATING STORAGE VARIABLES

for year_temp in np.arange(0, 30):

    income_net_of_commuting_costs = np.load(
        path_precalc_transp
        + 'GRID_incomeNetOfCommuting_' + str(year_temp) + '.npy')
    cal_average_income = np.load(
        path_precalc_transp + 'GRID_averageIncome_' + str(year_temp) + '.npy')
    # modal_shares = np.load(
    #     path_precalc_transp + 'GRID_modalShares' + str(year_temp) + '.npy')
    # od_flows = np.load(
    #     path_precalc_transp + 'GRID_ODflows' + str(year_temp) + '.npy')

    # All that changes
    (average_income, households_per_income_class
     ) = eqdyn.compute_average_income(
         spline_population_income_distribution,
         spline_income_distribution, param, year_temp)
    (param["subsidized_structure_value"]
     ) = (param["subsidized_structure_value_ref"]
          * (spline_inflation(year_temp) / spline_inflation(0)))
    (param["informal_structure_value"]
     ) = (param["informal_structure_value_ref"]
          * (spline_inflation(year_temp) / spline_inflation(0)))
    mean_income = spline_income(year_temp)
    interest_rate = eqdyn.interpolate_interest_rate(
        spline_interest_rate, year_temp)
    population = spline_population(year_temp)
    total_RDP = spline_RDP(year_temp)
    minimum_housing_supply = spline_minimum_housing_supply(year_temp)
    # income_mult = average_income / mean_income
    number_properties_RDP = spline_estimate_RDP(year_temp)
    # Scale factor needs to move to create monetary illusion in the
    # model, e.g. housing supply should not change when currency
    # changes and all prices move: this is where the formula comes
    # from
    construction_param = (
        (mean_income / param["income_year_reference"])
        ** (- param["coeff_b"]) * param["coeff_A"]
    )
    coeff_land = inpdt.import_coeff_land(
        spline_land_constraints, spline_land_backyard,
        spline_land_informal, spline_land_RDP, param, year_temp)

    (param, minimum_housing_supply, agricultural_rent
     ) = inpprm.import_construction_parameters(
        param, grid, housing_types_sp, data_sp["dwelling_size"],
        mitchells_plain_grid_baseline, grid_formal_density_HFA, coeff_land,
        interest_rate, options
        )

    agricultural_rent = inpprm.compute_agricultural_rent(
        spline_agricultural_price(year_temp), construction_param,
        interest_rate, param, options)

    # POPULATION OUTPUTS

    path_plots_temp = path_plots + str(year_temp) + '/'
    path_tables_temp = path_tables + str(year_temp) + '/'

    try:
        os.mkdir(path_plots_temp)
    except OSError as error:
        print(error)

    try:
        os.mkdir(path_tables_temp)
    except OSError as error:
        print(error)

    #  IN ONE DIMENSION

    # Now, we validate overall households density across space

    # We do the same for total number of households across space,
    # housing types and income groups
    dist_HH_per_housing_1d = outexp.simul_pop_housing_types(
        grid, simulation_households_housing_type[year_temp, :, :],
        path_plots_temp, path_tables_temp
    )

    dist_HH_per_income_1d = outexp.simul_pop_income_groups(
        grid, simulation_households_center[year_temp, :, :],
        path_plots_temp, path_tables_temp
    )

    # We also plot income groups across space (in 1D) for each housing type,
    # even if we cannot validate such output
    (dist_HH_per_housing_and_income_1d
     ) = outexp.simul_pop_htype_income(
         grid, simulation_households[year_temp, :, :, :],
         path_plots_temp, path_tables_temp)

    #  IN TWO DIMENSIONS

    #  For overall households
    sim_nb_households_tot = np.nansum(
        simulation_households_housing_type[year_temp, :, :], 0)
    total_sim = outexp.export_map(
        sim_nb_households_tot, grid, geo_grid, path_plots_temp, 'total_sim',
        "Total number of households (simulation)",
        path_tables_temp,
        ubnd=np.nanquantile(sim_nb_households_tot, 0.9999))

    #  Per housing type
    sim_nb_households_formal = simulation_households_housing_type[
        year_temp, 0, :]
    formal_sim = outexp.export_map(
        sim_nb_households_formal, grid, geo_grid,
        path_plots_temp, 'formal_sim',
        "Number of households in formal private (simulation)",
        path_tables_temp,
        ubnd=np.nanquantile(sim_nb_households_formal, 0.9999))

    sim_nb_households_backyard = simulation_households_housing_type[
        year_temp, 1, :]
    backyard_sim = outexp.export_map(
        sim_nb_households_backyard, grid, geo_grid, path_plots_temp,
        'backyard_sim',
        "Number of households in informal backyard (simulation)",
        path_tables_temp,
        ubnd=np.nanquantile(sim_nb_households_backyard, 0.9999))

    sim_nb_households_informal = simulation_households_housing_type[
        year_temp, 2, :]
    informal_sim = outexp.export_map(
        sim_nb_households_informal, grid, geo_grid, path_plots_temp,
        'informal_sim',
        "Number of households in informal settlements (simulation)",
        path_tables_temp,
        ubnd=np.nanquantile(sim_nb_households_informal, 0.9999))

    data_nb_households_rdp = simulation_households_housing_type[
        year_temp, 3, :]
    rdp_sim = outexp.export_map(
        data_nb_households_rdp, grid, geo_grid, path_plots_temp, 'rdp_sim',
        "Number of households in formal subsidized (data)",
        path_tables_temp,
        ubnd=np.nanquantile(data_nb_households_rdp, 0.9999))

    #  Per income group
    sim_nb_households_poor = simulation_households_center[year_temp, 0, :]
    poor_sim = outexp.export_map(
        sim_nb_households_poor, grid, geo_grid, path_plots_temp, 'poor_sim',
        "Number of poor households (simulation)",
        path_tables_temp,
        ubnd=np.nanquantile(sim_nb_households_poor, 0.9999))
    sim_nb_households_midpoor = simulation_households_center[year_temp, 1, :]
    midpoor_sim = outexp.export_map(
        sim_nb_households_midpoor, grid, geo_grid,
        path_plots_temp, 'midpoor_sim',
        "Number of mid-poor households (simulation)",
        path_tables_temp,
        ubnd=np.nanquantile(sim_nb_households_midpoor, 0.9999))
    sim_nb_households_midrich = simulation_households_center[year_temp, 2, :]
    midrich_sim = outexp.export_map(
        sim_nb_households_midrich, grid, geo_grid,
        path_plots_temp, 'midrich_sim',
        "Number of mid-rich households (simulation)",
        path_tables_temp,
        ubnd=np.nanquantile(sim_nb_households_midrich, 0.9999))
    sim_nb_households_rich = simulation_households_center[year_temp, 3, :]
    rich_sim = outexp.export_map(
        sim_nb_households_rich, grid, geo_grid, path_plots_temp, 'rich_sim',
        "Number of rich households (simulation)",
        path_tables_temp,
        ubnd=np.nanquantile(sim_nb_households_rich, 0.9999))

    # HOUSING SUPPLY OUTPUTS

    # By plotting the housing supply per unit of available land, we may check
    # whether the bell-shaped curve of urban development holds
    avg_hsupply_1d = outexp.valid_housing_supply(
        grid, simulation_housing_supply[year_temp, :, :],
        path_plots_temp, path_tables_temp)

    # We now consider overall land to recover building density
    housing_supply = simulation_housing_supply[year_temp,
                                               :, :] * coeff_land * 0.25
    hsupply_noland_1d = outexp.valid_housing_supply_noland(
        grid, housing_supply, path_plots_temp, path_tables_temp)

    hsupply_tot = np.nansum(housing_supply, 0)
    hsupply_2d_sim = outexp.export_map(
        hsupply_tot, grid, geo_grid, path_plots_temp, 'hsupply_2d_sim',
        "Total housing supply (in m²)",
        path_tables_temp,
        ubnd=np.nanquantile(hsupply_tot, 0.9999))
    # FAR = np.nansum(housing_supply, 0) / (0.25 * 1000000)
    # FAR_2d_sim = outexp.export_map(
    #     FAR, grid, geo_grid, path_plots_temp,  'FAR_2d_sim',
    #     "Overall floor-area ratio",
    #     path_tables_temp,
    #     ubnd=0.3)

    hsupply_formal = housing_supply[0, :]
    hsupply_formal_2d_sim = outexp.export_map(
        hsupply_formal, grid, geo_grid,
        path_plots_temp, 'hsupply_formal_2d_sim',
        "Total housing supply in private formal (in m²)",
        path_tables_temp,
        ubnd=np.nanquantile(hsupply_formal, 0.9999))
    # FAR_formal = housing_supply[0, :] / (0.25 * 1000000)
    # FAR_formal_2d_sim = outexp.export_map(
    #     FAR_formal, grid, geo_grid, path_plots_temp, 'FAR_formal_2d_sim',
    #     "Floor-area ratio in formal private",
    #     path_tables_temp,
    #     ubnd=0.15)

    # Pb of validation in hyper-centre is also reflected in price
    # sim_HFA_dens_formal = (simulation_housing_supply[year_temp, 0, :]
    #                        / 1000000)
    # HFA_dens_formal_2d_sim = outexp.export_map(
    #     sim_HFA_dens_formal, grid, geo_grid, path_plots_temp,
    #     'HFA_dens_formal_2d_sim',
    #     "Households density in formal private HFA (simulation)",
    #     path_tables_temp,
    #     ubnd=1)

    hsupply_backyard = housing_supply[1, :]
    hsupply_backyard_2d_sim = outexp.export_map(
        hsupply_backyard, grid, geo_grid, path_plots_temp,
        'hsupply_backyard_2d_sim',
        "Total housing supply in informal backyards (in m²)",
        path_tables_temp,
        ubnd=np.nanquantile(hsupply_backyard, 0.9999))
    # FAR_backyard = housing_supply[1, :] / (0.25 * 1000000)
    # FAR_backyard_2d_sim = outexp.export_map(
    #     FAR_backyard, grid, geo_grid, path_plots_temp, 'FAR_backyard_2d_sim',
    #     "Floor-area ratio in informal backyards",
    #     path_tables_temp,
    #     ubnd=0.10)

    hsupply_informal = housing_supply[2, :]
    hsupply_informal_2d_sim = outexp.export_map(
        hsupply_informal, grid, geo_grid, path_plots_temp,
        'hsupply_informal_2d_sim',
        "Total housing supply in informal settlements (in m²)",
        path_tables_temp,
        ubnd=np.nanquantile(hsupply_informal, 0.9999))
    # FAR_informal = housing_supply[2, :] / (0.25 * 1000000)
    # FAR_informal_2d_sim = outexp.export_map(
    #     FAR_informal, grid, geo_grid, path_plots_temp,
    #     'FAR_informal_2d_sim',
    #     "Floor-area ratio in informal settlements",
    #     path_tables_temp,
    #     ubnd=0.30)

    hsupply_rdp = housing_supply[3, :]
    hsupply_rdp_2d_sim = outexp.export_map(
        hsupply_rdp, grid, geo_grid, path_plots_temp, 'hsupply_rdp_2d_sim',
        "Total housing supply in formal subsidized (in m²)",
        path_tables_temp,
        ubnd=np.nanquantile(hsupply_rdp, 0.9999))
    # FAR_rdp = housing_supply[3, :] / (0.25 * 1000000)
    # FAR_rdp_2d_sim = outexp.export_map(
    #     FAR_rdp, grid, geo_grid, path_plots_temp, 'FAR_rdp_2d_sim',
    #     "Floor-area ratio in formal subsidized",
    #     path_tables_temp,
    #     ubnd=0.10)

    # As we do not know surface of built land (just of available land),
    # we need to rely on dwelling size to compute build heigth in
    # formal private

    # HOUSING PRICE OUTPUTS

    # First in one dimension

    land_price_1d = outexp.simulation_housing_price(
        grid, simulation_rent[year_temp, :, :],
        interest_rate, param, center,
        housing_types_sp, path_plots_temp, path_tables_temp,
        land_price=1)
    housing_price_1d = outexp.simulation_housing_price(
        grid, simulation_rent[year_temp, :, :],
        interest_rate, param, center,
        housing_types_sp, path_plots_temp, path_tables_temp,
        land_price=0)

    # Then in two dimensions

    rent_formal_simul = simulation_rent[year_temp, 0, :].copy()
    housing_price_formal_2d_sim = outexp.export_map(
        rent_formal_simul, grid, geo_grid,
        path_plots_temp, 'rent_formal_2d_sim',
        "Simulated average housing rents per location (private formal)",
        path_tables_temp,
        ubnd=np.nanquantile(rent_formal_simul, 0.9999))
    rent_backyard_simul = simulation_rent[year_temp, 1, :].copy()
    housing_price_backyard_2d_sim = outexp.export_map(
        rent_backyard_simul, grid, geo_grid, path_plots_temp,
        'rent_backyard_2d_sim',
        "Simulated average housing rents per location (informal backyards)",
        path_tables_temp,
        ubnd=np.nanquantile(rent_backyard_simul, 0.9999))
    rent_informal_simul = simulation_rent[year_temp, 2, :].copy()
    housing_price_informal_2d_sim = outexp.export_map(
        rent_informal_simul, grid, geo_grid, path_plots_temp,
        'rent_informal_2d_sim',
        "Simulated average housing rents per location (informal settlements)",
        path_tables_temp,
        ubnd=np.nanquantile(rent_informal_simul, 0.9999))

    land_price = (
        (simulation_rent[year_temp, 0:3, :] * param["coeff_A"])
        ** (1 / param["coeff_a"])
        * param["coeff_a"]
        * (param["coeff_b"] / (interest_rate + param["depreciation_rate"]))
        ** (param["coeff_b"] / param["coeff_a"])
        / interest_rate
    )
    landprice_formal_simul = land_price[0, :].copy()
    land_price_formal_2d_sim = outexp.export_map(
        landprice_formal_simul, grid, geo_grid,
        path_plots_temp, 'landprice_formal_2d_sim',
        "Simulated average land prices per location (private formal)",
        path_tables_temp,
        ubnd=np.nanquantile(landprice_formal_simul, 0.9999))
    landprice_backyard_simul = land_price[1, :].copy()
    land_price_backyard_2d_sim = outexp.export_map(
        landprice_backyard_simul, grid, geo_grid,
        path_plots_temp, 'landprice_backyard_2d_sim',
        "Simulated average land prices per location (informal backyards)",
        path_tables_temp,
        ubnd=np.nanquantile(landprice_backyard_simul, 0.9999))
    landprice_informal_simul = land_price[2, :].copy()
    land_price_informal_2d_sim = outexp.export_map(
        landprice_informal_simul, grid, geo_grid,
        path_plots_temp, 'landprice_informal_2d_sim',
        "Simulated average land prices per location (informal settlements)",
        path_tables_temp,
        ubnd=np.nanquantile(landprice_informal_simul, 0.9999))

    # DWELLING SIZE OUTPUTS

    dwelling_size_1d = outexp.simul_housing_demand(
        grid, center, simulation_dwelling_size[year_temp, :, :],
        simulation_households_housing_type[year_temp, :, :],
        path_plots_temp, path_tables_temp)

    formal_dwelling_size = simulation_dwelling_size[year_temp, 0, :]
    dwelling_size_2d = outexp.export_map(
        formal_dwelling_size, grid, geo_grid,
        path_plots_temp, 'formal_dwellingsize_2d_sim',
        "Simulated average dwelling size per location (formal private)",
        path_tables_temp,
        ubnd=np.nanquantile(formal_dwelling_size, 0.9999))

    # TRANSPORT OUTPUTS

    #  Income net of commuting costs
    netincome_poor = income_net_of_commuting_costs[0, :]
    netincome_poor_2d_sim = outexp.export_map(
        netincome_poor, grid, geo_grid,
        path_plots_temp, 'netincome_poor_2d_sim',
        "Estimated income net of commuting costs (poor)",
        path_tables_temp,
        ubnd=np.nanquantile(netincome_poor, 0.9999))
    netincome_midpoor = income_net_of_commuting_costs[1, :]
    netincome_midpoor_2d_sim = outexp.export_map(
        netincome_midpoor, grid, geo_grid, path_plots_temp,
        'netincome_midpoor_2d_sim',
        "Estimated income net of commuting costs (mid-poor)",
        path_tables_temp,
        ubnd=np.nanquantile(netincome_midpoor, 0.9999))
    netincome_midrich = income_net_of_commuting_costs[2, :]
    netincome_midrich_2d_sim = outexp.export_map(
        netincome_midrich, grid, geo_grid, path_plots_temp,
        'netincome_midrich_2d_sim',
        "Estimated income net of commuting costs (mid-rich)",
        path_tables_temp,
        ubnd=np.nanquantile(netincome_midrich, 0.9999))
    netincome_rich = income_net_of_commuting_costs[3, :]
    netincome_rich_2d_sim = outexp.export_map(
        netincome_rich, grid, geo_grid,
        path_plots_temp, 'netincome_rich_2d_sim',
        "Estimated income net of commuting costs (rich)",
        path_tables_temp,
        ubnd=np.nanquantile(netincome_rich, 0.9999))

    (avg_income_net_of_commuting_1d
     ) = outexp.plot_income_net_of_commuting_costs(
         grid, income_net_of_commuting_costs,
         path_plots_temp, path_tables_temp)

    #  Average income

    # avgincome_poor = cal_average_income[0, :]
    # avgincome_poor_2d_sim = outexp.export_map(
    #     avgincome_poor, grid, geo_grid,
    #     path_plots_temp, 'avgincome_poor_2d_sim',
    #     "Estimated average income (poor)",
    #     path_tables_temp,
    #     ubnd=25000, lbnd=10000)
    # avgincome_midpoor = cal_average_income[1, :]
    # avgincome_midpoor_2d_sim = outexp.export_map(
    #     avgincome_midpoor, grid, geo_grid, path_plots_temp,
    #     'avgincome_midpoor_2d_sim',
    #     "Estimated average income (mid-poor)",
    #     path_tables_temp,
    #     ubnd=70000, lbnd=25000)
    # avgincome_midrich = cal_average_income[2, :]
    # avgincome_midrich_2d_sim = outexp.export_map(
    #     avgincome_midrich, grid, geo_grid, path_plots_temp,
    #     'avgincome_midrich_2d_sim',
    #     "Estimated average income (mid-rich)",
    #     path_tables_temp,
    #     ubnd=200000, lbnd=100000)
    # avgincome_rich = cal_average_income[3, :]
    # avgincome_rich_2d_sim = outexp.export_map(
    #     avgincome_rich, grid, geo_grid,
    #     path_plots_temp, 'avgincome_rich_2d_sim',
    #     "Estimated average income (rich)",
    #     path_tables_temp,
    #     ubnd=850000, lbnd=550000)

    # (avg_income_1d
    #  ) = outexp.plot_average_income(
    #      grid, cal_average_income, path_plots_temp, path_tables_temp)

    # We also conduct validation with overall average income
    # Also do fit in 1D
    # np.seterr(divide='ignore', invalid='ignore')
    # overall_avg_income = (
    #     cal_average_income
    #     * simulation_households_center[year_temp, :, :]
    #     / np.nansum(simulation_households_center[year_temp, :, :], 0))
    # overall_avg_income[np.isnan(overall_avg_income)] = 0
    # overall_avg_income = np.nansum(overall_avg_income, 0)

    # avgincome_all_2d_sim = outexp.export_map(
    #     overall_avg_income, grid, geo_grid, path_plots_temp,
    #     'avgincome_all_2d_sim',
    #     "Estimated average income (all income groups)",
    #     path_tables_temp,
    #     ubnd=850000)

    # FLOOD DAMAGES

    fluvialu_floods = ['FU_5yr', 'FU_10yr', 'FU_20yr', 'FU_50yr', 'FU_75yr',
                       'FU_100yr', 'FU_200yr', 'FU_250yr', 'FU_500yr',
                       'FU_1000yr']
    pluvial_floods = ['P_5yr', 'P_10yr', 'P_20yr', 'P_50yr', 'P_75yr',
                      'P_100yr', 'P_200yr', 'P_250yr', 'P_500yr', 'P_1000yr']
    coastal_floods = ['C_MERITDEM_' + str(options['climate_change']) + '_0000',
                      'C_MERITDEM_' + str(options['climate_change']) + '_0002',
                      'C_MERITDEM_' + str(options['climate_change']) + '_0005',
                      'C_MERITDEM_' + str(options['climate_change']) + '_0010',
                      'C_MERITDEM_' + str(options['climate_change']) + '_0025',
                      'C_MERITDEM_' + str(options['climate_change']) + '_0050',
                      'C_MERITDEM_' + str(options['climate_change']) + '_0100',
                      'C_MERITDEM_' + str(options['climate_change']) + '_0250']

    # We get damages per housing type for one representative household!
    content_cost = outfld.compute_content_cost(
        simulation_households_center[year_temp, :, :],
        simulation_housing_supply[year_temp, :, :],
        income_net_of_commuting_costs, param,
        fraction_capital_destroyed, simulation_rent[year_temp, :, :],
        simulation_dwelling_size[year_temp, :, :], interest_rate)

    # Note that capital is in monetary values
    formal_structure_cost = outfld.compute_formal_structure_cost(
        simulation_rent[year_temp, :, :], param, interest_rate, coeff_land,
        simulation_households_housing_type[year_temp, :, :], param["coeff_A"])

    # We re-import flood data to compute ex-post damages in case agents do not
    # anticipate floods

    (fraction_capital_destroyed, structural_damages_small_houses,
     structural_damages_medium_houses, structural_damages_large_houses,
     content_damages, structural_damages_type1, structural_damages_type2,
     structural_damages_type3a, structural_damages_type3b,
     structural_damages_type4a, structural_damages_type4b
     ) = inpdt.import_full_floods_data(options, param, path_folder)

    # Then we run the aggregate tables

    fluvialu_damages_sim = outfld.compute_damages(
        fluvialu_floods, path_floods, param, content_cost,
        sim_nb_households_formal, data_nb_households_rdp,
        sim_nb_households_informal, sim_nb_households_backyard,
        simulation_dwelling_size[year_temp, :, :],
        formal_structure_cost, content_damages,
        structural_damages_type4b, structural_damages_type4a,
        structural_damages_type2, structural_damages_type3a, options,
        spline_inflation, year_temp, path_tables_temp, 'fluvialu_sim')

    pluvial_damages_sim = outfld.compute_damages(
        pluvial_floods, path_floods, param, content_cost,
        sim_nb_households_formal, data_nb_households_rdp,
        sim_nb_households_informal, sim_nb_households_backyard,
        simulation_dwelling_size[year_temp, :, :],
        formal_structure_cost, content_damages,
        structural_damages_type4b, structural_damages_type4a,
        structural_damages_type2, structural_damages_type3a, options,
        spline_inflation, year_temp, path_tables_temp, 'pluvial_sim')

    coastal_damages_sim = outfld.compute_damages(
        coastal_floods, path_floods, param, content_cost,
        sim_nb_households_formal, data_nb_households_rdp,
        sim_nb_households_informal, sim_nb_households_backyard,
        simulation_dwelling_size[year_temp, :, :],
        formal_structure_cost, content_damages,
        structural_damages_type4b, structural_damages_type4a,
        structural_damages_type2, structural_damages_type3a, options,
        spline_inflation, year_temp, path_tables_temp, 'coastal_sim')

    # We get aggregate graphs

    # outval.simul_damages(
    #     fluvialu_damages_sim,
    #     path_plots_temp, 'fluvialu', options)
    # outval.simul_damages(
    #     pluvial_damages_sim,
    #     path_plots_temp, 'pluvial', options)
    # outval.simul_damages(
    #     coastal_damages_sim,
    #     path_plots_temp, 'coastal', options)

    # Now in two dimensions

    fluvialu_damages_2d_sim = outfld.compute_damages_2d(
        fluvialu_floods, path_floods, param, content_cost,
        sim_nb_households_formal, data_nb_households_rdp,
        sim_nb_households_informal, sim_nb_households_backyard,
        simulation_dwelling_size[year_temp, :, :],
        formal_structure_cost, content_damages,
        structural_damages_type4b, structural_damages_type4a,
        structural_damages_type2, structural_damages_type3a, options,
        spline_inflation, year_temp, path_tables_temp, 'fluvialu_sim')

    pluvial_damages_2d_sim = outfld.compute_damages_2d(
        pluvial_floods, path_floods, param, content_cost,
        sim_nb_households_formal, data_nb_households_rdp,
        sim_nb_households_informal, sim_nb_households_backyard,
        simulation_dwelling_size[year_temp, :, :],
        formal_structure_cost, content_damages,
        structural_damages_type4b, structural_damages_type4a,
        structural_damages_type2, structural_damages_type3a, options,
        spline_inflation, year_temp, path_tables_temp, 'pluvial_sim')

    coastal_damages_2d_sim = outfld.compute_damages_2d(
        coastal_floods, path_floods, param, content_cost,
        sim_nb_households_formal, data_nb_households_rdp,
        sim_nb_households_informal, sim_nb_households_backyard,
        simulation_dwelling_size[year_temp, :, :],
        formal_structure_cost, content_damages,
        structural_damages_type4b, structural_damages_type4a,
        structural_damages_type2, structural_damages_type3a, options,
        spline_inflation, year_temp, path_tables_temp, 'coastal_sim')

    # Hence the maps and shapefiles

    fluvialu_damages_2d_sim_stacked = np.stack(
        [df for df in fluvialu_damages_2d_sim.values()])
    fluvialu_formal_structure_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        fluvialu_formal_structure_2d_sim[j] = outfld.annualize_damages(
            fluvialu_damages_2d_sim_stacked[:, j, 0],
            'fluvialu', 'formal', options)
    fluvialu_subsidized_structure_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        fluvialu_subsidized_structure_2d_sim[j] = outfld.annualize_damages(
            fluvialu_damages_2d_sim_stacked[:, j, 1],
            'fluvialu', 'subsidized', options)
    fluvialu_informal_structure_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        fluvialu_informal_structure_2d_sim[j] = outfld.annualize_damages(
            fluvialu_damages_2d_sim_stacked[:, j, 2],
            'fluvialu', 'informal', options)
    fluvialu_backyard_structure_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        fluvialu_backyard_structure_2d_sim[j] = outfld.annualize_damages(
            fluvialu_damages_2d_sim_stacked[:, j, 3],
            'fluvialu', 'backyard', options)
    fluvialu_formal_content_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        fluvialu_formal_content_2d_sim[j] = outfld.annualize_damages(
            fluvialu_damages_2d_sim_stacked[:, j, 4],
            'fluvialu', 'formal', options)
    fluvialu_informal_content_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        fluvialu_informal_content_2d_sim[j] = outfld.annualize_damages(
            fluvialu_damages_2d_sim_stacked[:, j, 5],
            'fluvialu', 'informal', options)
    fluvialu_backyard_content_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        fluvialu_backyard_content_2d_sim[j] = outfld.annualize_damages(
            fluvialu_damages_2d_sim_stacked[:, j, 6],
            'fluvialu', 'backyard', options)
    fluvialu_subsidized_content_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        fluvialu_subsidized_content_2d_sim[j] = outfld.annualize_damages(
            fluvialu_damages_2d_sim_stacked[:, j, 7],
            'fluvialu', 'subsidized', options)

    pluvial_damages_2d_sim_stacked = np.stack(
        [df for df in pluvial_damages_2d_sim.values()])
    pluvial_formal_structure_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        pluvial_formal_structure_2d_sim[j] = outfld.annualize_damages(
            pluvial_damages_2d_sim_stacked[:, j, 0],
            'pluvial', 'formal', options)
    pluvial_subsidized_structure_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        pluvial_subsidized_structure_2d_sim[j] = outfld.annualize_damages(
            pluvial_damages_2d_sim_stacked[:, j, 1],
            'pluvial', 'subsidized', options)
    pluvial_informal_structure_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        pluvial_informal_structure_2d_sim[j] = outfld.annualize_damages(
            pluvial_damages_2d_sim_stacked[:, j, 2],
            'pluvial', 'informal', options)
    pluvial_backyard_structure_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        pluvial_backyard_structure_2d_sim[j] = outfld.annualize_damages(
            pluvial_damages_2d_sim_stacked[:, j, 3],
            'pluvial', 'backyard', options)
    pluvial_formal_content_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        pluvial_formal_content_2d_sim[j] = outfld.annualize_damages(
            pluvial_damages_2d_sim_stacked[:, j, 4],
            'pluvial', 'formal', options)
    pluvial_informal_content_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        pluvial_informal_content_2d_sim[j] = outfld.annualize_damages(
            pluvial_damages_2d_sim_stacked[:, j, 5],
            'pluvial', 'informal', options)
    pluvial_backyard_content_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        pluvial_backyard_content_2d_sim[j] = outfld.annualize_damages(
            pluvial_damages_2d_sim_stacked[:, j, 6],
            'pluvial', 'backyard', options)
    pluvial_subsidized_content_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        pluvial_subsidized_content_2d_sim[j] = outfld.annualize_damages(
            pluvial_damages_2d_sim_stacked[:, j, 7],
            'pluvial', 'subsidized', options)

    coastal_damages_2d_sim_stacked = np.stack(
        [df for df in coastal_damages_2d_sim.values()])
    coastal_formal_structure_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        coastal_formal_structure_2d_sim[j] = outfld.annualize_damages(
            coastal_damages_2d_sim_stacked[:, j, 0],
            'coastal', 'formal', options)
    coastal_subsidized_structure_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        coastal_subsidized_structure_2d_sim[j] = outfld.annualize_damages(
            coastal_damages_2d_sim_stacked[:, j, 1],
            'coastal', 'subsidized', options)
    coastal_informal_structure_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        coastal_informal_structure_2d_sim[j] = outfld.annualize_damages(
            coastal_damages_2d_sim_stacked[:, j, 2],
            'coastal', 'informal', options)
    coastal_backyard_structure_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        coastal_backyard_structure_2d_sim[j] = outfld.annualize_damages(
            coastal_damages_2d_sim_stacked[:, j, 3],
            'coastal', 'backyard', options)
    coastal_formal_content_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        coastal_formal_content_2d_sim[j] = outfld.annualize_damages(
            coastal_damages_2d_sim_stacked[:, j, 4],
            'coastal', 'formal', options)
    coastal_informal_content_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        coastal_informal_content_2d_sim[j] = outfld.annualize_damages(
            coastal_damages_2d_sim_stacked[:, j, 5],
            'coastal', 'informal', options)
    coastal_backyard_content_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        coastal_backyard_content_2d_sim[j] = outfld.annualize_damages(
            coastal_damages_2d_sim_stacked[:, j, 6],
            'coastal', 'backyard', options)
    coastal_subsidized_content_2d_sim = np.zeros(24014)
    for j in np.arange(24014):
        coastal_subsidized_content_2d_sim[j] = outfld.annualize_damages(
            coastal_damages_2d_sim_stacked[:, j, 7],
            'coastal', 'subsidized', options)

    list_annualized_2d_damages = [
        fluvialu_backyard_structure_2d_sim,
        fluvialu_backyard_content_2d_sim,
        fluvialu_subsidized_structure_2d_sim,
        fluvialu_subsidized_content_2d_sim,
        fluvialu_informal_structure_2d_sim,
        fluvialu_informal_content_2d_sim,
        fluvialu_formal_structure_2d_sim,
        fluvialu_formal_content_2d_sim,
        pluvial_backyard_structure_2d_sim,
        pluvial_backyard_content_2d_sim,
        pluvial_subsidized_structure_2d_sim,
        pluvial_subsidized_content_2d_sim,
        pluvial_informal_structure_2d_sim,
        pluvial_informal_content_2d_sim,
        pluvial_formal_structure_2d_sim,
        pluvial_formal_content_2d_sim,
        coastal_backyard_structure_2d_sim,
        coastal_backyard_content_2d_sim,
        coastal_subsidized_structure_2d_sim,
        coastal_subsidized_content_2d_sim,
        coastal_informal_structure_2d_sim,
        coastal_informal_content_2d_sim,
        coastal_formal_structure_2d_sim,
        coastal_formal_content_2d_sim]

    list_annualized_2d_damages_formal = [
        fluvialu_formal_structure_2d_sim,
        fluvialu_formal_content_2d_sim,
        pluvial_formal_structure_2d_sim,
        pluvial_formal_content_2d_sim,
        coastal_formal_structure_2d_sim,
        coastal_formal_content_2d_sim]

    list_annualized_2d_damages_informal = [
        fluvialu_informal_structure_2d_sim,
        fluvialu_informal_content_2d_sim,
        pluvial_informal_structure_2d_sim,
        pluvial_informal_content_2d_sim,
        coastal_informal_structure_2d_sim,
        coastal_informal_content_2d_sim]

    list_annualized_2d_damages_backyard = [
        fluvialu_backyard_structure_2d_sim,
        fluvialu_backyard_content_2d_sim,
        pluvial_backyard_structure_2d_sim,
        pluvial_backyard_content_2d_sim,
        coastal_backyard_structure_2d_sim,
        coastal_backyard_content_2d_sim]

    list_annualized_2d_damages_subsidized = [
        fluvialu_subsidized_structure_2d_sim,
        fluvialu_subsidized_content_2d_sim,
        pluvial_subsidized_structure_2d_sim,
        pluvial_subsidized_content_2d_sim,
        coastal_subsidized_structure_2d_sim,
        coastal_subsidized_content_2d_sim]

    # NB: need to retrieve name of item
    for item in list_annualized_2d_damages:
        try:
            outexp.export_map(item, grid, geo_grid,
                              path_plots_temp, outexp.retrieve_name(
                                  item, -1), "",
                              path_tables_temp,
                              ubnd=np.quantile(item[item > 0], 0.9),
                              lbnd=np.min(item[item > 0]))
        except IndexError:
            pass

    # We could plot fraction of capital destroyed separately for each
    # flood type, but it would be similar due to bath-tub model
    # NB: note that content damage function is the same for all housing types
    for col in fraction_capital_destroyed.columns:
        value = fraction_capital_destroyed[col]
        outexp.export_map(value, grid, geo_grid,
                          path_plots_temp, col + '_fract_K_destroyed', "",
                          path_tables_temp,
                          ubnd=np.nanquantile(value, 0.9999))

    # Graphs with share of annual income destroyed

    selected_net_income_formal = np.empty(24014)
    cond = np.argmax(simulation_households[year_temp, 0, :, :], axis=0)
    for j in np.arange(24014):
        selected_net_income_formal[j] = (
            income_net_of_commuting_costs[cond[j], j])

    for item in list_annualized_2d_damages_formal:
        new_item = item / selected_net_income_formal
        try:
            outexp.export_map(new_item, grid, geo_grid,
                              path_plots_temp,
                              outexp.retrieve_name(item, -1) + '_shareinc', "",
                              path_tables_temp,
                              ubnd=np.nanquantile(new_item, 0.9999))
        except IndexError:
            pass

    selected_net_income_rdp = np.empty(24014)
    cond = np.argmax(simulation_households[year_temp, 3, :, :], axis=0)
    for j in np.arange(24014):
        selected_net_income_rdp[j] = (
            income_net_of_commuting_costs[cond[j], j])

    for item in list_annualized_2d_damages_subsidized:
        new_item = item / selected_net_income_rdp
        try:
            outexp.export_map(new_item, grid, geo_grid,
                              path_plots_temp,
                              outexp.retrieve_name(item, -1) + '_shareinc', "",
                              path_tables_temp,
                              ubnd=np.nanquantile(new_item, 0.9999))
        except IndexError:
            pass

    selected_net_income_backyard = np.empty(24014)
    cond = np.argmax(simulation_households[year_temp, 1, :, :], axis=0)
    for j in np.arange(24014):
        selected_net_income_backyard[j] = (
            income_net_of_commuting_costs[cond[j], j])

    for item in list_annualized_2d_damages_backyard:
        new_item = item / selected_net_income_backyard
        try:
            outexp.export_map(new_item, grid, geo_grid,
                              path_plots_temp,
                              outexp.retrieve_name(item, -1) + '_shareinc', "",
                              path_tables_temp,
                              ubnd=np.nanquantile(new_item, 0.9999))
        except IndexError:
            pass

    selected_net_income_informal = np.empty(24014)
    cond = np.argmax(simulation_households[year_temp, 2, :, :], axis=0)
    for j in np.arange(24014):
        selected_net_income_informal[j] = (
            income_net_of_commuting_costs[cond[j], j])

    for item in list_annualized_2d_damages_informal:
        new_item = item / selected_net_income_informal
        try:
            outexp.export_map(new_item, grid, geo_grid,
                              path_plots_temp,
                              outexp.retrieve_name(item, -1) + '_shareinc', "",
                              path_tables_temp,
                              ubnd=np.nanquantile(new_item, 0.9999))
        except IndexError:
            pass

    fluvialu_damages_2d_dyn.append(fluvialu_damages_2d_sim)
    pluvial_damages_2d_dyn.append(pluvial_damages_2d_sim)
    coastal_damages_2d_dyn.append(coastal_damages_2d_sim)


# %% DYNAMICS: GET OUT OF LOOP AFTER STORING WHAT'S NEEDED

# We get the evolution of aggregate flood damages over time (due to spatial
# sorting)
outval.simul_damages_time(
    fluvialu_damages_2d_dyn, path_plots, path_tables, 'fluvialu', options)
outval.simul_damages_time(
    pluvial_damages_2d_dyn, path_plots, path_tables, 'pluvial', options)
outval.simul_damages_time(
    coastal_damages_2d_dyn, path_plots, path_tables, 'coastal', options)
