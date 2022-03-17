# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 16:40:37 2022.

@author: vincentviguie
"""

# %% Preamble


# IMPORT PACKAGES

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import inputs.data as inpdt
import inputs.parameters_and_options as inpprm


# DEFINE FILE PATHS

path_code = '..'
path_folder = path_code + '/2. Data/'
path_precalc_inp = path_folder + '0. Precalculated inputs/'
path_outputs = path_code + '/4. Sorties/'

name = 'carbon_tax_car_bus_taxi_20211103_basile'


# LOAD SIMULATION DATA (from main.py)

simulation_households_center = np.load(
    path_outputs + name + '/simulation_households_center.npy')
simulation_dwelling_size = np.load(
    path_outputs + name + '/simulation_dwelling_size.npy')
simulation_rent = np.load(path_outputs + name + '/simulation_rent.npy')
simulation_households_housing_type = np.load(
    path_outputs + name + '/simulation_households_housing_type.npy')
simulation_households = np.load(
    path_outputs + name + '/simulation_households.npy')
simulation_utility = np.load(path_outputs + name + '/simulation_utility.npy')

income = np.load(
    path_folder+"precalculated_transport/incomeNetofCommuting_29.npy")

options = inpprm.import_options()
param = inpprm.import_param(path_precalc_inp, path_outputs)

if options["agents_anticipate_floods"] == 0:
    fraction_capital_destroyed, *_ = inpdt.import_floods_data(
        options, param, path_folder
        )  # need to add parameters?


# %% Plot output scenarios

# 0. Housing types

data = pd.DataFrame(
    {'2011': np.nansum(simulation_households_housing_type[0, :, :], 1),
     '2020': np.nansum(simulation_households_housing_type[9, :, :], 1),
     '2030': np.nansum(simulation_households_housing_type[19, :, :], 1),
     '2040': np.nansum(simulation_households_housing_type[28, :, :], 1)},
    index=["Formal private", "Informal in \n backyards",
           "Informal \n settlements", "Formal subsidized"]
    )
data.plot(kind="bar")
plt.tick_params(labelbottom=True)
plt.xticks(rotation='horizontal')
plt.ylabel("Number of households")
plt.ylim(0, 880000)


hh_2011 = simulation_households[0, :, :, :]
hh_2020 = simulation_households[9, :, :, :]
hh_2030 = simulation_households[19, :, :, :]
hh_2040 = simulation_households[29, :, :, :]

np.nansum(hh_2011, 2)
np.nansum(hh_2020, 2)
np.nansum(hh_2030, 2)
np.nansum(hh_2040, 2)

#  Formal
class_income = 3
income_class_2011 = np.argmax(simulation_households[29, :, :, :], 1)
subset = income_class_2011[0, :] == class_income
q = simulation_dwelling_size[29, 0, :][subset]
r = simulation_rent[29, 0, :][subset]
Y = income[class_income, :][subset]
B = 1
z = (Y - q*r) / (1 + (fraction_capital_destroyed.contents_formal[subset] * param["fraction_z_dwellings"]))
U = (z ** (1 - param["beta"])) * (q ** param["beta"]) * amenities[subset] * B

#class 3 (1/0) = 80300
#class 2 (1/0) = 17200
#class 1 (1/0) = 5000
#class 0 (1/0) = 500


#informal

class_income = 0
income_class_2011 = np.argmax(simulation_households[0, :, :, :], 1) 
subset = income_class_2011[2, :] == class_income
q = simulation_dwelling_size[0, 2, :][subset]
r = simulation_rent[0, 2, :][subset]
Y = income[class_income, :][subset]
B = np.load(path_outputs+'fluvial_and_pluvial/param_pockets.npy')

z = (Y - (q*r) - (fraction_capital_destroyed.structure_informal_settlements[subset] * (param["informal_structure_value"] * (spline_inflation(29) / spline_inflation(0))))) / (1 + (fraction_capital_destroyed.contents_informal[subset] * param["fraction_z_dwellings"]))
U = (z ** (1 - param["beta"])) * (q ** param["beta"]) * amenities[subset] * B[subset]

#class 3 (1/0) = ))))à#class 2 (1
#class 1 (1/0) = 5295
#class 0 (1/0) = 1885

#backyards

class_income = 1
income_class_2011 = np.argmax(simulation_households[0, :, :, :], 1) 
subset = income_class_2011[1, :] == class_income
q = simulation_dwelling_size[0, 1, :][subset]
r = simulation_rent[0, 1, :][subset]
Y = income[class_income, :][subset]
B = np.load(path_outputs+'fluvial_and_pluvial/param_backyards.npy')

z = (Y - (q*r) - (fraction_capital_destroyed.structure_backyards[subset] * (param["informal_structure_value"] * (spline_inflation(29) / spline_inflation(0))))) / (1 + (fraction_capital_destroyed.contents_backyard[subset] * param["fraction_z_dwellings"]))
U = (z ** (1 - param["beta"])) * (q ** param["beta"]) * amenities[subset] * B[subset]

#class 3 (1/0) = 
#class 2 (1/0) = 
#class 1 (1/0) = 5295
#class 0 (1/0) = 1885



sns.distplot(U, hist = True, kde = False)
np.nanmedian(U)
#### 1. Exposition des pop vulnérables

##  A. Evolution des prop flood prone areas

## FLUVIAL
stats_per_housing_type_2011_fluvial = compute_stats_per_housing_type(fluvial_floods, path_data, simulation_households_housing_type[0, 0, :], simulation_households_housing_type[0, 3, :], simulation_households_housing_type[0, 2, :], simulation_households_housing_type[0, 1, :], options, param, 0.01)
stats_per_housing_type_2040_fluvial = compute_stats_per_housing_type(fluvial_floods, path_data, simulation_households_housing_type[28, 0, :], simulation_households_housing_type[28, 3, :], simulation_households_housing_type[28, 2, :], simulation_households_housing_type[28, 1, :], options, param, 0.01)

label = ["Formal private", "Formal subsidized", "Informal \n settlements", "Informal \n in backyards"]
stats_2011_1 = [stats_per_housing_type_2011_fluvial.fraction_formal_in_flood_prone_area[2], stats_per_housing_type_2011_fluvial.fraction_subsidized_in_flood_prone_area[2], stats_per_housing_type_2011_fluvial.fraction_informal_in_flood_prone_area[2], stats_per_housing_type_2011_fluvial.fraction_backyard_in_flood_prone_area[2]]
stats_2011_2 = [stats_per_housing_type_2011_fluvial.fraction_formal_in_flood_prone_area[3], stats_per_housing_type_2011_fluvial.fraction_subsidized_in_flood_prone_area[3], stats_per_housing_type_2011_fluvial.fraction_informal_in_flood_prone_area[3], stats_per_housing_type_2011_fluvial.fraction_backyard_in_flood_prone_area[3]]
stats_2011_3 = [stats_per_housing_type_2011_fluvial.fraction_formal_in_flood_prone_area[5], stats_per_housing_type_2011_fluvial.fraction_subsidized_in_flood_prone_area[5], stats_per_housing_type_2011_fluvial.fraction_informal_in_flood_prone_area[5], stats_per_housing_type_2011_fluvial.fraction_backyard_in_flood_prone_area[5]]
stats_2040_1 = [stats_per_housing_type_2040_fluvial.fraction_formal_in_flood_prone_area[2], stats_per_housing_type_2040_fluvial.fraction_subsidized_in_flood_prone_area[2], stats_per_housing_type_2040_fluvial.fraction_informal_in_flood_prone_area[2], stats_per_housing_type_2040_fluvial.fraction_backyard_in_flood_prone_area[2]]
stats_2040_2 = [stats_per_housing_type_2040_fluvial.fraction_formal_in_flood_prone_area[3], stats_per_housing_type_2040_fluvial.fraction_subsidized_in_flood_prone_area[3], stats_per_housing_type_2040_fluvial.fraction_informal_in_flood_prone_area[3], stats_per_housing_type_2040_fluvial.fraction_backyard_in_flood_prone_area[3]]
stats_2040_3 = [stats_per_housing_type_2040_fluvial.fraction_formal_in_flood_prone_area[5], stats_per_housing_type_2040_fluvial.fraction_subsidized_in_flood_prone_area[5], stats_per_housing_type_2040_fluvial.fraction_informal_in_flood_prone_area[5], stats_per_housing_type_2040_fluvial.fraction_backyard_in_flood_prone_area[5]]
colors = ['#FF9999', '#00BFFF','#C1FFC1','#CAE1FF','#FFDEAD']
r = np.arange(4)
barWidth = 0.25
plt.figure(figsize=(10,7))
plt.bar(r, stats_2011_1, color=colors[0], edgecolor='white', width=barWidth, label="20 years")
plt.bar(r, np.array(stats_2011_2) - np.array(stats_2011_1), bottom=np.array(stats_2011_1), color=colors[1], edgecolor='white', width=barWidth, label='50 years')
plt.bar(r, np.array(stats_2011_3) - (np.array(stats_2011_2)), bottom=(np.array(stats_2011_2)), color=colors[2], edgecolor='white', width=barWidth, label='100 years')
plt.bar(r + 0.25, np.array(stats_2040_1), color=colors[0], edgecolor='white', width=barWidth)
plt.bar(r + 0.25, np.array(stats_2040_2) - np.array(stats_2040_1), bottom=np.array(stats_2040_1), color=colors[1], edgecolor='white', width=barWidth)
plt.bar(r + 0.25, np.array(stats_2040_3) - np.array(stats_2040_2), bottom=np.array(stats_2040_2), color=colors[2], edgecolor='white', width=barWidth)
plt.legend(loc = 'upper right')
plt.xticks(r, label)
plt.ylim(0, 75000)
plt.text(r[0] - 0.1, stats_per_housing_type_2011_fluvial.fraction_formal_in_flood_prone_area[5] + 0.005, "2011")
plt.text(r[1] - 0.1, stats_per_housing_type_2011_fluvial.fraction_subsidized_in_flood_prone_area[5] + 0.005, "2011") 
plt.text(r[2] - 0.1, stats_per_housing_type_2011_fluvial.fraction_informal_in_flood_prone_area[5] + 0.005, "2011") 
plt.text(r[3] - 0.1, stats_per_housing_type_2011_fluvial.fraction_backyard_in_flood_prone_area[5] + 0.005, "2011")
plt.text(r[0] + 0.15, stats_per_housing_type_2040_fluvial.fraction_formal_in_flood_prone_area[5] + 0.005, '2040')
plt.text(r[1] + 0.15, stats_per_housing_type_2040_fluvial.fraction_subsidized_in_flood_prone_area[5] + 0.005, '2040') 
plt.text(r[2] + 0.15, stats_per_housing_type_2040_fluvial.fraction_informal_in_flood_prone_area[5] + 0.005, '2040') 
plt.text(r[3] + 0.15, stats_per_housing_type_2040_fluvial.fraction_backyard_in_flood_prone_area[5] + 0.005, '2040') 
plt.tick_params(labelbottom=True)
plt.ylabel("Dwellings in flood-prone areas")
plt.show()

##v2 HARRIS
stats_per_housing_type_2011_fluvial["tot"] = stats_per_housing_type_2011_fluvial.fraction_formal_in_flood_prone_area + np.array(stats_per_housing_type_2011_fluvial.fraction_backyard_in_flood_prone_area) + stats_per_housing_type_2011_fluvial.fraction_subsidized_in_flood_prone_area + stats_per_housing_type_2011_fluvial.fraction_informal_in_flood_prone_area
stats_per_housing_type_2040_fluvial["tot"] = stats_per_housing_type_2040_fluvial.fraction_formal_in_flood_prone_area + np.array(stats_per_housing_type_2040_fluvial.fraction_backyard_in_flood_prone_area) + stats_per_housing_type_2040_fluvial.fraction_subsidized_in_flood_prone_area + stats_per_housing_type_2040_fluvial.fraction_informal_in_flood_prone_area

plt.figure(figsize=(10,7))
barWidth = 0.25
vec_2011_formal = stats_per_housing_type_2011_fluvial.fraction_formal_in_flood_prone_area / stats_per_housing_type_2011_fluvial["tot"]
vec_2011_formal = [np.nansum(simulation_households_housing_type[0, :, :], 1)[0] / sum(np.nansum(simulation_households_housing_type[0, :, :], 1)), vec_2011_formal[2], vec_2011_formal[3],vec_2011_formal[5]]
vec_2011_subsidized = stats_per_housing_type_2011_fluvial.fraction_subsidized_in_flood_prone_area / stats_per_housing_type_2011_fluvial["tot"]
vec_2011_subsidized = [np.nansum(simulation_households_housing_type[0, :, :], 1)[3]/ sum(np.nansum(simulation_households_housing_type[0, :, :], 1)), vec_2011_subsidized[2], vec_2011_subsidized[3],vec_2011_subsidized[5]]
vec_2011_informal = stats_per_housing_type_2011_fluvial.fraction_informal_in_flood_prone_area / stats_per_housing_type_2011_fluvial["tot"]
vec_2011_informal = [np.nansum(simulation_households_housing_type[0, :, :], 1)[2]/ sum(np.nansum(simulation_households_housing_type[0, :, :], 1)),vec_2011_informal[2], vec_2011_informal[3],vec_2011_informal[5]]
vec_2011_backyard = stats_per_housing_type_2011_fluvial.fraction_backyard_in_flood_prone_area / stats_per_housing_type_2011_fluvial["tot"]
vec_2011_backyard = [np.nansum(simulation_households_housing_type[0, :, :], 1)[1]/ sum(np.nansum(simulation_households_housing_type[0, :, :], 1)),vec_2011_backyard[2], vec_2011_backyard[3],vec_2011_backyard[5]]
vec_2040_formal = stats_per_housing_type_2040_fluvial.fraction_formal_in_flood_prone_area / stats_per_housing_type_2040_fluvial["tot"]
vec_2040_formal = [np.nansum(simulation_households_housing_type[28, :, :], 1)[0]/ sum(np.nansum(simulation_households_housing_type[28, :, :], 1)), vec_2040_formal[2], vec_2040_formal[3],vec_2040_formal[5]]
vec_2040_subsidized = stats_per_housing_type_2040_fluvial.fraction_subsidized_in_flood_prone_area / stats_per_housing_type_2040_fluvial["tot"]
vec_2040_subsidized = [np.nansum(simulation_households_housing_type[28, :, :], 1)[3]/ sum(np.nansum(simulation_households_housing_type[28, :, :], 1)), vec_2040_subsidized[2], vec_2040_subsidized[3],vec_2040_subsidized[5]]
vec_2040_informal = stats_per_housing_type_2040_fluvial.fraction_informal_in_flood_prone_area / stats_per_housing_type_2040_fluvial["tot"]
vec_2040_informal = [np.nansum(simulation_households_housing_type[28, :, :], 1)[2]/ sum(np.nansum(simulation_households_housing_type[28, :, :], 1)),vec_2040_informal[2], vec_2040_informal[3],vec_2040_informal[5]]
vec_2040_backyard = stats_per_housing_type_2040_fluvial.fraction_backyard_in_flood_prone_area / stats_per_housing_type_2040_fluvial["tot"]
vec_2040_backyard = [np.nansum(simulation_households_housing_type[28, :, :], 1)[1]/ sum(np.nansum(simulation_households_housing_type[28, :, :], 1)),vec_2040_backyard[2], vec_2040_backyard[3],vec_2040_backyard[5]]
plt.ylim(0, 1.3)
plt.ylabel("Fraction of dwellings of each housing type")
label = ["Over the city", "In 20-year \n return period \n flood zones","In 50-year \n return period \n flood zones","In 100-year \n return period \n flood zones"]
plt.bar(np.arange(4), vec_2011_formal, color=colors[0], edgecolor='white', width=barWidth, label="Formal private")
plt.bar(np.arange(4), vec_2011_subsidized, bottom=vec_2011_formal, color=colors[1], edgecolor='white', width=barWidth, label="Formal subsidized")
plt.bar(np.arange(4), vec_2011_informal, bottom=np.array(vec_2011_formal) + np.array(vec_2011_subsidized), color=colors[2], edgecolor='white', width=barWidth, label="Informal settlements")
plt.bar(np.arange(4), vec_2011_backyard, bottom=np.array(vec_2011_formal) + np.array(vec_2011_subsidized) + np.array(vec_2011_informal), color=colors[3], edgecolor='white', width=barWidth, label="Informal in backyards")
plt.bar(np.arange(4) + 0.25, vec_2040_formal, color=colors[0], edgecolor='white', width=barWidth)
plt.bar(np.arange(4) + 0.25, vec_2040_subsidized, bottom=vec_2040_formal, color=colors[1], edgecolor='white', width=barWidth)
plt.bar(np.arange(4) + 0.25, vec_2040_informal, bottom=np.array(vec_2040_formal) + np.array(vec_2040_subsidized), color=colors[2], edgecolor='white', width=barWidth)
plt.bar(np.arange(4) + 0.25, vec_2040_backyard, bottom=np.array(vec_2040_formal) + np.array(vec_2040_subsidized) + np.array(vec_2040_informal), color=colors[3], edgecolor='white', width=barWidth)
plt.legend(loc = 'upper left')
plt.xticks(np.arange(4), label)
plt.text(r[0] - 0.1, 1.005, "2011")
plt.text(r[1] - 0.1, 1.005, "2011") 
plt.text(r[2] - 0.1, 1.005, "2011") 
plt.text(r[3] - 0.1, 1.005, "2011")
plt.text(r[0] + 0.15, 1.005, '2040')
plt.text(r[1] + 0.15, 1.005, '2040') 
plt.text(r[2] + 0.15, 1.005, '2040') 
plt.text(r[3] + 0.15, 1.005, '2040') 

##PLUVIAL
stats_per_housing_type_2011_pluvial = compute_stats_per_housing_type(pluvial_floods, path_data, simulation_households_housing_type[0, 0, :], simulation_households_housing_type[0, 3, :], simulation_households_housing_type[0, 2, :], simulation_households_housing_type[0, 1, :], options, param, 0.01)
stats_per_housing_type_2040_pluvial = compute_stats_per_housing_type(pluvial_floods, path_data, simulation_households_housing_type[28, 0, :], simulation_households_housing_type[28, 3, :], simulation_households_housing_type[28, 2, :], simulation_households_housing_type[28, 1, :], options, param, 0.01)

label = ["Formal private", "Formal subsidized", "Informal \n settlements", "Informal \n in backyards"]
stats_2011_1 = [stats_per_housing_type_2011_pluvial.fraction_formal_in_flood_prone_area[2], stats_per_housing_type_2011_pluvial.fraction_subsidized_in_flood_prone_area[2], stats_per_housing_type_2011_pluvial.fraction_informal_in_flood_prone_area[2], stats_per_housing_type_2011_pluvial.fraction_backyard_in_flood_prone_area[2]]
stats_2011_2 = [stats_per_housing_type_2011_pluvial.fraction_formal_in_flood_prone_area[3], stats_per_housing_type_2011_pluvial.fraction_subsidized_in_flood_prone_area[3], stats_per_housing_type_2011_pluvial.fraction_informal_in_flood_prone_area[3], stats_per_housing_type_2011_pluvial.fraction_backyard_in_flood_prone_area[3]]
stats_2011_3 = [stats_per_housing_type_2011_pluvial.fraction_formal_in_flood_prone_area[5], stats_per_housing_type_2011_pluvial.fraction_subsidized_in_flood_prone_area[5], stats_per_housing_type_2011_pluvial.fraction_informal_in_flood_prone_area[5], stats_per_housing_type_2011_pluvial.fraction_backyard_in_flood_prone_area[5]]
stats_2040_1 = [stats_per_housing_type_2040_pluvial.fraction_formal_in_flood_prone_area[2], stats_per_housing_type_2040_pluvial.fraction_subsidized_in_flood_prone_area[2], stats_per_housing_type_2040_pluvial.fraction_informal_in_flood_prone_area[2], stats_per_housing_type_2040_pluvial.fraction_backyard_in_flood_prone_area[2]]
stats_2040_2 = [stats_per_housing_type_2040_pluvial.fraction_formal_in_flood_prone_area[3], stats_per_housing_type_2040_pluvial.fraction_subsidized_in_flood_prone_area[3], stats_per_housing_type_2040_pluvial.fraction_informal_in_flood_prone_area[3], stats_per_housing_type_2040_pluvial.fraction_backyard_in_flood_prone_area[3]]
stats_2040_3 = [stats_per_housing_type_2040_pluvial.fraction_formal_in_flood_prone_area[5], stats_per_housing_type_2040_pluvial.fraction_subsidized_in_flood_prone_area[5], stats_per_housing_type_2040_pluvial.fraction_informal_in_flood_prone_area[5], stats_per_housing_type_2040_pluvial.fraction_backyard_in_flood_prone_area[5]]
colors = ['#FF9999', '#00BFFF','#C1FFC1','#CAE1FF','#FFDEAD']
r = np.arange(4)
barWidth = 0.25
plt.figure(figsize=(10,7))
plt.bar(r, stats_2011_1, color=colors[0], edgecolor='white', width=barWidth, label="20 years")
plt.bar(r, np.array(stats_2011_2) - np.array(stats_2011_1), bottom=np.array(stats_2011_1), color=colors[1], edgecolor='white', width=barWidth, label='50 years')
plt.bar(r, np.array(stats_2011_3) - (np.array(stats_2011_2)), bottom=(np.array(stats_2011_2)), color=colors[2], edgecolor='white', width=barWidth, label='100 years')
plt.bar(r + 0.25, np.array(stats_2040_1), color=colors[0], edgecolor='white', width=barWidth)
plt.bar(r + 0.25, np.array(stats_2040_2) - np.array(stats_2040_1), bottom=np.array(stats_2040_1), color=colors[1], edgecolor='white', width=barWidth)
plt.bar(r + 0.25, np.array(stats_2040_3) - np.array(stats_2040_2), bottom=np.array(stats_2040_2), color=colors[2], edgecolor='white', width=barWidth)
plt.legend(loc = 'upper right')
plt.xticks(r, label)
plt.ylim(0, 290000)
plt.text(r[0] - 0.1, stats_per_housing_type_2011_pluvial.fraction_formal_in_flood_prone_area[5] + 0.005, "2011")
plt.text(r[1] - 0.1, stats_per_housing_type_2011_pluvial.fraction_subsidized_in_flood_prone_area[5] + 0.005, "2011") 
plt.text(r[2] - 0.1, stats_per_housing_type_2011_pluvial.fraction_informal_in_flood_prone_area[5] + 0.005, "2011") 
plt.text(r[3] - 0.1, stats_per_housing_type_2011_pluvial.fraction_backyard_in_flood_prone_area[5] + 0.005, "2011")
plt.text(r[0] + 0.15, stats_per_housing_type_2040_pluvial.fraction_formal_in_flood_prone_area[5] + 0.005, '2040')
plt.text(r[1] + 0.15, stats_per_housing_type_2040_pluvial.fraction_subsidized_in_flood_prone_area[5] + 0.005, '2040') 
plt.text(r[2] + 0.15, stats_per_housing_type_2040_pluvial.fraction_informal_in_flood_prone_area[5] + 0.005, '2040') 
plt.text(r[3] + 0.15, stats_per_housing_type_2040_pluvial.fraction_backyard_in_flood_prone_area[5] + 0.005, '2040') 
plt.tick_params(labelbottom=True)
plt.ylabel("Dwellings in flood-prone areas")
plt.show()

##v2 HARRIS
stats_per_housing_type_2011_pluvial["tot"] = stats_per_housing_type_2011_pluvial.fraction_formal_in_flood_prone_area + np.array(stats_per_housing_type_2011_pluvial.fraction_backyard_in_flood_prone_area) + stats_per_housing_type_2011_pluvial.fraction_subsidized_in_flood_prone_area + stats_per_housing_type_2011_pluvial.fraction_informal_in_flood_prone_area
stats_per_housing_type_2040_pluvial["tot"] = stats_per_housing_type_2040_pluvial.fraction_formal_in_flood_prone_area + np.array(stats_per_housing_type_2040_pluvial.fraction_backyard_in_flood_prone_area) + stats_per_housing_type_2040_pluvial.fraction_subsidized_in_flood_prone_area + stats_per_housing_type_2040_pluvial.fraction_informal_in_flood_prone_area

plt.figure(figsize=(10,7))
barWidth = 0.25
vec_2011_formal = stats_per_housing_type_2011_pluvial.fraction_formal_in_flood_prone_area / stats_per_housing_type_2011_pluvial["tot"]
vec_2011_formal = [np.nansum(simulation_households_housing_type[0, :, :], 1)[0] / sum(np.nansum(simulation_households_housing_type[0, :, :], 1)), vec_2011_formal[2], vec_2011_formal[3],vec_2011_formal[5]]
vec_2011_subsidized = stats_per_housing_type_2011_pluvial.fraction_subsidized_in_flood_prone_area / stats_per_housing_type_2011_pluvial["tot"]
vec_2011_subsidized = [np.nansum(simulation_households_housing_type[0, :, :], 1)[3]/ sum(np.nansum(simulation_households_housing_type[0, :, :], 1)), vec_2011_subsidized[2], vec_2011_subsidized[3],vec_2011_subsidized[5]]
vec_2011_informal = stats_per_housing_type_2011_pluvial.fraction_informal_in_flood_prone_area / stats_per_housing_type_2011_pluvial["tot"]
vec_2011_informal = [np.nansum(simulation_households_housing_type[0, :, :], 1)[2]/ sum(np.nansum(simulation_households_housing_type[0, :, :], 1)),vec_2011_informal[2], vec_2011_informal[3],vec_2011_informal[5]]
vec_2011_backyard = stats_per_housing_type_2011_pluvial.fraction_backyard_in_flood_prone_area / stats_per_housing_type_2011_pluvial["tot"]
vec_2011_backyard = [np.nansum(simulation_households_housing_type[0, :, :], 1)[1]/ sum(np.nansum(simulation_households_housing_type[0, :, :], 1)),vec_2011_backyard[2], vec_2011_backyard[3],vec_2011_backyard[5]]
vec_2040_formal = stats_per_housing_type_2040_pluvial.fraction_formal_in_flood_prone_area / stats_per_housing_type_2040_pluvial["tot"]
vec_2040_formal = [np.nansum(simulation_households_housing_type[28, :, :], 1)[0]/ sum(np.nansum(simulation_households_housing_type[28, :, :], 1)), vec_2040_formal[2], vec_2040_formal[3],vec_2040_formal[5]]
vec_2040_subsidized = stats_per_housing_type_2040_pluvial.fraction_subsidized_in_flood_prone_area / stats_per_housing_type_2040_pluvial["tot"]
vec_2040_subsidized = [np.nansum(simulation_households_housing_type[28, :, :], 1)[3]/ sum(np.nansum(simulation_households_housing_type[28, :, :], 1)), vec_2040_subsidized[2], vec_2040_subsidized[3],vec_2040_subsidized[5]]
vec_2040_informal = stats_per_housing_type_2040_pluvial.fraction_informal_in_flood_prone_area / stats_per_housing_type_2040_pluvial["tot"]
vec_2040_informal = [np.nansum(simulation_households_housing_type[28, :, :], 1)[2]/ sum(np.nansum(simulation_households_housing_type[28, :, :], 1)),vec_2040_informal[2], vec_2040_informal[3],vec_2040_informal[5]]
vec_2040_backyard = stats_per_housing_type_2040_pluvial.fraction_backyard_in_flood_prone_area / stats_per_housing_type_2040_pluvial["tot"]
vec_2040_backyard = [np.nansum(simulation_households_housing_type[28, :, :], 1)[1]/ sum(np.nansum(simulation_households_housing_type[28, :, :], 1)),vec_2040_backyard[2], vec_2040_backyard[3],vec_2040_backyard[5]]
plt.ylim(0, 1.3)
label = ["Over the city", "In 20-year \n return period \n flood zones","In 50-year \n return period \n flood zones","In 100-year \n return period \n flood zones"]
plt.bar(np.arange(4), vec_2011_formal, color=colors[0], edgecolor='white', width=barWidth, label="Formal private")
plt.bar(np.arange(4), vec_2011_subsidized, bottom=vec_2011_formal, color=colors[1], edgecolor='white', width=barWidth, label="Formal subsidized")
plt.bar(np.arange(4), vec_2011_informal, bottom=np.array(vec_2011_formal) + np.array(vec_2011_subsidized), color=colors[2], edgecolor='white', width=barWidth, label="Informal settlements")
plt.bar(np.arange(4), vec_2011_backyard, bottom=np.array(vec_2011_formal) + np.array(vec_2011_subsidized) + np.array(vec_2011_informal), color=colors[3], edgecolor='white', width=barWidth, label="Informal in backyards")
plt.bar(np.arange(4) + 0.25, vec_2040_formal, color=colors[0], edgecolor='white', width=barWidth)
plt.bar(np.arange(4) + 0.25, vec_2040_subsidized, bottom=vec_2040_formal, color=colors[1], edgecolor='white', width=barWidth)
plt.bar(np.arange(4) + 0.25, vec_2040_informal, bottom=np.array(vec_2040_formal) + np.array(vec_2040_subsidized), color=colors[2], edgecolor='white', width=barWidth)
plt.bar(np.arange(4) + 0.25, vec_2040_backyard, bottom=np.array(vec_2040_formal) + np.array(vec_2040_subsidized) + np.array(vec_2040_informal), color=colors[3], edgecolor='white', width=barWidth)
plt.legend(loc = 'upper left')
plt.ylabel("Fraction of dwellings of each housing type")
plt.xticks(np.arange(4), label)
plt.text(r[0] - 0.1, 1.005, "2011")
plt.text(r[1] - 0.1, 1.005, "2011") 
plt.text(r[2] - 0.1, 1.005, "2011") 
plt.text(r[3] - 0.1, 1.005, "2011")
plt.text(r[0] + 0.15, 1.005, '2040')
plt.text(r[1] + 0.15, 1.005, '2040') 
plt.text(r[2] + 0.15, 1.005, '2040') 
plt.text(r[3] + 0.15, 1.005, '2040') 

##  B. Evolution des dégâts par quartile de revenu

#spline_agricultural_rent, spline_interest_rate, spline_RDP, spline_population_income_distribution, spline_inflation, spline_income_distribution, spline_population, spline_interest_rate, spline_income, spline_minimum_housing_supply, spline_fuel = import_scenarios(income_2011, param, grid, path_scenarios) #we add required argument
formal_structure_cost_2011 = compute_formal_structure_cost_method2(simulation_rent[0, :, :], param, interpolate_interest_rate(spline_interest_rate, 0), import_coeff_land(spline_land_constraints, spline_land_backyard, spline_land_informal, spline_land_RDP, param, 0), simulation_households_housing_type[0, :, :], (spline_income(0) / param["income_year_reference"])**(-param["coeff_b"]) * param["coeff_A"])    
content_cost_2011 = compute_content_cost(simulation_households_center[0, :, :], np.load(path_folder+"precalculated_transport/incomeNetofCommuting_0.npy"), param, fraction_capital_destroyed, simulation_rent[0, :, :], simulation_dwelling_size[0, :, :], interpolate_interest_rate(spline_interest_rate, 0))
formal_structure_cost_2040 = compute_formal_structure_cost_method2(simulation_rent[28, :, :], param, interpolate_interest_rate(spline_interest_rate, 28), import_coeff_land(spline_land_constraints, spline_land_backyard, spline_land_informal, spline_land_RDP, param, 28), simulation_households_housing_type[28, :, :], (spline_income(28) / param["income_year_reference"])**(-param["coeff_b"]) * param["coeff_A"])    
content_cost_2040 = compute_content_cost(simulation_households_center[28, :, :], np.load(path_folder+"precalculated_transport/incomeNetofCommuting_29.npy"), param, fraction_capital_destroyed, simulation_rent[28, :, :], simulation_dwelling_size[28, :, :], interpolate_interest_rate(spline_interest_rate, 28))

item = 'FD_100yr'
path_data = path_folder+"/FATHOM/"
option = "percent" #"absolu"

df2011 = pd.DataFrame()
df2040 = pd.DataFrame()
type_flood = copy.deepcopy(item)
data_flood = np.squeeze(pd.read_excel(path_data + item + ".xlsx"))
    
formal_damages = structural_damages_type4a(data_flood['flood_depth'])
formal_damages[simulation_dwelling_size[0, 0, :] > param["threshold"]] = structural_damages_type4b(data_flood.flood_depth[simulation_dwelling_size[0, 0, :] > param["threshold"]])
subsidized_damages = structural_damages_type4a(data_flood['flood_depth'])
subsidized_damages[simulation_dwelling_size[0, 3, :] > param["threshold"]] = structural_damages_type4b(data_flood.flood_depth[simulation_dwelling_size[0, 3, :] > param["threshold"]])
        
df2011['formal_structure_damages'] = formal_structure_cost_2011 * formal_damages
df2011['subsidized_structure_damages'] = param["subsidized_structure_value_ref"] * subsidized_damages
df2011['informal_structure_damages'] = param["informal_structure_value_ref"] * structural_damages_type2(data_flood['flood_depth'])
df2011['backyard_structure_damages'] = ((16216 * (param["informal_structure_value_ref"] * structural_damages_type2(data_flood['flood_depth']))) + (74916 * (param["informal_structure_value_ref"] * structural_damages_type3a(data_flood['flood_depth'])))) / (74916 + 16216)
            
df2011['formal_content_damages'] =  content_cost_2011.formal * content_damages(data_flood['flood_depth'])
df2011['subsidized_content_damages'] = content_cost_2011.subsidized * content_damages(data_flood['flood_depth'])
df2011['informal_content_damages'] = content_cost_2011.informal * content_damages(data_flood['flood_depth'])
df2011['backyard_content_damages'] = content_cost_2011.backyard * content_damages(data_flood['flood_depth'])
    
formal_damages = structural_damages_type4a(data_flood['flood_depth'])
formal_damages[simulation_dwelling_size[28, 0, :] > param["threshold"]] = structural_damages_type4b(data_flood.flood_depth[simulation_dwelling_size[28, 0, :] > param["threshold"]])
subsidized_damages = structural_damages_type4a(data_flood['flood_depth'])
subsidized_damages[simulation_dwelling_size[28, 3, :] > param["threshold"]] = structural_damages_type4b(data_flood.flood_depth[simulation_dwelling_size[28, 3, :] > param["threshold"]])
   
df2040['formal_structure_damages'] = formal_structure_cost_2040 * formal_damages
df2040['subsidized_structure_damages'] = param["subsidized_structure_value_ref"] * (spline_inflation(28) / spline_inflation(0)) * subsidized_damages
df2040['informal_structure_damages'] = param["informal_structure_value_ref"] * (spline_inflation(28) / spline_inflation(0)) * structural_damages_type2(data_flood['flood_depth'])
df2040['backyard_structure_damages'] = ((16216 * (param["informal_structure_value_ref"] * (spline_inflation(28) / spline_inflation(0)) * structural_damages_type2(data_flood['flood_depth']))) + (74916 * (param["informal_structure_value_ref"] * (spline_inflation(28) / spline_inflation(0)) * structural_damages_type3a(data_flood['flood_depth'])))) / (74916 + 16216)
    
df2040['formal_content_damages'] =  content_cost_2040.formal * content_damages(data_flood['flood_depth'])
df2040['subsidized_content_damages'] = content_cost_2040.subsidized * content_damages(data_flood['flood_depth'])
df2040['informal_content_damages'] = content_cost_2040.informal * content_damages(data_flood['flood_depth'])
df2040['backyard_content_damages'] = content_cost_2040.backyard * content_damages(data_flood['flood_depth'])

df2011["formal_pop_flood_prone"] = simulation_households_housing_type[0, 0, :] * data_flood["prop_flood_prone"]
df2011["backyard_pop_flood_prone"] = simulation_households_housing_type[0, 1, :] * data_flood["prop_flood_prone"]
df2011["informal_pop_flood_prone"] = simulation_households_housing_type[0, 2, :] * data_flood["prop_flood_prone"]
df2011["subsidized_pop_flood_prone"] = simulation_households_housing_type[0, 3, :] * data_flood["prop_flood_prone"]
    
df2040["formal_pop_flood_prone"] = simulation_households_housing_type[28, 0, :] * data_flood["prop_flood_prone"]
df2040["backyard_pop_flood_prone"] = simulation_households_housing_type[28, 1, :] * data_flood["prop_flood_prone"]
df2040["informal_pop_flood_prone"] = simulation_households_housing_type[28, 2, :] * data_flood["prop_flood_prone"]
df2040["subsidized_pop_flood_prone"] = simulation_households_housing_type[28, 3, :] * data_flood["prop_flood_prone"]
    
df2011["formal_damages"] = df2011['formal_structure_damages'] + df2011['formal_content_damages']
df2011["informal_damages"] = df2011['informal_structure_damages'] + df2011['informal_content_damages']
df2011["subsidized_damages"] = df2011['subsidized_structure_damages'] + df2011['subsidized_content_damages']
df2011["backyard_damages"] = df2011['backyard_structure_damages'] + df2011['backyard_content_damages']
    
df2040["formal_damages"] = df2040['formal_structure_damages'] + df2040['formal_content_damages']    
df2040["informal_damages"] = df2040['informal_structure_damages'] + df2040['informal_content_damages']
df2040["subsidized_damages"] = df2040['subsidized_structure_damages'] + df2040['subsidized_content_damages']
df2040["backyard_damages"] = df2040['backyard_structure_damages'] + df2040['backyard_content_damages']
    
#import seaborn as sns
subset = df2011[(~np.isnan(df2011.formal_damages)) & (df2011.formal_pop_flood_prone > 0)]
sns.distplot(subset.formal_damages, hist = True, kde = False, hist_kws={'weights': subset.formal_pop_flood_prone})


######### 2. EAD

#spline_agricultural_rent, spline_interest_rate, spline_RDP, spline_population_income_distribution, spline_inflation, spline_income_distribution, spline_population, spline_interest_rate, spline_income, spline_minimum_housing_supply, spline_fuel = import_scenarios(income_2011, param, grid, path_scenarios) #we add required parameter

formal_structure_cost_2011 = compute_formal_structure_cost_method2(simulation_rent[0, :, :], param, interpolate_interest_rate(spline_interest_rate, 0), import_coeff_land(spline_land_constraints, spline_land_backyard, spline_land_informal, spline_land_RDP, param, 0), simulation_households_housing_type[0, :, :], (spline_income(0) / param["income_year_reference"])**(-param["coeff_b"]) * param["coeff_A"])    
content_cost_2011 = compute_content_cost(simulation_households_center[0, :, :], np.load(path_folder+"precalculated_transport/incomeNetofCommuting_0.npy"), param, fraction_capital_destroyed, simulation_rent[0, :, :], simulation_dwelling_size[0, :, :], interpolate_interest_rate(spline_interest_rate, 0))
damages_fluvial_2011 = compute_damages(fluvial_floods, path_data, param, content_cost_2011,
                    simulation_households_housing_type[0, 0, :], simulation_households_housing_type[0, 3, :], simulation_households_housing_type[0, 2, :], simulation_households_housing_type[0, 1, :], simulation_dwelling_size[0, :, :],
                    formal_structure_cost_2011, content_damages, structural_damages_type4b, structural_damages_type4a, structural_damages_type2, structural_damages_type3a, options, spline_inflation, 0)

formal_structure_cost_2040 = compute_formal_structure_cost_method2(simulation_rent[28, :, :], param, interpolate_interest_rate(spline_interest_rate, 28), import_coeff_land(spline_land_constraints, spline_land_backyard, spline_land_informal, spline_land_RDP, param, 28), simulation_households_housing_type[28, :, :], (spline_income(28) / param["income_year_reference"])**(-param["coeff_b"]) * param["coeff_A"])    
content_cost_2040 = compute_content_cost(simulation_households_center[28, :, :], np.load(path_folder+"precalculated_transport/incomeNetofCommuting_29.npy"), param, fraction_capital_destroyed, simulation_rent[28, :, :], simulation_dwelling_size[28, :, :], interpolate_interest_rate(spline_interest_rate, 28))
damages_fluvial_2040 = compute_damages(fluvial_floods, path_data, param, content_cost_2040,
                    simulation_households_housing_type[28, 0, :], simulation_households_housing_type[28, 3, :], simulation_households_housing_type[28, 2, :], simulation_households_housing_type[28, 1, :], simulation_dwelling_size[28, :, :],
                    formal_structure_cost_2040, content_damages, structural_damages_type4b, structural_damages_type4a, structural_damages_type2, structural_damages_type3a, options, spline_inflation, 28)

damages_fluvial_2011.backyard_damages = damages_fluvial_2011.backyard_content_damages + damages_fluvial_2011.backyard_structure_damages
damages_fluvial_2011.informal_damages = damages_fluvial_2011.informal_content_damages + damages_fluvial_2011.informal_structure_damages
damages_fluvial_2011.subsidized_damages = damages_fluvial_2011.subsidized_content_damages + damages_fluvial_2011.subsidized_structure_damages
damages_fluvial_2011.formal_damages = damages_fluvial_2011.formal_content_damages + damages_fluvial_2011.formal_structure_damages

damages_fluvial_2040.backyard_damages = damages_fluvial_2040.backyard_content_damages + damages_fluvial_2040.backyard_structure_damages
damages_fluvial_2040.informal_damages = damages_fluvial_2040.informal_content_damages + damages_fluvial_2040.informal_structure_damages
damages_fluvial_2040.subsidized_damages = damages_fluvial_2040.subsidized_content_damages + damages_fluvial_2040.subsidized_structure_damages
damages_fluvial_2040.formal_damages = damages_fluvial_2040.formal_content_damages + damages_fluvial_2040.formal_structure_damages


formal_structure_cost_2011 = compute_formal_structure_cost_method2(simulation_rent[0, :, :], param, interpolate_interest_rate(spline_interest_rate, 0), import_coeff_land(spline_land_constraints, spline_land_backyard, spline_land_informal, spline_land_RDP, param, 0), simulation_households_housing_type[0, :, :], (spline_income(0) / param["income_year_reference"])**(-param["coeff_b"]) * param["coeff_A"])    
content_cost_2011 = compute_content_cost(simulation_households_center[0, :, :], np.load(path_folder+"precalculated_transport/incomeNetofCommuting_0.npy"), param, fraction_capital_destroyed, simulation_rent[0, :, :], simulation_dwelling_size[0, :, :], interpolate_interest_rate(spline_interest_rate, 0))
damages_pluvial_2011 = compute_damages(pluvial_floods, path_data, param, content_cost_2011,
                    simulation_households_housing_type[0, 0, :], simulation_households_housing_type[0, 3, :], simulation_households_housing_type[0, 2, :], simulation_households_housing_type[0, 1, :], simulation_dwelling_size[0, :, :],
                    formal_structure_cost_2011, content_damages, structural_damages_type4b, structural_damages_type4a, structural_damages_type2, structural_damages_type3a, options, spline_inflation, 0)

formal_structure_cost_2040 = compute_formal_structure_cost_method2(simulation_rent[28, :, :], param, interpolate_interest_rate(spline_interest_rate, 28), import_coeff_land(spline_land_constraints, spline_land_backyard, spline_land_informal, spline_land_RDP, param, 28), simulation_households_housing_type[28, :, :], (spline_income(28) / param["income_year_reference"])**(-param["coeff_b"]) * param["coeff_A"])    
content_cost_2040 = compute_content_cost(simulation_households_center[28, :, :], np.load(path_folder+"precalculated_transport/incomeNetofCommuting_29.npy"), param, fraction_capital_destroyed, simulation_rent[28, :, :], simulation_dwelling_size[28, :, :], interpolate_interest_rate(spline_interest_rate, 28))
damages_pluvial_2040 = compute_damages(pluvial_floods, path_data, param, content_cost_2040,
                    simulation_households_housing_type[28, 0, :], simulation_households_housing_type[28, 3, :], simulation_households_housing_type[28, 2, :], simulation_households_housing_type[28, 1, :], simulation_dwelling_size[28, :, :],
                    formal_structure_cost_2040, content_damages, structural_damages_type4b, structural_damages_type4a, structural_damages_type2, structural_damages_type3a, options, spline_inflation, 28)

damages_pluvial_2011.backyard_damages = damages_pluvial_2011.backyard_content_damages + damages_pluvial_2011.backyard_structure_damages
damages_pluvial_2011.informal_damages = damages_pluvial_2011.informal_content_damages + damages_pluvial_2011.informal_structure_damages
damages_pluvial_2011.subsidized_damages = damages_pluvial_2011.subsidized_content_damages + damages_pluvial_2011.subsidized_structure_damages
damages_pluvial_2011.formal_damages = damages_pluvial_2011.formal_content_damages + damages_pluvial_2011.formal_structure_damages

damages_pluvial_2040.backyard_damages = damages_pluvial_2040.backyard_content_damages + damages_pluvial_2040.backyard_structure_damages
damages_pluvial_2040.informal_damages = damages_pluvial_2040.informal_content_damages + damages_pluvial_2040.informal_structure_damages
damages_pluvial_2040.subsidized_damages = damages_pluvial_2040.subsidized_content_damages + damages_pluvial_2040.subsidized_structure_damages
damages_pluvial_2040.formal_damages = damages_pluvial_2040.formal_content_damages + damages_pluvial_2040.formal_structure_damages

damages_pluvial_2011.formal_damages[0:3] = 0
damages_pluvial_2040.formal_damages[0:3] = 0
damages_pluvial_2011.backyard_damages[0:2] = 0
damages_pluvial_2040.backyard_damages[0:2] = 0
damages_pluvial_2011.subsidized_damages[0:2] = 0
damages_pluvial_2040.subsidized_damages[0:2] = 0

damages_2011 = pd.DataFrame()
damages_2011["backyard_damages"] = damages_pluvial_2011.backyard_damages + damages_fluvial_2011.backyard_damages
damages_2011["informal_damages"] = damages_pluvial_2011.informal_damages + damages_fluvial_2011.informal_damages
damages_2011["subsidized_damages"] = damages_pluvial_2011.subsidized_damages + damages_fluvial_2011.subsidized_damages
damages_2011["formal_damages"] = damages_pluvial_2011.formal_damages + damages_fluvial_2011.formal_damages


damages_2040 = pd.DataFrame()
damages_2040["backyard_damages"] = damages_pluvial_2040.backyard_damages + damages_fluvial_2040.backyard_damages
damages_2040["informal_damages"] = damages_pluvial_2040.informal_damages + damages_fluvial_2040.informal_damages
damages_2040["subsidized_damages"] = damages_pluvial_2040.subsidized_damages + damages_fluvial_2040.subsidized_damages
damages_2040["formal_damages"] = damages_pluvial_2040.formal_damages + damages_fluvial_2040.formal_damages


inflation = spline_inflation(28) / spline_inflation(0)
label = ["2011", "2040", "2040 (deflated)"]
stats_2011_formal = [annualize_damages(damages_2011.formal_damages),annualize_damages(damages_2040.formal_damages),annualize_damages(damages_2040.formal_damages)/inflation]
stats_2011_subsidized = [annualize_damages(damages_2011.subsidized_damages),annualize_damages(damages_2040.subsidized_damages),annualize_damages(damages_2040.subsidized_damages)/inflation]
stats_2011_informal = [annualize_damages(damages_2011.informal_damages), annualize_damages(damages_2040.informal_damages), annualize_damages(damages_2040.informal_damages)/inflation]
stats_2011_backyard = [annualize_damages(damages_2011.backyard_damages),annualize_damages(damages_2040.backyard_damages),annualize_damages(damages_2040.backyard_damages)/inflation]
colors = ['#FF9999', '#00BFFF','#C1FFC1','#CAE1FF','#FFDEAD']
r = np.arange(3)
barWidth = 0.5
plt.figure(figsize=(10,7))
plt.bar(r, stats_2011_formal, color=colors[0], edgecolor='white', width=barWidth, label="formal")
plt.bar(r, np.array(stats_2011_subsidized), bottom=np.array(stats_2011_formal), color=colors[1], edgecolor='white', width=barWidth, label='subsidized')
plt.bar(r, np.array(stats_2011_informal), bottom=(np.array(stats_2011_subsidized) + np.array(stats_2011_formal)), color=colors[2], edgecolor='white', width=barWidth, label='informal')
plt.bar(r, np.array(stats_2011_backyard), bottom=(np.array(stats_2011_informal) + np.array(stats_2011_subsidized) + np.array(stats_2011_formal)), color=colors[3], edgecolor='white', width=barWidth, label='backyard')
plt.legend()
plt.xticks(r, label)
plt.ylim(0, 600000000)
plt.tick_params(labelbottom=True)
plt.ylabel("Estimated annual damages (R)")
plt.show()

######### 2. Data for maps

pluvial_100yr = pd.read_excel(path_folder+"FATHOM/" + "P_100yr" + ".xlsx")
pluvial_100yr["flood_prone"] = ((pluvial_100yr.flood_depth > 0.05) & (pluvial_100yr.prop_flood_prone > 0.8)) #6695
#30% of the grid cell prone to floods of at least 1 cm every 100 yr
pop_2011 = np.nansum(simulation_households_housing_type[0, :, :], 0)
pop_2040 = np.nansum(simulation_households_housing_type[28, :, :], 0)
red_2011 = ((pop_2011 > 10) & (pluvial_100yr["flood_prone"] == 1))
grey_2011 = ((pop_2011 > 10) & (pluvial_100yr["flood_prone"] == 0))
blue_2011 = ((pop_2011 < 10) & (pluvial_100yr["flood_prone"] == 1))
red_2040 = ((pop_2040 > 10) & (pluvial_100yr["flood_prone"] == 1))
grey_2040 = ((pop_2040 > 10) & (pluvial_100yr["flood_prone"] == 0))
blue_2040 = ((pop_2040 < 10) & (pluvial_100yr["flood_prone"] == 1))
df = pd.DataFrame()
df["pop_2011"] = pop_2011
df["pop_2040"] = pop_2040
df["fraction_capital_destroyed"] = fraction_capital_destroyed["structure_formal_2"]
df.to_excel(path_outputs + name + "/map_data.xlsx")