# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 17:00:06 2020

@author: Charlotte Liotta
"""

import scipy
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import copy

from inputs.data import *

def compute_stats_per_housing_type(floods, path_data, nb_households_formal, nb_households_subsidized, nb_households_informal, nb_households_backyard, options, param, threshold):
    stats_per_housing_type = pd.DataFrame(columns = ['flood',
                                                     'fraction_formal_in_flood_prone_area', 'fraction_subsidized_in_flood_prone_area', 'fraction_informal_in_flood_prone_area', 'fraction_backyard_in_flood_prone_area',
                                                     'flood_depth_formal', 'flood_depth_subsidized', 'flood_depth_informal', 'flood_depth_backyard'])
    for flood in floods:
        type_flood = copy.deepcopy(flood)
        flood = np.squeeze(pd.read_excel(path_data + flood + ".xlsx"))
        
        flood.prop_flood_prone[flood.flood_depth < threshold] = 0
        flood.flood_depth[flood.flood_depth < threshold] = 0
        #if (type_flood == 'FD_20yr') & (options["WBUS2"] == 1):
        #   WBUS2_20yr = pd.read_excel("C:/Users/Charlotte Liotta/Desktop/cape_town/2. Data/Flood plains - from Claus/WBUS2_20yr.xlsx")
        #    flood.prop_flood_prone = np.fmax(flood.prop_flood_prone, WBUS2_20yr.prop_flood_prone)
        #    flood.flood_depth = np.maximum(flood.prop_flood_prone, param["depth_WBUS2_20yr"])
        #if (type_flood == 'FD_50yr') & (options["WBUS2"] == 1):
        #    WBUS2_50yr = pd.read_excel("C:/Users/Charlotte Liotta/Desktop/cape_town/2. Data/Flood plains - from Claus/WBUS2_50yr.xlsx")
        #    flood.prop_flood_prone = np.fmax(flood.prop_flood_prone, WBUS2_50yr.prop_flood_prone)
        #    flood.flood_depth = np.maximum(flood.prop_flood_prone, param["depth_WBUS2_50yr"])
        #if (type_flood == 'FD_100yr') & (options["WBUS2"] == 1):
        #    WBUS2_100yr = pd.read_excel("C:/Users/Charlotte Liotta/Desktop/cape_town/2. Data/Flood plains - from Claus/WBUS2_100yr.xlsx")
        #    flood.prop_flood_prone = np.fmax(flood.prop_flood_prone, WBUS2_100yr.prop_flood_prone)
        #    flood.flood_depth = np.maximum(flood.prop_flood_prone, param["depth_WBUS2_100yr"])
        print(type_flood)
        if ((type_flood == 'P_5yr') | (type_flood == 'P_10yr')):
            stats_per_housing_type = stats_per_housing_type.append({'flood': type_flood, 
                                                                    #'fraction_formal_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_formal) / sum(nb_households_formal), 
                                                                    #'fraction_subsidized_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_subsidized) / sum(nb_households_subsidized),
                                                                    #'fraction_informal_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_informal) / sum(nb_households_informal), 
                                                                    #'fraction_backyard_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_backyard) / sum(nb_households_backyard),
                                                                    'fraction_formal_in_flood_prone_area': 0, 
                                                                    'fraction_subsidized_in_flood_prone_area': 0,
                                                                    'fraction_informal_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_informal), 
                                                                    'fraction_backyard_in_flood_prone_area': 0,                                                              
                                                                    'flood_depth_formal': 0,
                                                                    'flood_depth_subsidized': 0,
                                                                    'flood_depth_informal': sum(flood['flood_depth'] * (flood['prop_flood_prone'] * nb_households_informal)  / sum(flood['prop_flood_prone'] * nb_households_informal)),
                                                                    'flood_depth_backyard': 0}, ignore_index = True)  
        elif (type_flood == 'P_20yr'):
            stats_per_housing_type = stats_per_housing_type.append({'flood': type_flood, 
                                                                    #'fraction_formal_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_formal) / sum(nb_households_formal), 
                                                                    #'fraction_subsidized_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_subsidized) / sum(nb_households_subsidized),
                                                                    #'fraction_informal_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_informal) / sum(nb_households_informal), 
                                                                    #'fraction_backyard_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_backyard) / sum(nb_households_backyard),
                                                                    'fraction_formal_in_flood_prone_area': 0, 
                                                                    'fraction_subsidized_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_subsidized),
                                                                    'fraction_informal_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_informal), 
                                                                    'fraction_backyard_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_backyard),                                                              
                                                                    'flood_depth_formal': 0,
                                                                    'flood_depth_subsidized': sum(flood['flood_depth'] * (flood['prop_flood_prone'] * nb_households_subsidized)  / sum(flood['prop_flood_prone'] * nb_households_subsidized)),
                                                                    'flood_depth_informal': sum(flood['flood_depth'] * (flood['prop_flood_prone'] * nb_households_informal)  / sum(flood['prop_flood_prone'] * nb_households_informal)),
                                                                    'flood_depth_backyard': sum(flood['flood_depth'] * (flood['prop_flood_prone'] * nb_households_backyard)  / sum(flood['prop_flood_prone'] * nb_households_backyard))}, ignore_index = True)  
        else:
            stats_per_housing_type = stats_per_housing_type.append({'flood': type_flood, 
                                                                    #'fraction_formal_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_formal) / sum(nb_households_formal), 
                                                                    #'fraction_subsidized_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_subsidized) / sum(nb_households_subsidized),
                                                                    #'fraction_informal_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_informal) / sum(nb_households_informal), 
                                                                    #'fraction_backyard_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_backyard) / sum(nb_households_backyard),
                                                                    'fraction_formal_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_formal), 
                                                                    'fraction_subsidized_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_subsidized),
                                                                    'fraction_informal_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_informal), 
                                                                    'fraction_backyard_in_flood_prone_area': np.sum(flood['prop_flood_prone'] * nb_households_backyard),                                                              
                                                                    'flood_depth_formal': sum(flood['flood_depth'] * (flood['prop_flood_prone'] * nb_households_formal)  / sum(flood['prop_flood_prone'] * nb_households_formal)),
                                                                    'flood_depth_subsidized': sum(flood['flood_depth'] * (flood['prop_flood_prone'] * nb_households_subsidized)  / sum(flood['prop_flood_prone'] * nb_households_subsidized)),
                                                                    'flood_depth_informal': sum(flood['flood_depth'] * (flood['prop_flood_prone'] * nb_households_informal)  / sum(flood['prop_flood_prone'] * nb_households_informal)),
                                                                    'flood_depth_backyard': sum(flood['flood_depth'] * (flood['prop_flood_prone'] * nb_households_backyard)  / sum(flood['prop_flood_prone'] * nb_households_backyard))}, ignore_index = True)   
        

            
    return stats_per_housing_type

def compute_damages(floods, path_data, param, content_cost,
                    nb_households_formal, nb_households_subsidized, nb_households_informal, nb_households_backyard, dwelling_size,
                    formal_structure_cost, content_damages, structural_damages_type4b, structural_damages_type4a, structural_damages_type2, structural_damages_type3a, options, spline_inflation, year_temp):
    
    #floods = ['FD_5yr', 'FD_10yr', 'FD_20yr', 'FD_50yr', 'FD_75yr', 'FD_100yr', 'FD_200yr', 'FD_250yr', 'FD_500yr', 'FD_1000yr']
    #path_data = "C:/Users/Charlotte Liotta/Desktop/cape_town/2. Data/FATHOM/"
    
    damages = pd.DataFrame(columns = ['flood',
                                      'formal_structure_damages',
                                      'subsidized_structure_damages',
                                      'informal_structure_damages',
                                      'backyard_structure_damages',
                                      'formal_content_damages',
                                      'subsidized_content_damages',
                                      'informal_content_damages',
                                      'backyard_content_damages'])
    
    for item in floods:
        type_flood = copy.deepcopy(item)
        data_flood = np.squeeze(pd.read_excel(path_data + item + ".xlsx"))
        
        #if (type_flood == 'FD_20yr') & (options["WBUS2"] == 1):
        #    WBUS2_20yr = pd.read_excel("C:/Users/Charlotte Liotta/Desktop/cape_town/2. Data/Flood plains - from Claus/WBUS2_20yr.xlsx")
        #    data_flood.prop_flood_prone = np.fmax(data_flood.prop_flood_prone, WBUS2_20yr.prop_flood_prone)
        #    data_flood.flood_depth = np.maximum(data_flood.prop_flood_prone, param["depth_WBUS2_20yr"])
        #if (type_flood == 'FD_50yr') & (options["WBUS2"] == 1):
        #    WBUS2_50yr = pd.read_excel("C:/Users/Charlotte Liotta/Desktop/cape_town/2. Data/Flood plains - from Claus/WBUS2_50yr.xlsx")
        #    data_flood.prop_flood_prone = np.fmax(data_flood.prop_flood_prone, WBUS2_50yr.prop_flood_prone)
        #    data_flood.flood_depth = np.maximum(data_flood.prop_flood_prone, param["depth_WBUS2_50yr"])
        #if (type_flood == 'FD_100yr') & (options["WBUS2"] == 1):
        #    WBUS2_100yr = pd.read_excel("C:/Users/Charlotte Liotta/Desktop/cape_town/2. Data/Flood plains - from Claus/WBUS2_100yr.xlsx")
        #    data_flood.prop_flood_prone = np.fmax(data_flood.prop_flood_prone, WBUS2_100yr.prop_flood_prone)
        #    data_flood.flood_depth = np.maximum(data_flood.prop_flood_prone, param["depth_WBUS2_100yr"])
        
        formal_damages = structural_damages_type4a(data_flood['flood_depth'])
        formal_damages[dwelling_size[0, :] > param["threshold"]] = structural_damages_type4b(data_flood.flood_depth[dwelling_size[0, :] > param["threshold"]])
        subsidized_damages = structural_damages_type4a(data_flood['flood_depth'])
        subsidized_damages[dwelling_size[3, :] > param["threshold"]] = structural_damages_type4b(data_flood.flood_depth[dwelling_size[3, :] > param["threshold"]])
        
        formal_structure_damages = np.nansum(nb_households_formal * data_flood["prop_flood_prone"] * formal_structure_cost * formal_damages)
        subsidized_structure_damages = np.nansum(nb_households_subsidized * data_flood["prop_flood_prone"] * param["subsidized_structure_value_ref"] * (spline_inflation(year_temp) / spline_inflation(0)) * subsidized_damages)

        informal_structure_damages = np.nansum(nb_households_informal * data_flood["prop_flood_prone"] * param["informal_structure_value_ref"] * (spline_inflation(year_temp) / spline_inflation(0)) * structural_damages_type2(data_flood['flood_depth']))
        backyard_structure_damages = ((16216 * (np.nansum(nb_households_backyard * data_flood["prop_flood_prone"] * param["informal_structure_value_ref"] * (spline_inflation(year_temp) / spline_inflation(0)) * structural_damages_type2(data_flood['flood_depth'])))) + (74916 * (np.nansum(nb_households_backyard * data_flood["prop_flood_prone"] * param["informal_structure_value_ref"] * (spline_inflation(year_temp) / spline_inflation(0)) * structural_damages_type3a(data_flood['flood_depth']))))) / (74916 + 16216)
        
        formal_content_damages = np.nansum(nb_households_formal * data_flood["prop_flood_prone"] * content_cost.formal * content_damages(data_flood['flood_depth']))
        subsidized_content_damages = np.nansum(nb_households_subsidized * data_flood["prop_flood_prone"] * content_cost.subsidized * content_damages(data_flood['flood_depth']))
        informal_content_damages = np.nansum(nb_households_informal * data_flood["prop_flood_prone"] * content_cost.informal * content_damages(data_flood['flood_depth']))
        backyard_content_damages = np.nansum(nb_households_backyard * data_flood["prop_flood_prone"] * content_cost.backyard * content_damages(data_flood['flood_depth']))
        
        damages = damages.append({'flood': type_flood,
                                  'formal_structure_damages': formal_structure_damages,
                                  'subsidized_structure_damages': subsidized_structure_damages,
                                  'informal_structure_damages': informal_structure_damages,
                                  'backyard_structure_damages': backyard_structure_damages,
                                  'formal_content_damages': formal_content_damages,
                                  'informal_content_damages': informal_content_damages,
                                  'backyard_content_damages': backyard_content_damages,
                                  'subsidized_content_damages': subsidized_content_damages}, ignore_index = True)
    
    return damages

def annualize_damages(array):
    interval0 = 1 - (1/5)    
    interval1 = (1/5) - (1/10)
    interval2 = (1/10) - (1/20)
    interval3 = (1/20) - (1/50)
    interval4 = (1/50) - (1/75)
    interval5 = (1/75) - (1/100)
    interval6 = (1/100) - (1/200)
    interval7 = (1/200) - (1/250)
    interval8 = (1/250) - (1/500)
    interval9 = (1/500) - (1/1000)
    interval10 = (1/1000)
    return 0.5 * ((interval0 * 0) + (interval1 * array[0]) + (interval2 * array[1]) + (interval3 * array[2]) + (interval4 * array[3]) + (interval5 * array[4]) + (interval6 * array[5]) + (interval7 * array[6]) + (interval8 * array[7]) + (interval9 * array[8]) + (interval10 * array[9]))

def compute_formal_structure_cost_method1(sp_price, dwelling_size_sp, SP_code, grid):
    
    formal_structure_cost = sp_price * dwelling_size_sp
    formal_structure_cost = SP_to_grid_2011_1(formal_structure_cost, SP_code, grid)
    formal_structure_cost[np.isinf(formal_structure_cost)] = np.nan
    formal_structure_cost[(formal_structure_cost) > 2000000] = 2000000
    return formal_structure_cost

def compute_formal_structure_cost_method2(initial_state_rent, param, interest_rate, coeff_land, initial_state_households_housing_types, construction_coeff):
    
    price_simul = (initial_state_rent[0, :] * construction_coeff * param["coeff_b"]/ (interest_rate + param["depreciation_rate"])) ** (1/param["coeff_a"])
    formal_structure_cost  = price_simul * (250000)  * coeff_land[0, :] / initial_state_households_housing_types[0, :]    
    formal_structure_cost[np.isinf(formal_structure_cost)] = np.nan
    formal_structure_cost[(formal_structure_cost) > 2000000] = 2000000
    return formal_structure_cost

def compute_content_cost(initial_state_household_centers, income_net_of_commuting_costs, param, fraction_capital_destroyed, initial_state_rent, initial_state_dwelling_size, interest_rate):
    
    content_cost = pd.DataFrame()
    
    income_class = np.nanargmax(initial_state_household_centers, 0)
    income_temp = np.empty(24014)
    income_temp[:] = np.nan
    for i in range(0, 24014):
        income_temp[i] = income_net_of_commuting_costs[int(income_class[i]), i]
    income_temp[income_temp < 0] = 0
    
    capital_destroyed = np.ones(len(fraction_capital_destroyed.structure_formal_2))
    capital_destroyed[initial_state_dwelling_size[3, :] > param["threshold"]] = fraction_capital_destroyed.structure_subsidized_2[initial_state_dwelling_size[3, :] > param["threshold"]]
    capital_destroyed[initial_state_dwelling_size[3, :] <= param["threshold"]] = fraction_capital_destroyed.structure_subsidized_1[initial_state_dwelling_size[3, :] <= param["threshold"]]
   
    
    fraction_backyard = (param["alpha"] * (param["RDP_size"] + param["backyard_size"] - param["q0"]) / (param["backyard_size"])) - (param["beta"] * (income_net_of_commuting_costs[0,:] - (capital_destroyed * param["subsidized_structure_value"])) / (param["backyard_size"] * initial_state_rent[1, :]))   
    fraction_backyard[initial_state_rent[1, :] == 0] = 0
    fraction_backyard = np.minimum(fraction_backyard, 1)
    fraction_backyard = np.maximum(fraction_backyard, 0) 
    
        
    content_cost["formal"] = param["fraction_z_dwellings"] * (1/(1 + (param["fraction_z_dwellings"] * fraction_capital_destroyed.contents_formal))) * (income_temp - (initial_state_rent[0, :] * initial_state_dwelling_size[0, :]))
    content_cost.formal[income_temp - (initial_state_rent[0, :] * initial_state_dwelling_size[0, :]) < 0] = np.nan
    content_cost.formal[content_cost.formal < (0.2 * income_temp)] = np.nan
    content_cost["informal"] = param["fraction_z_dwellings"] * (1/(1 + (param["fraction_z_dwellings"] * fraction_capital_destroyed.contents_informal))) * (income_temp - (initial_state_rent[2, :] * initial_state_dwelling_size[2, :]) - (fraction_capital_destroyed.structure_informal_settlements * param["informal_structure_value"]) - ((interest_rate + param["depreciation_rate"]) * param["informal_structure_value"]))
    content_cost["subsidized"] = param["fraction_z_dwellings"] * (1/(1 + (param["fraction_z_dwellings"] * fraction_capital_destroyed.contents_subsidized))) * (income_temp + (param["backyard_size"] * initial_state_rent[1, :] * fraction_backyard) - (capital_destroyed * param["subsidized_structure_value"]))
    content_cost["backyard"] = param["fraction_z_dwellings"] * (1/(1 + (param["fraction_z_dwellings"] * fraction_capital_destroyed.contents_backyard))) * (income_temp - (initial_state_rent[1, :] * initial_state_dwelling_size[1, :]) - (fraction_capital_destroyed.structure_backyards * param["informal_structure_value"]) - ((interest_rate + param["depreciation_rate"]) * param["informal_structure_value"]))    
    content_cost[content_cost < 0] = np.nan
    
    return content_cost
    