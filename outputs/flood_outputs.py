import pandas as pd
import numpy as np
import copy

# import inputs.data as inpdt


def compute_stats_per_housing_type(
        floods, path_floods, nb_households_formal, nb_households_subsidized,
        nb_households_informal, nb_households_backyard, path_tables,
        flood_categ, threshold=0.1):
    """
    Compute aggregate flood exposure statistics for a given flood type.

    More specifically, the function returns for each available return period
    and each housing type, an estimated total number of households exposed and
    associated average maximum flood depth level. The output is used as an
    argument by the validation_flood function (export_outputs_floods module)

    Parameters
    ----------
    floods : list
        List of file names for available flood maps per return period and flood
        type
    path_floods : str
        Path towards flood maps directory
    nb_households_formal : Series
        Number of households living in formal private housing, per grid cell
    nb_households_subsidized : Series
        Number of households living in formal subsidized housing, per grid cell
    nb_households_informal : Series
        Number of households living in informal settlements, per grid cell
    nb_households_backyard : Series
        Number of households living in informal backyards, per grid cell
    path_tables : str
        Path for saving output plots
    flood_categ : str
        Category of flood risks considered, used in name of output file
    threshold : float64, optional
        Maximum flood depth level (in m) below which we choose to discard flood
        risks. The default is 0.1, but is not used in benchmark version of the
        function.

    Returns
    -------
    stats_per_housing_type : DataFrame
        Table summarizing, for a given flood type and each associated return
        period, the estimated total number of households per housing type
        living in flood-prone areas, and the associated average maximum flood
        depth level

    """
    stats_per_housing_type = pd.DataFrame(
        columns=['flood', 'fraction_formal_in_flood_prone_area',
                 'fraction_subsidized_in_flood_prone_area',
                 'fraction_informal_in_flood_prone_area',
                 'fraction_backyard_in_flood_prone_area',
                 'flood_depth_formal', 'flood_depth_subsidized',
                 'flood_depth_informal', 'flood_depth_backyard']
        )

    for flood in floods:
        type_flood = copy.deepcopy(flood)
        flood = np.squeeze(pd.read_excel(path_floods + flood + ".xlsx"))

        # flood.prop_flood_prone[flood.flood_depth < threshold] = 0
        # flood.flood_depth[flood.flood_depth < threshold] = 0
        print(type_flood)

        if ((type_flood == 'P_5yr') | (type_flood == 'P_10yr')):
            stats_per_housing_type = pd.concat([
                stats_per_housing_type,
                pd.DataFrame(
                    {'flood': [type_flood],
                     'fraction_formal_in_flood_prone_area': [0],
                     'fraction_subsidized_in_flood_prone_area': [0],
                     'fraction_informal_in_flood_prone_area': [np.sum(
                         flood['prop_flood_prone'] * nb_households_informal)],
                     'fraction_backyard_in_flood_prone_area': [0],
                     'flood_depth_formal': [0],
                     'flood_depth_subsidized': [0],
                     'flood_depth_informal': [sum(
                         flood['flood_depth']
                         * (flood['prop_flood_prone'] * nb_households_informal)
                         / sum(
                             flood['prop_flood_prone']
                             * nb_households_informal)
                         )],
                     'flood_depth_backyard': [0]}
                    )],
                ignore_index=True)

        elif (type_flood == 'P_20yr'):
            stats_per_housing_type = pd.concat([
                stats_per_housing_type,
                pd.DataFrame(
                    {'flood': [type_flood],
                     'fraction_formal_in_flood_prone_area': [0],
                     'fraction_subsidized_in_flood_prone_area': [np.sum(
                         flood['prop_flood_prone']
                         * nb_households_subsidized)],
                     'fraction_informal_in_flood_prone_area': [np.sum(
                         flood['prop_flood_prone'] * nb_households_informal)],
                     'fraction_backyard_in_flood_prone_area': [np.sum(
                         flood['prop_flood_prone'] * nb_households_backyard)],
                     'flood_depth_formal': [0],
                     'flood_depth_subsidized': [sum(
                         flood['flood_depth']
                         * (flood['prop_flood_prone']
                            * nb_households_subsidized)
                         / sum(
                             flood['prop_flood_prone']
                             * nb_households_subsidized))],
                     'flood_depth_informal': [sum(
                         flood['flood_depth']
                         * (flood['prop_flood_prone'] * nb_households_informal)
                         / sum(
                             flood['prop_flood_prone']
                             * nb_households_informal)
                         )],
                     'flood_depth_backyard': [sum(
                         flood['flood_depth']
                         * (flood['prop_flood_prone'] * nb_households_backyard)
                         / sum(
                             flood['prop_flood_prone']
                             * nb_households_backyard)
                     )]}
                    )],
                ignore_index=True)

        else:
            stats_per_housing_type = pd.concat([
                stats_per_housing_type,
                pd.DataFrame(
                    {'flood': [type_flood],
                     'fraction_formal_in_flood_prone_area': [np.sum(
                         flood['prop_flood_prone'] * nb_households_formal)],
                     'fraction_subsidized_in_flood_prone_area': [np.sum(
                         flood['prop_flood_prone']
                         * nb_households_subsidized)],
                     'fraction_informal_in_flood_prone_area': [np.sum(
                         flood['prop_flood_prone'] * nb_households_informal)],
                     'fraction_backyard_in_flood_prone_area': [np.sum(
                         flood['prop_flood_prone'] * nb_households_backyard)],
                     'flood_depth_formal': [sum(
                         flood['flood_depth']
                         * (flood['prop_flood_prone'] * nb_households_formal)
                         / sum(
                             flood['prop_flood_prone']
                             * nb_households_formal)
                         )],
                     'flood_depth_subsidized': [sum(
                         flood['flood_depth']
                         * (flood['prop_flood_prone']
                            * nb_households_subsidized)
                         / sum(
                             flood['prop_flood_prone']
                             * nb_households_subsidized))],
                     'flood_depth_informal': [sum(
                         flood['flood_depth']
                         * (flood['prop_flood_prone']
                            * nb_households_informal)
                         / sum(
                             flood['prop_flood_prone']
                             * nb_households_informal)
                         )],
                     'flood_depth_backyard': [sum(
                         flood['flood_depth']
                         * (flood['prop_flood_prone'] * nb_households_backyard)
                         / sum(
                             flood['prop_flood_prone']
                             * nb_households_backyard)
                         )]}
                    )],
                ignore_index=True)

    stats_per_housing_type = stats_per_housing_type.fillna(value=0)
    stats_per_housing_type.to_csv(
        path_tables + flood_categ + '_stats_per_housing_type.csv')

    return stats_per_housing_type


def compute_stats_per_income_group(
        floods, path_floods, nb_households_rich, nb_households_midrich,
        nb_households_midpoor, nb_households_poor, path_tables,
        flood_categ, threshold=0.1):
    """
    Compute aggregate flood exposure statistics for a given flood type.

    More specifically, the function returns for each available return period
    and each income group, an estimated total number of households exposed and
    associated average maximum flood depth level.

    Parameters
    ----------
    floods : list
        List of file names for available flood maps per return period and flood
        type
    path_floods : str
        Path towards flood maps directory
    nb_households_rich : Series
        Number of rich households, per grid cell
    nb_households_midrich : Series
        Number of mid-rich households, per grid cell
    nb_households_midpoor : Series
        Number of mid-poor households, per grid cell
    nb_households_poor : Series
        Number of poor households, per grid cell
    path_tables : str
        Path for saving output plots
    flood_categ : str
        Category of flood risks considered, used in name of output file
    threshold : float64, optional
        Maximum flood depth level (in m) below which we choose to discard flood
        risks. The default is 0.1, but is not used in benchmark version of the
        function.

    Returns
    -------
    stats_per_income_group : DataFrame
        Table summatizing, for a given flood type and each associated return
        period, the estimated total number of households per income group
        living in flood-prone areas, and the associated average maximum flood
        depth level

    """
    stats_per_income_group = pd.DataFrame(
        columns=['flood', 'fraction_rich_in_flood_prone_area',
                 'fraction_midrich_in_flood_prone_area',
                 'fraction_midpoor_in_flood_prone_area',
                 'fraction_poor_in_flood_prone_area',
                 'flood_depth_rich', 'flood_depth_midrich',
                 'flood_depth_midpoor', 'flood_depth_poor']
        )

    for flood in floods:
        type_flood = copy.deepcopy(flood)
        flood = np.squeeze(pd.read_excel(path_floods + flood + ".xlsx"))

        # flood.prop_flood_prone[flood.flood_depth < threshold] = 0
        # flood.flood_depth[flood.flood_depth < threshold] = 0
        print(type_flood)

        stats_per_income_group = pd.concat([
            stats_per_income_group,
            pd.DataFrame(
                {'flood': [type_flood],
                 'fraction_rich_in_flood_prone_area': [np.sum(
                     flood['prop_flood_prone'] * nb_households_rich)],
                 'fraction_midrich_in_flood_prone_area': [np.sum(
                     flood['prop_flood_prone'] * nb_households_midrich)],
                 'fraction_midpoor_in_flood_prone_area': [np.sum(
                     flood['prop_flood_prone'] * nb_households_midpoor)],
                 'fraction_poor_in_flood_prone_area': [np.sum(
                     flood['prop_flood_prone'] * nb_households_poor)],
                 'flood_depth_rich': [sum(
                     flood['flood_depth']
                     * (flood['prop_flood_prone'] * nb_households_rich)
                     / sum(flood['prop_flood_prone'] * nb_households_rich)
                     )],
                 'flood_depth_midrich': [sum(
                     flood['flood_depth']
                     * (flood['prop_flood_prone'] * nb_households_midrich)
                     / sum(flood['prop_flood_prone'] * nb_households_midrich
                           ))],
                 'flood_depth_midpoor': [sum(
                     flood['flood_depth']
                     * (flood['prop_flood_prone'] * nb_households_midpoor)
                     / sum(flood['prop_flood_prone'] * nb_households_midpoor)
                     )],
                 'flood_depth_poor': [sum(
                     flood['flood_depth']
                     * (flood['prop_flood_prone'] * nb_households_poor)
                     / sum(flood['prop_flood_prone'] * nb_households_poor)
                     )]}
                )],
            ignore_index=True)

    stats_per_income_group = stats_per_income_group.fillna(value=0)
    stats_per_income_group.to_csv(
        path_tables + flood_categ + '_stats_per_income_group.csv')

    return stats_per_income_group


def compute_damages(floods, path_data, param, content_cost,
                    nb_households_formal, nb_households_subsidized,
                    nb_households_informal, nb_households_backyard,
                    dwelling_size, formal_structure_cost, content_damages,
                    structural_damages_type4b, structural_damages_type4a,
                    structural_damages_type2, structural_damages_type3a,
                    options, spline_inflation, year_temp,
                    path_tables, flood_categ):
    """
    Compute total structure and content damages per housing type.

    This function leverages the depth-damage functions from the literature
    to estimate the monetary value lost to floods based on the estimated total
    value of the underlying asset, per available return period. The logic is
    the same as in inputs.data.import_full_floods_data.

    Parameters
    ----------
    floods : list
        List of file names for available flood maps per return period and flood
        type
    path_data : str
        Path towards data used in the model
    param : dict
        Dictionary of default parameters
    content_cost : DataFrame
        Estimated value of composite good consumption that is considered as
        flood-prone, for each grid cell (24,014) and each housing type (4)
    nb_households_formal : Series
        Number of households living in formal private housing, per grid cell
    nb_households_subsidized : Series
        Number of households living in formal subsidized housing, per grid cell
    nb_households_informal : Series
        Number of households living in informal settlements, per grid cell
    nb_households_backyard : Series
        Number of households living in informal backyards, per grid cell
        DESCRIPTION.
    dwelling_size : ndarray(float64, ndim=2)
        Average dwelling size (in m²) per grid cell in each housing
        type (4)
    formal_structure_cost : ndarray(float64)
        Estimated construction cost of formal private housing structures, based
        on their market capital values, per grid cell (24,014)
    content_damages : interp1d
        Linear interpolation for fraction of capital destroyed (house
        contents) over maximum flood depth (in m) in a given area,
        from de Villiers et al., 2007
    structural_damages_type4b : interp1d
        Linear interpolation for fraction of capital destroyed (type-4b house
        structures) over maximum flood depth (in m) in a given area,
        from de Englhardt et al., 2019 (two-floor reinforced masonry/concrete
        and steel buildings)
    structural_damages_type4a : interp1d
        Linear interpolation for fraction of capital destroyed (type-4a house
        structures) over maximum flood depth (in m) in a given area,
        from de Englhardt et al., 2019 (one-floor reinforced masonry/concrete
        and steel buildings)
    structural_damages_type2 : interp1d
        Linear interpolation for fraction of capital destroyed (type-2 house
        structures) over maximum flood depth (in m) in a given area,
        from de Englhardt et al., 2019 (wooden buildings)
    structural_damages_type3a : interp1d
        Linear interpolation for fraction of capital destroyed (type-3a house
        structures) over maximum flood depth (in m) in a given area,
        from de Englhardt et al., 2019 (one-floor unreinforced masonry/concrete
        buildings)
    options : dict
        Dictionary of default options
    spline_inflation : interp1d
        Linear interpolation for inflation rate (in base 100 relative to
        baseline year) over the years (baseline year set at 0)
    year_temp : int
        Year (relative to baseline year set at 0) for which we want to
        run the function
    path_tables : str
        Path for saving output plots
    flood_categ : str
        Category of flood risks considered, used in name of output file

    Returns
    -------
    damages : DataFrame
        Table yielding, for each return period and housing types, the estimated
        total damages in terms of housing structures and contents

    """
    damages = pd.DataFrame(columns=['flood',
                                    'formal_structure_damages',
                                    'subsidized_structure_damages',
                                    'informal_structure_damages',
                                    'backyard_structure_damages',
                                    'formal_content_damages',
                                    'subsidized_content_damages',
                                    'informal_content_damages',
                                    'backyard_content_damages'])

    for item in floods:

        print(item)

        type_flood = copy.deepcopy(item)
        data_flood = np.squeeze(pd.read_excel(path_data + item + ".xlsx"))

        formal_damages = structural_damages_type4a(data_flood['flood_depth'])
        formal_damages[dwelling_size[0, :] > param["threshold"]
                       ] = structural_damages_type4b(
            data_flood.flood_depth[dwelling_size[0, :] > param["threshold"]])
        subsidized_damages = structural_damages_type4a(
            data_flood['flood_depth'])
        subsidized_damages[dwelling_size[3, :] > param["threshold"]
                           ] = structural_damages_type4b(
            data_flood.flood_depth[dwelling_size[3, :] > param["threshold"]])

        formal_structure_damages = np.nansum(
            nb_households_formal * data_flood["prop_flood_prone"]
            * formal_structure_cost * formal_damages)
        subsidized_structure_damages = np.nansum(
            nb_households_subsidized * data_flood["prop_flood_prone"]
            * param["subsidized_structure_value_ref"]
            * (spline_inflation(year_temp) / spline_inflation(0))
            * subsidized_damages)

        informal_structure_damages = np.nansum(
            nb_households_informal * data_flood["prop_flood_prone"]
            * param["informal_structure_value_ref"]
            * (spline_inflation(year_temp) / spline_inflation(0))
            * structural_damages_type2(data_flood['flood_depth']))

        # backyard_structure_damages = (
        #     16216 * np.nansum(
        #         nb_households_backyard * data_flood["prop_flood_prone"]
        #         * param["informal_structure_value_ref"]
        #         * (spline_inflation(year_temp) / spline_inflation(0))
        #         * structural_damages_type2(data_flood['flood_depth']))
        #     + 74916 * np.nansum(
        #         nb_households_backyard * data_flood["prop_flood_prone"]
        #         * param["informal_structure_value_ref"]
        #         * (spline_inflation(year_temp) / spline_inflation(0))
        #         * structural_damages_type3a(data_flood['flood_depth']))
        #     ) / (74916 + 16216)

        # In our benchmark, we only consider informal backyards
        backyard_structure_damages = np.nansum(
            nb_households_backyard * data_flood["prop_flood_prone"]
            * param["informal_structure_value_ref"]
            * (spline_inflation(year_temp) / spline_inflation(0))
            * structural_damages_type3a(data_flood['flood_depth']))

        formal_content_damages = np.nansum(
            nb_households_formal * data_flood["prop_flood_prone"]
            * content_cost.formal * content_damages(data_flood['flood_depth']))
        subsidized_content_damages = np.nansum(
            nb_households_subsidized * data_flood["prop_flood_prone"]
            * content_cost.subsidized
            * content_damages(data_flood['flood_depth']))
        informal_content_damages = np.nansum(
            nb_households_informal * data_flood["prop_flood_prone"]
            * content_cost.informal
            * content_damages(data_flood['flood_depth']))
        backyard_content_damages = np.nansum(
            nb_households_backyard * data_flood["prop_flood_prone"]
            * content_cost.backyard
            * content_damages(data_flood['flood_depth']))

        damages = pd.concat([
            damages,
            pd.DataFrame(
                {'flood': [type_flood],
                 'formal_structure_damages': [formal_structure_damages],
                 'subsidized_structure_damages': [
                     subsidized_structure_damages],
                 'informal_structure_damages': [informal_structure_damages],
                 'backyard_structure_damages': [backyard_structure_damages],
                 'formal_content_damages': [formal_content_damages],
                 'informal_content_damages': [informal_content_damages],
                 'backyard_content_damages': [backyard_content_damages],
                 'subsidized_content_damages': [subsidized_content_damages]}
                )],
            ignore_index=True)

    damages = damages.fillna(value=0)
    # damages[damages < 0] = 0
    damages.to_csv(
        path_tables + flood_categ + '_damages.csv')

    return damages


def compute_damages_2d(floods, path_data, param, content_cost,
                       nb_households_formal, nb_households_subsidized,
                       nb_households_informal, nb_households_backyard,
                       dwelling_size, formal_structure_cost, content_damages,
                       structural_damages_type4b, structural_damages_type4a,
                       structural_damages_type2, structural_damages_type3a,
                       options, spline_inflation, year_temp,
                       path_tables, flood_categ):
    """
    Compute structure and content damages per housing type across space.

    This function leverages the depth-damage functions from the literature
    to estimate the monetary value lost to floods based on the estimated total
    value of the underlying asset, per available return period. Here, we get
    spatial, and not aggregate, data. The use of this function instead of its
    1D counterpart depends on the plots we are interested in.

    Parameters
    ----------
    floods : list
        List of file names for available flood maps per return period and flood
        type
    path_data : str
        Path towards data used in the model
    param : dict
        Dictionary of default parameters
    content_cost : DataFrame
        Estimated value of composite good consumption that is considered as
        flood-prone, for each grid cell (24,014) and each housing type (4)
    nb_households_formal : Series
        Number of households living in formal private housing, per grid cell
    nb_households_subsidized : Series
        Number of households living in formal subsidized housing, per grid cell
    nb_households_informal : Series
        Number of households living in informal settlements, per grid cell
    nb_households_backyard : Series
        Number of households living in informal backyards, per grid cell
        DESCRIPTION.
    dwelling_size : ndarray(float64, ndim=2)
        Average dwelling size (in m²) per grid cell in each housing
        type (4)
    formal_structure_cost : ndarray(float64)
        Estimated construction cost of formal private housing structures, based
        on their market capital values, per grid cell (24,014)
    content_damages : interp1d
        Linear interpolation for fraction of capital destroyed (house
        contents) over maximum flood depth (in m) in a given area,
        from de Villiers et al., 2007
    structural_damages_type4b : interp1d
        Linear interpolation for fraction of capital destroyed (type-4b house
        structures) over maximum flood depth (in m) in a given area,
        from de Englhardt et al., 2019 (two-floor reinforced masonry/concrete
        and steel buildings)
    structural_damages_type4a : interp1d
        Linear interpolation for fraction of capital destroyed (type-4a house
        structures) over maximum flood depth (in m) in a given area,
        from de Englhardt et al., 2019 (one-floor reinforced masonry/concrete
        and steel buildings)
    structural_damages_type2 : interp1d
        Linear interpolation for fraction of capital destroyed (type-2 house
        structures) over maximum flood depth (in m) in a given area,
        from de Englhardt et al., 2019 (wooden buildings)
    structural_damages_type3a : interp1d
        Linear interpolation for fraction of capital destroyed (type-3a house
        structures) over maximum flood depth (in m) in a given area,
        from de Englhardt et al., 2019 (one-floor unreinforced masonry/concrete
        buildings)
    options : dict
        Dictionary of default options
    spline_inflation : interp1d
        Linear interpolation for inflation rate (in base 100 relative to
        baseline year) over the years (baseline year set at 0)
    year_temp : int
        Year (relative to baseline year set at 0) for which we want to
        run the function
    path_tables : str
        Path for saving output plots
    flood_categ : str
        Category of flood risks considered, used in name of output file

    Returns
    -------
    dict_damages : dict
        Dictionary yielding, for each return period, the estimated damages per
        grid cell (24,014) and housing type (4) in terms of housing structures
        and contents

    """
    dict_damages = {}

    for item in floods:

        print(item)

        data_flood = np.squeeze(pd.read_excel(path_data + item + ".xlsx"))

        formal_damages = structural_damages_type4a(data_flood['flood_depth'])
        formal_damages[dwelling_size[0, :] > param["threshold"]
                       ] = structural_damages_type4b(
            data_flood.flood_depth[dwelling_size[0, :] > param["threshold"]])
        subsidized_damages = structural_damages_type4a(
            data_flood['flood_depth'])
        subsidized_damages[dwelling_size[3, :] > param["threshold"]
                           ] = structural_damages_type4b(
            data_flood.flood_depth[dwelling_size[3, :] > param["threshold"]])

        formal_structure_damages = (
            nb_households_formal * data_flood["prop_flood_prone"]
            * formal_structure_cost * formal_damages)

        subsidized_structure_damages = (
            nb_households_subsidized * data_flood["prop_flood_prone"]
            * param["subsidized_structure_value_ref"]
            * (spline_inflation(year_temp) / spline_inflation(0))
            * subsidized_damages)

        informal_structure_damages = (
            nb_households_informal * data_flood["prop_flood_prone"]
            * param["informal_structure_value_ref"]
            * (spline_inflation(year_temp) / spline_inflation(0))
            * structural_damages_type2(data_flood['flood_depth']))

        # backyard_structure_damages = (
        #     16216 * np.nansum(
        #         nb_households_backyard * data_flood["prop_flood_prone"]
        #         * param["informal_structure_value_ref"]
        #         * (spline_inflation(year_temp) / spline_inflation(0))
        #         * structural_damages_type2(data_flood['flood_depth']))
        #     + 74916 * np.nansum(
        #         nb_households_backyard * data_flood["prop_flood_prone"]
        #         * param["informal_structure_value_ref"]
        #         * (spline_inflation(year_temp) / spline_inflation(0))
        #         * structural_damages_type3a(data_flood['flood_depth']))
        #     ) / (74916 + 16216)

        # In our benchmark, we only consider informal backyards
        backyard_structure_damages = (
            nb_households_backyard * data_flood["prop_flood_prone"]
            * param["informal_structure_value_ref"]
            * (spline_inflation(year_temp) / spline_inflation(0))
            * structural_damages_type3a(data_flood['flood_depth']))

        formal_content_damages = (
            nb_households_formal * data_flood["prop_flood_prone"]
            * content_cost.formal * content_damages(data_flood['flood_depth']))
        subsidized_content_damages = (
            nb_households_subsidized * data_flood["prop_flood_prone"]
            * content_cost.subsidized
            * content_damages(data_flood['flood_depth']))
        informal_content_damages = (
            nb_households_informal * data_flood["prop_flood_prone"]
            * content_cost.informal
            * content_damages(data_flood['flood_depth']))
        backyard_content_damages = (
            nb_households_backyard * data_flood["prop_flood_prone"]
            * content_cost.backyard
            * content_damages(data_flood['flood_depth']))

        damages = pd.DataFrame(
            {'formal_structure_damages': formal_structure_damages,
             'subsidized_structure_damages': subsidized_structure_damages,
             'informal_structure_damages': informal_structure_damages,
             'backyard_structure_damages': backyard_structure_damages,
             'formal_content_damages': formal_content_damages,
             'informal_content_damages': informal_content_damages,
             'backyard_content_damages': backyard_content_damages,
             'subsidized_content_damages': subsidized_content_damages})

        damages = damages.fillna(value=0)
        damages[damages < 0] = 0
        # damages.to_csv(
        #     path_tables + flood_categ + '_' + item + '_damages_2d.csv')
        dict_damages[item] = damages

    return dict_damages


def annualize_damages(array_init, type_flood, housing_type, options):
    """
    Return expected value of flood damages for given location and housing type.

    The logic is the same as in inputs.data.compute_fraction_capital_destroyed.

    Parameters
    ----------
    array_init : ndarray(float64)
        Array containing estimated damage values per available return period,
        for a given flood type, housing type, and grid cell
    type_flood : str
        Type of flood risk considered, used to determine the number of return
        periods available (depends on FATHOM/DELTARES data source)
    housing_type : str
        Housing type considered, used to determine which corrections to apply
        for pluvial flood risks
    options : dict
        Dictionary of default options

    Returns
    -------
    float64
        Annualized / expected value of future damage flows for a given flood
        type, housing type, and grid cell

    """
    array = array_init.copy()
    if type_flood == 'pluvial' and options["correct_pluvial"] == 1:
        if housing_type == 'formal':
            array[0] = 0
            array[1] = 0
            array[2] = 0
        elif housing_type == 'subsidized' or housing_type == 'backyard':
            array[0] = 0
            array[1] = 0

    if (type_flood == 'pluvial' or type_flood == 'fluviald'
       or type_flood == 'fluvialu'):
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
        damages0 = array[0]
        damages1 = array[0] + array[1]
        damages2 = array[1] + array[2]
        damages3 = array[2] + array[3]
        damages4 = array[3] + array[4]
        damages5 = array[4] + array[5]
        damages6 = array[5] + array[6]
        damages7 = array[6] + array[7]
        damages8 = array[7] + array[8]
        damages9 = array[8] + array[9]
        damages10 = array[9] + array[9]

        return (0.5
                * ((interval0 * damages0) + (interval1 * damages1)
                    + (interval2 * damages2) + (interval3 * damages3)
                    + (interval4 * damages4) + (interval5 * damages5)
                    + (interval6 * damages6) + (interval7 * damages7)
                    + (interval8 * damages8) + (interval9 * damages9)
                    + (interval10 * damages10)))

    elif type_flood == 'coastal':
        interval0 = 1 - (1/2)
        interval1 = (1/2) - (1/5)
        interval2 = (1/5) - (1/10)
        interval3 = (1/10) - (1/25)
        interval4 = (1/25) - (1/50)
        interval5 = (1/50) - (1/100)
        interval6 = (1/100) - (1/250)
        interval7 = (1/250)
        damages0 = array[0] + array[1]
        damages1 = array[1] + array[2]
        damages2 = array[2] + array[3]
        damages3 = array[3] + array[4]
        damages4 = array[4] + array[5]
        damages5 = array[5] + array[6]
        damages6 = array[6] + array[7]
        damages7 = array[7] + array[7]

        return (0.5
                * ((interval0 * damages0) + (interval1 * damages1)
                    + (interval2 * damages2) + (interval3 * damages3)
                    + (interval4 * damages4) + (interval5 * damages5)
                    + (interval6 * damages6) + (interval7 * damages7)))


def compute_formal_structure_cost(
        initial_state_rent, param, interest_rate, coeff_land,
        initial_state_households_housing_types, construction_coeff):
    """
    Estimate construction costs of formal private housing structures in space.

    Note that the estimation process relies on a theoretical relation linking
    market prices to capital values. The value obtained is therefore an
    outcome of the model, and may not correspond to accounting estimates that
    are common in the impact evaluation literature.

    Parameters
    ----------
    initial_state_rent : ndarray(float64, ndim=2)
        Average annual rent (in rands) per grid cell for each housing type (4)
    param : dict
        Dictionary of default parameters
    interest_rate : float64
        Real interest rate for the overall economy, corresponding to an average
        over past years
    coeff_land : ndarray(float64, ndim=2)
        Updated land availability for each grid cell (24,014) and each
        housing type (4: formal private, informal backyards, informal
        settlements, formal subsidized)
    initial_state_households_housing_types : ndarray(float64, ndim=2)
        Number of households per grid cell in each housing type (4)
    construction_coeff : ndarray(float64)
        (Calibrated) scale factor for the construction function of formal
        private developers

    Returns
    -------
    formal_structure_cost : ndarray(float64)
        Estimated construction cost of formal private housing structures, based
        on their market capital values, per grid cell (24,014)

    """
    # We convert price to capital per unit of land
    price_simul = (
        initial_state_rent[0, :] * construction_coeff * param["coeff_b"]
        / (interest_rate + param["depreciation_rate"])
        ) ** (1/param["coeff_a"])
    # TODO: Difference with variable below?
    # price_simul = initial_state_capital_land[0, :]

    # We multiply by available land area, and average the output across
    # households
    np.seterr(divide='ignore', invalid='ignore')
    formal_structure_cost = (
        price_simul * (250000) * coeff_land[0, :]
        / initial_state_households_housing_types[0, :])
    formal_structure_cost[np.isinf(formal_structure_cost)] = np.nan

    return formal_structure_cost


def compute_content_cost(
        initial_state_household_centers, initial_state_housing_supply,
        income_net_of_commuting_costs, param,
        fraction_capital_destroyed, initial_state_rent,
        initial_state_dwelling_size, interest_rate):
    """
    Estimate value of flood-prone composite good consumption in space.

    Again, this is based on model outcomes. Since cost estimates are specific
    to housing type, we rely on an estimation of average income per housing
    type that we plug into households' budget constraint.

    Parameters
    ----------
    initial_state_household_centers : ndarray(float64, ndim=2)
        Number of households per grid cell in each income group (4)
    initial_state_housing_supply : ndarray(float64, ndim=2)
        Housing supply per unit of available land (in m² per km²)
        for each housing type (4) in each grid cell
    income_net_of_commuting_costs : ndarray(float64, ndim=2)
        Expected annual income net of commuting costs (in rands, for
        one household), for each geographic unit, by income group (4)
    param : dict
        Dictionary of default parameters
    fraction_capital_destroyed : DataFrame
        Data frame of expected fractions of capital destroyed, for housing
        structures and contents in different housing types, in each
        grid cell (24,014)
    initial_state_rent : ndarray(float64, ndim=2)
        Average annual rent (in rands) per grid cell for each housing type (4)
    initial_state_dwelling_size : ndarray(float64, ndim=2)
        Average dwelling size (in m²) per grid cell in each housing type (4)
    interest_rate : float64
        Real interest rate for the overall economy, corresponding to an average
        over past years

    Returns
    -------
    content_cost : DataFrame
        Estimated value of composite good consumption that is considered as
        flood-prone, for each grid cell (24,014) and each housing type (4)

    """
    content_cost = pd.DataFrame()

    # We estimate average incomes net of commuting costs in space, per housing
    # type
    np.seterr(divide='ignore', invalid='ignore')
    income_formal = np.nansum(
        income_net_of_commuting_costs * initial_state_household_centers
        / np.nansum(initial_state_household_centers, 0), 0)
    income_formal[income_formal < 0] = np.nan
    np.seterr(divide='ignore', invalid='ignore')
    income_informal = np.nansum(
        income_net_of_commuting_costs[0:2, :]
        * initial_state_household_centers[0:2, :]
        / np.nansum(initial_state_household_centers[0:2, :], 0), 0)
    income_informal[income_informal < 0] = np.nan
    income_subsidized = income_net_of_commuting_costs[0, :]
    income_subsidized[income_subsidized < 0] = np.nan

    # We also define other useful variables

    # First, the fraction of capital destroyed in formal subsidized housing
    # (that pops up in those households' budget constraint)
    capital_destroyed = np.zeros(
        len(fraction_capital_destroyed.structure_formal_2))
    # capital_destroyed[:] = np.nan
    (capital_destroyed[initial_state_dwelling_size[3, :] > param["threshold"]]
     ) = fraction_capital_destroyed.structure_subsidized_2[
         initial_state_dwelling_size[3, :] > param["threshold"]]
    (capital_destroyed[initial_state_dwelling_size[3, :] <= param["threshold"]]
     ) = fraction_capital_destroyed.structure_subsidized_1[
         initial_state_dwelling_size[3, :] <= param["threshold"]]

    # Then, the fraction of RDP backyard (per grid cell) that is rented out
    # in equilibrium (in m²/m²). We remind that this corresponds to the
    # housing supply per unit of available land for informal backyards.
    fraction_backyard = initial_state_housing_supply[1, :] / 1000000

    # We just multiply the amount of composite good from the budget constraint
    # by the share parameter to obtain cost estimates

    content_cost["formal"] = (
        param["fraction_z_dwellings"]
        / (1 + param["fraction_z_dwellings"]
           * fraction_capital_destroyed.contents_formal)
        * (income_formal
           - initial_state_rent[0, :] * initial_state_dwelling_size[0, :])
        )

    content_cost["informal"] = (
        param["fraction_z_dwellings"]
        / (1 + param["fraction_z_dwellings"]
           * fraction_capital_destroyed.contents_informal)
        * (income_informal
           - initial_state_rent[2, :] * initial_state_dwelling_size[2, :]
           - fraction_capital_destroyed.structure_informal_settlements
           * param["informal_structure_value"]
           - (interest_rate + param["depreciation_rate"])
           * param["informal_structure_value"])
        )

    content_cost["subsidized"] = (
        param["fraction_z_dwellings"]
        / (1 + param["fraction_z_dwellings"]
           * fraction_capital_destroyed.contents_subsidized)
        * (income_subsidized
           + param["backyard_size"] * initial_state_rent[1, :]
           * fraction_backyard
           - (capital_destroyed + param["depreciation_rate"])
           * param["subsidized_structure_value"])
        )

    content_cost["backyard"] = (
        param["fraction_z_dwellings"] /
        (1 + param["fraction_z_dwellings"]
         * fraction_capital_destroyed.contents_backyard)
        * (income_informal
           - initial_state_rent[1, :] * initial_state_dwelling_size[1, :]
           - fraction_capital_destroyed.structure_backyards
           * param["informal_structure_value"]
           - (interest_rate + param["depreciation_rate"])
           * param["informal_structure_value"])
        )

    content_cost[content_cost < 0] = np.nan

    return content_cost


def create_flood_dict(floods, path_floods, path_tables,
                      sim_nb_households_poor, sim_nb_households_midpoor,
                      sim_nb_households_midrich, sim_nb_households_rich):
    """
    Create dictionary for household distribution in a given set of flood maps.

    The spatial distribution is broken into income groups, as this is the
    relevant dimension in the plot_flood_severity_distrib function
    (export_outputs_floods module) that takes the output of this function as an
    argument.

    Parameters
    ----------
    floods : list
        List of file names for available flood maps per return period and flood
        type
    path_floods : str
        Path towards flood maps directory
    path_tables : str
        Path for saving output plots
    sim_nb_households_poor : Series
        Number of poor households, per grid cell
    sim_nb_households_midpoor : Series
        Number of mid-poor households, per grid cell
    sim_nb_households_midrich : Series
        Number of mid-rich households, per grid cell
    sim_nb_households_rich : Series
        Number of rich households, per grid cell

    Returns
    -------
    dictio : dict
        Dictionary yielding, for each return period of a given flood type, the
        spatial distribution of households broken into income groups, along
        with the maximum flood depth level and fraction of flood-prone area in
        their residential location

    """
    dictio = {}
    for flood in floods:
        print(flood)
        flood_data = np.squeeze(pd.read_excel(path_floods + flood + ".xlsx"))
        sim_poor_index = pd.DataFrame(sim_nb_households_poor)
        sim_midpoor_index = pd.DataFrame(sim_nb_households_midpoor)
        sim_midrich_index = pd.DataFrame(sim_nb_households_midrich)
        sim_rich_index = pd.DataFrame(sim_nb_households_rich)
        sim_poor_index = sim_poor_index.rename(
            columns={sim_poor_index.columns[0]: 'sim_poor'})
        sim_midpoor_index = sim_midpoor_index.rename(
            columns={sim_midpoor_index.columns[0]: 'sim_midpoor'})
        sim_midrich_index = sim_midrich_index.rename(
            columns={sim_midrich_index.columns[0]: 'sim_midrich'})
        sim_rich_index = sim_rich_index.rename(
            columns={sim_rich_index.columns[0]: 'sim_rich'})
        flood_df = pd.merge(flood_data, sim_poor_index,
                            left_index=True, right_index=True)
        flood_df = pd.merge(flood_df, sim_midpoor_index,
                            left_index=True, right_index=True)
        flood_df = pd.merge(flood_df, sim_midrich_index,
                            left_index=True, right_index=True)
        flood_df = pd.merge(flood_df, sim_rich_index,
                            left_index=True, right_index=True)
        flood_df.to_csv(path_tables + flood + 'distrib_households.csv')
        dictio[flood] = flood_df
    return dictio
