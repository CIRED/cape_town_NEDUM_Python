import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

import outputs.flood_outputs as outfld


# %% Floods

def validation_flood(stats1, stats2, legend1, legend2, type_flood,
                     path_plots):
    """
    Validation bar plot for household spatial distribution across flood zones.

    The validation is done across some selected return periods for illustrative
    purposes. Breakdown is given by housing type. The output plots show the
    average maximum flood depth level to which households are exposed, and the
    total number of households living in flood-prone areas.

    Parameters
    ----------
    stats1 : DataFrame
        Table summarizing, for a given flood type and each associated return
        period, the estimated total number of households per housing type
        living in flood-prone areas, and the associated average maximum flood
        depth level. We expect validation data here.
    stats2 : DataFrame
        Table summarizing, for a given flood type and each associated return
        period, the estimated total number of households per housing type
        living in flood-prone areas, and the associated average maximum flood
        depth level. We expect simulation results here.
    legend1 : str
        Legend for first data frame
    legend2 : str
        Legend for second data frame
    type_flood : str
        Type of flood risk considered, used to define output file name
    path_plots : str
        Path for saving output plots

    Returns
    -------
    None.

    """
    label = ["Formal private", "Formal subsidized",
             "Informal \n settlements", "Informal \n backyards"]

    # FLOOD DEPTH

    # First for validation data
    # RP = 20 yrs
    flood_depth_20_1 = [stats1.flood_depth_formal[2],
                        stats1.flood_depth_subsidized[2],
                        stats1.flood_depth_informal[2],
                        stats1.flood_depth_backyard[2]]
    # RP = 50 yrs
    flood_depth_50_1 = [stats1.flood_depth_formal[3],
                        stats1.flood_depth_subsidized[3],
                        stats1.flood_depth_informal[3],
                        stats1.flood_depth_backyard[3]]
    # RP = 100 yrs
    flood_depth_100_1 = [
        stats1.flood_depth_formal[5], stats1.flood_depth_subsidized[5],
        stats1.flood_depth_informal[5], stats1.flood_depth_backyard[5]]

    # Then for simulation
    # RP = 20 yrs
    flood_depth_20_2 = [stats2.flood_depth_formal[2],
                        stats2.flood_depth_subsidized[2],
                        stats2.flood_depth_informal[2],
                        stats2.flood_depth_backyard[2]]
    # RP = 50 yrs
    flood_depth_50_2 = [stats2.flood_depth_formal[3],
                        stats2.flood_depth_subsidized[3],
                        stats2.flood_depth_informal[3],
                        stats2.flood_depth_backyard[3]]
    # RP = 100 yrs
    flood_depth_100_2 = [stats2.flood_depth_formal[5],
                         stats2.flood_depth_subsidized[5],
                         stats2.flood_depth_informal[5],
                         stats2.flood_depth_backyard[5]]

    colors = ['#FF9999', '#00BFFF', '#C1FFC1', '#CAE1FF', '#FFDEAD']
    r = np.arange(len(label))
    barWidth = 0.25
    plt.figure(figsize=(10, 7))
    plt.bar(r, np.array(flood_depth_20_1),
            color=colors[1], edgecolor='white', width=barWidth,
            label='20 years')
    # valueb = np.maximum(np.array(flood_depth_50_1) - np.array(flood_depth_20_1),
    # np.full(4, 0.003))
    # floorb = np.maximum(np.array(flood_depth_50_), np.array(flood_depth_20_))
    plt.bar(r, np.array(flood_depth_50_1) - np.array(flood_depth_20_1),
            bottom=np.array(flood_depth_20_1), color=colors[2],
            edgecolor='white', width=barWidth, label='50 years')
    # valuec = np.maximum(np.array(flood_depth_100_1) - floorb, np.full(4, 0.003))
    plt.bar(r, np.array(flood_depth_100_1) - np.array(flood_depth_50_1),
            bottom=np.array(flood_depth_50_1), color=colors[3],
            edgecolor='white', width=barWidth, label='100 years')
    plt.bar(r + barWidth, np.array(flood_depth_20_2),
            color=colors[1], edgecolor='white', width=barWidth)
    # valueb2 = np.maximum(np.array(flood_depth_50_2) - np.array(flood_depth_20_2),
    # np.full(4, 0.003))
    # floorb2 = np.maximum(np.array(flood_depth_50_2), np.array(flood_depth_20_2))
    plt.bar(r + barWidth,
            np.array(flood_depth_50_2) - np.array(flood_depth_20_2),
            bottom=np.array(flood_depth_20_2), color=colors[2],
            edgecolor='white', width=barWidth)
    # valuec2 = np.maximum(np.array(flood_depth_100_2) - floorb2,
    # np.full(4, 0.003))
    plt.bar(r + barWidth,
            np.array(flood_depth_100_2) - np.array(flood_depth_50_2),
            bottom=np.array(flood_depth_50_2), color=colors[3],
            edgecolor='white', width=barWidth)
    plt.legend()
    plt.xticks(r + barWidth/2, label)
    # plt.ylim(0, 1)

    plt.text(r[0] - 0.1,
             np.maximum(
                 stats1.flood_depth_formal[5], stats1.flood_depth_formal[2]
                 ) + 0.002,
             legend1)
    plt.text(r[1] - 0.1,
             np.maximum(
                 stats1.flood_depth_subsidized[5],
                 stats1.flood_depth_subsidized[2]
                 ) + 0.002,
             legend1)
    plt.text(r[2] - 0.1,
             np.maximum(
                 stats1.flood_depth_informal[5], stats1.flood_depth_informal[2]
                 ) + 0.002,
             legend1)
    plt.text(r[3] - 0.1,
             np.maximum(
                 stats1.flood_depth_backyard[5], stats1.flood_depth_backyard[2]
                 ) + 0.002,
             legend1)
    plt.text(r[0] + 0.15,
             np.maximum(
                 stats2.flood_depth_formal[5], stats2.flood_depth_formal[2]
                 ) + 0.002,
             legend2)
    plt.text(r[1] + 0.15,
             np.maximum(
                 stats2.flood_depth_subsidized[5],
                 stats2.flood_depth_subsidized[2]
                 ) + 0.002,
             legend2)
    plt.text(r[2] + 0.15,
             np.maximum(
                 stats2.flood_depth_informal[5], stats2.flood_depth_informal[2]
                 ) + 0.002,
             legend2)
    plt.text(r[3] + 0.15,
             np.maximum(
                 stats2.flood_depth_backyard[5], stats2.flood_depth_backyard[2]
                 ) + 0.002,
             legend2)
    plt.ylabel("Average flood depth (m)", labelpad=15)
    plt.tick_params(labelbottom=True)
    plt.savefig(path_plots + 'validation_flood_depth_' + type_flood + '.png')
    # plt.show()
    plt.close()

    # FLOOD-PRONE AREA

    flood_area_20_1 = [
        stats1.fraction_formal_in_flood_prone_area[2],
        stats1.fraction_subsidized_in_flood_prone_area[2],
        stats1.fraction_informal_in_flood_prone_area[2],
        stats1.fraction_backyard_in_flood_prone_area[2]]
    flood_area_50_1 = [
        stats1.fraction_formal_in_flood_prone_area[3],
        stats1.fraction_subsidized_in_flood_prone_area[3],
        stats1.fraction_informal_in_flood_prone_area[3],
        stats1.fraction_backyard_in_flood_prone_area[3]]
    flood_area_100_1 = [
        stats1.fraction_formal_in_flood_prone_area[5],
        stats1.fraction_subsidized_in_flood_prone_area[5],
        stats1.fraction_informal_in_flood_prone_area[5],
        stats1.fraction_backyard_in_flood_prone_area[5]]

    flood_area_20_2 = [
        stats2.fraction_formal_in_flood_prone_area[2],
        stats2.fraction_subsidized_in_flood_prone_area[2],
        stats2.fraction_informal_in_flood_prone_area[2],
        stats2.fraction_backyard_in_flood_prone_area[2]]
    flood_area_50_2 = [
        stats2.fraction_formal_in_flood_prone_area[3],
        stats2.fraction_subsidized_in_flood_prone_area[3],
        stats2.fraction_informal_in_flood_prone_area[3],
        stats2.fraction_backyard_in_flood_prone_area[3]]
    flood_area_100_2 = [
        stats2.fraction_formal_in_flood_prone_area[5],
        stats2.fraction_subsidized_in_flood_prone_area[5],
        stats2.fraction_informal_in_flood_prone_area[5],
        stats2.fraction_backyard_in_flood_prone_area[5]]

    colors = ['#FF9999', '#00BFFF', '#C1FFC1', '#CAE1FF', '#FFDEAD']
    r = np.arange(len(label))
    barWidth = 0.25
    plt.figure(figsize=(10, 7))
    plt.bar(r, flood_area_20_1, color=colors[0], edgecolor='white',
            width=barWidth, label="20 years")
    plt.bar(r, np.array(flood_area_50_1) - np.array(flood_area_20_1),
            bottom=np.array(flood_area_20_1), color=colors[1],
            edgecolor='white', width=barWidth, label='50 years')
    plt.bar(r, np.array(flood_area_100_1) - np.array(flood_area_50_1),
            bottom=np.array(flood_area_50_1), color=colors[2],
            edgecolor='white', width=barWidth, label='100 years')
    plt.bar(r + barWidth, np.array(flood_area_20_2),
            color=colors[0], edgecolor='white', width=barWidth)
    plt.bar(r + barWidth,
            np.array(flood_area_50_2) - np.array(flood_area_20_2),
            bottom=np.array(flood_area_20_2), color=colors[1],
            edgecolor='white', width=barWidth)
    plt.bar(r + barWidth,
            np.array(flood_area_100_2) - np.array(flood_area_50_2),
            bottom=np.array(flood_area_50_2), color=colors[2],
            edgecolor='white', width=barWidth)

    # plt.legend(loc='upper right')
    plt.legend()
    plt.xticks(r + barWidth/2, label)
    plt.text(
        r[0] - 0.1,
        stats1.fraction_formal_in_flood_prone_area[5] + 1000,
        legend1)
    plt.text(
        r[1] - 0.1,
        stats1.fraction_subsidized_in_flood_prone_area[5] + 1000,
        legend1)
    plt.text(
        r[2] - 0.1,
        stats1.fraction_informal_in_flood_prone_area[5] + 1000,
        legend1)
    plt.text(
        r[3] - 0.1,
        stats1.fraction_backyard_in_flood_prone_area[5] + 1000,
        legend1)
    plt.text(
        r[0] + 0.15,
        stats2.fraction_formal_in_flood_prone_area[5] + 1000,
        legend2)
    plt.text(
        r[1] + 0.15,
        stats2.fraction_subsidized_in_flood_prone_area[5] + 1000,
        legend2)
    plt.text(
        r[2] + 0.15,
        stats2.fraction_informal_in_flood_prone_area[5] + 1000,
        legend2)
    plt.text(
        r[3] + 0.15,
        stats2.fraction_backyard_in_flood_prone_area[5] + 1000,
        legend2)
    plt.tick_params(labelbottom=True)
    plt.ylabel("Dwellings in flood-prone areas", labelpad=15)
    plt.savefig(path_plots + 'validation_flood_proportion_'
                + type_flood + '.png')
    # plt.show()
    plt.close()


def validation_flood_coastal(stats1, stats2, legend1, legend2, type_flood,
                             path_plots):
    """
    Validation bar plot for household spatial distribution across flood zones.

    The validation is done across some selected return periods for illustrative
    purposes. Breakdown is given by housing type. The output plots show the
    average maximum flood depth level to which households are exposed, and the
    total number of households living in flood-prone areas. This function is
    specific to coastal floods as the underlying flood maps are not available
    for the same return periods as fluvial / pluvial flood maps.

    Parameters
    ----------
    stats1 : DataFrame
        Table summarizing, for a given flood type and each associated return
        period, the estimated total number of households per housing type
        living in flood-prone areas, and the associated average maximum flood
        depth level. We expect validation data here.
    stats2 : DataFrame
        Table summarizing, for a given flood type and each associated return
        period, the estimated total number of households per housing type
        living in flood-prone areas, and the associated average maximum flood
        depth level. We expect simulation results here.
    legend1 : str
        Legend for first data frame
    legend2 : str
        Legend for second data frame
    type_flood : str
        Type of flood risk considered, used to define output file name
    path_plots : str
        Path for saving output plots

    Returns
    -------
    None.

    """
    label = ["Formal private", "Formal subsidized",
             "Informal \n settlements", "Informal \n backyards"]

    # FLOOD DEPTH

    # First for validation data
    # RP = 25 yrs
    flood_depth_25_1 = [stats1.flood_depth_formal[4],
                        stats1.flood_depth_subsidized[4],
                        stats1.flood_depth_informal[4],
                        stats1.flood_depth_backyard[4]]
    # RP = 50 yrs
    flood_depth_50_1 = [stats1.flood_depth_formal[5],
                        stats1.flood_depth_subsidized[5],
                        stats1.flood_depth_informal[5],
                        stats1.flood_depth_backyard[5]]
    # RP = 100 yrs
    flood_depth_100_1 = [stats1.flood_depth_formal[6],
                         stats1.flood_depth_subsidized[6],
                         stats1.flood_depth_informal[6],
                         stats1.flood_depth_backyard[6]]

    # Then for simulation
    # RP = 20 yrs
    flood_depth_25_2 = [stats2.flood_depth_formal[4],
                        stats2.flood_depth_subsidized[4],
                        stats2.flood_depth_informal[4],
                        stats2.flood_depth_backyard[4]]
    # RP = 50 yrs
    flood_depth_50_2 = [stats2.flood_depth_formal[5],
                        stats2.flood_depth_subsidized[5],
                        stats2.flood_depth_informal[5],
                        stats2.flood_depth_backyard[5]]
    # RP = 100 yrs
    flood_depth_100_2 = [stats2.flood_depth_formal[6],
                         stats2.flood_depth_subsidized[6],
                         stats2.flood_depth_informal[6],
                         stats2.flood_depth_backyard[6]]

    colors = ['#FF9999', '#00BFFF', '#C1FFC1', '#CAE1FF', '#FFDEAD']
    r = np.arange(len(label))
    barWidth = 0.25
    plt.figure(figsize=(10, 7))
    plt.bar(r, np.array(flood_depth_25_1),
            color=colors[1], edgecolor='white', width=barWidth,
            label='25 years')
    # valueb = np.maximum(np.array(flood_depth_50_1) - np.array(flood_depth_25_1),
    # np.full(4, 0.003))
    # floorb = np.maximum(np.array(flood_depth_50_1), np.array(flood_depth_25_1))
    plt.bar(r, np.array(flood_depth_50_1) - np.array(flood_depth_25_1),
            bottom=np.array(flood_depth_25_1), color=colors[2],
            edgecolor='white', width=barWidth, label='50 years')
    # valuec = np.maximum(np.array(flood_depth_100_1) - floorb, np.full(4, 0.003))
    plt.bar(r, np.array(flood_depth_100_1) - np.array(flood_depth_50_1),
            bottom=np.array(flood_depth_50_1), color=colors[3],
            edgecolor='white', width=barWidth, label='100 years')

    plt.bar(r + barWidth, np.array(flood_depth_25_2),
            color=colors[1], edgecolor='white', width=barWidth)
    # valueb2 = np.maximum(np.array(flood_depth_50_2) - np.array(flood_depth_25_2),
    # np.full(4, 0.003))
    # floorb2 = np.maximum(np.array(flood_depth_50_2), np.array(flood_depth_25_2))
    plt.bar(r + barWidth,
            np.array(flood_depth_50_2) - np.array(flood_depth_25_2),
            bottom=np.array(flood_depth_25_2), color=colors[2],
            edgecolor='white', width=barWidth)
    # valuec2 = np.maximum(np.array(flood_depth_100_2) - floorb2,
    # np.full(4, 0.003))
    plt.bar(r + barWidth,
            np.array(flood_depth_100_2) - np.array(flood_depth_50_2),
            bottom=np.array(flood_depth_50_2), color=colors[3],
            edgecolor='white', width=barWidth)
    plt.legend()
    plt.xticks(r + barWidth/2, label)
    # plt.ylim(0, 1)

    plt.text(r[0] - 0.1,
             np.maximum(
                 stats1.flood_depth_formal[6], stats1.flood_depth_formal[4]
                 ) + 0.01,
             legend1)
    plt.text(r[1] - 0.1,
             np.maximum(
                 stats1.flood_depth_subsidized[6],
                 stats1.flood_depth_subsidized[4]
                 ) + 0.01,
             legend1)
    plt.text(r[2] - 0.1,
             np.maximum(
                 stats1.flood_depth_informal[6], stats1.flood_depth_informal[4]
                 ) + 0.01,
             legend1)
    plt.text(r[3] - 0.1,
             np.maximum(
                 stats1.flood_depth_backyard[6], stats1.flood_depth_backyard[4]
                 ) + 0.01,
             legend1)
    plt.text(r[0] + 0.15,
             np.maximum(
                 stats2.flood_depth_formal[6], stats2.flood_depth_formal[4]
                 ) + 0.01,
             legend2)
    plt.text(r[1] + 0.15,
             np.maximum(
                 stats2.flood_depth_subsidized[6],
                 stats2.flood_depth_subsidized[4]
                 ) + 0.01,
             legend2)
    plt.text(r[2] + 0.15,
             np.maximum(
                 stats2.flood_depth_informal[6], stats2.flood_depth_informal[4]
                 ) + 0.01,
             legend2)
    plt.text(r[3] + 0.15,
             np.maximum(
                 stats2.flood_depth_backyard[6], stats2.flood_depth_backyard[4]
                 ) + 0.01,
             legend2)
    plt.ylabel("Average flood depth (m)", labelpad=15)
    plt.tick_params(labelbottom=True)
    plt.savefig(path_plots + 'validation_flood_depth_' + type_flood + '.png')
    # plt.show()
    plt.close()

    # FLOOD-PRONE AREA

    flood_area_25_1 = [
        stats1.fraction_formal_in_flood_prone_area[4],
        stats1.fraction_subsidized_in_flood_prone_area[4],
        stats1.fraction_informal_in_flood_prone_area[4],
        stats1.fraction_backyard_in_flood_prone_area[4]]
    flood_area_50_1 = [
        stats1.fraction_formal_in_flood_prone_area[5],
        stats1.fraction_subsidized_in_flood_prone_area[5],
        stats1.fraction_informal_in_flood_prone_area[5],
        stats1.fraction_backyard_in_flood_prone_area[5]]
    flood_area_100_1 = [
        stats1.fraction_formal_in_flood_prone_area[6],
        stats1.fraction_subsidized_in_flood_prone_area[6],
        stats1.fraction_informal_in_flood_prone_area[6],
        stats1.fraction_backyard_in_flood_prone_area[6]]
    flood_area_25_2 = [
        stats2.fraction_formal_in_flood_prone_area[4],
        stats2.fraction_subsidized_in_flood_prone_area[4],
        stats2.fraction_informal_in_flood_prone_area[4],
        stats2.fraction_backyard_in_flood_prone_area[4]]
    flood_area_50_2 = [
        stats2.fraction_formal_in_flood_prone_area[5],
        stats2.fraction_subsidized_in_flood_prone_area[5],
        stats2.fraction_informal_in_flood_prone_area[5],
        stats2.fraction_backyard_in_flood_prone_area[5]]
    flood_area_100_2 = [
        stats2.fraction_formal_in_flood_prone_area[6],
        stats2.fraction_subsidized_in_flood_prone_area[6],
        stats2.fraction_informal_in_flood_prone_area[6],
        stats2.fraction_backyard_in_flood_prone_area[6]]

    colors = ['#FF9999', '#00BFFF', '#C1FFC1', '#CAE1FF', '#FFDEAD']
    r = np.arange(len(label))
    barWidth = 0.25
    plt.figure(figsize=(10, 7))
    plt.bar(r, flood_area_25_1, color=colors[0], edgecolor='white',
            width=barWidth, label="25 years")
    plt.bar(r, np.array(flood_area_50_1) - np.array(flood_area_25_1),
            bottom=np.array(flood_area_25_1), color=colors[1],
            edgecolor='white', width=barWidth, label='50 years')
    plt.bar(r, np.array(flood_area_100_1) - np.array(flood_area_50_1),
            bottom=np.array(flood_area_50_1), color=colors[2],
            edgecolor='white', width=barWidth, label='100 years')
    plt.bar(r + barWidth, np.array(flood_area_25_2),
            color=colors[0], edgecolor='white', width=barWidth)
    plt.bar(r + barWidth,
            np.array(flood_area_50_2) - np.array(flood_area_25_2),
            bottom=np.array(flood_area_25_2), color=colors[1],
            edgecolor='white', width=barWidth)
    plt.bar(r + barWidth,
            np.array(flood_area_100_2) - np.array(flood_area_50_2),
            bottom=np.array(flood_area_50_2), color=colors[2],
            edgecolor='white', width=barWidth)

    # plt.legend(loc='upper right')
    plt.legend()
    plt.xticks(r + barWidth/2, label)
    plt.text(
        r[0] - 0.1,
        stats1.fraction_formal_in_flood_prone_area[6] + 20,
        legend1)
    plt.text(
        r[1] - 0.1,
        stats1.fraction_subsidized_in_flood_prone_area[6] + 20,
        legend1)
    plt.text(
        r[2] - 0.1,
        stats1.fraction_informal_in_flood_prone_area[6] + 20,
        legend1)
    plt.text(
        r[3] - 0.1,
        stats1.fraction_backyard_in_flood_prone_area[6] + 20,
        legend1)
    plt.text(
        r[0] + 0.15,
        stats2.fraction_formal_in_flood_prone_area[6] + 20,
        legend2)
    plt.text(
        r[1] + 0.15,
        stats2.fraction_subsidized_in_flood_prone_area[6] + 20,
        legend2)
    plt.text(
        r[2] + 0.15,
        stats2.fraction_informal_in_flood_prone_area[6] + 20,
        legend2)
    plt.text(
        r[3] + 0.15,
        stats2.fraction_backyard_in_flood_prone_area[6] + 20,
        legend2)
    plt.tick_params(labelbottom=True)
    plt.ylabel("Dwellings in flood-prone areas", labelpad=15)
    plt.savefig(path_plots + 'validation_flood_proportion_'
                + type_flood + '.png')
    # plt.show()
    plt.close()


def valid_damages(damages1, damages2, path_plots, flood_categ, options):
    """
    Validation bar plot for estimated damages for a given flood type.

    Breakdown is given by housing type as depth-damage functions used to
    estimate costs are specific to building materials. Here, we just annualize
    previously estimated total damages across available return periods. The
    function returns two separate plots for damages done to housing structures
    and contents.

    Parameters
    ----------
    damages1 : DataFrame
        Table yielding, for each return period and housing types, the estimated
        total damages in terms of housing structures and contents. Here, we
        expect simulation results.
    damages2 : DataFrame
        Table yielding, for each return period and housing types, the estimated
        total damages in terms of housing structures and contents. Here, we
        expect validation data.
    path_plots : str
        Path for saving output plots
    flood_categ : str
        Type of flood risk considered, used to define output file name
    options : dict
        Dictionary of default options

    Returns
    -------
    None.

    """
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    plt.rcParams['xtick.bottom'] = False
    plt.rcParams['xtick.labelbottom'] = False
    plt.rcParams['ytick.left'] = True
    plt.rcParams['ytick.labelleft'] = True

    data1 = [
        [outfld.annualize_damages(
            damages1.formal_structure_damages, flood_categ, 'formal', options)
            / 1000000,
         outfld.annualize_damages(
             damages1.subsidized_structure_damages, flood_categ, 'subsidized',
             options) / 1000000,
         outfld.annualize_damages(
             damages1.informal_structure_damages, flood_categ, 'informal',
             options) / 1000000,
         outfld.annualize_damages(
             damages1.backyard_structure_damages, flood_categ, 'backyard',
             options) / 1000000],
        [outfld.annualize_damages(
            damages1.formal_content_damages, flood_categ, 'formal', options)
            / 1000000,
         outfld.annualize_damages(
             damages1.subsidized_content_damages, flood_categ, 'subsidized',
             options) / 1000000,
         outfld.annualize_damages(
             damages1.informal_content_damages, flood_categ, 'informal',
             options) / 1000000,
         outfld.annualize_damages(
             damages1.backyard_content_damages, flood_categ, 'backyard',
             options) / 1000000]
        ]
    data2 = [
        [outfld.annualize_damages(
            damages2.formal_structure_damages, flood_categ, 'formal', options)
            / 1000000,
         outfld.annualize_damages(
             damages2.subsidized_structure_damages, flood_categ, 'subsidized',
             options) / 1000000,
         outfld.annualize_damages(
             damages2.informal_structure_damages, flood_categ, 'informal',
             options) / 1000000,
         outfld.annualize_damages(
             damages2.backyard_structure_damages, flood_categ, 'backyard',
             options) / 1000000],
        [outfld.annualize_damages(
            damages2.formal_content_damages, flood_categ, 'formal', options)
            / 1000000,
         outfld.annualize_damages(
             damages2.subsidized_content_damages, flood_categ, 'subsidized',
             options) / 1000000,
         outfld.annualize_damages(
             damages2.informal_content_damages, flood_categ, 'informal',
             options) / 1000000,
         outfld.annualize_damages(
             damages2.backyard_content_damages, flood_categ, 'backyard',
             options) / 1000000]
        ]
    barWidth = 0.25
    X = np.arange(4)
    colors = ['#FF9999', '#00BFFF', '#C1FFC1', '#CAE1FF', '#FFDEAD']
    plt.figure(figsize=(10, 7))
    plt.bar(
        X - barWidth/2, data1[0], color=colors[1], width=barWidth,
        label="Structures (sim)")
    plt.bar(
        X + barWidth/2, data2[0], color=colors[2], width=barWidth,
        label="Structures (data)")
    plt.legend()
    # plt.ylim(0, 25)
    quarter = ["Formal private", "Formal subsidized",
               "Informal \n settlements", "Informal \n in backyards"]
    plt.xticks(X, quarter)
    plt.tick_params(labelbottom=True)
    plt.ylabel("Million R per year")
    plt.savefig(path_plots + 'valid_' + flood_categ
                + '_structures_damages.png')
    # plt.show()
    plt.close()

    plt.figure(figsize=(10, 7))
    plt.bar(
        X - barWidth/2, data1[1], color=colors[1], width=barWidth,
        label="Contents (sim)")
    plt.bar(
        X + barWidth/2, data2[1], color=colors[2], width=barWidth,
        label="Contents (data)")
    plt.legend()
    # plt.ylim(0, 25)
    quarter = ["Formal private", "Formal subsidized",
               "Informal \n settlements", "Informal \n in backyards"]
    plt.xticks(X, quarter)
    plt.tick_params(labelbottom=True)
    plt.ylabel("Million R per year")
    plt.savefig(path_plots + 'valid_' + flood_categ + '_contents_damages.png')
    # plt.show()
    plt.close()


def simul_damages(damages, path_plots, flood_categ, options):
    """
    Bar plot for estimated damages for a given flood type.

    Breakdown is given by housing type as depth-damage functions used to
    estimate costs are specific to building materials. Here, we just annualize
    previously estimated total damages across available return periods. The
    function returns two separate plots for damages done to housing structures
    and contents. This function (as opposed to plot_damages) is specifically
    used for subsequent periods when validation data is not available.

    Parameters
    ----------
    damages : DataFrame
        Table yielding, for each return period and housing types, the estimated
        total damages in terms of housing structures and contents
    path_plots : str
        Path for saving output plots
    flood_categ : str
        Type of flood risk considered, used to define output file name
    options : dict
        Dictionary of default options

    Returns
    -------
    None.

    """
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    plt.rcParams['xtick.bottom'] = False
    plt.rcParams['xtick.labelbottom'] = False
    plt.rcParams['ytick.left'] = True
    plt.rcParams['ytick.labelleft'] = True

    data = [
        [outfld.annualize_damages(
            damages.formal_structure_damages, flood_categ, 'formal', options)
            / 1000000,
         outfld.annualize_damages(
             damages.subsidized_structure_damages, flood_categ, 'subsidized',
             options) / 1000000,
         outfld.annualize_damages(
             damages.informal_structure_damages, flood_categ, 'informal',
             options) / 1000000,
         outfld.annualize_damages(
             damages.backyard_structure_damages, flood_categ, 'backyard',
             options) / 1000000],
        [outfld.annualize_damages(
            damages.formal_content_damages, flood_categ, 'formal', options)
            / 1000000,
         outfld.annualize_damages(
             damages.subsidized_content_damages, flood_categ, 'subsidized',
             options) / 1000000,
         outfld.annualize_damages(
             damages.informal_content_damages, flood_categ, 'informal',
             options) / 1000000,
         outfld.annualize_damages(
             damages.backyard_content_damages, flood_categ, 'backyard',
             options) / 1000000]
        ]

    barWidth = 0.25
    X = np.arange(4)
    colors = ['#FF9999', '#00BFFF', '#C1FFC1', '#CAE1FF', '#FFDEAD']
    plt.figure(figsize=(10, 7))
    plt.bar(
        X, data[0], color=colors[1], width=barWidth,
        label="Structures")
    plt.legend()
    # plt.ylim(0, 25)
    quarter = ["Formal private", "Formal subsidized",
               "Informal \n settlements", "Informal \n in backyards"]
    plt.xticks(X, quarter)
    plt.tick_params(labelbottom=True)
    plt.ylabel("Million R per year")
    plt.savefig(path_plots + flood_categ + '_structures_damages.png')
    # plt.show()
    plt.close()

    plt.figure(figsize=(10, 7))
    plt.bar(
        X, data[1], color=colors[1], width=barWidth,
        label="Contents")
    plt.legend()
    # plt.ylim(0, 25)
    quarter = ["Formal private", "Formal subsidized",
               "Informal \n settlements", "Informal \n in backyards"]
    plt.xticks(X, quarter)
    plt.tick_params(labelbottom=True)
    plt.ylabel("Million R per year")
    plt.savefig(path_plots + flood_categ + '_contents_damages.png')
    # plt.show()
    plt.close()


def simul_damages_time(list_damages, path_plots, path_tables,
                       flood_categ, options):
    """
    Plot evolution of aggregate annualized damages per housing type over time.

    The function returns different plots for each housing type. Again, they
    are broken down across structures and contents damages. The value obtained
    for each year corresponds to the annualized value of the damages summed
    across all locations.

    Parameters
    ----------
    list_damages : list
        List, for each simulation year, of the output of the
        outputs.flood_outputs.compute_damages_2d function: a dictionary
        yielding, for each return period, the estimated damages per grid cell
        (24,014) and housing type (4) in terms of housing structures and
        contents
    path_plots : str
        Path for saving output plots
    path_tables : str
        Path for saving output plots
    flood_categ : str
        Type of flood risk considered, used to define output file name
    options : dict
        Dictionary of default options

    Returns
    -------
    None.

    """
    list_data = []
    for damages in list_damages:
        new_damages = {}
        for key in damages:
            new_damages[key] = damages[key].sum(axis=0)
        clean_damages = pd.DataFrame.from_dict(new_damages, orient='index')
        data = [
            [outfld.annualize_damages(
                clean_damages.formal_structure_damages, flood_categ,
                'formal', options) / 1000000,
             outfld.annualize_damages(
                 clean_damages.subsidized_structure_damages, flood_categ,
                 'subsidized', options) / 1000000,
             outfld.annualize_damages(
                 clean_damages.informal_structure_damages,
                 flood_categ, 'informal', options) / 1000000,
             outfld.annualize_damages(
                 clean_damages.backyard_structure_damages,
                 flood_categ, 'backyard', options) / 1000000],
            [outfld.annualize_damages(
                clean_damages.formal_content_damages,
                flood_categ, 'formal', options) / 1000000,
             outfld.annualize_damages(
                 clean_damages.subsidized_content_damages,
                 flood_categ, 'subsidized', options) / 1000000,
             outfld.annualize_damages(
                 clean_damages.informal_content_damages,
                 flood_categ, 'informal', options) / 1000000,
             outfld.annualize_damages(
                 clean_damages.backyard_content_damages,
                 flood_categ, 'backyard', options) / 1000000]
            ]
        list_data.append(data)

    list_data = np.array(list_data)
    np.save(path_tables + flood_categ + 'evol_damages', list_data)
    print(flood_categ + 'evol_damages table saved')

    years_simul = np.arange(2011, 2011 + 30)
    colors = ['#FF9999', '#00BFFF', '#C1FFC1', '#CAE1FF', '#FFDEAD']

    # It is best to separate housing types for visualisation

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(years_simul, list_data[:, 0, 0],
            color=colors[1], label="Structure")
    ax.plot(years_simul, list_data[:, 1, 0],
            color=colors[2], label="Contents")
    ax.set_ylim(0)
    ax.yaxis.set_major_formatter(
        mpl.ticker.StrMethodFormatter('{x:,.0f}'))
    plt.legend()
    plt.tick_params(labelbottom=True)
    plt.ylabel("Million R per year", labelpad=15)
    plt.savefig(path_plots + flood_categ + '_evol_FP_damages.png')
    plt.close()
    print(flood_categ + '_evol_FP_damages done')

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(years_simul, list_data[:, 0, 1],
            color=colors[1], label="Structure")
    ax.plot(years_simul, list_data[:, 1, 1],
            color=colors[2], label="Contents")
    ax.set_ylim(0)
    ax.yaxis.set_major_formatter(
        mpl.ticker.StrMethodFormatter('{x:,.0f}'))
    plt.legend()
    plt.tick_params(labelbottom=True)
    plt.ylabel("Million R per year)", labelpad=15)
    plt.savefig(path_plots + flood_categ + '_evol_FS_damages.png')
    plt.close()
    print(flood_categ + '_evol_FS_damages done')

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(years_simul, list_data[:, 0, 2],
            color=colors[1], label="Structure")
    ax.plot(years_simul, list_data[:, 1, 2],
            color=colors[2], label="Contents")
    ax.set_ylim(0)
    ax.yaxis.set_major_formatter(
        mpl.ticker.StrMethodFormatter('{x:,.0f}'))
    plt.legend()
    plt.tick_params(labelbottom=True)
    plt.ylabel("Million R per year", labelpad=15)
    plt.savefig(path_plots + flood_categ + '_evol_IS_damages.png')
    plt.close()
    print(flood_categ + '_evol_IS_damages done')

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(years_simul, list_data[:, 0, 3],
            color=colors[1], label="Structure")
    ax.plot(years_simul, list_data[:, 1, 3],
            color=colors[2], label="Contents")
    ax.set_ylim(0)
    ax.yaxis.set_major_formatter(
        mpl.ticker.StrMethodFormatter('{x:,.0f}'))
    plt.legend()
    plt.tick_params(labelbottom=True)
    plt.ylabel("Million R per year", labelpad=15)
    plt.savefig(path_plots + flood_categ + '_evol_IB_damages.png')
    plt.close()
    print(flood_categ + '_evol_IB_damages dones')


def round_nearest(x, a):
    """
    Return rounded value to nearest decimal number above some threshold.

    This is a technical calculation that is called by the
    plot_flood_severity_distrib function: it will allow to split values across
    bins, while discarding small values (smaller than half bar width) to avoid
    a peak at the origin of the plot.

    Parameters
    ----------
    x : Series
        Any series of numbers with many decimals
    a : float
        Level below which we want to discard below which we want to
        discard small values

    Returns
    -------
    float
        Value rounded to the nearest decimal

    """
    return round(round(x / a) * a, 2)


def plot_flood_severity_distrib(barWidth, transparency, dictio, flood_type,
                                path_plots, ylim):
    """
    Bar plot distribution of households across flood zones of varying severity.

    Only selected return periods are shown for illustrative purposes. Output
    plot is specific to a given flood type and is broken down into income
    groups as a proxy of households' vulnerability. Households are grouped into
    bins for different maximum flood depth levels in their residential location
    and the resulting population distribution is displayed without the majority
    living in low-risk zones (below the flood depth level given by barWidth
    argument), to preserve plot readability.

    Parameters
    ----------
    barWidth : float
        Bar width argument expected by pyplot: ocrresponds to the step used to
        group maximum flood depth (in m) levels into bins
    transparency : list
        List of transparency indices to distinguish across plotted return
        periods
    dictio : dict
        Dictionary yielding, for each return period of a given flood type, the
        spatial distribution of households broken into income groups, along
        with the maximum flood depth level and fraction of flood-prone area in
        their residential location
    flood_type : TYPE
        Code for flood type used as the first component of dictionary keys
    path_plots : str
        Path for saving output plots
    ylim : int
        Maximum value on the y-axis

    Returns
    -------
    None.

    """
    if flood_type == "FD":
        df_1 = dictio[flood_type + '_20yr']
        df_2 = dictio[flood_type + '_50yr']
        df_3 = dictio[flood_type + '_100yr']
    elif flood_type == "FU":
        df_1 = dictio[flood_type + '_20yr']
        df_2 = dictio[flood_type + '_50yr']
        df_3 = dictio[flood_type + '_100yr']
    if flood_type == "P":
        df_1 = dictio[flood_type + '_20yr']
        df_2 = dictio[flood_type + '_50yr']
        df_3 = dictio[flood_type + '_100yr']
    if flood_type == "C_MERITDEM_1":
        df_1 = dictio[flood_type + '_0025']
        df_2 = dictio[flood_type + '_0050']
        df_3 = dictio[flood_type + '_0100']

    df = pd.DataFrame(data=np.transpose(
        np.array([df_1.flood_depth, df_1.sim_poor, df_1.sim_midpoor,
                  df_1.sim_midrich, df_1.sim_rich,
                  df_2.flood_depth, df_2.sim_poor, df_2.sim_midpoor,
                  df_2.sim_midrich, df_2.sim_rich,
                  df_3.flood_depth, df_3.sim_poor, df_3.sim_midpoor,
                  df_3.sim_midrich, df_3.sim_rich])),
        columns=["x_1", "ypoor_1", "ymidpoor_1", "ymidrich_1", "yrich_1",
                 "x_2", "ypoor_2", "ymidpoor_2", "ymidrich_2", "yrich_2",
                 "x_3", "ypoor_3", "ymidpoor_3", "ymidrich_3", "yrich_3"])
    df["round_1"] = round_nearest(df.x_1, barWidth)
    df["round_2"] = round_nearest(df.x_2, barWidth)
    df["round_3"] = round_nearest(df.x_3, barWidth)
    new_df_1 = df[["round_1", "ypoor_1", "ymidpoor_1", "ymidrich_1", "yrich_1"]
                  ].groupby(['round_1']).sum()
    new_df_1["rounded"] = new_df_1.index
    new_df_2 = df[["round_2", "ypoor_2", "ymidpoor_2", "ymidrich_2", "yrich_2"]
                  ].groupby(['round_2']).sum()
    new_df_2["rounded"] = new_df_2.index
    new_df_3 = df[["round_3", "ypoor_3", "ymidpoor_3", "ymidrich_3", "yrich_3"]
                  ].groupby(['round_3']).sum()
    new_df_3["rounded"] = new_df_3.index

    fig, ax = plt.subplots(2, 2, figsize=(12, 12))

    ax[0, 0].bar(new_df_1.rounded[barWidth:3], new_df_1.ypoor_1[barWidth:3],
                 width=barWidth, color='royalblue', alpha=transparency[0],
                 label='20 years')
    ax[0, 0].bar(new_df_2.rounded[barWidth:3], new_df_2.ypoor_2[barWidth:3],
                 width=barWidth, color='cornflowerblue', alpha=transparency[1],
                 label='50 years')
    ax[0, 0].bar(new_df_3.rounded[barWidth:3], new_df_3.ypoor_3[barWidth:3],
                 width=barWidth, color='lightsteelblue', alpha=transparency[2],
                 label='100 years')
    ax[0, 0].set_ylabel("Households (nb)")
    ax[0, 0].set_xlabel("Severity of floods (m)")
    ax[0, 0].yaxis.set_major_formatter(
        mpl.ticker.StrMethodFormatter('{x:,.0f}'))
    ax[0, 0].set_title("Poor")
    ax[0, 0].set_ylim([0, ylim])
    ax[0, 0].legend()

    ax[0, 1].bar(new_df_1.rounded[barWidth:3], new_df_1.ymidpoor_1[barWidth:3],
                 width=barWidth, color='royalblue', alpha=transparency[0],
                 label='25 years')
    ax[0, 1].bar(new_df_2.rounded[barWidth:3], new_df_2.ymidpoor_2[barWidth:3],
                 width=barWidth, color='cornflowerblue', alpha=transparency[1],
                 label='50 years')
    ax[0, 1].bar(new_df_3.rounded[barWidth:3], new_df_3.ymidpoor_3[barWidth:3],
                 width=barWidth, color='lightsteelblue', alpha=transparency[2],
                 label='100 years')
    ax[0, 1].set_ylabel("Households (nb)")
    ax[0, 1].set_xlabel("Severity of floods (m)")
    ax[0, 1].yaxis.set_major_formatter(
        mpl.ticker.StrMethodFormatter('{x:,.0f}'))
    ax[0, 1].set_title("Mid-poor")
    ax[0, 1].set_ylim([0, ylim])
    ax[0, 1].legend()
    ax[1, 0].bar(new_df_1.rounded[barWidth:3], new_df_1.ymidrich_1[barWidth:3],
                 width=barWidth, color='royalblue', alpha=transparency[0],
                 label='25 years')

    ax[1, 0].bar(new_df_2.rounded[barWidth:3], new_df_2.ymidrich_2[barWidth:3],
                 width=barWidth, color='cornflowerblue', alpha=transparency[1],
                 label='50 years')
    ax[1, 0].bar(new_df_3.rounded[barWidth:3], new_df_3.ymidrich_3[barWidth:3],
                 width=barWidth, color='lightsteelblue', alpha=transparency[2],
                 label='100 years')
    ax[1, 0].set_ylabel("Households (nb)")
    ax[1, 0].set_xlabel("Severity of floods (m)")
    ax[1, 0].yaxis.set_major_formatter(
        mpl.ticker.StrMethodFormatter('{x:,.0f}'))
    ax[1, 0].set_title("Mid-rich")
    ax[1, 0].set_ylim([0, ylim])
    ax[1, 0].legend()
    ax[1, 1].bar(new_df_1.rounded[barWidth:3], new_df_1.yrich_1[barWidth:3],
                 width=barWidth, color='royalblue', alpha=transparency[0],
                 label='25 years')

    ax[1, 1].bar(new_df_2.rounded[barWidth:3], new_df_2.yrich_2[barWidth:3],
                 width=barWidth, color='cornflowerblue', alpha=transparency[1],
                 label='50 years')
    ax[1, 1].bar(new_df_3.rounded[barWidth:3], new_df_3.yrich_3[barWidth:3],
                 width=barWidth, color='lightsteelblue', alpha=transparency[2],
                 label='100 years')
    ax[1, 1].set_ylabel("Households (nb)")
    ax[1, 1].set_xlabel("Severity of floods (m)")
    ax[1, 1].yaxis.set_major_formatter(
        mpl.ticker.StrMethodFormatter('{x:,.0f}'))
    ax[1, 1].set_title("Rich")
    ax[1, 1].set_ylim([0, ylim])
    ax[1, 1].legend()

    plt.savefig(path_plots + flood_type + '_severity_distrib.png')
    plt.show()
    plt.close()
