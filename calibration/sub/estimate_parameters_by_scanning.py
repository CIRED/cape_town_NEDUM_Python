# -*- coding: utf-8 -*-
"""



"""

import numpy as np
import math
import statsmodels.api as sm

import calibration.sub.loglikelihood as callog


def EstimateParametersByScanning(incomeNetOfCommuting, dataRent,
                                 dataDwellingSize, dataIncomeGroup,
                                 selectedSP, tableAmenities,
                                 variablesRegression, initRho, listBeta,
                                 listBasicQ, initUti2, listUti3, listUti4,
                                 options):
    """
    Estimate parameters by maximizing log likelihood over scanned values.

    This function scans over values of endogenous parameters entering
    households' utility function, and returns those which maximize a
    composite log-likelihood. The latter is defined as the sum of
    log-likelihoods measuring the fit of the model on several observed data
    moments, namely dwelling sizes, exogenous amenities, and income sorting.
    To compute those separate log-likelihoods, this function calls on the
    calibration.sub.loglikelihood module after providing sample and variable
    selection for the estimation, as well as theoretical partial relations
    from the structure of the model, used in regressions.
    The purpose of this function is to provide an informed guess on initial
    values for the estimate_parameters_by_optimization module to run a proper
    gradient descent optimization, by starting within the set of interior
    feasible solutions.
    By default, we only fit the data moment on exogenous amenities (and forget
    about the other dimensions of the composite likelihood) when beta and q0
    are pinned down, as it is the only relevant moment for the amenity score
    calibration. The other dimensions can be recovered by commenting out the
    cancelling out terms at the end of the script, if we want to recalibrate
    beta and q0 internally (and potentially improve the fit of the model at
    the price of empirical validity).

    Parameters
    ----------
    incomeNetOfCommuting : ndarray(float64, ndim=2)
        Expected annual income net of commuting costs (in rands, for
        one household), for SP (1,046), by income group (4)
    dataRent : Series
        Theoretical average annual rent for formal private housing, computed
        from data on average land prices, for each SP (1,046)
    dataDwellingSize : Series
        Average dwelling size (in m²) for each SP (1,046), from SP data
    dataIncomeGroup : ndarray(float64)
        Categorical variable indicating, for each Small Place (1,046), the
        dominant income group (from 0 to 3)
    selectedSP : Series
        Dummy variable used to select SPs (1,046) with enough formal private
        housing to identify the regressions used in the function (less
        stringent than selectedDensity)
    tableAmenities : DataFrame
        Table yielding, for each selected geographic unit, a set of dummy
        variables corresponding to available exogenous amenites
    variablesRegression : list
        List containing labels for exogenous amenity variables used in the
        final regressions
    initRho : int
        Spatial autocorrelation parameter (not used in practice, as is set to
        zero)
    listBeta : ndarray(float64)
        List of values over which to scan for the surplus housing elasticity
        parameter from households' utility function
    listBasicQ : ndarray(float64)
        List of values over which to scan for the basic need in housing (in m²)
        from households' utility function
    initUti2 : int32
        Target utility level of the second poorest income group
    listUti3 : ndarray(float64)
        List of values over which to scan for the utility level of the second
        richest income group
    listUti4 : ndarray(float64)
        List of values over which to scan for the utility level of the richest
        income group
    options : dict
        Dictionary of default options

    Returns
    -------
    parameters : ndarray(float64)
        Vector of "calibrated parameters". In the order: surplus housing
        elasticity, basic need in housing, and utility levels for income groups
        3 and 4
    scoreTot : float64
        Maximum value of composite log-likelihood for the fit on observed
        amenities, dwelling sizes, and population sorting by income: provides
        a performance metric for the calibration process
    parametersAmenities : ndarray(float64)
        List of estimates for the impact of exogenous (dummy) amenities on the
        calibrated amenity index (in log-form). In the order: intercept,
        distance to the ocean <2km, distance to the ocean between 2 and 4km,
        slope between 1 and 5%, slope >5%, being located within the airport
        cone, distance to district parks <2km, distance to biosphere reserve
        <2km, distance to train station <2km, distance to urban heritage site
        <2km
    modelAmenities : regression.linear_model.RegressionResultsWrapper
        Object summarizing the results of the log-regressions of the
        theoretical amenity index over observed exogenous amenity dummies
    parametersHousing : int
        List of estimates related to the fit of the model on building density
        / housing supply: this is not included in this version of the model
        (vector is set equal to zero) as we already exploit this relation to
        estimate parameters of the construction function (compared to other
        versions where we use construction costs)

    """
    # SAMPLE AND VARIABLE SELECTION

    # We remove poorest income group as it is crowded out of formal private
    # housing sector in practice
    net_income = incomeNetOfCommuting[1:4, :]
    net_income[net_income < 0] = np.nan

    # We generate a matrix of dummies for dominant income group (with positive
    # net income) in each SP
    groupLivingSpMatrix = (net_income > 0)
    for i in range(0, 3):
        groupLivingSpMatrix[i, dataIncomeGroup != i+1] = np.zeros(1, 'bool')

    # We define a set of selection arrays that will serve in regressions for
    # model fit over different variables

    # For the fit on observed income sorting (highest bid-rents), we select SPs
    # with enough formal private housing (selectedSP), where values are well
    # defined (~np.isnan(dataRent))
    selectedRents = selectedSP & ~np.isnan(dataRent)

    # For the fit on dwelling size, we select SPs with enough formal private
    # housing (selectedSP), where all values (entering the regression) are well
    # defined (~np.isnan(dataDwellingSize) and ~np.isnan(dataRent))
    selectedDwellingSize = (selectedSP
                            & ~np.isnan(dataDwellingSize)
                            & ~np.isnan(dataRent))

    # We also select variables for the fit on observed amenities.
    # Note that we also use sample selection from selectedRents vector as
    # rents do enter the theoretical formula for the amenity index
    tableRegression = tableAmenities.loc[selectedRents, :]
    predictorsAmenitiesMatrix = tableRegression.loc[:, variablesRegression]
    predictorsAmenitiesMatrix = sm.add_constant(predictorsAmenitiesMatrix)
    # predictorsAmenitiesMatrix = np.vstack(
    #     [np.ones(predictorsAmenitiesMatrix.shape[0]),
    #      predictorsAmenitiesMatrix.T]
    #     ).T

    # %% FUNCTIONS USED FOR THE FIT ON DWELLING SIZE AND AMENITIES

    # Relationship between rents and dwelling sizes (see technical
    # documentation)
    CalculateDwellingSize = (
        lambda beta, basic_q, incomeTemp, rentTemp:
            beta * incomeTemp / rentTemp + (1 - beta) * basic_q
            )

    # Log likelihood for a lognormal law of mean 0: we will assume that
    # dwelling size and amenity residuals follow such a law (see math
    # appendix)
    ComputeLogLikelihood = (
        lambda sigma, error:
            np.nansum(- np.log(2 * math.pi * sigma ** 2) / 2
                      - 1 / (2 * sigma ** 2) * (error) ** 2)
            )

    # %% OPTIMIZATION ALGORITHM

    # We decide whether we want to use GLM or OLS estimation for the fit on
    # exogenous amenities.
    # Note that GLM is not stable in this version of the code (default option
    # set as zero)
    optionRegression = options["glm"]

    # Initial value of parameters (all possible combinations)
    # Note that we do not consider spatial autocorrelation
    combinationInputs = np.array(
        np.meshgrid(listBeta, listBasicQ, listUti3, listUti4)).T.reshape(-1, 4)

    # Initial score (log-likelihood) values for each regression
    # Note that fit on housing supply / building density is included, as the
    # estimation (where it will be cancelled out) is done in the
    # calibration.sub.loglikelihood module.
    scoreAmenities = - 10000 * np.ones(combinationInputs.shape[0])
    scoreDwellingSize = - 10000 * np.ones(combinationInputs.shape[0])
    scoreIncomeSorting = - 10000 * np.ones(combinationInputs.shape[0])
    scoreHousing = - 10000 * np.ones(combinationInputs.shape[0])
    scoreTotal = - 10000 * np.ones(combinationInputs.shape[0])

    print('\nStart scanning')
    print('\n')

    # We update the scores for each set of parameters
    # Note that
    for index in range(0, combinationInputs.shape[0]):
        # print(index)
        (scoreTotal[index], scoreAmenities[index], scoreDwellingSize[index],
         scoreIncomeSorting[index], scoreHousing[index], parametersAmenities,
         modelAmenities, parametersHousing) = callog.LogLikelihoodModel(
             combinationInputs[index, :], initUti2, net_income,
             groupLivingSpMatrix, dataDwellingSize, selectedDwellingSize,
             dataRent, selectedRents, predictorsAmenitiesMatrix,
             tableRegression, variablesRegression, CalculateDwellingSize,
             ComputeLogLikelihood, optionRegression, options)

    print('\nScanning complete')
    print('\n')

    # We just pick the parameters associated to the maximum score
    # Note that scoreHousing is set as zero and does not impact parameter
    # selection.

    # By default, we cancel out the scores associated with observed dwelling
    # sizes and income sorting: this is because the fit on exogenous amenities
    # is the only relevant dimension when beta and q0 are pinned down and the
    # amenity score is the only parameter left to calibrate
    scoreDwellingSize = 0
    scoreIncomeSorting = 0

    scoreVect = (scoreAmenities + scoreDwellingSize + scoreIncomeSorting
                 + scoreHousing)
    scoreTot = np.amax(scoreVect)
    which = np.argmax(scoreVect)
    parameters = combinationInputs[which, :]

    return (parameters, scoreTot, parametersAmenities, modelAmenities,
            parametersHousing)
