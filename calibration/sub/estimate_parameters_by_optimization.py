# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 10:49:58 2020.

@author: Charlotte Liotta
"""

import numpy as np
import math
import scipy
import statsmodels.api as sm
import warnings

import calibration.sub.loglikelihood as callog


def EstimateParametersByOptimization(
        incomeNetOfCommuting, dataRent, dataDwellingSize, dataIncomeGroup,
        selectedSP, tableAmenities, variablesRegression, initRho, initBeta,
        initBasicQ, initUti2, initUti3, initUti4, options):
    """
    Estimate parameters by maximizing log likelihood via gradient descent.

    This function runs an interior-point algoithm over values of endogenous
    parameters entering households' utility function, and returns those which
    maximize a composite log-likelihood. The latter is defined as the sum of
    log-likelihoods measuring the fit of the model on several observed data
    moments, namely dwelling sizes, exogenous amenities, and income sorting.
    To compute those separate log-likelihoods, this function calls on the
    calibration.sub.loglikelihood module after providing sample and variable
    selection for the estimation, as well as theoretical partial relations
    from the structure of the model, used in regressions.
    It leverages initial guess from estimate_parameters_by_scanning module on
    parameter values to start within the set of feasible solutions, and to
    converge towards an interior (as opposed to corner) solution.
    The algorithm only converges when considering the fit along all dimensions
    (not only exogenous amenities). This function is therefore not appropriate
    when we want to focus calibration on the amenity score only.

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
    initBeta : float64
        Initial guess for the surplus housing elasticity from households'
        utility function: comes from the estimate_parameters_by_scanning module
    initBasicQ : float64
        Initial guess for the basic need in housing (in m²) from households'
        utility function: comes from the estimate_parameters_by_scanning module
    initUti2 : int32
        Target utility level of the second poorest income group
    initUti3 : float64
        Initial guess for the utility level of the second richest income group:
        comes from the estimate_parameters_by_scanning module
    initUti4 : float64
        Initial guess for the utility level of the richest income group:
        comes from the estimate_parameters_by_scanning module
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

    # Initial value of parameters (from pre-scan)
    # This allows to start the optimization within the feasible allocation
    # area, and not to get stuck in potential corner solutions
    initialVector = np.array(
        [initBeta, initBasicQ, initUti3, initUti4])

    # Determines function that will be minimized (we keep only the total
    # composite score, indexed by 0)
    minusLogLikelihoodModel = (
        lambda X0:
            - callog.LogLikelihoodModel(
                X0, initUti2, net_income, groupLivingSpMatrix,
                dataDwellingSize, selectedDwellingSize, dataRent,
                selectedRents,
                predictorsAmenitiesMatrix, tableRegression,
                variablesRegression, CalculateDwellingSize,
                ComputeLogLikelihood, optionRegression, options)[0]
            )

    # Now, we optimize using interior-point minimization algorithm

    # By default, we pin down values for beta and q0. Then, we allow utility
    # levels for the two richest income groups (making the most part of
    # formal private housing) to vary widely: this should converge towards an
    # interior solution for a good initial guess on those values
    # NB: again, this can be changed if needed
    bnds = ((initBeta, initBeta),
            (initBasicQ, initBasicQ),
            (0, 10**4),
            (0, 10**5))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # We run the algorithm and store the output
        res = scipy.optimize.minimize(
            minusLogLikelihoodModel, initialVector, bounds=bnds
            )
    parameters = res.x
    scoreTot = res.fun
    exitFlag = res.success
    print(exitFlag)

    # Estimate the function to get the associated amenity parameters
    if options["glm"] == 1:
        optionRegression = 1
        (*_, parametersAmenities, modelAmenity, parametersHousing
         ) = callog.LogLikelihoodModel(
             parameters, initUti2, net_income, groupLivingSpMatrix,
             dataDwellingSize, selectedDwellingSize, dataRent,
             selectedRents,
             predictorsAmenitiesMatrix, tableRegression, variablesRegression,
             CalculateDwellingSize, ComputeLogLikelihood, optionRegression,
             options)
    elif options["glm"] == 0:
        optionRegression = 0
        (*_, parametersAmenities, modelAmenity, parametersHousing
         ) = callog.LogLikelihoodModel(
             parameters, initUti2, net_income, groupLivingSpMatrix,
             dataDwellingSize, selectedDwellingSize, dataRent,
             selectedRents,
             predictorsAmenitiesMatrix, tableRegression, variablesRegression,
             CalculateDwellingSize, ComputeLogLikelihood, optionRegression,
             options)

    print('*** Estimation of beta and q0 done ***')

    return (parameters, scoreTot, parametersAmenities, modelAmenity,
            parametersHousing)
