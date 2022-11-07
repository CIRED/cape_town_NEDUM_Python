# -*- coding: utf-8 -*-
"""



"""

import numpy as np
import statsmodels.api as sm


def LogLikelihoodModel(X0, Uo2, net_income, groupLivingSpMatrix,
                       dataDwellingSize,
                       selectedDwellingSize, dataRent, selectedRents,
                       predictorsAmenitiesMatrix,
                       tableRegression, variables_regression,
                       CalculateDwellingSize, ComputeLogLikelihood,
                       optionRegression, options):
    """
    Compute composite log-likelihood for model fit given scanned parameters.

    This function computes three separate log-likelihoods that capture the fit
    of the model along distinct data moments. Then, it sums them to give a
    composite log-likelihood that will be maximized to calibrate the utility
    function parameters (as part of scanning or smooth optimization process).
    More precisely, it starts by regressing theoretical values of the log
    amenity index on observed exogenous dummy amenity variables. The first
    log-likelihood is computed based on the value of residuals. Then, it
    defines a second log-likelihood for the fit on income sorting (matching
    the observed dominant income group with the highest bidding income group).
    Finally, it gets the residuals from the log-difference between theoretical
    and observed dwelling sizes, and computes the third log-likelihood from
    there.

    Parameters
    ----------
    X0 : ndarray(float64)
        Set of inputs (namely, surplus housing elasticity, basic need in
        housing, utility levels for income groups 3 and 4) tested for a given
        iteration
    Uo2 : int32
        Target utility level for income group 2
    net_income : ndarray(float64, ndim=2)
        Expected annual income net of commuting costs (in rands, for
        one household), for SP (1,046), by income group excluding the poorest
        (3)
    groupLivingSpMatrix : ndarray(bool, ndim=2)
        Dummy variable indicating, for each income group excluding the poorest
        (3) whether it is dominant in each SP (1,046)
    dataDwellingSize : Series
        Average dwelling size (in m²) for each SP (1,046), from SP data
    selectedDwellingSize : Series
        Dummy variable indicating, for each SP (1,046), whether it is selected
        into the sample used when regressing observed dwelling sizes on their
        theoretical equivalent
    dataRent : Series
        Theoretical average annual rent for formal private housing, computed
        from data on average land prices, for each SP (1,046)
    selectedRents : Series
        Dummy variable indicating, for each SP (1,046), whether it is selected
        into the sample used when estimating the discrete choice logit model
        associated with income sorting (identifying observed rents with highest
        bid-rents from the dominant income group)
    predictorsAmenitiesMatrix : ndarray(float64, ndim=2)
        Values of selected exogenous dummy amenity variables (10, including the
        intercept) in each selected SP (according to selectedRents)
    tableRegression : DataFrame
        Values of all exogenous dummy amenity variables (16, excluding the
        intercept) in each selected SP (according to selectedRents)
    variables_regression : list
        List of labels for selected exogenous dummy amenity variables (9,
        excluding the intercept)
    CalculateDwellingSize : function
        Function defining the relationship between rents and dwelling sizes
        in the formal private sector (see technical documentation)
    ComputeLogLikelihood : function
        Log-likelihood function for a lognormal law of mean 0: we will assume
        that dwelling size and amenity residuals follow such a law
    optionRegression : int
        Option to run GLM (instead of OLS) regression for the estimation of
        exogenous amenity estimates: default is set as zero as GLM is unstable
    options : dict
        Dictionary of default options

    Returns
    -------
    scoreTotal : float64
        Value of the composite log-likelihood for the set of parameters scanned
    scoreAmenities : float64
        Value of the log-likelihood for the fit on exogenous amenities
    scoreDwellingSize : float64
        Value of the log-likelihood for the fit on observed dwelling sizes
    scoreIncomeSorting : float64
        Value of the log-likelihood for the fit on observed income sorting
        (matching observed rents to highest bid-rents from dominant income
        group)
    scoreHousing : float64
        Value of the log-likelihood for the fit on observed housing supply
        / building density: this is not used in this version of the model as
        the relation is already used for the calibration of construction
        function parameters (hence is set equal to zero)
    parametersAmenities : ndarray(float64)
        List of estimates for the impact of exogenous (dummy) amenities on the
        calibrated amenity index (in log-form). In the order: distance to the
        ocean <2km, distance to the ocean between 2 and 4km, slope between
        1 and 5%, slope >5%, being located within the airport cone, distance
        to district parks <2km, distance to biosphere reserve <2km, distance
        to train station <2km, distance to urban heritage site <2km
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
    # We extract scanned parameters from input vector
    beta = X0[0]
    basicQ = X0[1]
    Uo = np.array([Uo2, X0[2], X0[3]])

    # %% Error on the amenity index

    # Theoretical formula for log amenity index (from technical documentation)
    # Note that we stay under sample selection from selectedRents variable
    logAmenityIndex = (
        np.log(np.array(Uo)[:, None])
        - np.log(
            (1 - beta) ** (1 - beta) * beta ** beta
            * (net_income[:, selectedRents]
               - basicQ * np.array(dataRent)[None, selectedRents])
            / (np.array(dataRent)[None, selectedRents] ** beta))
        )
    # We select values for dominant income groups only and flatten the array:
    # this allows to select the appropriate net income and utility to identify
    # the regression
    logAmenityIndex = np.nansum(
        logAmenityIndex * groupLivingSpMatrix[:, selectedRents], 0)
    # logAmenityIndex[np.abs(logAmenityIndex.imag) > 0] = np.nan
    # logAmenityIndex[logAmenityIndex == -np.inf] = np.nan

    # We get residuals for the regression of log amenity index on exogenous
    # dummy variables.
    # Note that we identify dummies with logs of original amenity values
    # (normalized to take values between 1 and e): see technical documentation

    # OLS estimation
    if (optionRegression == 0):
        A = predictorsAmenitiesMatrix  # [~np.isnan(logAmenityIndex), :]
        y = (logAmenityIndex[~np.isnan(logAmenityIndex)]).real
        # parametersAmenities, residuals, rank, s = np.linalg.lstsq(A, y,
        #                                                           rcond=None)
        # res = scipy.optimize.lsq_linear(A, y)
        # parametersAmenities = res.x
        # residuals = res.fun
        modelSpecification = sm.OLS(y, A, missing='drop')
        modelAmenities = modelSpecification.fit()
        parametersAmenities = modelAmenities.params
        errorAmenities = modelAmenities.resid_pearson
        # errorAmenities = y - np.nansum(A * parametersAmenities, 1)
        # modelAmenities = 0

    # GLM estimation (unstable)
    elif (optionRegression == 1):
        residu = logAmenityIndex.real
        A = predictorsAmenitiesMatrix[~np.isnan(logAmenityIndex), :]
        y = (logAmenityIndex[~np.isnan(logAmenityIndex)]).real
        parametersAmenities, residuals, rank, s = np.linalg.lstsq(A, y,
                                                                  rcond=None)
        modelSpecification = sm.GLM(
            residu, tableRegression.loc[:, variables_regression])
        modelAmenities = modelSpecification.fit()
        errorAmenities = modelAmenities.resid_pearson

    # Residuals follow a log-normal law, hence the associated log-likelihood
    # from ComputeLogLikelihood function (see math appendix)
    scoreAmenities = ComputeLogLikelihood(
        np.sqrt(np.nansum(errorAmenities ** 2)
                / np.nansum(~np.isnan(errorAmenities))),
        errorAmenities)

    # %% Error on income sorting

    # We start by defining input variables
    utility_over_amenity = Uo[:, None] / np.exp(logAmenityIndex[None, :])
    net_income_select = net_income[:, selectedRents]

    # We then apply formula from technical documentation
    bidRents = (
        (net_income_select * (1 - beta) ** (1 - beta) * beta ** beta)
        / (utility_over_amenity
            - basicQ * (1 - beta) ** (1 - beta) * beta ** beta)
        )

    # We select SPs where at least some income group makes a positive bid
    selectedBidRents = (np.nansum(bidRents, 0) > 0)
    # We again select the dominant income group
    incomeGroupSelectedRents = groupLivingSpMatrix[:, selectedRents]

    # We create a function for the minus log-likelihood of income sorting:
    # see technical documentation for math formula. We consider minus form
    # to be able to use scipy.optimize.minimize() function. We will then return
    # minus form for the minimum obtained to recover the maximum value of the
    # log-likelihood (across scale factors).

    minusLogLikIncomeSorting = (
        lambda scaleParam:
            - np.nansum(np.nansum(
                bidRents[:, selectedBidRents] / scaleParam
                * incomeGroupSelectedRents[:, selectedBidRents], 0))
            + np.nansum(np.log(np.nansum(
                np.exp(bidRents[:, selectedBidRents] / scaleParam),
                0)))
            )

    # Actually, since we cannot identify the gravity parameter, we pin it at
    # 10, which allows to get a score of the same order of magnitude as for
    # the fit on amenities. It is also a realistic and conservative guess:
    # the order of magnitude is the same as for the gravity parameter from the
    # commuting choice model, and it is low enough (given the function
    # definition, the higher the parameter, the smaller the output).
    # Also note that the improvement is marginal for higher orders of magnitude

    # We keep the optimization in comments for reference, in case we want to be
    # less conservative about the value of the log-likelihood

    scoreIncomeSorting = - minusLogLikIncomeSorting(10)

    # bnds = {(0, 10**10)}
    # initScale = 10**5
    # res = scipy.optimize.minimize(
    #     minusLogLikIncomeSorting, initScale, bounds=bnds,
    #     options={'maxiter': 100, 'disp': False})
    # optiminusLogLikIncomeSorting = res.fun
    # exitFlag = res.success
    # print(exitFlag)

    # %% Error on dwelling sizes

    # We get theoretical values for market rents that we identify to highest
    # bid-rents from dominant income group
    # NB: As bid rents are already defined for sample of SPs selected under
    # selectedRents, we need to take the subsection
    # selectedDwellingSize[selectedRents] to operate another selection based on
    # selectedDwellingSize while respecting array dimensions
    simulatedRents = np.nansum(
        bidRents[:, selectedDwellingSize[selectedRents]]
        * groupLivingSpMatrix[:, selectedDwellingSize],
        0)

    # We get theoretical values for dwelling size based on pre-defined function
    # (see technical documentation for math formula), leveraging the above
    # definition of simulatedRents
    dwellingSize = CalculateDwellingSize(
        beta,
        basicQ,
        np.nansum(net_income[:, selectedDwellingSize]
                  * groupLivingSpMatrix[:, selectedDwellingSize], 0),
        simulatedRents)

    # The (log) error on observed dwelling sizes is directly obtain by taking
    # the log-difference with the theoretical counterpart (staying under the
    # same sample selection)
    errorDwellingSize = (
        np.log(dwellingSize)
        - np.log(dataDwellingSize[selectedDwellingSize])
        )

    # Residuals follow a log-normal law, hence the associated log-likelihood
    # from ComputeLogLikelihood function (see math appendix)
    scoreDwellingSize = ComputeLogLikelihood(
        np.sqrt(np.nansum(errorDwellingSize ** 2)
                / np.nansum(~np.isnan(errorDwellingSize))),
        errorDwellingSize)

    # %% Total

    # The sum of logs is the same as the log of a product, hence we can define
    # our composite log-likelihood function as the sum of our separate
    # log-likelihoods

    scoreTotal = scoreAmenities + scoreDwellingSize + scoreIncomeSorting
    # scoreTotal = scoreAmenities

    # We may also include a measure of the fit for housing supply / household
    # density, which has not been retained in this version of the model, as
    # the underlying relation is already used to calibrate parameters of the
    # construction function.
    # NB: We still define the variables that we set equal to zero to serve as
    # placeholders
    scoreHousing = 0
    parametersHousing = 0

    return (scoreTotal, scoreAmenities, scoreDwellingSize, scoreIncomeSorting,
            scoreHousing, parametersAmenities, modelAmenities,
            parametersHousing)
