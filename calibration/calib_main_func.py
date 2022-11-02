# -*- coding: utf-8 -*-

import numpy as np
import numpy.matlib
import statsmodels.api as sm
import os
import math
import scipy

import calibration.sub.compute_income as calcmp
import calibration.sub.import_employment_data as calemp
import calibration.sub.estimate_parameters_by_scanning as calscan
import calibration.sub.estimate_parameters_by_optimization as calopt
import calibration.sub.import_amenities as calam


def estim_construct_func_param(options, param, data_sp,
                               threshold_income_distribution,
                               income_distribution, data_rdp, housing_types_sp,
                               data_number_formal, data_income_group,
                               selected_density,
                               path_data, path_precalc_inp, path_folder):
    """
    Estimate coefficients of housing production function (Cobb-Douglas).

    This function leverages a partial relation from our general equilibrium
    model, that is estimated on SP data which does not enter the simulations
    as an input. More precisely, it combines the expression of the optimal
    housing supply in the formal private sector with the highest bid condition
    (see technical documentation for math formulas).

    Parameters
    ----------
    options : dict
        Dictionary of default options
    param : dict
        Dictionary of default parameters
    data_sp : DataFrame
        Table yielding, for each Small Place (1,046), the average dwelling size
        (in m²), the average land price and annual income level (in rands),
        the size of unconstrained area for construction (in m²), the total area
        (in km²), the distance to the city centre (in km), whether or not the
        location belongs to Mitchells Plain, and the SP code
    threshold_income_distribution : ndarray(int32)
        Annual income level (in rands) above which a household is taken as
        being part of one of the 4 income groups in the model
    income_distribution : ndarray(uint16, ndim=2)
        Exogenous number of households in each Small Place (1,046) for each
        income group in the model (4)
    data_rdp : DataFrame
        Table yielding, for each grid cell (24,014), the associated cumulative
        count of cells with some formal subsidized housing, and the associated
        area (in m²) dedicated to such housing
    housing_types_sp : DataFrame
        Table yielding, for each Small Place (1,046), the number of informal
        backyards, of informal settlements, and total dwelling units, as well
        as their (centroid) x and y coordinates
    data_number_formal : Series
        Number of formal private housing units considered for each Small Place
        (1,046)
    data_income_group : ndarray(float64)
        Categorical variable indicating, for each Small Place (1,046), the
        dominant income group (from 0 to 3)
    selected_density : Series
        Dummy variable allowing for sample selection across Small Places
        (1,046) for regressions that are only valid in the formal private
        housing sector
    path_data : str
        Path towards data used in the model
    path_precalc_inp : str
        Path for precalcuted input data (calibrated parameters)
    path_folder : str
        Path towards the root data folder

    Returns
    -------
    coeff_b : float64
        Calibrated capital elasticity in housing production function
    coeff_a : float64
        Calibrated land elasticity in housing production function
    coeffKappa : float64
        Calibrated scale factor in housing production function

    """
    # We define our outcome variable
    y = np.log(data_number_formal[selected_density])
    # We define our independent variables
    # Note that we use data_sp["unconstrained_area"] (which is accurate data at
    # SP level) rather than coeff_land (which is an estimate at grid level)
    X = np.transpose(
        np.array([np.ones(len(data_sp["price"][selected_density])),
                  np.log(data_sp["price"][selected_density]),
                  np.log(data_sp["dwelling_size"][selected_density]),
                  np.log(param["max_land_use"]
                         * data_sp["unconstrained_area"][selected_density])])
        )
    # NB: Our data set for dwelling sizes only provides the average dwelling
    # size at the Sub-Place level, aggregating formal and informal housing

    # We run the linear regression
    modelSpecification = sm.OLS(y, X, missing='drop')
    model_construction = modelSpecification.fit()
    print(model_construction.summary())
    parametersConstruction = model_construction.params

    # We export outputs of the model
    coeff_b = parametersConstruction["x1"]
    coeff_a = 1 - coeff_b
    # Scale factor formula comes from zero profit condition combined with
    # footnote 16 from Pfeiffer et al. (typo in original paper)
    if options["correct_kappa"] == 1:
        coeffKappa = ((1 / (coeff_b / coeff_a) ** coeff_b)
                      * np.exp(parametersConstruction["const"]))
    elif options["correct_kappa"] == 0:
        coeffKappa = ((1 / (coeff_b) ** coeff_b)
                      * np.exp(parametersConstruction["const"]))

    try:
        os.mkdir(path_precalc_inp)
    except OSError as error:
        print(error)

    np.save(path_precalc_inp + 'calibratedHousing_b.npy', coeff_b)
    np.save(path_precalc_inp + 'calibratedHousing_kappa.npy', coeffKappa)

    return coeff_b, coeff_a, coeffKappa


def estim_incomes_and_gravity(param, grid, list_lambda,
                              households_per_income_class,
                              average_income, income_distribution,
                              spline_inflation, spline_fuel,
                              spline_population_income_distribution,
                              spline_income_distribution,
                              path_data, path_precalc_inp,
                              path_precalc_transp, options):
    """
    Estimate incomes per job center and income group, with gravity parameter.

    This function leverages theoretical formulas from
    calibration.sub.compute_income and data imported through the
    calibration.sub.import_employment_data module. It first import transport
    costs and observed number of commuters per selected job centre and income
    group, then estimates the associated incomes for a given gravity parameter
    by minimizing the error over the simulated number of commuters. Then, it
    selects among a list of scanned values the final value of the gravity
    parameter (and the associated incomes) by minimizing the error over the
    distribution of commuters along their residence-workplace distances.

    Parameters
    ----------
    param : dict
        Dictionary of default parameters
    grid : DataFrame
        Table yielding, for each grid cell (24,014), its x and y
        (centroid) coordinates, and its distance (in km) to the city centre
    list_lambda : ndarray(float64)
        List of values over which to scan for the gravity parameter used in
        the commuting choice model
    households_per_income_class : ndarray(float64)
        Exogenous total number of households per income group (excluding people
        out of employment, for 4 groups)
    average_income : ndarray(float64)
        Average median income for each income group in the model (4)
    income_distribution : ndarray(uint16, ndim=2)
        Exogenous number of households in each Small Place (1,046) for each
        income group in the model (4)
    spline_inflation : interp1d
        Linear interpolation for inflation rate (in base 100 relative to
        baseline year) over the years (baseline year set at 0)
    spline_fuel : interp1d
        Linear interpolation for fuel price (in rands per km)
        over the years (baseline year set at 0)
    spline_population_income_distribution : interp1d
        Linear interpolation for total population per income group in the data
        (12) over the years (baseline year set at 0)
    spline_income_distribution : interp1d
        Linear interpolation for median annual income (in rands) per income
        group in the data (12) over the years (baseline year set at 0)
    path_data : str
        Path towards data used in the model
    path_precalc_inp : str
        Path for precalcuted input data (calibrated parameters)
    path_precalc_transp : str
        Path for precalcuted transport inputs (intermediate outputs from
        commuting choice model)
    options : dict
        Dictionary of default options

    Returns
    -------
    incomeCentersKeep : ndarray(float64, ndim=2)
        Calibrated average annual household income (including unemployment)
        for each income group (4), per grid cell (24,014)
    lambdaKeep : float64
        Calibrated gravity parameter from the commuting choice model
    cal_avg_income : ndarray(float64)
        Overall calibrated average income across income groups (4), for
        validation only
    scoreKeep : ndarray(float64)
        Mean ratio of simulated over observed number of commuters per job
        center (185), collapsed at the income-group level (4): captures the
        error over our calibrated parameters
    bhattacharyyaDistances : ndarray(float64)
        Bhattacharyya distances (measure the similarity of two probability
        distributions) between the calculated distribution of commuting
        distances and aggregates from the Transport Survey, for each scanned
        value of the gravity parameter. This is used as an auxiliary measure
        to pin down a unique gravity parameter (and associated matrix of
        incomes), and is only given as an output of the function for reference.

    """
    # We import number of workers in each selected job center.
    # Note that it is rescaled to match aggregate income distribution in census
    job_centers = calemp.import_employment_data(
        households_per_income_class, param, path_data)

    # We import transport cost data.
    # Note that we reason at the SP level here. Also note that we are
    # considering round trips and households made up of two people.
    (timeOutput, distanceOutput, monetaryCost, costTime
     ) = calcmp.import_transport_costs(
         grid, param, 0, households_per_income_class,
         spline_inflation, spline_fuel, spline_population_income_distribution,
         spline_income_distribution,
         path_precalc_inp, path_precalc_transp, 'SP', options)

    # Note that this is long to run.
    # Here again, we are considering rescaled income data.
    (incomeCenters, distanceDistribution, scoreMatrix
     ) = calcmp.EstimateIncome(
         param, timeOutput, distanceOutput[:, :, 0], monetaryCost, costTime,
         job_centers, average_income, income_distribution, list_lambda,
         options)

    # Gives aggregate statistics for % of commuters per distance bracket
    # NB: bracketsDistance = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 200])
    # with floor residence-workplace distances in km
    # (see calibration.sub.compute_income)

    # NB: we could check differential mobility patterns as a robustness check
    # by estimating a distinct gravity parameter for each income group if we
    # get access to a breakdown of residence-workplace distance distribution
    # per income group

    # Residence-workplace distance distribution comes from the CoCT's 2013
    # Transport Survey
    data_distance_distribution = np.array(
        [45.6174222, 18.9010734, 14.9972971, 9.6725616, 5.9425438, 2.5368754,
         0.9267125, 0.3591011, 1.0464129])

    # Compute accessibility index
    # NB1: Bhattacharyya distance measures the similarity of two probability
    # distributions (here, data vs. simulated % of commuters)
    # NB2: Mahalanobis distance is a particular case of the Bhattacharyya
    # distance when the standard deviations of the two classes are the same
    bhattacharyyaDistances = (
        - np.log(np.nansum(np.sqrt(data_distance_distribution[:, None]
                                   / 100 * distanceDistribution), 0))
        )
    whichLambda = np.argmin(bhattacharyyaDistances)

    # Hence, we keep the lambda that minimizes the distance and the associated
    # income vector
    lambdaKeep = list_lambda[whichLambda]
    incomeCentersKeep = incomeCenters[:, :, whichLambda]

    # We also keep the associated error metric
    # scoreKeep = scoreMatrix[whichLambda, :]
    scoreKeep = scoreMatrix

    # Note that income is set to -inf for job centers and income groups in
    # which it could not be calibrated

    np.save(path_precalc_inp + 'incomeCentersKeep.npy', incomeCentersKeep)
    np.save(path_precalc_inp + 'lambdaKeep.npy', lambdaKeep)

    # Note that it is unclear whether "average" income from data includes
    # unemployment or not: a priori, it does for short spells (less than one
    # year) and should therefore be slightly bigger than calibrated income
    # (which should reflect all unemployment): this is what we observe in
    # practice
    incomeCentersKeep[incomeCentersKeep < 0] = math.nan
    cal_avg_income = np.nanmean(incomeCentersKeep, 0)

    return (incomeCentersKeep, lambdaKeep, cal_avg_income, scoreKeep,
            bhattacharyyaDistances)


def estim_util_func_param(data_number_formal, data_income_group,
                          housing_types_sp, data_sp,
                          coeff_a, coeff_b, coeffKappa, interest_rate,
                          incomeNetOfCommuting,
                          path_data, path_precalc_inp,
                          options, param):
    """
    Calibrate utility function parameters.

    This function leverages the following modules: import_amenities,
    estimate_parameters_by_scanning, and estimate_parameters_by_optimization.
    As before, we use partial relations coming from our general equilibrium
    structure (see technical documentation for math formulas). This time, we
    look at the utility function parameters that maximize a composite
    likelihood function for the fit on observed amenities, dwelling sizes, and
    population sorting by income (see calibration.sub.loglikelihood module).
    We proceed first by scanning over a discrete range of parameter values,
    then by running a smooth solver taking outputs from scanning as initial
    values.

    Parameters
    ----------
    data_number_formal : Series
        Number of formal private housing units considered for each Small Place
        (1,046)
    data_income_group : ndarray(float64)
        Categorical variable indicating, for each Small Place (1,046), the
        dominant income group (from 0 to 3)
    housing_types_sp : DataFrame
        Table yielding, for each Small Place (1,046), the number of informal
        backyards, of informal settlements, and total dwelling units, as well
        as their (centroid) x and y coordinates
    data_sp : DataFrame
        Table yielding, for each Small Place (1,046), the average dwelling size
        (in m²), the average land price and annual income level (in rands),
        the size of unconstrained area for construction (in m²), the total area
        (in km²), the distance to the city centre (in km), whether or not the
        location belongs to Mitchells Plain, and the SP code
    coeff_a : float64
        Calibrated land elasticity in housing production function
    coeff_b : float64
        Calibrated capital elasticity in housing production function
    coeffKappa : float64
        Calibrated scale factor in housing production function
    interest_rate : float64
        Interest rate for the overall economy, corresponding to an average
        over past years
    incomeNetOfCommuting : ndarray(float64, ndim=2)
        Expected annual income net of commuting costs (in rands, for
        one household), for each geographic unit, by income group (4)
    path_data : str
        Path towards data used in the model
    path_precalc_inp : str
        Path for precalcuted input data (calibrated parameters)
    options : dict
        Dictionary of default options
    param : dict
        Dictionary of default parameters

    Returns
    -------
    calibratedUtility_beta : float64
        Calibrated surplus housing elasticity in households' utility function
    calibratedUtility_q0 : float64
        Parametric basic need in housing (in m²). Note that this parameter is
        not an output of the calibration per se, as it is exogenously set. It
        is included here for reference as it enters the households' utility
        function and enters the optimization programme as an input. Note that
        this could be optimized over (as in Pfeiffer et al.), but only within
        a narrow range of values to preserve feasibilty of allocations.
    cal_amenities : ndarray(float64)
        Calibrated amenity index for each grid cell (24,014): this is not
        normalized yet.

    """
    # We select in which areas we actually measure the likelihood.
    # We use a less stringent definition as for the fit on building density
    # (estimation of construction function parameters), as we are going to
    # separately add conditions more specific to each regression used in the
    # process.
    # NB: There is a trade-off between empirical validity and statistical power
    selectedSP = (
        (data_number_formal > 0.90 * housing_types_sp.total_dwellings_SP_2011)
        & (data_income_group > 0)
        )

    # We are going to scan over values on which to optimize

    # As explained in the parameters_and_options module, we consider that
    # benchmark value from Finlay and Williams (2022) is more robust than
    # our own estimation for the elasticity of housing demand. We therefore
    # pin this value down as default. Note that we can define a range of
    # acceptable values for sensitivity checks.

    listBeta = scipy.io.loadmat(
                path_precalc_inp + 'calibratedUtility_beta.mat'
                )["calibratedUtility_beta"].squeeze()
    # listBeta = 0.25
    # listBeta = np.arange(0.1, 0.41, 0.02)

    # Again, we pin the value of basic need in housing to the one estimated
    # in Pfeiffer et al., to serve as a placeholder for later empirical
    # calibration (see parameters_and_options module) and to reduce the
    # numerical complexity of the optimization process (that may converge
    # towards multiple optima)

    listBasicQ = scipy.io.loadmat(
                path_precalc_inp + 'calibratedUtility_q0.mat'
                )["calibratedUtility_q0"].squeeze()
    # listBasicQ = 4

    # Coefficient for spatial autocorrelation (not used)
    listRho = 0

    # Target utility levels used to define amenity score: we take levels close
    # to what we expect in equilibrium
    # utilityTarget = np.array([400, 900, 5400, 26000])
    utilityTarget = np.array([1200, 4800, 16000, 77000])

    # Then, we only allow utilities for the two richest income groups to vary,
    # for the sake of numerical simplicity and as they will drive most of the
    # changes in absolute values. Actually, we even exclude the poorest income
    # group from the analysis, again for the sake of numerical simplicity and
    # as it will in practice be crowded out of the formal private sector
    # (which drives most of our identification).

    # Compared to Pfeiffer et al., we only fit the data moment on exogenous
    # amenities when pinning down the values of beta and q0 (see definition
    # of the composite log-likelihood in the estimate_parameters modules).
    # This is because the amenity score is the only parameter left to
    # calibrate, and we want it to be predicted in the most acurrate way
    # by our explanatory covariates. Consequently, the selected values for the
    # utility levels will not necessarily match the endogenous values we get as
    # an equilibrium outcome: this is because they are selected so as to
    # minimize the component of the theoretical amenity score that is left
    # unexplained by chosen exogenous covariates. Indeed, we will keep the
    # explained component for the definition of the amenity score used in our
    # simulations.

    # Alternatively, we could obtain this score from model inversion (as we do
    # later with the disamenity factor from living in informal backyards
    # / settlements), implicitly defining it as a residual that allows the
    # model to fit the data better. In more intuitive terms, this boils down to
    # optimizing over the value of the amenity score at the same time as
    # utility levels when solving the equilibrium (remember that we consider
    # a closed-city model where total population per income group is exogenous
    # and utility levels are endogenous). This approach is standard in
    # quantitative spatial economic models (see Redding and Rossi-Hansberg,
    # 2017 for more details). However, we choose to model the score explicitly,
    # defining it a priori as a calibrated parameter, and only solving for
    # utility levels (and disamenity factors) in equilibrium.
    # NB: The approach taken in Pfeiffer et al., where the score is calibrated
    # a priori, but along other data moments, can be seen as an intermediate
    # between the two approaches presented above.

    # This approach comes with the risk of a poorer overall fit. However, if
    # validation shows that model predictions are good enough for use in spite
    # of this limitation, the approach also brings clear benefits. First,
    # it prevents the model from overfitting, which increases its external
    # validity when assessing counterfactuals. Then, the amenity score has a
    # clear empirical interpretation, which allows to assess more specifically
    # policies directed towards amenities (translating into a change of values
    # for some exogenous covariates, that can be direcly mapped to a change in
    # equilibrium outcomes).

    listVariation = np.arange(0.5, 1.5, 0.01)
    initUti2 = utilityTarget[1]
    listUti3 = utilityTarget[2] * listVariation
    listUti4 = utilityTarget[3] * listVariation

    # We define our equation on formal rents: see technical documentation
    # for math formulas
    if options["correct_kappa"] == 1 and options["deprec_land"] == 1:
        dataRent = (
            data_sp["price"] ** (coeff_a)
            * (param["depreciation_rate"]
               + interest_rate)
            / (coeffKappa * coeff_b ** coeff_b * coeff_a ** coeff_a)
            )
    elif options["correct_kappa"] == 1 and options["deprec_land"] == 0:
        dataRent = (
            (interest_rate * data_sp["price"]) ** coeff_a
            * (param["depreciation_rate"]
               + interest_rate) ** coeff_b
            / (coeffKappa * coeff_b ** coeff_b * coeff_a ** coeff_a)
            )
    elif options["correct_kappa"] == 0 and options["deprec_land"] == 1:
        dataRent = (
            data_sp["price"] ** (coeff_a)
            * (param["depreciation_rate"]
               + interest_rate)
            / (coeffKappa * coeff_b ** coeff_b)
            )
    elif options["correct_kappa"] == 0 and options["deprec_land"] == 0:
        dataRent = (
            (data_sp["price"] * interest_rate) ** coeff_a
            * (param["depreciation_rate"]
               + interest_rate) ** coeff_b
            / (coeffKappa * coeff_b ** coeff_b)
            )

    # Note that the above formulas consider that we observe land prices in the
    # data, hence the conversion to housing rents through the zero profit
    # condition for formal private developers. If we observe housing prices
    # instead, we should use to below formula, expressing rents as the annual
    # coupon associated with an infinite bond valued at housing price
    dataRent = data_sp["price"] * interest_rate

    # We import amenity data at the SP level
    amenities_sp = calam.import_exog_amenities(
        path_data, path_precalc_inp, 'SP')
    # We select amenity variables to be used in regressions
    # NB: choice has to do with relevance and exogeneity of variables
    # Choice set is relevant but may only refer to a small subset of locations
    # (train stations are often dysfunctional, for instance)
    variables_regression = [
        'distance_distr_parks', 'distance_ocean', 'distance_ocean_2_4',
        'distance_urban_herit', 'airport_cone2', 'slope_1_5', 'slope_5',
        'distance_biosphere_reserve', 'distance_train']

    # We run the parameter scanning.
    # Note that this may be long to run as it depends on the combination of all
    # inputs: we need to make sure that the estimated parameters fall within
    # the predefined value ranges to exclude corner solutions.
    (parametersScan, scoreScan, parametersAmenitiesScan, modelAmenityScan,
     parametersHousing) = calscan.EstimateParametersByScanning(
         incomeNetOfCommuting, dataRent, data_sp["dwelling_size"],
         data_income_group, selectedSP, amenities_sp, variables_regression,
         listRho, listBeta, listBasicQ, initUti2, listUti3, listUti4, options)

    print(modelAmenityScan.summary())

    # Now we run the optimization algorithm with identified value of the
    # parameters: this corresponds to an interior-point algorithm.
    # Note that this require above scanning process to be well-defined in
    # order to converge to an interior solution: we provide the option of not
    # running it in case this leads to absurd results.

    if options["param_optim"] == 1:

        initBeta = parametersScan[0]
        initBasicQ = parametersScan[1]
        initUti3 = parametersScan[2]
        initUti4 = parametersScan[3]

        (parameters, scoreTot, parametersAmenities, modelAmenity,
         parametersHousing) = calopt.EstimateParametersByOptimization(
             incomeNetOfCommuting, dataRent, data_sp["dwelling_size"],
             data_income_group, selectedSP, amenities_sp, variables_regression,
             listRho, initBeta, initBasicQ, initUti2, initUti3, initUti4,
             options)

        print(modelAmenity.summary())

    # Exporting and saving outputs

    amenities_grid = calam.import_exog_amenities(
        path_data, path_precalc_inp, 'grid')
    predictors_grid = amenities_grid.loc[:, variables_regression]
    predictors_grid = sm.add_constant(predictors_grid)
    # predictors_grid = np.vstack(
    #     [np.ones(predictors_grid.shape[0]),
    #      predictors_grid.T]
    #     ).T

    # NB: we only retain the explained component of the theoretical amenity
    # score
    if options["param_optim"] == 1:
        cal_amenities = np.exp(
            np.nansum(predictors_grid * parametersAmenities, 1))
        calibratedUtility_beta = parameters[0]
        calibratedUtility_q0 = parameters[1]
        calibratedUtility_u3 = parameters[2]
        calibratedUtility_u4 = parameters[3]
    elif options["param_optim"] == 0:
        cal_amenities = np.exp(
            np.nansum(predictors_grid * parametersAmenitiesScan, 1))
        calibratedUtility_beta = parametersScan[0]
        calibratedUtility_q0 = parametersScan[1]
        calibratedUtility_u3 = parametersScan[2]
        calibratedUtility_u4 = parametersScan[3]

    np.save(path_precalc_inp + 'calibratedUtility_beta',
            calibratedUtility_beta)
    np.save(path_precalc_inp + 'calibratedUtility_q0', calibratedUtility_q0)
    np.save(path_precalc_inp + 'calibratedAmenities', cal_amenities)

    return (calibratedUtility_beta, calibratedUtility_q0, cal_amenities,
            calibratedUtility_u3, calibratedUtility_u4)
