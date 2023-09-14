# -*- coding: utf-8 -*-

# Run this app with `python3 sens_matrix_dashboard.py` and
# view the plots at http://127.0.0.1:8050/ in your web browser.
# (To open a web browser on a larson-group computer,
# login to malan with `ssh -X` and then type `firefox &`.)

def main():

    import dash
    #from dash import dcc
    import dash_core_components as dcc
    #from dash import html
    import dash_html_components as html
    import plotly.express as px
    import plotly.graph_objects as go
    import pandas as pd

    import numpy as np
    import pdb
    import sklearn
    import plotly.figure_factory as ff
    from plotly.subplots import make_subplots
    #from plotly.figure_factory import create_quiver
    from itertools import chain

    from analyze_sensitivity_matrix import \
            analyzeSensMatrix, \
            findOutliers, findParamsUsingElastic
    from test_analyzeSensMatrix import write_test_netcdf_files
    from set_up_dashboard_inputs import setUpInputs, setUpPreliminaries, \
                                        setupDefaultMetricValsCol

    #print("New run --------------------------------------")

    # The user should input all tuning data into file set_up_dashboard_inputs.py
    metricsNames, metricsWeights, metricsNorms, \
    obsMetricValsDict, \
    paramsNames, paramsScales, \
    transformedParamsNames, \
    prescribedParamsNames, prescribedParamsScales, \
    prescribedTransformedParamsNames, \
    prescribedParamValsRow, \
    prescribedSensNcFilenames, prescribedSensNcFilenamesExt, \
    sensNcFilenames, sensNcFilenamesExt, \
    defaultNcFilename, linSolnNcFilename, \
    reglrCoef \
    = \
        setUpInputs()


    # Number of metrics
    numMetrics = len(metricsNames)

    obsMetricValsCol, normMetricValsCol, \
    defaultBiasesCol, \
    defaultParamValsOrigRow,\
    magParamValsRow, \
    dnormlzdPrescribedParams, \
    magPrescribedParamValsRow, \
    = setUpPreliminaries(metricsNames, metricsNorms,
                           obsMetricValsDict,
                           paramsNames, transformedParamsNames,
                           prescribedParamsNames, prescribedParamValsRow,
                           prescribedTransformedParamsNames,
                           sensNcFilenames,
                           defaultNcFilename
                          )

    #print("Orig defaultBiasesCol=", defaultBiasesCol)
    #print("metricsNames=", metricsNames)
    #print("obsMetricValsCol=", obsMetricValsCol)

    # Construct numMetrics x numParams matrix of second derivatives, d2metrics/dparams2.
    # The derivatives are normalized by observed metric values and max param values.
    normlzdCurvMatrix, normlzdSensMatrixPoly, normlzdConstMatrix, \
    normlzdOrdDparamsMin, normlzdOrdDparamsMax = \
        constructNormlzdCurvMatrix(metricsNames, paramsNames, transformedParamsNames, \
                                   metricsWeights, obsMetricValsCol, normMetricValsCol, magParamValsRow, \
                                   sensNcFilenames, sensNcFilenamesExt, defaultNcFilename)

    #print("normlzdSensMatrixPoly=", normlzdSensMatrixPoly)
    #print("normlzdCurvMatrix=", normlzdCurvMatrix)

    # In order to weight certain metrics, multiply each row of normlzdSensMatrixPoly
    # by metricsWeights
    normlzdWeightedSensMatrixPoly = np.diag(np.transpose(metricsWeights)[0]) @ normlzdSensMatrixPoly

    # For prescribed parameters, construct numMetrics x numParams matrix of second derivatives, d2metrics/dparams2.
    # The derivatives are normalized by observed metric values and max param values.
    normlzdPrescribedCurvMatrix, normlzdPrescribedSensMatrixPoly, normlzdPrescribedConstMatrix, \
    normlzdPrescribedOrdDparamsMin, normlzdPrescribedOrdDparamsMax = \
        constructNormlzdCurvMatrix(metricsNames, prescribedParamsNames, prescribedTransformedParamsNames, \
                                   metricsWeights, obsMetricValsCol, normMetricValsCol, magPrescribedParamValsRow, \
                                   prescribedSensNcFilenames, prescribedSensNcFilenamesExt, defaultNcFilename)

    #print("normlzdPrecribedSensMatrixPoly=", normlzdPrescribedSensMatrixPoly)
    #print("normlzdPrecribedCurvMatrix=", normlzdPrescribedCurvMatrix)

    # This is the prescribed correction to the metrics that appears on the left-hand side of the Taylor equation.
    #   It is not a bias from the obs.  It is a correction to the simulated default metric values
    #   based on prescribed param values.
    normlzdPrescribedBiasesCol = \
         fwdFnc( dnormlzdPrescribedParams, normlzdPrescribedSensMatrixPoly, normlzdPrescribedCurvMatrix, numMetrics )

    #print("normlzdPrecribedBiasesCol=", normlzdPrescribedBiasesCol)

    prescribedBiasesCol = normlzdPrescribedBiasesCol * np.abs(normMetricValsCol)

    #print("precribedBiasesCol=", prescribedBiasesCol)

    # defaultBiasesCol + prescribedBiasesCol = -fwdFnc_tuned_params  (see objFnc).
    # This lumps the prescribed-parameter adjustment into defaultBiasesCol
    defaultBiasesCol = defaultBiasesCol + prescribedBiasesCol

    #print("Orig+Prescribed defaultBiasesCol=", defaultBiasesCol)

    #pdb.set_trace()

    print("Optimizing parameter values . . . ")

    defaultBiasesApproxNonlin, \
    dnormlzdParamsSolnNonlin, paramsSolnNonlin, \
    dnormlzdParamsSolnLin, paramsSolnLin, \
    defaultBiasesApproxNonlin2x, \
    defaultBiasesApproxNonlinNoCurv, defaultBiasesApproxNonlin2xCurv = \
        solveUsingNonlin(metricsNames, paramsNames, transformedParamsNames,
                         metricsWeights, normMetricValsCol, magParamValsRow,
                         sensNcFilenames, sensNcFilenamesExt, defaultNcFilename,
                         defaultParamValsOrigRow,
                         normlzdSensMatrixPoly, defaultBiasesCol,
                         #normlzdSensMatrixPoly, defaultBiasesCol-prescribedBiasesCol,
                         normlzdCurvMatrix,
                         reglrCoef,
                         beVerbose=False)

    print("Tuned parameter values (paramsSolnNonlin)")
    for idx in range(0,len(paramsNames)): print("{:33s} {:7.7g}".format( paramsNames[idx], paramsSolnNonlin[idx][0] ) )

    normlzdLinplusSensMatrixPoly = normlzdSemiLinMatrixFnc(
                                        dnormlzdParamsSolnNonlin, normlzdSensMatrixPoly, 
                                        normlzdCurvMatrix, numMetrics)

    #normlzdLinplusSensMatrixPoly = normlzdSensMatrixPoly \
    #           + 0.5 * normlzdCurvMatrix * ( np.ones((len(metricsNames),1)) @ dnormlzdParamsSolnNonlin.T )

    #print("dnormlzdPrecribedParams=", dnormlzdPrescribedParams)
    #print("dnormlzdParamsSolnNonlin=", dnormlzdParamsSolnNonlin)
    #print("paramsSolnNonlin=", paramsSolnNonlin)
    #print("magParamValsRow=", magParamValsRow)
    #print("magPrescribedParamValsRow=", magPrescribedParamValsRow)

    #pdb.set_trace()

    # Calculate changes in parameter values needed to match metrics.
    #defaultMetricValsCol, defaultBiasesCol, \
    #defaultBiasesApprox, defaultBiasesApproxLowVals, defaultBiasesApproxHiVals, \
    #defaultBiasesApproxPC, defaultBiasesApproxLowValsPC, defaultBiasesApproxHiValsPC, \
    #normlzdWeightedDefaultBiasesApprox, normlzdWeightedDefaultBiasesApproxPC, \
    #defaultBiasesOrigApprox, defaultBiasesOrigApproxPC, \
    #sensMatrixOrig, sensMatrix, normlzdSensMatrix, \
    #normlzdWeightedSensMatrix, biasNormlzdSensMatrix, svdInvrsNormlzdWeighted, \
    #vhNormlzd, uNormlzd, sNormlzd, \
    #vhNormlzdWeighted, uNormlzdWeighted, sNormlzdWeighted, \
    #magParamValsRow, \
    #defaultParamValsOrigRow, dparamsSoln, dnormlzdParamsSoln, \
    #dparamsSolnPC, dnormlzdParamsSolnPC, \
    #paramsSoln, paramsLowVals, paramsHiVals, \
    #paramsSolnPC, paramsLowValsPC, paramsHiValsPC = \
    #      analyzeSensMatrix(metricsNames, paramsNames, transformedParamsNames,
    #                        metricsWeights,
    #                        sensNcFilenames, sensNcFilenamesExt, defaultNcFilename,
    #                        obsMetricValsCol, normMetricValsCol)

    # For the remaining calculations, define linear sensitivity based on polynomial
    #    (quadratic) curve fit, rather than a finite difference between simulated regional values.
    #normlzdSensMatrix = normlzdSensMatrixPoly

    # Find best-fit params by use of the Elastic Net algorithm
    defaultBiasesApproxElastic, defaultBiasesApproxElasticNonlin, \
    dnormlzdParamsSolnElastic, paramsSolnElastic = \
        findParamsUsingElastic(normlzdSensMatrixPoly, normlzdWeightedSensMatrixPoly,
                     defaultBiasesCol, normMetricValsCol, metricsWeights,
                     magParamValsRow, defaultParamValsOrigRow,
                     normlzdCurvMatrix,
                     beVerbose=False)

    defaultBiasesApproxElasticCheck = ( normlzdWeightedSensMatrixPoly @ dnormlzdParamsSolnElastic ) \
                            * np.reciprocal(metricsWeights) * np.abs(normMetricValsCol)

    #print("defaultBiasesApproxElastic = ", defaultBiasesApproxElastic)
    #print("defaultBiasesApproxElasticCheck = ", defaultBiasesApproxElasticCheck)

    #pdb.set_trace()

    # Estimate non-linearity of the global model to perturbations in parameter values.
    # To do so, calculate radius of curvature of the three points from the default simulation
    #   and the two sensitivity simulations.
    #calcNormlzdRadiusCurv(metricsNames, paramsNames, transformedParamsNames, paramsScales,
    #                      metricsWeights, obsMetricValsCol,
    #                      sensNcFilenames, sensNcFilenamesExt, defaultNcFilename)

    # Create scatterplot to look at outliers
    #createPcaBiplot(normlzdLinplusSensMatrixPoly, defaultBiasesCol, normMetricValsCol, metricsNames, paramsNames)

    ## Find outliers by use of the ransac algorithm
    #outlier_mask, defaultBiasesApproxRansac, normlzdWeightedDefaultBiasesApproxRansac, \
    #dnormlzdParamsSolnRansac, paramsSolnRansac = \
    #    findOutliers(normlzdSensMatrix, normlzdWeightedSensMatrix, \
    #                 defaultBiasesCol, normMetricValsCol, magParamValsRow, defaultParamValsOrigRow)
    #print( "ransac_outliers = ", metricsNames[outlier_mask] )
    #print( "ransac_inliers = ", metricsNames[~outlier_mask] )
    ##pdb.set_trace()



    ##############################################
    #
    #    Create plots
    #
    ##############################################


    print("Creating plots . . .")


    # Calculate symmetric error bars on fitted parameter values,
    #    based on difference in sensitivity matrix, i.e., based on size of nonlinear terms.
    #    For use in figures such as paramsBarChart.
    paramsLowValsPCBound, paramsHiValsPCBound = \
        calcParamsBounds(metricsNames, paramsNames, transformedParamsNames,
                     metricsWeights, obsMetricValsCol, normMetricValsCol,
                     magParamValsRow,
                     sensNcFilenames, sensNcFilenamesExt, defaultNcFilename)

    threeDotFig = \
        createThreeDotFig(metricsNames, paramsNames, transformedParamsNames,
                               metricsWeights, obsMetricValsCol, normMetricValsCol, magParamValsRow,
                               normlzdCurvMatrix, normlzdSensMatrixPoly, normlzdConstMatrix,
                               normlzdOrdDparamsMin, normlzdOrdDparamsMax,
                               sensNcFilenames, sensNcFilenamesExt, defaultNcFilename)


    # Check whether the minimizer actually reduces chisqd
    # Initial value of chisqd, which assumes parameter perturbations are zero
    normlzdDefaultBiasesCol = defaultBiasesCol/np.abs(normMetricValsCol)
    chisqdZero = objFnc(np.zeros_like(defaultParamValsOrigRow), \
                        normlzdSensMatrixPoly, normlzdDefaultBiasesCol, metricsWeights, \
                        normlzdCurvMatrix, reglrCoef, numMetrics)
    # Optimized value of chisqd, which uses optimal values of parameter perturbations
    chisqdMin = objFnc(dnormlzdParamsSolnNonlin.T, \
                        normlzdSensMatrixPoly, normlzdDefaultBiasesCol, metricsWeights, \
                        normlzdCurvMatrix, reglrCoef, numMetrics)

    #print("chisqdZero =", chisqdZero)
    #print("chisqdMin =", chisqdMin)
    print("chisqdMinRatio =", chisqdMin/chisqdZero)

    chisqdUnweightedZero = objFnc(np.zeros_like(defaultParamValsOrigRow), \
                        normlzdSensMatrixPoly, normlzdDefaultBiasesCol, np.ones_like(metricsWeights), \
                        normlzdCurvMatrix, reglrCoef, numMetrics)
    # Optimized value of chisqd, which uses optimal values of parameter perturbations
    chisqdUnweightedMin = objFnc(dnormlzdParamsSolnNonlin.T, \
                        normlzdSensMatrixPoly, normlzdDefaultBiasesCol, np.ones_like(metricsWeights), \
                        normlzdCurvMatrix, reglrCoef, numMetrics)

    #print("chisqdUnweightedZero =", chisqdUnweightedZero)
    #print("chisqdUnweightedMin =", chisqdUnweightedMin)
    print("chisqdUnweightedMinRatio =", chisqdUnweightedMin/chisqdUnweightedZero)

    # Set up a column vector of metric values from the global simulation based on optimized
    #     parameter values.
    linSolnMetricValsCol = setupDefaultMetricValsCol(metricsNames, linSolnNcFilename)

    # Store biases in default simulation, ( global_model - obs )
    linSolnBiasesCol = np.subtract(linSolnMetricValsCol, obsMetricValsCol)
    #linSolnBiasesCol = linSolnBiasesCol + prescribedBiasesCol

    # Check whether the minimizer actually reduces chisqd
    # Initial value of chisqd, which assumes parameter perturbations are zero
    normlzdLinSolnBiasesCol = linSolnBiasesCol/np.abs(normMetricValsCol)
    chisqdLinSolnMin = objFnc(np.zeros_like(defaultParamValsOrigRow), \
                        normlzdSensMatrixPoly, normlzdLinSolnBiasesCol, metricsWeights, \
                        normlzdCurvMatrix, reglrCoef, numMetrics)

    #print("chisqdZero =", chisqdZero)
    #print("chisqdMin =", chisqdMin)
    print("chisqdLinSolnMinRatio =", chisqdLinSolnMin/chisqdZero)

    chisqdUnweightedLinSolnMin = objFnc(np.zeros_like(defaultParamValsOrigRow), \
                        normlzdSensMatrixPoly, normlzdLinSolnBiasesCol, np.ones_like(metricsWeights), \
                        normlzdCurvMatrix, reglrCoef, numMetrics)

    #print("chisqdUnweightedZero =", chisqdUnweightedZero)
    #print("chisqdUnweightedMin =", chisqdUnweightedMin)
    print("chisqdUnweightedLinSolnMinRatio =", chisqdUnweightedLinSolnMin/chisqdUnweightedZero)
    print("-----------------------------------------------------")

    # Calculate the fraction of the default-sim bias that remains after tuning.
    # This is unweighted and hence is not necessarily less than one.
    # defaultBiasesApprox = J*delta_p = ( fwd - def )
    # numerator = ( fwd - def ) + ( def - obs ) = ( fwd - obs )
#    Bias = ( defaultBiasesApprox + defaultBiasesCol )
    # defaultBiasesCol = delta_b = ( default - obs ) = denominator
#    BiasMagRatio = np.linalg.norm(Bias/np.abs(normMetricValsCol))**2 / \
#                   np.linalg.norm(defaultBiasesCol/np.abs(normMetricValsCol))**2

    # Calculate the fraction of the default-sim bias that remains after tuning,
    #    but using a truncated PC observation.
    # This is unweighted and hence is not necessarily less than one.
    # defaultBiasesApproxPC = J*delta_p = ( fwd - def )
    # numerator = ( fwd - def ) + ( def - obs ) = ( fwd - obs )
#    BiasPC = ( defaultBiasesApproxPC + defaultBiasesCol )
    # defaultBiasesCol = delta_b = ( default - obs ) = denominator
#    BiasPCMagRatio = np.linalg.norm(BiasPC/np.abs(normMetricValsCol))**2 / \
#                     np.linalg.norm(defaultBiasesCol/np.abs(normMetricValsCol))**2

    ## Calculate the fraction of the default-sim bias that remains after tuning,
    ##    but using a truncated PC observation.
    ## This is unweighted and hence is not necessarily less than one.
    ## defaultBiasesApproxRansac = J*delta_p = ( fwd - def )
    ## numerator = ( fwd - def ) + ( def - obs ) = ( fwd - obs )
    #BiasRansac = ( defaultBiasesApproxRansac + defaultBiasesCol )
    ## defaultBiasesCol = delta_b = ( default - obs ) = denominator
    #BiasRansacMagRatio = np.linalg.norm(BiasRansac/np.abs(normMetricValsCol))**2 / \
    #                 np.linalg.norm(defaultBiasesCol/np.abs(normMetricValsCol))**2

    # Calculate the fraction of the default-sim bias that remains after tuning,
    #    but using a truncated PC observation.
    # This is unweighted and hence is not necessarily less than one.
    # defaultBiasesApproxElastic = J*delta_p = ( fwd - def )
    # numerator = ( fwd - def ) + ( def - obs ) = ( fwd - obs )
    BiasElastic = ( defaultBiasesApproxElastic + defaultBiasesCol )
    # defaultBiasesCol = delta_b = ( default - obs ) = denominator
    BiasElasticMagRatio = np.linalg.norm(BiasElastic/np.abs(normMetricValsCol))**2 / \
                     np.linalg.norm(defaultBiasesCol/np.abs(normMetricValsCol))**2

    # Calculate the global-model bias relative to the default-sim bias.
    # This is unweighted and hence is not necessarily less than one.
    # defaultBiasesApprox = J*delta_p = ( fwd - def )
    # numerator = ( linSoln - def ) + ( def - obs ) = ( linSoln - obs )
    #linSolnBias = ( linSolnBiasesCol + defaultBiasesCol )
    # defaultBiasesCol = delta_b = ( default - obs ) = denominator
    #linSolnBiasMagRatio = np.linalg.norm(linSolnBias/np.abs(normMetricValsCol))**2 / \
    #                      np.linalg.norm(defaultBiasesCol/np.abs(normMetricValsCol))**2

    # Calculate the fraction of bias removed by the non-PC soln, but normalized and weighted,
    # like the equations that the SVD actually solves, so that according to theory,
    # the value should be < 1.
    # But I'm not sure if it will be < 1 if the parameters are transformed to log space.
    normlzdMDeltaB = metricsWeights * defaultBiasesCol / np.abs(normMetricValsCol) # right-hand side
#    weightedBiasNumer = normlzdWeightedDefaultBiasesApprox + normlzdMDeltaB
#    weightedBiasDenom = normlzdMDeltaB
#    weightedBiasMagRatio = np.linalg.norm(weightedBiasNumer)**2 / np.linalg.norm(weightedBiasDenom)**2

    # Calculate the fraction of bias removed by PC soln, but normalized and weighted,
    # like the equations that the SVD actually solves, so that according to theory,
    # the value should be < 1.
    # But I'm not sure if it will be < 1 if the parameters are transformed to log space.
#    weightedBiasPCNumer = normlzdWeightedDefaultBiasesApproxPC + normlzdMDeltaB
#    weightedBiasPCDenom = normlzdMDeltaB
#    weightedBiasPCMagRatio = np.linalg.norm(weightedBiasPCNumer)**2 / np.linalg.norm(weightedBiasPCDenom)**2

    ## Calculate the fraction of bias removed by Ransac soln, but normalized and weighted,
    ## like the equations that the SVD actually solves, so that according to theory,
    ## the value should be < 1.
    ## But I'm not sure if it will be < 1 if the parameters are transformed to log space.
    #weightedBiasRansacNumer = normlzdWeightedDefaultBiasesApproxRansac + normlzdMDeltaB
    #weightedBiasRansacDenom = normlzdMDeltaB
    #weightedBiasRansacMagRatio = np.linalg.norm(weightedBiasRansacNumer)**2 / np.linalg.norm(weightedBiasRansacDenom)**2

    # Calculate the fraction of bias removed by Elastic soln, but normalized and weighted,
    # like the equations that the SVD actually solves, so that according to theory,
    # the value should be < 1.
    # But I'm not sure if it will be < 1 if the parameters are transformed to log space.
    weightedBiasElasticNumer = defaultBiasesApproxElasticNonlin + normlzdMDeltaB
    weightedBiasElasticDenom = normlzdMDeltaB
    weightedBiasElasticMagRatio = np.linalg.norm(weightedBiasElasticNumer)**2 / np.linalg.norm(weightedBiasElasticDenom)**2

# Fraction of bias that is removed by the tuned, non-linear forward solution
# weightedBiasLin = metricsWeights * ( fwd - obs ) / (def - obs )
#weightedBiasLin = metricsWeights * ( 1.0 - linSolnBiasesCol / defaultBiasesCol )
#weightedBiasLinMag = np.dot( np.transpose(weightedBiasLin), weightedBiasLin )
#weightedBiasLinMagRatio = weightedBiasLinMag / np.dot( np.transpose(metricsWeights), metricsWeights )
# weightedBiasLinMag = numerator
#weightedBiasLinMag = np.dot( np.transpose(weightedBiasLin), weightedBiasLin )
#weightedBiasLinMagRatio = weightedBiasLinMag / \
#                          np.dot( np.transpose(metricsWeights*defaultBiasesCol),
#                                  metricsWeights*defaultBiasesCol )

    # weightedBiasLin = metricsWeights * ( lin - obs ) = numerator
    #weightedBiasLinSoln = metricsWeights * ( linSolnBiasesCol + defaultBiasesCol ) / np.abs(normMetricValsCol)
    #weightedBiasLinSolnMagRatio = np.linalg.norm(weightedBiasLinSoln)**2 / np.linalg.norm(normlzdMDeltaB)**2


    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

    sensMatrixDashboard = dash.Dash(__name__, external_stylesheets=external_stylesheets)

    # Plot the biases of the default simulation and the SVD approximation of that
#    biasesMatrix = np.dstack((-defaultBiasesCol,
#                          defaultBiasesApprox,
#                          defaultBiasesApproxPC,
#                          #defaultBiasesApproxElastic,
#                          defaultBiasesApproxElasticNonlin,
#                          defaultBiasesApproxNonlin,
#                          linSolnBiasesCol
#                         )).squeeze()
#    fracBiasesMatrix = np.diagflat(np.reciprocal(np.abs(normMetricValsCol))) @ biasesMatrix
#    df = pd.DataFrame(fracBiasesMatrix,
#                    index=metricsNames,
#                  columns= ['fracDefBias',
#                            'fracDefBiasesApprox',
#                            'fracDefBiasesApproxPC',
#                            'fracDefBiasesApproxElasticNonlin',
#                            'fracDefBiasesApproxNonlin',
#                            'fracLinSolnBiasesCol'
#                           ])
#    biasesFig = px.line(df, x=df.index, y=df.columns,
#              title = """Fractional biases of default simulation and approximations thereof.<br>
#                    Plotted quantities have the structure -(def-obs), -(def-fwd), -(def-lin)""")
#    biasesFig.update_yaxes(title="-(Def-Sim) / abs(obs metric value)")
#    biasesFig.update_xaxes(title="Regional metric")
#    biasesFig.layout.legend.title = "Default or which approximation"
#    biasesFig.update_layout(hovermode="x")
#    biasesFig.data[1].name = "fracDefBiasesApprox, " \
#                         + "{:.2f}".format(weightedBiasMagRatio) \
#                        + ", {:.2f}".format(BiasMagRatio)
#    biasesFig.data[2].name = "fracDefBiasesApproxPC, " \
#                         + "{:.2f}".format(weightedBiasPCMagRatio) \
#                         + ", {:.2f}".format(BiasPCMagRatio)
#    biasesFig.data[3].name = "fracDefBiasesApproxElastic, " \
#                         + "{:.2f}".format(weightedBiasElasticMagRatio) \
#                         + ", {:.2f}".format(BiasElasticMagRatio)
#    biasesFig.data[4].name = "fracDefBiasesApproxNonlin, " \
#                         + "{:.2f}".format(-99) \
#                         + ", {:.2f}".format(-99)
#    biasesFig.data[5].name = "fracLinSolnBiasesCol, " \
#                          + "{:.2f}".format(weightedBiasLinSolnMagRatio) \
#                          + ", {:.2f}".format(linSolnBiasMagRatio)


    # Create plot showing how well the regional biases are actually removed
    metricsSens = np.linalg.norm(normlzdWeightedSensMatrixPoly, axis=1) # measure of sensitivity of each metric
    # metricsSensOrdered = (rankdata(metricsSens) - 1).astype(int)  # this ordering doesn't work as an index
    metricsSensOrdered = metricsSens.argsort()
    metricsNamesOrdered = metricsNames[metricsSensOrdered]  # list of metrics names, ordered from least to most sensitive
    # Plot a black dot for each default-run bias
    biasesOrderMatrix = np.dstack(( -defaultBiasesCol[metricsSensOrdered] )).squeeze()
    fracBiasesOrderMatrix = np.diagflat(np.reciprocal(np.abs(normMetricValsCol[metricsSensOrdered]))) @ biasesOrderMatrix
    df = pd.DataFrame(fracBiasesOrderMatrix,
                      index=metricsNamesOrdered,
                      columns= ['fracDefBias'])
    biasesOrderedArrowFig = px.line(df, x=df.index, y=df.columns,
              title = """<span style='color:blue'>Predicted</span> and <span style='color:red'>actual</span> removal of regional biases""")
    biasesOrderedArrowFig.update_yaxes(title="-Normalized bias")
    biasesOrderedArrowFig.update_xaxes(title="Regional metric")
    biasesOrderedArrowFig.update_layout(hovermode="x")
    biasesOrderedArrowFig.update_layout(showlegend=False)
    biasesOrderedArrowFig.update_traces(mode='markers', line_color='black')  # Plot default biases as black dots
    biasesOrderedArrowFig.update_yaxes(visible=True,zeroline=True,zerolinewidth=1,zerolinecolor='gray') # Plot x axis
    biasesOrderedArrowFig.update_layout( width=700, height=500  )
    # Now plot an arrow for each region that points from default-run bias to new bias after tuning
    xArrow = np.arange(len(metricsNamesOrdered)) # x-coordinate of arrows
    yArrow = -defaultBiasesCol[metricsSensOrdered,0]/np.abs(normMetricValsCol[metricsSensOrdered,0])
    gap = 0.2  # horizontal spacing between arrows
    # Plot error bar on prediction arrow.  Bar runs between 0- and 2x-curvature solns.
    for i, item in enumerate(metricsNamesOrdered):
        biasesOrderedArrowFig.add_annotation(
        x =  xArrow[i]-gap,  # ith arrow's head
        # ith arrow's head:
        y = (-defaultBiasesApproxNonlinNoCurv-defaultBiasesCol)[metricsSensOrdered[i],0]/np.abs(normMetricValsCol[metricsSensOrdered[i],0]),
        ax =  xArrow[i]-gap,  # ith arrow's head
        # ith arrow's head:
        ay = (-defaultBiasesApproxNonlin2xCurv-defaultBiasesCol)[metricsSensOrdered[i],0]/np.abs(normMetricValsCol[metricsSensOrdered[i],0]),
        font = dict(family="bold", color="blue", size=30),        
        showarrow=True,
        xref='x',
        yref='y',
        axref='x',
        ayref='y',
        text='',  # blank because we want only the arrow
        arrowhead=0,
        arrowsize=1,
        arrowwidth=6,
        arrowcolor='lightskyblue' # https://stackoverflow.com/questions/72496150/user-friendly-names-for-plotly-css-colors
	    )
    #biasesOrderedArrowFig.add_scatter(x=df.index, y=df.columns, line_color='pink')  # attempt to make black dot appear on top
    biasesOrderedArrowFig.update_traces(mode='markers', line_color='black')
    # Plot arrows showing the tuner's nonlinear predicted bias removal
    for i, item in enumerate(metricsNamesOrdered):
        biasesOrderedArrowFig.add_annotation(
        x=  xArrow[i] - gap,  # ith arrow's head
        # ith arrow's head:
        y= (-defaultBiasesApproxNonlin-defaultBiasesCol)[metricsSensOrdered[i],0]/np.abs(normMetricValsCol[metricsSensOrdered[i],0]),
        #y= (-defaultBiasesApproxNonlinNoCurv-defaultBiasesCol)[metricsSensOrdered[i],0]/np.abs(normMetricValsCol[metricsSensOrdered[i],0]),
        ax= xArrow[i] - gap,  # ith arrow's tail
        ay=  yArrow[i],  # ith arrow's tail
        xref='x',
        yref='y',
        axref='x',
        ayref='y',
        text='',  # blank because we want only the arrow
        showarrow=True,
        arrowhead=3,
        arrowsize=1,
        arrowwidth=2,
        arrowcolor='blue'
        )
    # Add a hand-made legend
    biasesOrderedArrowFig.add_annotation(text='Tuner prediction of bias removal',
                                  font=dict(color='blue'),
                                  align='left', xref='paper', yref='paper', 
                                  x=0.05, y=0.98, showarrow=False)
    biasesOrderedArrowFig.add_annotation(text='Error bar on tuner prediction',
                                  font=dict(color='skyblue'),
                                  align='left', xref='paper', yref='paper', 
                                  x=0.05, y=0.91, showarrow=False)
    biasesOrderedArrowFig.add_annotation(text='Realized E3SM bias removal',
                                  font=dict(color='red'),  # keep E3SM legend 
                                  #font=dict(color='rgba(255,0,0,0.0)'), # omit E3SM legend
                                  align='left', xref='paper', yref='paper', 
                                  x=0.05, y=0.84, showarrow=False)
    # Plot arrows showing the bias removal of E3SM's solution
    for i, item in enumerate(metricsNamesOrdered):
        biasesOrderedArrowFig.add_annotation(
        x=  xArrow[i]+gap,  # ith arrow's head
        # ith arrow's head:
        #y= (-linSolnBiasesCol-defaultBiasesCol)[metricsSensOrdered[i],0]/np.abs(normMetricValsCol[metricsSensOrdered[i],0]),
        y= (-linSolnBiasesCol)[metricsSensOrdered[i],0]/np.abs(normMetricValsCol[metricsSensOrdered[i],0]),
        ax= xArrow[i]+gap,  # ith arrow's tail
        ay=  yArrow[i],  # ith arrow's tail
        xref='x',
        yref='y',
        axref='x',
        ayref='y',
        text='',  # blank because we want only the arrow
        showarrow=True,
        arrowhead=3,
        arrowsize=1,
        arrowwidth=2,
        arrowcolor='red'  
        #,opacity=0.0  # omit arrow for E3SM's solution
	    )
#     # Plot 0-curvature error bars on prediction arrow
#     for i, item in enumerate(metricsNamesOrdered):
#         biasesOrderedArrowFig.add_annotation(
#         x =  xArrow[i]-gap,  # ith arrow's head
#         # ith arrow's head:
#         y = (-defaultBiasesApproxNonlinNoCurv-defaultBiasesCol)[metricsSensOrdered[i],0]/np.abs(normMetricValsCol[metricsSensOrdered[i],0]),
#         text ='-',  # plot horizontal line
#         font = dict(family="bold", color="blue", size=30),        
#         showarrow=False
# 	    )
#     # Plot 2x-curvature error bars on prediction arrow
#     for i, item in enumerate(metricsNamesOrdered):
#         biasesOrderedArrowFig.add_annotation(
#         x =  xArrow[i]-gap,  # ith arrow's head
#         # ith arrow's head:
#         y = (-defaultBiasesApproxNonlin2xCurv-defaultBiasesCol)[metricsSensOrdered[i],0]/np.abs(normMetricValsCol[metricsSensOrdered[i],0]),
#         text ='-',  # plot horizontal line
#         font = dict(family="bold", color="blue", size=30),        
#         showarrow=False
# 	    )               


    #biasesOrderedArrowFig.add_trace(go.Scatter(x=xArrow, y=yArrow,
    #                      name='Region of improvement', mode='markers',
    #                       marker=dict(color='green', size=14)))
    #pdb.set_trace()

    #biasesOrderedArrowFig.write_image('biasesOrderedArrowFig.png', scale=6)

    ## Create plot showing how well the regional biases are actually removed
    metricsSens = np.linalg.norm(normlzdWeightedSensMatrixPoly, axis=1) # measure of sensitivity of each metric
    ## metricsSensOrdered = (rankdata(metricsSens) - 1).astype(int)  # this ordering doesn't work as an index
    metricsSensOrdered = metricsSens.argsort()
    metricsNamesOrdered = metricsNames[metricsSensOrdered]  # list of metrics names, ordered from least to most sensitive
    normlzdSensMatrixOrdered = normlzdSensMatrixPoly[metricsSensOrdered,:]
    ## Form matrix of parameter perturbations, for later multiplication into the sensitivity matrix
    dnormlzdParamsSolnNonlinMatrix = np.ones((len(metricsNames),1)) @ dnormlzdParamsSolnNonlin.T
    normlzdSensParamsMatrixOrdered = normlzdSensMatrixOrdered * dnormlzdParamsSolnNonlinMatrix
    #df = pd.DataFrame(-1*normlzdSensParamsMatrixOrdered,
    #                  index=metricsNamesOrdered,
    #                  columns=paramsNames)
    #biasContrOrderFig = px.bar(df, x=df.index, y=df.columns,
    #          title = """Linear contributions to actual removal of regional biases""")
    #biasContrOrderFig.update_yaxes(title="-(Def-Sim) / abs(obs metric value)")
    #biasContrOrderFig.update_xaxes(title="Regional metric")
    #biasContrOrderFig.update_layout(hovermode="x")
    #biasContrOrderFig.update_layout(showlegend=True)
    ##biasContrOrderFig.update_traces(mode='markers', line_color='black')  # Plot default biases as black dots
    #biasContrOrderFig.update_yaxes(visible=True,zeroline=True,zerolinewidth=1,zerolinecolor='gray') # Plot x axis
    #biasContrOrderFig.update_layout( width=800, height=500  )
    # Now plot an arrow for each region that points from default-run bias to new bias after tuning
    # xArrow = np.arange(len(metricsNamesOrdered)) # x-coordinate of arrows
    # yArrow = -defaultBiasesCol[metricsSensOrdered,0]/np.abs(normMetricValsCol[metricsSensOrdered,0])
    # gap = 0.1  # horizontal spacing between arrows
    # Plot arrows showing the tuner's nonlinear predicted bias removal
    # for i, item in enumerate(metricsNamesOrdered):
    #     biasContrOrderFig.add_annotation(
    #     x=  xArrow[i] - gap,  # ith arrow's head
    #     # ith arrow's head:
    #     y= (-defaultBiasesApproxNonlin-defaultBiasesCol)[metricsSensOrdered[i],0]/np.abs(normMetricValsCol[metricsSensOrdered[i],0]),
    #     #y= (-defaultBiasesApproxNonlinNoCurv-defaultBiasesCol)[metricsSensOrdered[i],0]/np.abs(normMetricValsCol[metricsSensOrdered[i],0]),
    #     ax= xArrow[i] - gap,  # ith arrow's tail
    #     ay=  yArrow[i],  # ith arrow's tail
    #     xref='x',
    #     yref='y',
    #     axref='x',
    #     ayref='y',
    #     text='',  # blank because we want only the arrow
    #     showarrow=True,
    #     arrowhead=3,
    #     arrowsize=1,
    #     arrowwidth=2,
    #     arrowcolor='blue'
    #     )
    # # Add a hand-made legend
    # biasContrOrderFig.add_annotation(text='tuner prediction of bias removal',
    #                               font=dict(color='blue'),
    #                               align='left', xref='paper', yref='paper', x=0.05, y=0.9, showarrow=False)
    # biasContrOrderFig.add_annotation(text='realized E3SM bias removal',
    #                               font=dict(color='red'), #'rgba(255,0,0,0.0)'),
    #                               align='left', xref='paper', yref='paper', x=0.05, y=0.8, showarrow=False)



    # Create plot showing lumped linear+nonlinear contributions to each metric
    # Form matrix of parameter perturbations, for later multiplication into the sensitivity matrix
    dnormlzdParamsSolnNonlinMatrix = np.ones((len(metricsNames),1)) @ dnormlzdParamsSolnNonlin.T
    curvParamsMatrixOrdered = 0.5 * normlzdCurvMatrix[metricsSensOrdered,:] * dnormlzdParamsSolnNonlinMatrix**2
    #print("Sum rows=", np.sum(-normlzdSensParamsMatrixOrdered-curvParamsMatrixOrdered, axis=1))
    minusNonlinMatrixOrdered = -1*curvParamsMatrixOrdered + -1*normlzdSensParamsMatrixOrdered

    biasTotContrbBarFig = \
          createBarChart( minusNonlinMatrixOrdered, index=metricsNamesOrdered, columns=paramsNames,
                          orientation = 'v',
                          title="""Linear + nonlinear contributions of parameters to actual removal of regional biases""",
                          xlabel="Regional metric", ylabel="Contribution to bias removal", 
                          width=800, height=500 )

#    paramsTotContrbBarFig = \
#          createBarChart( minusNonlinMatrixOrdered.T, index=paramsNames, columns=metricsNamesOrdered,
#                          orientation = 'v',
#                          title="""Linear + nonlinear contributions of parameters to actual removal of regional biases""",
#                          xlabel="Parameter", ylabel="Contribution to bias removal",
#                          width=800, height=500 )

    sensMatrixBarFig = \
          createBarChart( normlzdSensMatrixOrdered.T, index=paramsNames, columns=metricsNamesOrdered,
                          orientation = 'v',
                          title="""Linear contributions of parameters to actual removal of regional biases""",
                          xlabel="Parameter", ylabel="Contribution to bias removal",
                          width=800, height=500 )

    linplusSensMatrixBarFig = \
          createBarChart( normlzdLinplusSensMatrixPoly[metricsSensOrdered,:].T, index=paramsNames, columns=metricsNamesOrdered,
                          orientation = 'v',
                          title="""Linplus contributions of parameters to actual removal of regional biases""",
                          xlabel="Parameter", ylabel="Contribution to bias removal", 
                          width=800, height=500 )

    #dnormlzdParamsSolnNonlinMatrix = np.ones((len(metricsNames),1)) @ dnormlzdParamsSolnNonlin.T
    #curvParamsMatrixOrdered = 0.5 * normlzdCurvMatrix[metricsSensOrdered,:] * dnormlzdParamsSolnNonlinMatrix**2
    #print("Sum rows=", np.sum(-normlzdSensParamsMatrixOrdered-curvParamsMatrixOrdered, axis=1))
    dfLin = pd.DataFrame(-1*normlzdSensParamsMatrixOrdered,
                      index=metricsNamesOrdered,
                      columns=paramsNames)
    # biasLinNlIndivContrbBarFig = go.Figure()
    # biasLinNlIndivContrbBarFig.add_trace( go.Bar(x=df.index, y=list(curvParamsMatrixOrdered[:,2:3])) )
    #df_long = pd.wide_to_long( df, i=df.index, j=df.columns, stubnames=[''] )
    dfLin = dfLin.reset_index()
    dfLin.rename(columns = {'index':'metricsNamesOrdered'}, inplace = True)
    #print("biasContrGrouped df=", df.to_string())
    #print("df.columns.values=", df.columns.values)
    #df.columns[1] = ['metricsNamesOrdered']

    dfLin_long = dfLin.melt( id_vars='metricsNamesOrdered', 
                            var_name='paramsNames', value_name='Contribution to bias removal')
    dfLin_long.insert(0, 'isNonlin', ['linear'] * len(paramsNames) * len(metricsNamesOrdered) )
    #print("df_long=", df_long.to_string())
    
    dfNonlin = pd.DataFrame(-1*curvParamsMatrixOrdered,
                      index=metricsNamesOrdered,
                      columns=paramsNames)
    # biasLinNlIndivContrbBarFig = go.Figure()
    # biasLinNlIndivContrbBarFig.add_trace( go.Bar(x=df.index, y=list(curvParamsMatrixOrdered[:,2:3])) )
    #df_long = pd.wide_to_long( df, i=df.index, j=df.columns, stubnames=[''] )
    dfNonlin = dfNonlin.reset_index()
    dfNonlin.rename(columns = {'index':'metricsNamesOrdered'}, inplace = True)
    #print("biasContrGrouped df=", df.to_string())
    #print("df.columns.values=", df.columns.values)
    #df.columns[1] = ['metricsNamesOrdered']

    dfNonlin_long = dfNonlin.melt( id_vars='metricsNamesOrdered', 
                            var_name='paramsNames', value_name='Contribution to bias removal')
    dfNonlin_long.insert(0, 'isNonlin', ['nonlinear'] * len(paramsNames) * len(metricsNamesOrdered) )
    #print("df_long=", df_long.to_string())    
    
    dfLinNonlin_long = pd.concat([dfLin_long, dfNonlin_long], ignore_index=True)
    
    biasLinNlIndivContrbBarFig = px.bar(dfLinNonlin_long,
                                 facet_col='metricsNamesOrdered', y='Contribution to bias removal', 
                                 x='isNonlin', color='paramsNames') #,
              #title = """Long: Linear ++ nonlinear contributions to actual removal of regional biases""")
    #biasLinNlIndivContrbBarFig.update_yaxes(title="-(Def-Sim) / abs(obs metric value)")
    #biasLinNlIndivContrbBarFig.update_xaxes(title="Regional metric")
    ##biasLinNlIndivContrbBarFig.update_xaxes(visible=False)
    ##biasLinNlIndivContrbBarFig.update_yaxes(visible=False)
    biasLinNlIndivContrbBarFig.update_layout(hovermode="x")
    biasLinNlIndivContrbBarFig.update_layout(showlegend=True)
    ##biasLinNlIndivContrbBarFig.update_xaxes(showticklabels=False).update_yaxes(showticklabels=False)
    #biasLinNlIndivContrbBarFig.for_each_annotation(lambda a: a.update(text=''))
    biasLinNlIndivContrbBarFig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    biasLinNlIndivContrbBarFig.update_annotations(textangle=-90)
    biasLinNlIndivContrbBarFig.update_layout(margin = dict(t = 160))
    for axis in biasLinNlIndivContrbBarFig.layout:
        #if type(biasLinNlIndivContrbBarFig.layout[axis]) == go.layout.YAxis:
        #    biasLinNlIndivContrbBarFig.layout[axis].title.text = 'Contribution to bias removal'
        if type(biasLinNlIndivContrbBarFig.layout[axis]) == go.layout.XAxis:
            biasLinNlIndivContrbBarFig.layout[axis].title.text = ''
    #biasLinNlIndivContrbBarFig.layout.title.text = ''
    #biasLinNlIndivContrbBarFig.update_yaxes(visible=True,zeroline=True,zerolinewidth=1,zerolinecolor='gray') # Plot x axis
    biasLinNlIndivContrbBarFig.update_layout( width=1000, height=450  )
    #print("curvParams =", -1*curvParamsMatrixOrdered)
    #print("normlzdSens =", -1*normlzdSensParamsMatrixOrdered)
    

    # Plot a scatterplot of default-simulation bias and SVD approximation of that bias.
    # Each column tells us how all metrics vary with a single parameter.
    biasSensDirMatrix = np.concatenate((defaultBiasesApproxNonlin/np.abs(normMetricValsCol),
                                       (-defaultBiasesCol/np.abs(normMetricValsCol))), axis=1)
    biasAndParamsNames = ["biasApproxNonlin", "bias"]
    #biasAndParamsNames = np.append(["bias", "bias_approx_pc"], paramsNames)
    df = pd.DataFrame(biasSensDirMatrix,
                  index=metricsNames,
                  columns=biasAndParamsNames)
    biasSensDirMatrixScatter = px.scatter(df, x="biasApproxNonlin", y="bias", 
                                          text=metricsNames, title="Bias approx vs bias" )
    biasSensDirMatrixOneOneLine = px.line(df, x="bias", y="bias")
    #biasSensDirMatrixOneMOneLine = px.line(df, x="bias", y=-df.loc[:,"bias"])
    biasVsBiasApproxScatterplot = go.Figure(data=biasSensDirMatrixScatter.data
                                              + biasSensDirMatrixOneOneLine.data)
    #biasRange = (max(df.loc[:,"bias"]), min(df.loc[:,"bias"]))
    #biasVsBiasApproxScatterplot.add_trace(go.Scatter(x=biasRange, y=biasRange, fill='tozeroy',
    #                           name='Region of improvement', mode='none',
    #                           fillcolor='rgba(253,253,150,0.7)'))
    biasVsBiasApproxScatterplot.update_xaxes(title="(defaultBiasesApproxNonlin)/obs")
    biasVsBiasApproxScatterplot.update_yaxes(title="-defaultBiasesCol/obs")
    biasVsBiasApproxScatterplot.update_traces(textposition='top center')
    biasVsBiasApproxScatterplot.update_yaxes(visible=True,zeroline=True,zerolinewidth=2,zerolinecolor='lightblue') # Plot x axis
    biasVsBiasApproxScatterplot.update_layout( width=800, height=500  )
    biasVsBiasApproxScatterplot.update_layout(title="Bias approx vs bias")


    # Plot a scatterplot of default-simulation bias and SVD approximation of that bias.
    # Each column tells us how all metrics vary with a single parameter.
    biasSensMatrix = np.concatenate((-defaultBiasesCol/np.abs(normMetricValsCol),
                                    (-defaultBiasesApproxElastic-defaultBiasesCol)/np.abs(normMetricValsCol)), axis=1)
                                     #defaultBiasesApproxElastic/np.abs(normMetricValsCol)), axis=1)
    biasAndParamsNames = ["bias", "bias_approx_pc"]
    #biasAndParamsNames = np.append(["bias", "bias_approx_pc"], paramsNames)
    df = pd.DataFrame(biasSensMatrix,
                  index=metricsNames,
                  columns=biasAndParamsNames)
    biasSensMatrixScatter = px.scatter(df, x="bias", y="bias_approx_pc", text=metricsNames,
    #biasSensMatrixScatter = px.scatter(df, x=np.append(["bias_approx_pc"], paramsNames), y="bias",
              title = """Columns of normalized sensitivity matrix.<br>
                       vs. bias vector.<br>
                       """ )
    biasSensMatrixOneOneLine = px.line(df, x="bias", y="bias")
    biasSensMatrixOneMOneLine = px.line(df, x="bias", y=-df.loc[:,"bias"])
    biasSensMatrixScatterFig = go.Figure(data=biasSensMatrixScatter.data
                                              + biasSensMatrixOneOneLine.data
                                              + biasSensMatrixOneMOneLine.data)
    biasRange = (max(df.loc[:,"bias"]), min(df.loc[:,"bias"]))
    biasSensMatrixScatterFig.add_trace(go.Scatter(x=biasRange, y=biasRange, fill='tozeroy',
                               name='Region of improvement', mode='none',
                               fillcolor='rgba(253,253,150,0.7)'))
    biasSensMatrixScatterFig.update_yaxes(title="(-defaultBiasesApproxElastic-defaultBiasesCol)/obs")
    biasSensMatrixScatterFig.update_xaxes(title="-defaultBiasesCol/obs")
    biasSensMatrixScatterFig.update_traces(textposition='top center')


    # Plot a scatterplot of minimum parameter perturbation vs. fractional default bias approximation
    # Calculate lower bound on normalized parameter perturbations
    #normlzdDefaultBiasesCol = ( metricsWeights * (-defaultBiasesCol) /
    normlzdDefaultBiasesCol = ( (defaultBiasesCol) /
                                np.abs(normMetricValsCol) )
    #sensMatrixRowMag = np.linalg.norm(normlzdWeightedSensMatrixPoly, axis=1)
    sensMatrixRowMag = np.linalg.norm(normlzdSensMatrixPoly, axis=1)
    #sensMatrixRowMag = np.amax(np.abs(normlzdSensMatrixPoly), axis=1)
    dpMin = np.abs(normlzdDefaultBiasesCol) / np.atleast_2d(sensMatrixRowMag).T
    #u_dot_b = np.atleast_2d(sensMatrixRowMag).T * -normlzdDefaultBiasesCol
    dpMinMatrix = np.dstack((np.reciprocal(dpMin),
    #dpMinMatrix = np.dstack((np.abs(u_dot_b),
    ##dpMinMatrix = np.dstack((np.atleast_2d(sensMatrixRowMag).T,
    #                      np.abs(defaultBiasesApproxElastic)/np.abs(normMetricValsCol)
                          np.abs(defaultBiasesCol)/np.abs(normMetricValsCol)
                         )).squeeze()
    biasAndParamsNames = ["dpMinInvrs", "bias_approx"]
    df = pd.DataFrame(dpMinMatrix,
                  index=metricsNames,
                  columns=biasAndParamsNames)
    dpMinMatrixScatter = px.scatter(df, x="dpMinInvrs", y="bias_approx", text=metricsNames,
              title = """dpMinInvrs  vs. |approx bias vector|.<br>
                       """ )
    dpMinMatrixScatterFig = go.Figure(data=dpMinMatrixScatter.data)
    dpMinMatrixScatterFig.update_yaxes(title="|defaultBiasesApproxElastic|")
    dpMinMatrixScatterFig.update_xaxes(title="dpMinInvrs")
    dpMinMatrixScatterFig.update_traces(textposition='top center')


    # Plot the sensitivity of each regional metric.
    #    More specifically, plot the maximum magnitude value of each row of the sensitivity matrix.
    df = pd.DataFrame(np.max(np.abs(normlzdSensMatrixPoly), axis=1), # max of absolute val of each row
                  index=metricsNames,
                  columns=['Max abs normlzd sensitivity'])
    maxSensMetricsFig = px.line(df, x=df.index, y=df.columns,
              title = """Maximum normalized sensitivity of each metric with respect to parameters.<br>
                       (Low sensitivity means that the metric is unbudgeable by these parameters.)""" )
    maxSensMetricsFig.update_yaxes(title="Max |sens row|")
    maxSensMetricsFig.update_xaxes(title="Regional metric")
    maxSensMetricsFig.update_layout(hovermode="x")
    maxSensMetricsFig.update_traces(mode='lines+markers')


    # Plot the biases versus sensitivity of each regional metric.
    #    More specifically, plot the maximum magnitude value of each row of the sensitivity matrix.
    #df = pd.DataFrame({'Max abs normlzd sensitivity': np.max(np.abs(normlzdSensMatrixPoly), axis=1), # max |row elements|
    df = pd.DataFrame({'Max abs normlzd sensitivity': np.sum(normlzdWeightedSensMatrixPoly, axis=1), # sum of row elements
    #df = pd.DataFrame({'Max abs normlzd sensitivity': np.linalg.norm(normlzdWeightedSensMatrixPoly, axis=1), # rms of row elements
    #df = pd.DataFrame({'Max abs normlzd sensitivity': np.linalg.norm(normlzdSensMatrixPoly, axis=1), # rms of row elements
    #df = pd.DataFrame({'Max abs normlzd sensitivity':
    #                    -defaultBiasesCol[:,0]/np.abs(normMetricValsCol[:,0])*np.linalg.norm(normlzdWeightedSensMatrixPoly, axis=1), # sum of row elements
                       'Max abs normlzd lin+ sensitivity': np.linalg.norm(normlzdLinplusSensMatrixPoly, axis=1), # sum of row elements
                       'Default biases': -defaultBiasesCol[:,0]/np.abs(normMetricValsCol[:,0]),
    #                   'revised tuning': (-defaultBiasesApproxElastic-defaultBiasesCol)[:,0]/np.abs(normMetricValsCol[:,0])
                      }, index=metricsNames )
    #biasesVsSensMagScatterplot = px.scatter(df, x='Max abs normlzd sensitivity', y=df.columns[1:2],
    #biasesVsSensMagScatterplot = px.scatter(df, x='Max abs normlzd sensitivity', y='default tuning',
    biasesVsSensMagScatterplot = px.scatter(df, x=['Max abs normlzd sensitivity','Max abs normlzd lin+ sensitivity'], y='Default biases',
                                 text=metricsNames, 
                                 title = """Regional biases vs. magnitude of sensitivity.""")
    biasesVsSensMagScatterplot.update_yaxes(title="Regional biases")
    biasesVsSensMagScatterplot.update_xaxes(title="Magnitude of sensitivity of regional metrics to parameter changes")
    biasesVsSensMagScatterplot.update_layout(showlegend=False)
    biasesVsSensMagScatterplot.update_layout(hovermode="x")
    biasesVsSensMagScatterplot.update_layout( width=700, height=500  )

    # Compute length of arrows between default and tuned biases
    #metricsNamesPadded = ",,".join(metricsNames).split(",")
    #metricsNamesPadded = ",,".join(metricsNamesPadded).split(",")
    #metricsNamesPadded = np.append(metricsNamesPadded, ["", "", ""], axis=0)
    xArrow = np.linalg.norm(normlzdWeightedSensMatrixPoly, axis=1)
    yArrow = -defaultBiasesCol[:,0]/np.abs(normMetricValsCol[:,0])
    #uArrow = np.zeros_like(xArrow)
    #vArrow = (-defaultBiasesApproxElasticNonlin)[:,0]/np.abs(normMetricValsCol[:,0])
    #arrowFig = create_quiver(xArrow, yArrow, uArrow, vArrow,
    #                         scale=1,text=metricsNamesPadded)
    #arrowFig.update_yaxes(title="Regional biases")
    #arrowFig.update_traces(mode='lines+text')  # make labels appear in plot, not just hovermode


    # Plot biases vs. sensitivity, but with arrows indicating the degree of bias reduction
    df = pd.DataFrame({'Max abs normlzd sensitivity': np.linalg.norm(normlzdWeightedSensMatrixPoly, axis=1), # sum of row elements

                       'default tuning': -defaultBiasesCol[:,0]/np.abs(normMetricValsCol[:,0]),
                       'revised tuning': (-defaultBiasesApproxNonlin-defaultBiasesCol)[:,0]/np.abs(normMetricValsCol[:,0])
                      }, index=metricsNames )
    biasesVsSensArrowFig = px.scatter(df, x='Max abs normlzd sensitivity', y=df.columns[1:2],
                                 text=metricsNames,
                                 title = """Regional biases with default and nonlin tuning versus sensitivity, with arrows.""" )
    biasesVsSensArrowFig.update_traces(textposition="middle right")
    for i, item in enumerate(metricsNames):
        biasesVsSensArrowFig.add_annotation(
        x=  xArrow[i],  # ith arrow's head
        y= (-defaultBiasesApproxNonlin-defaultBiasesCol)[i,0]/np.abs(normMetricValsCol[i,0]),  # ith arrow's length
        #y= (-defaultBiasesApproxNonlin2x-defaultBiasesCol)[i,0]/np.abs(normMetricValsCol[i,0]),  # ith arrow's length
        ax= xArrow[i],  # ith arrow's tail
        ay=  yArrow[i],  # ith arrow's tail
        xref='x',
        yref='y',
        axref='x',
        ayref='y',
        text='',  # if you want only the arrow
        showarrow=True,
        arrowhead=3,
        arrowsize=1,
        arrowwidth=2,
        arrowcolor='black'
        )
    biasesVsSensArrowFig.update_yaxes(title="Regional biases")
    biasesVsSensArrowFig.update_xaxes(title="Sensitivity of regional metrics to parameter changes")
    biasesVsSensArrowFig.update_layout(hovermode="x")
    biasesVsSensArrowFig.update_traces(cliponaxis=False)
    biasesVsSensArrowFig.update_yaxes(automargin=True)


    # Plot the residual*sensitivity vs. bias*sensitivity
    # The goal is to separate out which regional are amenable to tuning
    normlzdBias = -defaultBiasesCol[:,0]/np.abs(normMetricValsCol[:,0])
    normlzdResid = (-defaultBiasesApproxElastic-defaultBiasesCol)[:,0]/np.abs(normMetricValsCol[:,0])
    df = pd.DataFrame({'bias times sensitivity': normlzdBias*np.sum(normlzdSensMatrixPoly, axis=1), # sum of row elements
                       'residual times sensitivity': normlzdResid*np.sum(normlzdSensMatrixPoly, axis=1) # sum of row elements
                      }, index=metricsNames )
    residVsBiasScatter = px.scatter(df, x='bias times sensitivity', y=df.columns[1:],
                                 text=metricsNames, title = """Residual times sensitivity versus bias times sensitivity.""" )
    residVsBiasOneOneLine = px.line(df, x="bias times sensitivity", y="bias times sensitivity")
    residVsBiasFig = go.Figure(data=residVsBiasScatter.data
                                    + residVsBiasOneOneLine.data)
    residVsBiasFig.update_yaxes(title="Residual times sensitivity")
    residVsBiasFig.update_xaxes(title="Bias times sensitivity")
    residVsBiasFig.update_layout(hovermode="x")

    # Plot the relative biases versus sensitivity of each regional metric.
    #    More specifically, plot the maximum magnitude value of each row of the sensitivity matrix.
    #pdb.set_trace()
    relBiasNumerator = np.abs(-defaultBiasesApproxElastic-defaultBiasesCol)[:,0]/np.abs(normMetricValsCol[:,0])
    relBiasDenom = np.maximum(0.02, np.abs(-defaultBiasesCol[:,0]/np.abs(normMetricValsCol[:,0])) )
    df = pd.DataFrame({'Max abs normlzd sensitivity': np.max(np.abs(normlzdSensMatrixPoly), axis=1), # max |row elements|
                       'relative bias': relBiasNumerator / relBiasDenom
                      }, index=metricsNames )
    relBiasesVsSensFig = px.scatter(df, x='Max abs normlzd sensitivity', y=df.columns[1:],
                                 text=metricsNames, title = """Ratio of abs(approx_bias)/abs(default_bias) versus sensitivity.""" )
    relBiasesVsSensFig.update_yaxes(title="Relative biases")
    relBiasesVsSensFig.update_xaxes(title="Max abs normlzd sensitivity")
    relBiasesVsSensFig.update_layout(hovermode="x")

    # Plot the relative biases versus sensitivity of each regional metric.
    #    More specifically, plot the maximum magnitude value of each row of the sensitivity matrix.
    absBiasTuned = np.abs(-defaultBiasesApproxElastic-defaultBiasesCol)[:,0]/np.abs(normMetricValsCol[:,0])
    absBiasDefault = np.abs(-defaultBiasesCol[:,0]/np.abs(normMetricValsCol[:,0]))
    df = pd.DataFrame({'Max abs normlzd sensitivity': np.max(np.abs(normlzdSensMatrixPoly), axis=1), # max |row elements|
                       'bias difference': absBiasDefault - absBiasTuned
                      }, index=metricsNames )
    diffBiasesVsSensFig = px.scatter(df, x='Max abs normlzd sensitivity', y=df.columns[1:],
                                     text=metricsNames,
                                     title = """abs(default_bias) - abs(approx_bias) versus sensitivity.""" )
    diffBiasesVsSensFig.update_yaxes(title="Bias difference")
    diffBiasesVsSensFig.update_xaxes(title="Max abs normlzd sensitivity")
    diffBiasesVsSensFig.update_layout(hovermode="x")

    # Plot each column of normalized sensitivity matrix as a separate line.
    # Each column tells us how all metrics vary with a single parameter.
    df = pd.DataFrame( np.hstack( (-defaultBiasesCol/np.abs(normMetricValsCol),normlzdSensMatrixPoly) ),
                  index=metricsNames,
                  columns=np.append('Norm bias', paramsNames) )
    normlzdSensMatrixColsFig = px.line(df, x=df.index, y=df.columns,
              title =  """Columns of normalized, unweighted sensitivity matrix (plus the bias
                                       vector).<br>
                       Each column (line) shows how sensitive the metrics are to a change in a single parameter value.<br>
                       (A positive value means that an increase in parameter value brings the default simulation closer to obs.)""" )
    normlzdSensMatrixColsFig.update_yaxes(title="Norml sens, (|param|/|obsmetric|) * dmetric/dparam")
    normlzdSensMatrixColsFig.update_xaxes(title="Regional metric")
    normlzdSensMatrixColsFig.layout.legend.title = "Parameter"
    normlzdSensMatrixColsFig.update_layout(hovermode="x")
    #pdb.set_trace()

    # Plot each row of normalized sensitivity matrix as a separate line.
    # Each row tells us how a single metric varies with all parameters.
    df = pd.DataFrame(np.transpose(normlzdSensMatrixPoly),
                  index=paramsNames,
                  columns=metricsNames)
    normlzdSensMatrixRowsFig = px.line(df, x=df.index, y=df.columns,
              title = """Rows of normalized, unweighted sensitivity matrix.<br>
                       Each row (line) tells us how a single metric varies with all parameters.<br>
                       (A positive value means that an increase in parameter value brings the default simulation closer to obs.)""" )
    normlzdSensMatrixRowsFig.update_yaxes(title="Norml sens, (|param|/|obsmetric|) * dmetric/dparam")
    normlzdSensMatrixRowsFig.update_xaxes(title="Parameter")
    normlzdSensMatrixRowsFig.layout.legend.title = "Metric"
    normlzdSensMatrixRowsFig.update_layout(hovermode="x")

    # Plot each column of right-singular vector matrix, V.
#    rightSingVectorNums = (np.arange(paramsNames.shape[0])+1).astype(str)
#    df = pd.DataFrame(np.transpose(vhNormlzd),
#                  index=paramsNames,
#                  columns=rightSingVectorNums)
#    vhNormlzdColsFig = px.line(df, x=df.index, y=df.columns,
#              title = """Columns of normalized, unweighted right-singular vector matrix, V.<br>
#                        Each column (line) is a vector of parameter values associated with a singular value.<br>""" )
#    vhNormlzdColsFig.update_yaxes(title="Right-singular vector")
#    vhNormlzdColsFig.update_xaxes(title="Parameter")
#    vhNormlzdColsFig.layout.legend.title = "Singular vector"
#    vhNormlzdColsFig.update_layout(hovermode="x")
#    for idx, val in np.ndenumerate(sNormlzd):
#        vhNormlzdColsFig.data[idx[0]].name = "{}".format(idx[0]+1) + ", " + "{:.2e}".format(val)

    # Plot each column of normalized, unweighted left-singular vector matrix, U.
#    df = pd.DataFrame(uNormlzd,
#                  index=metricsNames,
#                  columns=rightSingVectorNums)
#    uNormlzdColsFig = px.line(df, x=df.index, y=df.columns,
#              title = """Columns of normalized, unweighted left-singular vector matrix, U.<br>
#                       Each column (line) is a vector of metric values associated with a singular value.<br>""" )
#    uNormlzdColsFig.update_yaxes(title="Left-singular vector")
#    uNormlzdColsFig.update_xaxes(title="Parameter")
#    uNormlzdColsFig.layout.legend.title = "Singular vector"
#    uNormlzdColsFig.update_layout(hovermode="x")
#    for idx, val in np.ndenumerate(sNormlzd):
#        uNormlzdColsFig.data[idx[0]].name = "{}".format(idx[0]+1) + ", " + "{:.2e}".format(val)

    # Plot each column of normalized, weighted left-singular vector matrix, U.
#    df = pd.DataFrame(uNormlzdWeighted,
#                  index=metricsNames,
#                  columns=rightSingVectorNums)
#    uNormlzdWeightedColsFig = px.line(df, x=df.index, y=df.columns,
#              title = """Columns of normalized, weighted left-singular vector matrix, U.<br>
#                       Each column (line) is a vector of metric values associated with a singular value.<br>""" )
#    uNormlzdWeightedColsFig.update_yaxes(title="Left-singular vector")
#    uNormlzdWeightedColsFig.update_xaxes(title="Parameter")
#    uNormlzdWeightedColsFig.layout.legend.title = "Singular vector"
#    uNormlzdWeightedColsFig.update_layout(hovermode="x")
#    for idx, val in np.ndenumerate(sNormlzdWeighted):
#        uNormlzdWeightedColsFig.data[idx[0]].name = "{}".format(idx[0]+1) + ", " + "{:.2e}".format(val)

    # Plot each column of left-singular vector matrix, U, multiplied by biases.
    # Plot each column of left-singular vector matrix, U, multiplied by biases.
#    df = pd.DataFrame(uNormlzd*defaultBiasesCol/np.abs(normMetricValsCol),
#
#                   index=metricsNames,
#                  columns=rightSingVectorNums)
#    uNormlzdBiasColsFig = px.line(df, x=df.index, y=df.columns,
#              title = """Columns of normalized, unweighted left-singular vector matrix, U, dotted with delta_b (defaultBiasesCol).<br>
#                       Each column (line) is a vector of metric values associated with a singular value.<br>""" )
#    uNormlzdBiasColsFig.update_yaxes(title="uNormlzd dot delta_b")
#    uNormlzdBiasColsFig.update_xaxes(title="Parameter")
#    uNormlzdBiasColsFig.layout.legend.title = "Singular vector"
#    uNormlzdBiasColsFig.update_layout(hovermode="x")
#    for idx, val in np.ndenumerate(sNormlzd):
#        uNormlzdBiasColsFig.data[idx[0]].name = "{}".format(idx[0]+1) + ", " + "{:.2e}".format(val)
    #pdb.set_trace()

    # Plot box and whiskers plot of optimal parameter values.
    # Multiply in the user-designated scale factors before plotting.
    paramsAbbrv = np.char.replace( paramsNames, 'clubb_', '' )
    paramsAbbrv = np.char.replace( paramsAbbrv, 'c_invrs_tau_', '' )
    paramsAbbrv = np.char.replace( paramsAbbrv, 'wpxp_n2', 'n2' )
    paramsAbbrv = np.char.replace( paramsAbbrv, 'wpxp_n2', 'n2' )
    paramsAbbrv = np.char.replace( paramsAbbrv, 'threshold', 'thresh' )
    df = pd.DataFrame( np.hstack( defaultParamValsOrigRow[0,:]*paramsScales ),
                  index=paramsAbbrv, columns=["Default plus error bars"] )
    df["err_minus"] = ( defaultParamValsOrigRow[0,:] -  paramsLowValsPCBound[:,0] ) * paramsScales
    df["err_plus"]  = ( paramsHiValsPCBound[:,0] - defaultParamValsOrigRow[0,:] ) * paramsScales
    paramsBarChart = px.scatter(df, x=df.index, y=df.columns,
              error_y="err_plus", error_y_minus="err_minus",
              title =  """Best-fit parameter values""" )
    paramsBarChart.update_traces( go.Scatter(
                            mode='markers',
                            marker=dict(color='black', size=14),
                            error_y=dict( color='black', thickness=2, width=10 )
                                    ))
    #paramsBarChart.add_trace(go.Scatter(x=paramsNames, y=paramsLowValsPCBound[:,0]*paramsScales,
    #                               name=r'$paramsSolnPC - \sigma$',
    #                               line=dict(color='white', width=0), mode='lines', showlegend=False))
    #paramsBarChart.add_trace(go.Scatter(x=paramsNames, y=paramsHiValsPCBound[:,0]*paramsScales, fill='tonexty',
    #                           name='Default Parameter Values +- sigma', mode='none',
    #                               fillcolor='rgba(253,253,150,1.0)'))
    #paramsBarChart.add_trace(go.Scatter(x=paramsNames, y=defaultParamValsOrigRow[0,:]*paramsScales,
    #                               name='Default Parameter Values', line=dict(color='black', width=6) ))
    paramsBarChart.add_trace(go.Scatter(x=paramsAbbrv, y=paramsSolnLin[:,0]*paramsScales,
                                   mode='markers',
                                   marker=dict(color='green', size=8),
                                   name='Linear regression, |dp|=' 
                                       + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSolnLin)) ))
    paramsBarChart.add_trace(go.Scatter(x=paramsAbbrv, y=paramsSolnNonlin[:,0]*paramsScales,
                                   mode='markers',
                                   marker_symbol='x',
                                   marker=dict(color='orange',  size=12),
                                   name='paramsSolnNonlin, |dpPC|='
                                       + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSolnNonlin)) ))
    paramsBarChart.add_trace(go.Scatter(x=paramsAbbrv, y=paramsSolnElastic[:,0]*paramsScales,
                                   mode='markers',
                                   marker_symbol='square',
                                   marker=dict(color='cyan', size=8),
                                   name='Lasso regression, |dpLasso|='
                                        + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSolnElastic)) ,
                                   line=dict(color='red', width=2)  ))
    #paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnPCBound[:,0]*paramsScales,
    #                                name='paramsSolnPCBound, |dpBound|='
    #                               + '{:.2e}'.format(0.0) ))
    paramsBarChart.update_yaxes(title="User-scaled parameter value")
    paramsBarChart.update_xaxes(title="Parameter Name")
    paramsBarChart.update_layout(hovermode="x")
    paramsBarChart.update_layout( width=1000, height=500  )


    # Plot the parameter values recommended by SVD.
    # Multiply in the user-designated scale factors before plotting.
    paramsFig = go.Figure()
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsLowValsPCBound[:,0]*paramsScales,
                                   name=r'$paramsSolnPC - \sigma$',
                                   line=dict(color='white', width=0), mode='lines', showlegend=False))
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsHiValsPCBound[:,0]*paramsScales, fill='tonexty',
                               name='Default Parameter Values +- sigma', mode='none',
                                   fillcolor='rgba(253,253,150,1.0)'))
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=defaultParamValsOrigRow[0,:]*paramsScales,
                                   name='Default Parameter Values', line=dict(color='black', width=6) ))
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnLin[:,0]*paramsScales,
                                   name='Linear regression, |dp|=' + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSolnLin)) ))
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnNonlin[:,0]*paramsScales,
                                   name='paramsSolnNonlin, |dpNonlin|='
                                   + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSolnNonlin)) ))
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnElastic[:,0]*paramsScales,
                                    name='Lasso regression, |dpLasso|='
                                   + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSolnElastic)) ,
                                   line=dict(color='red', width=2)  ))
    #paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnPCBound[:,0]*paramsScales,
    #                                name='paramsSolnPCBound, |dpBound|='
    #                               + '{:.2e}'.format(0.0) ))
    paramsFig.update_yaxes(title="User-scaled parameter value")
    paramsFig.update_xaxes(title="Parameter Name")
    paramsFig.update_layout(hovermode="x")
    paramsFig.update_layout( width=1000, height=500  )

    #pdb.set_trace()

    roundedNormlzdSensMatrix = np.around( normlzdLinplusSensMatrixPoly, decimals=2)
    df_sensmat = pd.DataFrame(roundedNormlzdSensMatrix,
                  index=metricsNames,
                  columns=paramsNames)
    matMaxAbs = np.max(np.abs(roundedNormlzdSensMatrix))
    normlzdSensMatrixFig = ff.create_annotated_heatmap(
                   z=df_sensmat.to_numpy(),
                   x=df_sensmat.columns.tolist(),
                   y=df_sensmat.index.tolist(),
                   #coloraxis="coloraxis",
                   #colorscale=px.colors.diverging.balance,
                   colorscale='balance',
                   zmin=-matMaxAbs, zmax=matMaxAbs,
                   showscale=False, ygap=1, xgap=1
                   )
    normlzdSensMatrixFig.update_xaxes(side="bottom")
    normlzdSensMatrixFig.update_layout(
    title_text='Color-coded normalized sensitivity matrix', 
    title_x=0.5, 
    #width=800, 
    #height=1400,
    xaxis_showgrid=False,
    yaxis_showgrid=False,
    xaxis_zeroline=False,
    yaxis_zeroline=False,
    #yaxis_autorange='reversed',
    template='plotly_white'
    )
    df_biasArray = pd.DataFrame( -np.around(defaultBiasesCol/np.abs(normMetricValsCol), decimals=2),
                   index=metricsNames,
                   columns= ['-Normalized Biases'])
    normlzdBiasArrayFig = ff.create_annotated_heatmap(
                   z=df_biasArray.to_numpy(),
                   x=df_biasArray.columns.tolist(),
                   y=df_biasArray.index.tolist(),
                   #colorscale=normlzdSensMatrixFig.colorscale,
                   #colorscale=px.colors.diverging.balance,
                   #coloraxis="coloraxis",
                   colorscale='balance',
                   zmin=-matMaxAbs, zmax=matMaxAbs,
                   showscale=True, ygap=1, xgap=1
                   )
    normlzdBiasArrayFig.update_layout(
    title_text='', 
    title_x=0.5, 
    #width=10, 
    #height=1400,
    xaxis_showgrid=False,
    yaxis_showgrid=False,
    xaxis_zeroline=False,
    yaxis_zeroline=False,
    #yaxis_autorange='reversed',
    #template='plotly_white'
    )
    sensMatrixBiasFig = make_subplots(
    rows=1, cols=2,
    column_widths=[0.9, 0.1],
    horizontal_spacing=0.2,
    )
    sensMatrixBiasFig.add_trace(normlzdSensMatrixFig.data[0], row=1, col=1)
    #sensMatrixBiasFig.add_trace(normlzdBiasArrayFig.data[0], 1, 1)
    sensMatrixBiasFig.add_trace(normlzdBiasArrayFig.data[0], row=1, col=2)
    sensMatrixBiasFig.update_layout(
        title_text='Color-coded normalized sensitivity matrix',
        height=700,
        template='plotly_white') 
    sensMatrixBiasFig.update_layout(coloraxis=dict(colorscale='RdBu',cmin=-1,cmax=1), showlegend=False)

    #pdb.set_trace()

    cosAnglesMatrix = calcMatrixAngles( normlzdLinplusSensMatrixPoly )
    roundedCosAnglesMatrix = np.around(cosAnglesMatrix, decimals=2)
    df = pd.DataFrame(roundedCosAnglesMatrix,
                  index=metricsNames,
                  columns=metricsNames)
    upTriMask = np.logical_not( np.tril(np.ones_like(roundedCosAnglesMatrix, dtype=bool)) )
    df_mask = df.mask(upTriMask)
    #maskedRoundedCosAnglesMatrix = np.ma.masked_array(roundedCosAnglesMatrix, mask=upTriMask)
    #maskedRoundedCosAnglesMatrix.filled(np.nan)
    #print("maskedAngles =", maskedRoundedCosAnglesMatrix)
    #print("cosAnglesMatrix =", cosAnglesMatrix)
    #print("upTriMask =", upTriMask)
    metricsCorrArrayFig = ff.create_annotated_heatmap(
                   z=df_mask.to_numpy(),
                   x=df_mask.columns.tolist(),
                   y=df_mask.columns.tolist(),
                   colorscale=px.colors.diverging.balance,
                   showscale=True, ygap=1, xgap=1
                   )
    #metricsCorrArrayFig = go.Figure(data=go.Heatmap(
    #                z=roundedCosAnglesMatrix,  
    ##                labels=dict(x="Metrics", y="Metrics")x=['SWCF_GLB', 'SWCF_DYCOMS', 'SWCF_HAWAII', 'SWCF_VOCAL', 'SWCF_VOCAL_near', 'SWCF_LBA', 'SWCF_WP', 'SWCF_EP', 'SWCF_NP', 'SWCF_SP', 'SWCF_CAF', 'SWCF_Nambian', 'SWCF_Nambian_near', 'LWCF_GLB', 'PRECT_GLB'])
    ##                 labels=dict(x="hullo")
    #                x=metricsNames.tolist(),
    #                y=metricsNames.tolist() )
    ##                text_auto=True  )
    #                )
    #    metricsCorrArrayFig = px.imshow(
    #                   img=roundedCosAnglesMatrix,
    #                   x=metricsNames.tolist(),
    #                   y=metricsNames.tolist(),
    #                   color_continuous_scale=px.colors.diverging.balance
    #                    )
    #    metricsCorrArrayFig.update_traces(text=roundedCosAnglesMatrix)
    metricsCorrArrayFig.update_xaxes(side="bottom")
    metricsCorrArrayFig.update_layout(
    title_text='cos(angle) among metrics (i.e., rows of sens matrix)', 
    title_x=0.5, 
    width=800, 
    height=700,
    xaxis_showgrid=False,
    yaxis_showgrid=False,
    xaxis_zeroline=False,
    yaxis_zeroline=False,
    yaxis_autorange='reversed',
    template='plotly_white'
    )

    # Create color-coded matrix that displays correlations among parameter vectors
    normlzdSensMatrixConcatBiases = np.hstack((normlzdLinplusSensMatrixPoly, normlzdDefaultBiasesCol))
    cosAnglesMatrix = calcMatrixAngles( normlzdSensMatrixConcatBiases.T )
    roundedCosAnglesMatrix = np.around(cosAnglesMatrix, decimals=2)
    df = pd.DataFrame(roundedCosAnglesMatrix,
                  index=np.concatenate((paramsNames,['bias vector'])),
                  columns=np.concatenate((paramsNames,['bias vector'])))
    upTriMask = np.logical_not( np.tril(np.ones_like(roundedCosAnglesMatrix, dtype=bool)) )
    df_mask = df.mask(upTriMask)
    paramsCorrArrayFig = ff.create_annotated_heatmap(
                   z=df_mask.to_numpy(),
                   x=df_mask.columns.tolist(),
                   y=df_mask.columns.tolist(),
                   colorscale=px.colors.diverging.balance,
                   showscale=True, ygap=1, xgap=1
                   )
    paramsCorrArrayFig.update_xaxes(side="bottom")
    paramsCorrArrayFig.update_layout(
    title_text='cos(angle) among parameters (i.e., columns of sens matrix)',
    title_x=0.5,
    width=800,
    height=700,
    xaxis_showgrid=False,
    yaxis_showgrid=False,
    xaxis_zeroline=False,
    yaxis_zeroline=False,
    yaxis_autorange='reversed',
    template='plotly_white'
    )

    cosAnglesMatrix = calcMatrixAngles( normlzdSensMatrixPoly )
    invrsCosFactorMinusMatrix = np.power( np.maximum( np.finfo(float).eps, 2. * ( 1. - cosAnglesMatrix ) ) , -0.5 )
    invrsCosFactorPlusMatrix = np.power( np.maximum( np.finfo(float).eps, 2. * ( 1. + cosAnglesMatrix ) ), -0.5 )
    dbOnAbsSensVector = \
        -defaultBiasesCol/np.abs(normMetricValsCol) \
            / np.linalg.norm(normlzdSensMatrixPoly, axis=1).reshape(-1, 1)
    dbOnAbsSensMatrix1 = np.ones((len(metricsNames),1)) @ dbOnAbsSensVector.T
    dbOnAbsSensMatrix2 = dbOnAbsSensVector @ np.ones((1,len(metricsNames)) )
    dpMin2PtMinusMatrix = invrsCosFactorMinusMatrix * \
        np.abs( dbOnAbsSensMatrix2 - dbOnAbsSensMatrix1 )
    dpMin2PtPlusMatrix = invrsCosFactorPlusMatrix * \
        np.abs( dbOnAbsSensMatrix2 + dbOnAbsSensMatrix1 )
    dpMin2PtMatrix = np.maximum( dpMin2PtMinusMatrix, dpMin2PtPlusMatrix )
    roundedDpMin2PtMatrix = np.around(dpMin2PtMatrix, decimals=2)
    dpMin2PtMatrix = np.fill_diagonal(roundedDpMin2PtMatrix, np.nan)
    df = pd.DataFrame(roundedDpMin2PtMatrix,
                  index=metricsNames,
                  columns=metricsNames)
    #upTriMask = np.logical_not( np.tril(np.ones_like(roundedDpMin2PtMatrix, dtype=bool)) )
    upTriMask = np.triu(np.ones_like(roundedDpMin2PtMatrix, dtype=bool))
    df_mask = df.mask(upTriMask)
    dpMin2PtFig = ff.create_annotated_heatmap(
                   z=df_mask.to_numpy(),
                   x=df_mask.columns.tolist(),
                   y=df_mask.columns.tolist(),
                   colorscale=px.colors.sequential.Bluered,
                   showscale=True, ygap=1, xgap=1
                   )
    dpMin2PtFig.update_xaxes(side="bottom")
    dpMin2PtFig.update_layout(
    title_text='dpMin between 2 metrics (i.e., rows of sens matrix)', 
    title_x=0.5, 
    width=800, 
    height=700,
    xaxis_showgrid=False,
    yaxis_showgrid=False,
    xaxis_zeroline=False,
    yaxis_zeroline=False,
    yaxis_autorange='reversed',
    template='plotly_white'
    )




    sensMatrixDashboard.layout = html.Div(children=[
        html.H1(children='Sensitivity matrix diagnostics'),

        html.Div(children=''' '''),
        
        dcc.Graph( id='paramsBarChart', figure=paramsBarChart ),
        dcc.Graph( id='biasesOrderedArrowFig', figure=biasesOrderedArrowFig ),
        dcc.Graph( id='biasTotContrbBarFig', figure=biasTotContrbBarFig ),
#        dcc.Graph( id='paramsTotContrbBarFig', figure=paramsTotContrbBarFig ),
        dcc.Graph( id='sensMatrixBarFig', figure=sensMatrixBarFig ),
        dcc.Graph( id='linplusSensMatrixBarFig', figure=linplusSensMatrixBarFig ),
        dcc.Graph( id='biasLinNlIndivContrbBarFig', figure=biasLinNlIndivContrbBarFig ),
        dcc.Graph( id='biasesVsSensMagScatterplot', figure=biasesVsSensMagScatterplot ),
        dcc.Graph( id='biasVsBiasApproxScatterplot', figure=biasVsBiasApproxScatterplot ),
#  config= { 'toImageButtonOptions': { 'scale': 6 } }
        dcc.Graph( id='sensMatrixBiasFig', figure=sensMatrixBiasFig ),   
        dcc.Graph( id='paramsCorrArrayFig', figure=paramsCorrArrayFig ),
        dcc.Graph( id='metricsCorrArrayFig', figure=metricsCorrArrayFig ),
        dcc.Graph( id='dpMin2PtFig', figure=dpMin2PtFig ),
        dcc.Graph( id='threeDotFig', figure=threeDotFig ),
        dcc.Graph( id='biasesSensScatterFig', figure=biasSensMatrixScatterFig ),
        dcc.Graph( id='dpMinScatterFig', figure=dpMinMatrixScatterFig ),
        dcc.Graph( id='maxSensMetricsFig', figure=maxSensMetricsFig ),
        dcc.Graph( id='residVsBiasFig', figure=residVsBiasFig ),
        dcc.Graph( id='diffBiasesVsSensFig', figure=diffBiasesVsSensFig ),
        dcc.Graph( id='relBiasesVsSensFig', figure=relBiasesVsSensFig ),
        dcc.Graph( id='normlzdSensMatrixColsFig', figure=normlzdSensMatrixColsFig ),
        dcc.Graph( id='normlzdSensMatrixRowsFig', figure=normlzdSensMatrixRowsFig ),
 
        #dcc.Graph( id='paramsFig', figure=paramsFig ),
        #dcc.Graph( id='biasesFig', figure=biasesFig ),
        #dcc.Graph( id='biasContrOrderFig', figure=biasContrOrderFig ),
        #dcc.Graph( id='biasesVsSensArrowFig', figure=biasesVsSensArrowFig ),

        #dcc.Graph( id='vhNormlzdColsFig', figure=vhNormlzdColsFig ),
        #dcc.Graph( id='uNormlzdWeightedColsFig', figure=uNormlzdWeightedColsFig ),
        #dcc.Graph( id='uNormlzdColsFig', figure=uNormlzdColsFig ),
        #dcc.Graph( id='uNormlzdBiasColsFig', figure=uNormlzdBiasColsFig )
    ])

    sensMatrixDashboard.run_server(debug=True)

    return

# Calculate semi-linear matrix, sensMatrix + curvMatrix*dp, for use in forward solution
def normlzdSemiLinMatrixFnc(dnormlzdParams, normlzdSensMatrix, normlzdCurvMatrix, numMetrics):
    import numpy as np

    normlzdSemiLinMatrix = \
        normlzdSensMatrix \
        + 0.5 * normlzdCurvMatrix * ( np.ones((numMetrics,1)) @ dnormlzdParams.T )

    return normlzdSemiLinMatrix

# Calculate forward nonlinear solution, normalized but not weighted
def fwdFnc(dnormlzdParams, normlzdSensMatrix, normlzdCurvMatrix, numMetrics):
    import numpy as np

    normlzdDefaultBiasesApproxNonlin = \
        normlzdSemiLinMatrixFnc(dnormlzdParams, normlzdSensMatrix, normlzdCurvMatrix, numMetrics) @ dnormlzdParams

            #normlzdSensMatrix @ dnormlzdParams \
            #+ 0.5 * normlzdCurvMatrix @ (dnormlzdParams * dnormlzdParams) 

    return normlzdDefaultBiasesApproxNonlin

# Define objective function that is to be minimized.
def objFnc(dnormlzdParams, normlzdSensMatrix, normlzdDefaultBiasesCol, metricsWeights,
           normlzdCurvMatrix, reglrCoef, numMetrics):

    import numpy as np
    import pdb

    dnormlzdParams = np.atleast_2d(dnormlzdParams).T # convert from 1d row array to 2d column array
    chisqd = np.linalg.norm( (-normlzdDefaultBiasesCol \
                              - fwdFnc(dnormlzdParams, normlzdSensMatrix, normlzdCurvMatrix, numMetrics) \
                             ) * metricsWeights \
                                , ord=2 \
                           )**1  \
                + reglrCoef * np.linalg.norm( dnormlzdParams, ord=1 )

    return chisqd

def solveUsingNonlin(metricsNames, paramsNames, transformedParamsNames,
                     metricsWeights, normMetricValsCol, magParamValsRow,
                     sensNcFilenames, sensNcFilenamesExt, defaultNcFilename,
                     defaultParamValsOrigRow,
                     normlzdSensMatrix, defaultBiasesCol,
                     normlzdCurvMatrix,
                     reglrCoef,
                     beVerbose):

    import numpy as np
    import pdb
    from scipy.optimize import minimize


    numMetrics = len(metricsNames)

    # Perform nonlinear optimization
    normlzdDefaultBiasesCol = defaultBiasesCol/np.abs(normMetricValsCol)
    #dnormlzdParamsSolnNonlin = minimize(objFnc,x0=np.ones_like(np.transpose(defaultParamValsOrigRow)), \
    dnormlzdParamsSolnNonlin = minimize(objFnc,x0=np.zeros_like(np.transpose(defaultParamValsOrigRow)), \
    #dnormlzdParamsSolnNonlin = minimize(objFnc,dnormlzdParamsSoln, \
                               args=(normlzdSensMatrix, normlzdDefaultBiasesCol, metricsWeights,
                               normlzdCurvMatrix, reglrCoef, numMetrics),\
                               method='Powell', tol=1e-12)
    dnormlzdParamsSolnNonlin = np.atleast_2d(dnormlzdParamsSolnNonlin.x).T


    dparamsSolnNonlin = dnormlzdParamsSolnNonlin * np.transpose(magParamValsRow)
    paramsSolnNonlin = np.transpose(defaultParamValsOrigRow) + dparamsSolnNonlin
    if beVerbose:
        print("paramsSoln.T=", paramsSoln.T)
        print("paramsSolnNonlin.T=", paramsSolnNonlin.T)
        print("normlzdSensMatrix@dnPS.x.T=", normlzdSensMatrix @ dnormlzdParamsSolnNonlin)
        print("normlzdDefaultBiasesCol.T=", normlzdDefaultBiasesCol.T)
        print("normlzdSensMatrix=", normlzdSensMatrix)

    normlzdWeightedDefaultBiasesApproxNonlin = \
             fwdFnc(dnormlzdParamsSolnNonlin, normlzdSensMatrix, normlzdCurvMatrix, numMetrics) \
             * metricsWeights

    scale = 2
    normlzdWeightedDefaultBiasesApproxNonlin2x = \
             fwdFnc(scale*dnormlzdParamsSolnNonlin, normlzdSensMatrix, 1*normlzdCurvMatrix, numMetrics) \
             * metricsWeights

    # defaultBiasesApprox = (forward model soln - default soln)
    defaultBiasesApproxNonlin = normlzdWeightedDefaultBiasesApproxNonlin \
                                * np.reciprocal(metricsWeights) * np.abs(normMetricValsCol)

    defaultBiasesApproxNonlin2x = normlzdWeightedDefaultBiasesApproxNonlin2x \
                                * np.reciprocal(metricsWeights) * np.abs(normMetricValsCol)

    # To provide error bars, calculate solution with no nonlinear term and double the nonlinear term
    defaultBiasesApproxNonlinNoCurv = \
             fwdFnc(dnormlzdParamsSolnNonlin, normlzdSensMatrix, 0*normlzdCurvMatrix, numMetrics) \
             * np.abs(normMetricValsCol)

    defaultBiasesApproxNonlin2xCurv = \
             fwdFnc(dnormlzdParamsSolnNonlin, normlzdSensMatrix, 2*normlzdCurvMatrix, numMetrics) \
             * np.abs(normMetricValsCol)



    dnormlzdParamsSolnLin = minimize(objFnc,x0=np.zeros_like(np.transpose(defaultParamValsOrigRow)), \
                               args=(normlzdSensMatrix, normlzdDefaultBiasesCol, metricsWeights,
                               0*normlzdCurvMatrix, reglrCoef, numMetrics),\
                               method='Powell')
    dnormlzdParamsSolnLin = np.atleast_2d(dnormlzdParamsSolnLin.x).T
    dparamsSolnLin = dnormlzdParamsSolnLin * np.transpose(magParamValsRow)
    paramsSolnLin = np.transpose(defaultParamValsOrigRow) + dparamsSolnLin


    return (defaultBiasesApproxNonlin, \
            dnormlzdParamsSolnNonlin, paramsSolnNonlin, \
            dnormlzdParamsSolnLin, paramsSolnLin, \
            defaultBiasesApproxNonlin2x, \
            defaultBiasesApproxNonlinNoCurv, defaultBiasesApproxNonlin2xCurv \
           )

def constructNormlzdCurvMatrix(metricsNames, paramsNames, transformedParamsNames,
                               metricsWeights, obsMetricValsCol, normMetricValsCol, magParamValsRow,
                               sens1NcFilenames, sens2NcFilenames, defaultNcFilename):

    """
    For nonlinear 2nd-order term of Taylor series: 0.5*dp^2*d2m/dp2+...,
    construct a numMetrics x numParams matrix of 2nd-order derivatives, d2m/dp2.
    Each row is a different metric.  Each column is a different parameter.
    The matrix is nondimensionalized by the observed values of metrics
    and maximum values of parameters.
    """
    import numpy as np
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    import sys
    import netCDF4
    #import matplotlib.pyplot as plt
    import pdb

    from analyze_sensitivity_matrix import setupSensArrays
    from set_up_dashboard_inputs import setupDefaultParamVectors, \
                                        setupDefaultMetricValsCol
    from scipy.interpolate import UnivariateSpline

    if ( len(paramsNames) != len(sens1NcFilenames)   ):
        print("Number of parameters does not equal number of netcdf files.")
        quit()

    # Number of tunable parameters
    numParams = len(paramsNames)

    # Number of metrics
    numMetrics = len(metricsNames)

    # For use in normalizing metrics matrices
    invrsObsMatrix = np.reciprocal(np.abs(normMetricValsCol)) @ np.ones((1,numParams))

    # Set up a column vector of metric values from the default simulation
    defaultMetricValsCol = \
        setupDefaultMetricValsCol(metricsNames, defaultNcFilename)
    defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1,numParams))
    normlzdDefaultMetricValsMatrix = defaultMetricValsMatrix * invrsObsMatrix

    # Based on the default simulation,
    #    set up a column vector of metrics and a row vector of parameter values.
    defaultParamValsRow, defaultParamValsOrigRow = \
        setupDefaultParamVectors(paramsNames, transformedParamsNames,
                                 numParams,
                                 defaultNcFilename)
    normlzdDefaultParamValsRow = defaultParamValsRow * np.reciprocal(magParamValsRow)

    # Based on the numParams sensitivity simulations,
    #    set up a row vector of modified parameter values.
    # Also set up numMetrics x numParams matrix,
    #    each column of which lists the metrics
    #    from one of the sensitivity simulations
    sens1MetricValsMatrix, sens1ParamValsRow, sens1ParamValsOrigRow = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sens1NcFilenames,
                        beVerbose=False)
    normlzdSens1ParamValsRow = sens1ParamValsRow * np.reciprocal(magParamValsRow)
    normlzdSens1MetricValsMatrix = sens1MetricValsMatrix * invrsObsMatrix

    # Set up sensitivity-simulation matrices from the extended sensitivity simulation
    sens2MetricValsMatrix, sens2ParamValsRow, sens2ParamValsOrigRow = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sens2NcFilenames,
                        beVerbose=False)
    normlzdSens2ParamValsRow = sens2ParamValsRow * np.reciprocal(magParamValsRow)
    normlzdSens2MetricValsMatrix = sens2MetricValsMatrix * invrsObsMatrix

    # Initialize matrix to store second derivatives of metrics w.r.t. parameters
    normlzdCurvMatrix = np.zeros_like(sens1MetricValsMatrix)
    normlzdCurvMatrix2 = np.zeros_like(sens1MetricValsMatrix)  # 2nd way of calculating derivs
    normlzdCurvMatrixPoly = np.zeros_like(sens1MetricValsMatrix)  # 3rd way of calculating derivs
    normlzdSensMatrixPoly = np.zeros_like(sens1MetricValsMatrix)  # Approx of linear sensitivity
    normlzdConstMatrixPoly = np.zeros_like(sens1MetricValsMatrix)  # Approx of linear sensitivity
    normlzdOrdDparamsMin = np.zeros_like(sens1MetricValsMatrix)
    normlzdOrdDparamsMax = np.zeros_like(sens1MetricValsMatrix)

    #pdb.set_trace()

    # Compute quadratic coefficients using a polynomial fit to metric and parameters
    for arrayCol in np.arange(numParams):
        for arrayRow in np.arange(numMetrics):

            # Set up three (x,y) points whose 2nd-order derivative we wish to calculate.
            # For the spline code below, the x points need to be ordered from least to greatest.
            if normlzdSens1ParamValsRow[0,arrayCol] < normlzdDefaultParamValsRow[0,arrayCol] < normlzdSens2ParamValsRow[0,arrayCol]:
                normlzdOrdParams = [ normlzdSens1ParamValsRow[0,arrayCol],
                       normlzdDefaultParamValsRow[0,arrayCol],
                       normlzdSens2ParamValsRow[0,arrayCol] ]
                normlzdOrdDparams = normlzdOrdParams - normlzdDefaultParamValsRow[0,arrayCol]
                normlzdOrdMetrics = [ normlzdSens1MetricValsMatrix[arrayRow,arrayCol],
                        normlzdDefaultMetricValsMatrix[arrayRow,arrayCol],
                        normlzdSens2MetricValsMatrix[arrayRow,arrayCol] ]
            elif normlzdSens2ParamValsRow[0,arrayCol] < normlzdDefaultParamValsRow[0,arrayCol] < normlzdSens1ParamValsRow[0,arrayCol]:
                normlzdOrdParams = [ normlzdSens2ParamValsRow[0,arrayCol],
                       normlzdDefaultParamValsRow[0,arrayCol],
                       normlzdSens1ParamValsRow[0,arrayCol] ]
                normlzdOrdDparams = normlzdOrdParams - normlzdDefaultParamValsRow[0,arrayCol]
                normlzdOrdMetrics = [ normlzdSens2MetricValsMatrix[arrayRow,arrayCol],
                        normlzdDefaultMetricValsMatrix[arrayRow,arrayCol],
                        normlzdSens1MetricValsMatrix[arrayRow,arrayCol] ]
            elif normlzdDefaultParamValsRow[0,arrayCol] < normlzdSens1ParamValsRow[0,arrayCol] < normlzdSens2ParamValsRow[0,arrayCol]:
                normlzdOrdParams = [ normlzdDefaultParamValsRow[0,arrayCol],
                       normlzdSens1ParamValsRow[0,arrayCol],
                       normlzdSens2ParamValsRow[0,arrayCol] ]
                normlzdOrdDparams = normlzdOrdParams - normlzdDefaultParamValsRow[0,arrayCol]
                normlzdOrdMetrics = [ normlzdDefaultMetricValsMatrix[arrayRow,arrayCol],
                        normlzdSens1MetricValsMatrix[arrayRow,arrayCol],
                        normlzdSens2MetricValsMatrix[arrayRow,arrayCol] ]
            elif normlzdDefaultParamValsRow[0,arrayCol] < normlzdSens2ParamValsRow[0,arrayCol] < normlzdSens1ParamValsRow[0,arrayCol]:
                normlzdOrdParams = [ normlzdDefaultParamValsRow[0,arrayCol],
                       normlzdSens2ParamValsRow[0,arrayCol],
                       normlzdSens1ParamValsRow[0,arrayCol] ]
                normlzdOrdDparams = normlzdOrdParams - normlzdDefaultParamValsRow[0,arrayCol]
                normlzdOrdMetrics = [ normlzdDefaultMetricValsMatrix[arrayRow,arrayCol],
                        normlzdSens2MetricValsMatrix[arrayRow,arrayCol],
                        normlzdSens1MetricValsMatrix[arrayRow,arrayCol] ]
            elif normlzdSens1ParamValsRow[0,arrayCol] < normlzdSens2ParamValsRow[0,arrayCol] < normlzdDefaultParamValsRow[0,arrayCol]:
                normlzdOrdParams = [ normlzdSens1ParamValsRow[0,arrayCol],
                       normlzdSens2ParamValsRow[0,arrayCol],
                       normlzdDefaultParamValsRow[0,arrayCol] ]
                normlzdOrdDparams = normlzdOrdParams - normlzdDefaultParamValsRow[0,arrayCol]
                normlzdOrdMetrics = [ normlzdSens1MetricValsMatrix[arrayRow,arrayCol],
                        normlzdSens2MetricValsMatrix[arrayRow,arrayCol],
                        normlzdDefaultMetricValsMatrix[arrayRow,arrayCol] ]
            elif normlzdSens2ParamValsRow[0,arrayCol] < normlzdSens1ParamValsRow[0,arrayCol] < normlzdDefaultParamValsRow[0,arrayCol]:
                normlzdOrdParams = [ normlzdSens2ParamValsRow[0,arrayCol],
                       normlzdSens1ParamValsRow[0,arrayCol],
                       normlzdDefaultParamValsRow[0,arrayCol] ]
                normlzdOrdDparams = normlzdOrdParams - normlzdDefaultParamValsRow[0,arrayCol]
                normlzdOrdMetrics = [ normlzdSens2MetricValsMatrix[arrayRow,arrayCol],
                        normlzdSens1MetricValsMatrix[arrayRow,arrayCol],
                        normlzdDefaultMetricValsMatrix[arrayRow,arrayCol] ]
            else:
                print("Sensitivity parameter values are equal to each other or the default value in constructNormlzdCurvMatrix.")
                print("normlzdSens1ParamValsRow=", normlzdSens1ParamValsRow)
                print("normlzdSens2ParamValsRow=", normlzdSens2ParamValsRow)
                print("normlzdDefaultParamValsRow=", normlzdDefaultParamValsRow)
                pdb.set_trace()
                quit()

            normlzdOrdDparamsMin[arrayRow,arrayCol] = np.min(normlzdOrdDparams)
            normlzdOrdDparamsMax[arrayRow,arrayCol]  = np.max(normlzdOrdDparams)

            # Calculate second-order spline based on three given (x,y) points.
            metricValsSpline = UnivariateSpline(normlzdOrdParams,normlzdOrdMetrics,s=0,k=2)
            # Based on spline, find 2nd derivative at arbitrary point (1).
            # I hope that the derivative has the same value at all points,
            #    since it is a parabola.
            normlzdCurvMatrix[arrayRow,arrayCol] = metricValsSpline.derivative(n=2)(1)

            # Check results using a 2nd calculation
            polyCoefs = np.polyfit(normlzdOrdParams, normlzdOrdMetrics, 2)
            normlzdCurvMatrix2[arrayRow,arrayCol] = 2*polyCoefs[0]
            #pdb.set_trace()

            # Check results using a 3rd calculation
            polyCoefs = np.polyfit(normlzdOrdDparams, normlzdOrdMetrics, 2)
            normlzdCurvMatrixPoly[arrayRow][arrayCol] = 2. * polyCoefs[0]
            normlzdSensMatrixPoly[arrayRow,arrayCol] = polyCoefs[1]
            normlzdConstMatrixPoly[arrayRow,arrayCol] = polyCoefs[2]


    return ( normlzdCurvMatrixPoly, normlzdSensMatrixPoly, normlzdConstMatrixPoly, \
             normlzdOrdDparamsMin, normlzdOrdDparamsMax )

def createBarChart( matrix, index, columns, 
                        orientation,
                        title, 
                        xlabel, ylabel,
                        width, height):

    import plotly.express as px
    import pandas as pd

    df = pd.DataFrame(matrix,
                      index=index,
                      columns=columns)
    barChart = px.bar(df, x=df.index, y=df.columns, orientation=orientation,
                          title = title)
    barChart.update_xaxes(title=xlabel)
    barChart.update_yaxes(title=ylabel)
    barChart.update_layout(hovermode="x")
    barChart.update_layout(showlegend=True)
    barChart.update_yaxes(visible=True,zeroline=True,zerolinewidth=1,zerolinecolor='gray') # Plot x axis
    barChart.update_layout( width=width, height=height  )

    return barChart

def createThreeDotFig(metricsNames, paramsNames, transformedParamsNames,
                               metricsWeights, obsMetricValsCol, normMetricValsCol, magParamValsRow,
                               normlzdCurvMatrix, normlzdSensMatrixPoly, normlzdConstMatrix,
                               normlzdOrdDparamsMin, normlzdOrdDparamsMax,
                               sens1NcFilenames, sens2NcFilenames, defaultNcFilename):

    """
    For nonlinear 2nd-order term of Taylor series: 0.5*dp^2*d2m/dp2+...,
    construct a numMetrics x numParams matrix of 2nd-order derivatives, d2m/dp2.
    Each row is a different metric.  Each column is a different parameter.
    The matrix is nondimensionalized by the observed values of metrics
    and maximum values of parameters.
    """
    import numpy as np
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    import sys
    import netCDF4
    #import matplotlib.pyplot as plt
    import pdb

    from analyze_sensitivity_matrix import setupSensArrays
    from set_up_dashboard_inputs import setupDefaultParamVectors, \
                                        setupDefaultMetricValsCol
    from scipy.interpolate import UnivariateSpline

    if ( len(paramsNames) != len(sens1NcFilenames)   ):
        print("Number of parameters must equal number of netcdf files.")
        quit()

    # Number of tunable parameters
    numParams = len(paramsNames)

    # Number of metrics
    numMetrics = len(metricsNames)

    # Based on the default simulation,
    #    set up a column vector of metrics and a row vector of parameter values.
    defaultParamValsRow, defaultParamValsOrigRow = \
        setupDefaultParamVectors(paramsNames, transformedParamsNames,
                                 numParams,
                                 defaultNcFilename)

    # Set up a column vector of metric values from the default simulation
    defaultMetricValsCol = \
        setupDefaultMetricValsCol(metricsNames, defaultNcFilename)
    defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1,numParams))

    # Based on the numParams sensitivity simulations,
    #    set up a row vector of modified parameter values.
    # Also set up numMetrics x numParams matrix,
    #    each column of which lists the metrics
    #    from one of the sensitivity simulations
    sens1MetricValsMatrix, sens1ParamValsRow, sens1ParamValsOrigRow = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sens1NcFilenames,
                        beVerbose=False)

    # Set up sensitivity-simulation matrices from the extended sensitivity simulation
    sens2MetricValsMatrix, sens2ParamValsRow, sens2ParamValsOrigRow = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sens2NcFilenames,
                        beVerbose=False)

    threeDotFig = make_subplots( rows=numMetrics, cols=numParams,
                                shared_xaxes=True
                                #horizontal_spacing = 0.1/numParams,
                                #vertical_spacing = 0.1/numMetrics
                                )
    for arrayCol in np.arange(numParams):
        for arrayRow in np.arange(numMetrics):

            paramVals = [defaultParamValsRow[0][arrayCol], 
                         sens1ParamValsRow[0][arrayCol], 
                         sens2ParamValsRow[0][arrayCol]]
            metricVals = [defaultMetricValsMatrix[arrayRow][arrayCol], 
                          sens1MetricValsMatrix[arrayRow][arrayCol],
                          sens2MetricValsMatrix[arrayRow][arrayCol]]

            # Plot 3 dots at metric values
            threeDotFig.add_trace(
                go.Scatter(x=paramVals, y=metricVals, 
                               mode='markers',
                               marker=dict(color='black', size=5)),
                           row=arrayRow+1, 
                           col=arrayCol+1
                                  )

            # Plot quadratic curve passing through 3 dots
            paramPts = np.linspace( np.min(paramVals), np.max(paramVals), num=60 )
            #print("paramPts =", paramPts)
            dnormlzdParamPts = np.linspace( normlzdOrdDparamsMin[arrayRow,arrayCol],
                                            normlzdOrdDparamsMax[arrayRow,arrayCol], num=60 )
            metricPts = np.zeros_like(dnormlzdParamPts)
            for idx, dnormlzdParamPt in enumerate(dnormlzdParamPts):
                metricPts[idx] = 0.5 * normlzdCurvMatrix[arrayRow,arrayCol] * dnormlzdParamPt**2 \
                                 + normlzdSensMatrixPoly[arrayRow,arrayCol] * dnormlzdParamPt \
                                 + normlzdConstMatrix[arrayRow,arrayCol]
                metricPts[idx] = metricPts[idx] * np.abs(normMetricValsCol[arrayRow])

            #print("metricPts =", metricPts)
            threeDotFig.add_trace(
                        go.Scatter(x=paramPts, y=metricPts,
                               mode='lines',
                               line=dict(color='blue', width=2)),
                        row=arrayRow+1,
                        col=arrayCol+1
                                  )

            # Calculate horizontal line at observed value
            threeObsMetricVals = np.squeeze(obsMetricValsCol[arrayRow][0]*np.ones((3,1)))
            #threeObsMetricValsList = threeObsMetricVals.tolist()
            #print("obsMetricVals=", threeObsMetricValsList)
            #print("paramVals=", paramVals)            
            threeDotFig.add_trace(
                go.Scatter(x=paramVals, y=threeObsMetricVals, 
                               mode='lines',
                               line=dict(color='red', width=2)),
                           row=arrayRow+1,
                           col=arrayCol+1
                                  )

            # Label the metric and parameter for the subplots
            if (arrayRow == numMetrics-1):  # Put params labels only along bottom of plot
                #threeDotFig.update_xaxes(title=dict(text=paramsNames[arrayCol], 
                #                         tickangle=45),
                threeDotFig.update_xaxes(title_text=paramsNames[arrayCol]\
                                         .replace('clubb_','').replace('c_invrs_tau_','')\
                                         .replace('wpxp_n2','n2').replace('threshold','thresh'),
                                         #title_font_size=8,
                                         tickangle=45,
                                         row=arrayRow+1, col=arrayCol+1
                                         )
                threeDotFig.update_xaxes(tickangle=45,row=arrayRow+1, col=arrayCol+1)
            if (arrayCol == 0): # Insert metrics label only along left edge of plot
                threeDotFig.update_yaxes(title_text=metricsNames[arrayRow], row=arrayRow+1, col=arrayCol+1)
            threeDotFig.update_layout(showlegend=False,
                    title_text="Simulated metric values vs. parameter values for each metric", 
                    height=2500)
            
    threeDotFig.update_xaxes(tickangle=45) # Put params label at 45-degree angle

    return ( threeDotFig )


def calcNormlzdRadiusCurv(metricsNames, paramsNames, transformedParamsNames, paramsScales,
                          metricsWeights, obsMetricValsCol,
                          sensNcFilenames, sensNcFilenamesExt, defaultNcFilename):

    """
    Calculate radius of curvature of output from 2 sensitivity simulations plus the default
    simulation.
    """
    import numpy as np
    import sys
    import netCDF4
    import matplotlib.pyplot as plt
    import pdb

    from analyze_sensitivity_matrix import setupDefaultParamVectors, \
                                           setupSensArrays
    from set_up_dashboard_inputs import setupDefaultMetricValsCol

    if ( len(paramsNames) != len(sensNcFilenames)   ):
        print("Number of parameters must equal number of netcdf files.")
        quit()

    # Number of tunable parameters
    numParams = len(paramsNames)

    # Number of metrics
    numMetrics = len(metricsNames)

    # Set up a column vector of metric values from the default simulation
    defaultMetricValsCol = \
        setupDefaultMetricValsCol(metricsNames, defaultNcFilename)

    # Based on the default simulation,
    #    set up a column vector of metrics and a row vector of parameter values.
    defaultParamValsRow, defaultParamValsOrigRow = \
        setupDefaultParamVectors(paramsNames, transformedParamsNames,
                                 numParams,
                                 defaultNcFilename)

    defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1,numParams))

    # Based on the numParams sensitivity simulations,
    #    set up a row vector of modified parameter values.
    # Also set up numMetrics x numParams matrix,
    #    each column of which lists the metrics
    #    from one of the sensitivity simulations
    sensMetricValsMatrix, sensParamValsRow, sensParamValsOrigRow = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sensNcFilenames,
                        beVerbose=False)

    # Set up sensitivity-simulation matrices from the extended sensitivity simulation
    sensMetricValsMatrixExt, sensParamValsRowExt, sensParamValsOrigRowExt = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sensNcFilenamesExt,
                        beVerbose=False)

    normlzd_radius_of_curv = np.full_like(sensMetricValsMatrix, 0.0)

    # Calculate differences in parameter values between default, sensitivity,
    #    and extended sensitivity runs.
    delta_params_def_sens = sensParamValsRow - defaultParamValsRow
    delta_params_def_sensExt = sensParamValsRowExt - defaultParamValsRow
    delta_params_sens_sensExt = sensParamValsRowExt - sensParamValsRow

    # Calculate numMetrics x numParams matrix of metric values.
    delta_metrics_def_sens = sensMetricValsMatrix - defaultMetricValsMatrix
    delta_metrics_def_sensExt = sensMetricValsMatrixExt - defaultMetricValsMatrix
    delta_metrics_sens_sensExt = sensMetricValsMatrixExt - sensMetricValsMatrix
    for col in np.arange(numParams):
        for row in np.arange(numMetrics):
            # Distance between points in simulations = sqrt(dparam**2 + dmetric**2)
            length_def_sens = np.linalg.norm([delta_params_def_sens[0][col],
                                              delta_metrics_def_sens[row][col]])
            length_def_sensExt = np.linalg.norm([delta_params_def_sensExt[0][col],
                                                 delta_metrics_def_sensExt[row][col]])
            length_sens_sensExt = np.linalg.norm([delta_params_sens_sensExt[0][col],
                                                  delta_metrics_sens_sensExt[row][col]])
            semi_perim = 0.5 * ( length_def_sens + length_def_sensExt + length_sens_sensExt )
            # area of triangle formed by points.  Use Heron's formula.
            area = np.sqrt( semi_perim *
                           (semi_perim-length_def_sens) *
                           (semi_perim-length_def_sensExt) *
                           (semi_perim-length_sens_sensExt)
                          )
            if (area == 0.0):
                print( '\nIn calcNormlzdRadiusCurv, area == 0.0 for param ', paramsNames[col],
                        'and metric ', metricsNames[row] )

            # Greatest distance between parameter values in the 3 simulations:
            max_params_width = \
            np.max(np.abs([delta_params_def_sens[0][col],
                        delta_params_def_sensExt[0][col],
                        delta_params_sens_sensExt[0][col]]))
            if (max_params_width == 0.0):
                print( '\nIn calcNormlzdRadiusCurv, max_params_width == 0.0 for param ', paramsNames[col],
                        'and metric ', metricsNames[row] )

            # Calculate Menger curvature from triangle area and distance between points:
            normlzd_radius_of_curv[row][col] = 0.25 * length_def_sens*length_def_sensExt*length_sens_sensExt \
                                                / area / max_params_width

    #pdb.set_trace()
    fig, axs = plt.subplots(numMetrics, numParams, figsize=(24,36))
    for col in np.arange(numParams):
        for row in np.arange(numMetrics):

            paramVals = [defaultParamValsRow[0][col], sensParamValsRow[0][col], sensParamValsRowExt[0][col]]
            metricVals = [defaultMetricValsMatrix[row][col], sensMetricValsMatrix[row][col],
                  sensMetricValsMatrixExt[row][col]]

            axs[row, col].plot( paramVals, metricVals, marker=".", ls="" )
            axs[row, col].plot( paramVals, obsMetricValsCol[row][0] * np.ones((3,1)), color="r" )
            axs[row, col].set_xlabel(paramsNames[col])
            axs[row, col].set_ylabel(metricsNames[row])
            #fig.show()

    plt.show()
    plt.savefig('param_metric_scatter.png')
    #pdb.set_trace()

    return


def calcParamsBounds(metricsNames, paramsNames, transformedParamsNames,
                          metricsWeights, obsMetricValsCol, normMetricValsCol,
                          magParamValsRow,
                          sensNcFilenames, sensNcFilenamesExt, defaultNcFilename):

    """
    Calculate the maximum parameter perturbations based on the non-linearity of the global model
    simulation.
    """
    import numpy as np
    import sys
    import netCDF4
    #import matplotlib.pyplot as plt
    import pdb

    from analyze_sensitivity_matrix import setupSensArrays, calcSvdInvrs, calcParamsSoln
    from set_up_dashboard_inputs import setupDefaultParamVectors, \
                                        setupDefaultMetricValsCol

    if ( len(paramsNames) != len(sensNcFilenames)   ):
        print("Number of parameters must equal number of netcdf files.")
        quit()

    # Number of tunable parameters
    numParams = len(paramsNames)

    # Number of metrics
    numMetrics = len(metricsNames)

    # Set up a column vector of metric values from the default simulation
    defaultMetricValsCol = \
        setupDefaultMetricValsCol(metricsNames, defaultNcFilename)

    # Based on the default simulation,
    #    set up a column vector of metrics and a row vector of parameter values.
    defaultParamValsRow, defaultParamValsOrigRow = \
        setupDefaultParamVectors(paramsNames, transformedParamsNames,
                                 numParams,
                                 defaultNcFilename)

    defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1,numParams))

    defaultBiasesCol, sensMatrix, normlzdWeightedSensMatrix = \
        calcParamsBoundsHelper(metricsNames, paramsNames, transformedParamsNames,
                           metricsWeights, obsMetricValsCol, normMetricValsCol,
                           numMetrics, numParams,
                           magParamValsRow,
                           defaultMetricValsCol, defaultParamValsRow, defaultParamValsOrigRow,
                           sensNcFilenames)

    defaultBiasesColExt, sensMatrixExt, normlzdWeightedSensMatrixExt = \
        calcParamsBoundsHelper(metricsNames, paramsNames, transformedParamsNames,
                           metricsWeights, obsMetricValsCol, normMetricValsCol,
                           numMetrics, numParams,
                           magParamValsRow,
                           defaultMetricValsCol, defaultParamValsRow, defaultParamValsOrigRow,
                           sensNcFilenamesExt)

    normlzdWeightedSensMatrixDiff = normlzdWeightedSensMatrixExt-normlzdWeightedSensMatrix
    sensMatrixDiff = sensMatrixExt-sensMatrix


    # sValsRatio = a threshold ratio of largest singular value to a smaller singular value.
    # If sValsRatio is large enough, then all singular vectors will be kept.
    # If sValsRatio is 1, then only the first singular vector will be kept.
    sValsRatio = 300.

    # Calculate inverse of the singular value decomposition.
    # This gives the recommended changes to parameter values.
    svdInvrsNormlzdWeighted, svdInvrsNormlzdWeightedPC, \
    sValsTruncInvNormlzdWeighted, sValsTruncInvNormlzdWeightedPC, \
    vhNormlzdWeighted, uNormlzdWeighted, sNormlzdWeighted = \
         calcSvdInvrs(normlzdWeightedSensMatrixDiff, sValsRatio, beVerbose=False)

    paramsSolnPC, paramsLowValsPC, paramsHiValsPC, dparamsSolnPCBound, dnormlzdParamsSolnPC, \
    defaultBiasesApproxPC, defaultBiasesApproxLowValsPC, \
    defaultBiasesApproxHiValsPC = \
             calcParamsSoln(svdInvrsNormlzdWeightedPC, metricsWeights, magParamValsRow, \
                            sensMatrixDiff, normlzdWeightedSensMatrixDiff, \
                            obsMetricValsCol, normMetricValsCol, defaultBiasesCol,
                            defaultParamValsOrigRow, \
                            sValsTruncInvNormlzdWeightedPC,
                            vhNormlzdWeighted, \
                            numParams, paramsNames,
                            transformedParamsNames )

    paramsLowValsPCBound = defaultParamValsOrigRow.T - 0.5*np.abs(dparamsSolnPCBound)
    paramsHiValsPCBound  = defaultParamValsOrigRow.T + 0.5*np.abs(dparamsSolnPCBound)

    #pdb.set_trace()

    return ( paramsLowValsPCBound, paramsHiValsPCBound )

def calcParamsBoundsHelper(metricsNames, paramsNames, transformedParamsNames,
                           metricsWeights, obsMetricValsCol, normMetricValsCol,
                           numMetrics, numParams,
                           magParamValsRow,
                           defaultMetricValsCol, defaultParamValsRow, defaultParamValsOrigRow,
                           sensNcFilenames):

    """
    Calculate the maximum parameter perturbations based on the non-linearity of the global model
    simulation.
    """
    import numpy as np
    import sys
    import netCDF4
    #import matplotlib.pyplot as plt
    import pdb

    from analyze_sensitivity_matrix import setupSensArrays, \
                                           constructSensMatrix, calcSvdInvrs, calcParamsSoln
    from set_up_dashboard_inputs import setupDefaultMetricValsCol

    # Based on the numParams sensitivity simulations,
    #    set up a row vector of modified parameter values.
    # Also set up numMetrics x numParams matrix,
    #    each column of which lists the metrics
    #    from one of the sensitivity simulations
    sensMetricValsMatrix, sensParamValsRow, sensParamValsOrigRow = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sensNcFilenames,
                        beVerbose=False)

    # Matrix of metric values from default simulation
    # Each column in the matrix is repeated numParams times, for later multiplication
    defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1,numParams))

    # Calculate the sensitivity matrix and the sensitivity matrix
    # normalized by the discrepancies from observations in default simulation.
    # Use transformed parameter values.
    defaultBiasesCol, sensMatrix, normlzdSensMatrix, biasNormlzdSensMatrix = \
             constructSensMatrix(sensMetricValsMatrix, sensParamValsRow,
                                 defaultMetricValsMatrix, defaultParamValsRow,
                                 defaultMetricValsCol,
                                 magParamValsRow,
                                 obsMetricValsCol, normMetricValsCol,
                                 numMetrics, numParams,
                                 beVerbose=False)

    # In order to weight certain metrics, multiply each row of normlzdSensMatrix
    # by metricsWeights
    normlzdWeightedSensMatrix = np.diag(np.transpose(metricsWeights)[0]) @ normlzdSensMatrix

    # Calculate inverse of the singular value decomposition.
    # This gives the recommended changes to parameter values.
    #svdInvrsNormlzdWeighted, svdInvrsNormlzdWeightedPC, \
    #sValsTruncInvNormlzdWeighted, sValsTruncInvNormlzdWeightedPC, \
    #vhNormlzdWeighted, uNormlzdWeighted, sNormlzdWeighted = \
    #     calcSvdInvrs(normlzdWeightedSensMatrix)

    #paramsSolnPC, paramsLowValsPC, paramsHiValsPC, dparamsSolnPC, dnormlzdParamsSolnPC, \
    #defaultBiasesApproxPC, defaultBiasesApproxLowValsPC, \
    #defaultBiasesApproxHiValsPC = \
    #         calcParamsSoln(svdInvrsNormlzdWeightedPC, metricsWeights, magParamValsRow, \
    #                        sensMatrix, normlzdWeightedSensMatrix, \
    #                        obsMetricValsCol, normMetricValsCol, defaultBiasesCol,
    #                        defaultParamValsOrigRow, \
    #                        sValsTruncInvNormlzdWeightedPC,
    #                        vhNormlzdWeighted, \
    #                        numParams, paramsNames,
    #                        transformedParamsNames )


    #pdb.set_trace()

    return ( defaultBiasesCol, sensMatrix, normlzdWeightedSensMatrix )

def calcMatrixVectorAngles( matrix, row ):
    '''Calculate cos(angle) between one row of a matrix and all rows of the same matrix.
       Returns a column vector, with length equal to the number of rows in the matrix.'''
    
    import sklearn
    
    normed_matrix = sklearn.preprocessing.normalize( matrix, axis=1, norm='l2' )
    
    cosAngles = normed_matrix @ normed_matrix[row, :].T
    
    return cosAngles

def calcMatrixAngles( matrix ):
    '''Calculate cos(angle) among all rows of the same matrix.'''
    
    import sklearn
    
    normed_matrix = sklearn.preprocessing.normalize( matrix, axis=1, norm='l2' )
    
    cosAnglesMatrix = normed_matrix @ normed_matrix.T
    
    return cosAnglesMatrix

def createPcaBiplot(normlzdSensMatrix, defaultBiasesCol, normMetricValsCol, metricsNames, paramsNames):

    import numpy as np
    from pca import pca
    import pdb

    # reduce the data towards 2 PCs
    model = pca(n_components=2, detect_outliers='ht2')

    # Augmented array with LHS and RHS
    augMatrix = np.concatenate((normlzdSensMatrix, -defaultBiasesCol / np.abs(normMetricValsCol) ), axis=1)

    paramsList = list(paramsNames)
    paramsList.append('dbias')
    augParamsNames = np.asarray(paramsList)

    # Fit transform
    results = model.fit_transform(augMatrix, row_labels=metricsNames, col_labels=augParamsNames)
    #results = model.fit_transform(normlzdSensMatrix, row_labels=metricsNames, col_labels=paramsNames)

    PC_test = model.transform(augMatrix)
    #PC_test = model.transform(normlzdSensMatrix)
    outliers, outliers_params = model.compute_outliers(PC=PC_test)
    print("PCA outliers = ", outliers)

    # Plot explained variance
    #fig, ax = model.plot()

    # Scatter first 2 PCs
    #fig, ax = model.scatter()

    # Make biplot with the number of features
    fig, ax = model.biplot(n_feat=paramsNames.size+1)
    #fig, ax = model.biplot(n_feat=paramsNames.size)

    #pdb.set_trace()

    return

if __name__ == '__main__':
    main()
#        sensMatrixDashboard.run_server(debug=True)
