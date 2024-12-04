# -*- coding: utf-8 -*-

# Run this app with `python3 sens_matrix_dashboard.py` and
# view the plots at http://127.0.0.1:8050/ in your web browser.
# (To open a web browser on a larson-group computer,
# login to malan with `ssh -X` and then type `firefox &`.)

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import plotly.colors as pc
#import pca
import dash
from dash import dcc
#import dash_core_components as dcc
from dash import html
#import dash_html_components as html
#import fnmatch
from sens_matrix_dashboard import lossFncMetrics, approxMatrixWithSvd

def createFigs(metricsNames,
               paramsNames, transformedParamsNames, paramsScales,
               metricsWeights, obsMetricValsCol, normMetricValsCol, magParamValsRow,
               defaultBiasesCol, defaultBiasesApproxNonlin, defaultBiasesApproxElastic,
               defaultBiasesApproxNonlinNoCurv, defaultBiasesApproxNonlin2xCurv,
               normlzdDefaultBiasesCol,
               normlzdCurvMatrix, normlzdSensMatrixPoly, normlzdConstMatrix,
               normlzdOrdDparamsMin, normlzdOrdDparamsMax,
               normlzdWeightedSensMatrixPoly,
               dnormlzdParamsSolnNonlin,
               defaultParamValsOrigRow,
               linSolnBiasesCol, normlzdLinplusSensMatrixPoly,
               paramsSolnLin, dnormlzdParamsSolnLin,
               paramsSolnNonlin,
               paramsSolnElastic, dnormlzdParamsSolnElastic,
               sensNcFilenames, sensNcFilenamesExt, defaultNcFilename,
               beVerbose):
    ##############################################
    #
    #    Create plots
    #
    ##############################################

    print("Creating plots . . .")

    downloadConfig = {
        'toImageButtonOptions': {
            'format': 'png',  # 'svg', 'jpeg', 'webp'
            'filename': 'plotly',
            #'width': 800,
            #'height': 600,
            'scale': 4  # Increase resolution
        }
    }

    normlzdResid = (-defaultBiasesApproxNonlin - defaultBiasesCol)[:, 0] \
                   / np.abs(normMetricValsCol[:, 0])

    # Use these flags to determine whether or not to create specific plots
    plot_paramsErrorBarsFig = True
    plot_biasesOrderedArrowFig = False  #True
    plot_threeDotFig = True
    plot_metricsBarChart = True
    plot_paramsIncrsBarChart = True
    plot_paramsAbsIncrsBarChart = True
    plot_paramsTotContrbBarChart = False
    plot_biasesVsDiagnosticScatterplot = False
    plot_dpMin2PtFig = True
    plot_dpMinMatrixScatterFig = True
    plot_projectionMatrixFigs = False #True
    plot_biasesVsSensMagScatterplot = True
    plot_biasesVsSvdScatterplot = True
    plot_paramsCorrArrayFig = True
    plot_sensMatrixAndBiasVecFig = True
    plot_PcaBiplot = False
    plot_PcSensMap = True
    plot_vhMatrixFig = True

    # Remove prefixes from CLUBB variable names in order to shorten them
    paramsAbbrv = abbreviateClubbParamsNames(paramsNames)
    paramsAbbrvPlusBias = np.append(paramsAbbrv, 'bias').tolist()

    metricsNamesNoprefix = np.char.replace(metricsNames, "SWCF_", "")

    # Create a way to order the metrics by sensitivity, for later use in plots
    metricsSens = np.linalg.norm(normlzdWeightedSensMatrixPoly, axis=1)  # measure of sensitivity of each metric
    # metricsSensOrder = (rankdata(metricsSens) - 1).astype(int)  # this ordering doesn't work as an index
    metricsSensOrder = metricsSens.argsort()
    metricsNamesOrdered = metricsNames[metricsSensOrder]  # list of metrics names, ordered from least to most sensitive

    # These are the metrics that we want to include in the plots (whitelisted variables).
    # I.e., set whitelistedMetricsMask=T for variables that we want to plot.
    whitelistedMetricsMask = ((metricsNames == 'SWCF_4_12') \
                              | (metricsNames == 'SWCF_4_1') \
                              | (metricsNames == 'SWCF_2_17') \
                              | (metricsNames == 'SWCF_6_14') \
                              #| (metricsNames == 'SWCF_6_18') \
                              #| (metricsNames == 'SWCF_6_1') \
                              #| (metricsNames == 'SWCF_6_10') \
                              #| (metricsNames == 'SWCF_6_13') \
                              #| (metricsNames == 'SWCF_6_15') \
                              #| (metricsNames == 'SWCF_2_14') \
                              #| (metricsNames == 'SWCF_2_15') \
                              #| (metricsNames == 'SWCF_2_10') \
                              #| (metricsNames == 'SWCF_7_14') \
                              #| (metricsNames == 'SWCF_7_2') \
                              #| (metricsNames == 'SWCF_5_8') \
                              #| (metricsNames == 'SWCF_1_14') \
                              #| (metricsNames == 'SWCF_1_15') \
                              #| (metricsNames == 'SWCF_4_15') \
                              # | (metricsNames == 'SWCF_5_16') \
                              # | (metricsNames == 'SWCF_1_13') \
                              #| (metricsNames == 'SWCF_6_18') \
                              #| (metricsNames == 'SWCF_4_16') \
                              #| (metricsNames == 'SWCF_6_7') \
                              #| (metricsNames == 'SWCF_9_13') \
                              #| (metricsNames == 'SWCF_9_6') \
                              #| (metricsNames == 'SWCF_9_5') \
                              #| (metricsNames == 'SWCF_6_4') \
                              #| (metricsNames == 'SWCF_6_6') \
                              #| (metricsNames == 'SWCF_5_5') \
                              #| (metricsNames == 'SWCF_1_8') \
                              #| (metricsNames == 'SWCF_1_9') \
                              #| (metricsNames == 'LWCF_6_18') \
                              #| (metricsNames == 'LWCF_9_18') \
                              #| (metricsNames == 'LWCF_1_7')  \
                              | (metricsNames == 'PRECT_9_4') \
                              | (metricsNames == 'PRECT_9_9') \
                              | (metricsNames == 'PRECT_9_8') \
                              | (metricsNames == 'PRECT_4_14') \
                              | (metricsNames == 'PRECT_6_7') \
                              | (metricsNames == 'PRECT_6_6') \
                              | (metricsNames == 'PRECT_6_13') \
                              | (metricsNames == 'PRECT_6_14') \
                              | (metricsNames == 'PRECT_4_3') \
                              )

    # Use this line if you want to exclude (blacklist) some variables:
    #whitelistedMetricsMask = (metricsNames != 'SWCF_4_1') & (metricsNames != 'SWCF_4_2')

    # Use this line if you want to include all metrics:
    #whitelistedMetricsMask = (metricsNames != 'noRealMetricWouldHaveThisName')

    # These are the params that we want to include in the plots (whitelisted variables).
    # I.e., set maskParamsNames=T for variables that we want to plot.
    #maskParamsNames = (paramsNames == 'clubb_c_k10') | (paramsNames == 'clubb_altitude_threshold') \
    #                  | (paramsNames == 'clubb_c_invrs_tau_shear') | (paramsNames == 'clubb_c_invrs_tau_bkgnd') \
    #                  | (paramsNames == 'clubb_c_invrs_tau_n2_wp2') | (paramsNames == 'clubb_c_invrs_tau_n2_wp2') \
    #                  | (paramsNames == 'clubb_c_invrs_tau_n2')

    # Use this line if you want to exclude (blacklist) some variables:
    #maskParamsNames = (paramsNames != 'clubb_c_invrs_tau_n2') & (paramsNames != 'clubb_c_k10')

    tunedLossCol = lossFncMetrics(dnormlzdParamsSolnNonlin, normlzdSensMatrixPoly,
                   normlzdDefaultBiasesCol, metricsWeights,
                   normlzdCurvMatrix, len(metricsNames))

    defaultLossCol = lossFncMetrics(np.zeros_like(dnormlzdParamsSolnNonlin), normlzdSensMatrixPoly,
                   normlzdDefaultBiasesCol, metricsWeights,
                   normlzdCurvMatrix, len(metricsNames))

    tunedLossChange = tunedLossCol - defaultLossCol

    numImprovedMetrics = 8
    mostImprovedIdxs = np.argpartition(tunedLossChange, numImprovedMetrics, axis=0)
    whitelistedMetricsMask = np.zeros_like(metricsNames, dtype=bool) # Initialize to False
    whitelistedMetricsMask[mostImprovedIdxs[:numImprovedMetrics, 0]] = True
    mostDegradedIdxs = np.argpartition(-tunedLossChange, numImprovedMetrics, axis=0)
    whitelistedMetricsMask[mostDegradedIdxs[:numImprovedMetrics, 0]] = True

    # Use this line if you want to include all params:
    maskParamsNames = (paramsNames != 'noRealParamWouldHaveThisName')

    #    #whitelistedMetricsMask = np.logical_not(whitelistedMetricsMask)  # get rid of named elements
    # Apply mask
    metricsNamesMasked = np.ma.masked_array(metricsNames, np.logical_not(whitelistedMetricsMask)).compressed()
    paramsNamesMasked = np.ma.masked_array(paramsNames, np.logical_not(maskParamsNames)).compressed()
    normMetricValsColMasked = normMetricValsCol[whitelistedMetricsMask]
    defaultBiasesColMasked = defaultBiasesCol[whitelistedMetricsMask]
    normlzdDefaultBiasesColMasked = normlzdDefaultBiasesCol[whitelistedMetricsMask]
    defaultBiasesApproxNonlinMasked = defaultBiasesApproxNonlin[whitelistedMetricsMask]
    normlzdSensMatrixPolyMasked = normlzdSensMatrixPoly[whitelistedMetricsMask].T[maskParamsNames].T
    normlzdLinplusSensMatrixPolyMetricsMasked = normlzdLinplusSensMatrixPoly[whitelistedMetricsMask]
    normlzdLinplusSensMatrixPolyParamsMetricsMasked = normlzdLinplusSensMatrixPolyMetricsMasked.T[maskParamsNames].T

    defaultBiasesApproxNonlinNoCurvMasked = defaultBiasesApproxNonlinNoCurv[whitelistedMetricsMask]
    defaultBiasesApproxNonlin2xCurvMasked = defaultBiasesApproxNonlin2xCurv[whitelistedMetricsMask]

    linSolnBiasesColMasked = linSolnBiasesCol[whitelistedMetricsMask]

    metricsWeightsMasked = metricsWeights[whitelistedMetricsMask]

    obsMetricValsColMasked = obsMetricValsCol[whitelistedMetricsMask]
    normlzdCurvMatrixMasked = normlzdCurvMatrix[whitelistedMetricsMask].T[maskParamsNames].T
    normlzdConstMatrixMasked = normlzdConstMatrix[whitelistedMetricsMask].T[maskParamsNames].T
    normlzdOrdDparamsMinMasked = normlzdOrdDparamsMin[whitelistedMetricsMask].T[maskParamsNames].T
    normlzdOrdDparamsMaxMasked = normlzdOrdDparamsMax[whitelistedMetricsMask].T[maskParamsNames].T

    magParamValsRowMasked = magParamValsRow.T[maskParamsNames].T

    sensNcFilenamesMasked = np.ma.masked_array(sensNcFilenames, np.logical_not(maskParamsNames)).compressed()
    sensNcFilenamesExtMasked = np.ma.masked_array(sensNcFilenamesExt, np.logical_not(maskParamsNames)).compressed()

    metricsSensMasked = metricsSens[whitelistedMetricsMask]
    metricsSensMaskedOrder = metricsSensMasked.argsort()

    metricsNamesMaskedOrdered = metricsNamesMasked[metricsSensMaskedOrder]
    normMetricValsColMaskedOrdered = normMetricValsColMasked[metricsSensMaskedOrder, 0]
    defaultBiasesColMaskedOrdered = defaultBiasesColMasked[metricsSensMaskedOrder, 0]
    defaultBiasesApproxNonlinMaskedOrdered = defaultBiasesApproxNonlinMasked[metricsSensMaskedOrder, 0]

    #    metricsSensOrderMasked = metricsSensOrder[whitelistedMetricsMask]

    normlzdSensMatrixOrdered = normlzdSensMatrixPoly[metricsSensOrder, :]
    ## Form matrix of parameter perturbations, for later multiplication into the sensitivity matrix
    dnormlzdParamsSolnNonlinMatrix = np.ones((len(metricsNames), 1)) @ dnormlzdParamsSolnNonlin.T
    normlzdSensParamsMatrixOrdered = normlzdSensMatrixOrdered * dnormlzdParamsSolnNonlinMatrix
    normlzdSensParamsMatrix= normlzdSensMatrixPoly * dnormlzdParamsSolnNonlinMatrix

    # Create plot showing lumped linear+nonlinear contributions to each metric
    # Form matrix of parameter perturbations, for later multiplication into the sensitivity matrix
    #dnormlzdParamsSolnNonlinMatrix = np.ones((len(metricsNames),1)) @ dnormlzdParamsSolnNonlin.T
    curvParamsMatrixOrdered = 0.5 * normlzdCurvMatrix[metricsSensOrder, :] * dnormlzdParamsSolnNonlinMatrix ** 2
    curvParamsMatrix = 0.5 * normlzdCurvMatrix * dnormlzdParamsSolnNonlinMatrix ** 2
    #print("Sum rows=", np.sum(-normlzdSensParamsMatrixOrdered-curvParamsMatrixOrdered, axis=1))
    minusNonlinMatrixDparamsOrdered = -1 * curvParamsMatrixOrdered + -1 * normlzdSensParamsMatrixOrdered
    minusNonlinMatrixDparamsOrderedMasked = \
        minusNonlinMatrixDparamsOrdered[whitelistedMetricsMask[metricsSensOrder]]
    nonlinMatrixDparams = curvParamsMatrix + normlzdSensParamsMatrix
    nonlinMatrixDparamsMasked = nonlinMatrixDparams[whitelistedMetricsMask]
    nonlinMatrixDparamsMasked = nonlinMatrixDparamsMasked \
                    * ( np.abs(normMetricValsCol[whitelistedMetricsMask]) @ np.ones_like(paramsSolnNonlin.T) ) \
                    + ( defaultBiasesColMasked + obsMetricValsColMasked ) @ np.ones_like(paramsSolnNonlin.T)
    #minusNonlinMatrixDparamsMasked = ( defaultBiasesColMasked + obsMetricValsColMasked ) @ np.ones_like(paramsSolnNonlin.T)


    if plot_paramsErrorBarsFig:
        print("Creating paramsErrorBarsFig . . .")
        # Calculate symmetric error bars on fitted parameter values,
        #    based on difference in sensitivity matrix, i.e., based on size of nonlinear terms.
        #    For use in figures such as paramsErrorBarsFig.
        paramsLowValsPCBound, paramsHiValsPCBound = \
            calcParamsBounds(metricsNames, paramsNames, transformedParamsNames,
                             metricsWeights, obsMetricValsCol, normMetricValsCol,
                             magParamValsRow,
                             sensNcFilenames, sensNcFilenamesExt, defaultNcFilename)
        paramsErrorBarsFig = \
            createParamsErrorBarsFig(paramsAbbrv, defaultParamValsOrigRow, paramsScales,
                                     paramsLowValsPCBound, paramsHiValsPCBound,
                                     paramsSolnLin, dnormlzdParamsSolnLin,
                                     paramsSolnNonlin, dnormlzdParamsSolnNonlin,
                                     paramsSolnElastic, dnormlzdParamsSolnElastic)

    if plot_threeDotFig:
        print("Creating threeDotFig . . .")
        #    threeDotFig = \
        #        createThreeDotFig(metricsNames, paramsNames, transformedParamsNames,
        #                          metricsWeights, obsMetricValsCol, normMetricValsCol, magParamValsRow,
        #                          normlzdCurvMatrix, normlzdSensMatrixPoly, normlzdConstMatrix,
        #                          normlzdOrdDparamsMin, normlzdOrdDparamsMax,
        #                          sensNcFilenames, sensNcFilenamesExt, defaultNcFilename)
        threeDotFig = \
            createThreeDotFig(metricsNamesMasked, paramsNamesMasked, transformedParamsNames,
                              paramsAbbrv,
                              paramsSolnNonlin,
                              nonlinMatrixDparamsMasked,
                              metricsWeightsMasked, obsMetricValsColMasked, normMetricValsColMasked,
                              magParamValsRowMasked,
                              normlzdCurvMatrixMasked, normlzdSensMatrixPolyMasked, normlzdConstMatrixMasked,
                              normlzdOrdDparamsMinMasked, normlzdOrdDparamsMaxMasked,
                              sensNcFilenamesMasked, sensNcFilenamesExtMasked, defaultNcFilename)

    if plot_biasesOrderedArrowFig:
        print("Creating biasesOrderedArrowFig . . .")
        biasesOrderedArrowFig = \
            createBiasesOrderedArrowFig(metricsSensMaskedOrder, metricsNamesMaskedOrdered,
                                        defaultBiasesColMasked, normMetricValsColMasked,
                                        defaultBiasesApproxNonlinNoCurvMasked, defaultBiasesApproxNonlin2xCurvMasked,
                                        defaultBiasesApproxNonlinMasked,
                                        linSolnBiasesColMasked)
        #createBiasesOrderedArrowFig(metricsSensOrder, metricsNamesOrdered,
        #                            defaultBiasesCol, normMetricValsCol,
        #                            defaultBiasesApproxNonlinNoCurv, defaultBiasesApproxNonlin2xCurv,
        #                            defaultBiasesApproxNonlin,
        #                            linSolnBiasesCol)


    # Nota bene on masked ordering:
    # if metricsNamesOrdered = metricsNames[metricsSensOrder], then
    # metricsNamesOrdered[whitelistedMetricsMask[metricsSensOrder]] gives the same list and order as
    # metricsNamesMasked[metricsSensMaskedOrder]
    # metricsNames[whitelistedMetricsMask] gives the same list, but in a different order

    if plot_paramsTotContrbBarChart:
        paramsTotContrbBarChart = \
            createBarChart(minusNonlinMatrixDparamsOrdered.T, index=paramsNames, columns=metricsNamesOrdered,
                           orientation='v',
                           title="""Linear + nonlinear contributions to removal of biases <br>
                                   (dMetrics=sensMatrix*dParams)""",
                           xlabel="Parameter", ylabel="Contribution to bias removal",
                           width=800, height=500)

    if False:
        linplusSensMatrixBarFig = \
            createBarChart(normlzdLinplusSensMatrixPoly[metricsSensOrder, :].T, index=paramsNames,
                           columns=metricsNamesOrdered,
                           #                barBase = np.zeros_like(paramsScales),
                           orientation='v',
                           title="""Contributions to columns of linplus sensitivity matrix""",
                           xlabel="Parameter", ylabel="Contribution to column",
                           width=800, height=500)

    if False:
        metricsCorrArrayFig = createCorrArrayFig(normlzdLinplusSensMatrixPoly, metricsNames,
                                                 title='cos(angle) among metrics (i.e., rows of sens matrix)')

    if plot_metricsBarChart:
        print("Creating metricsBarChart . . .")
        #    minusNormlzdDefaultBiasesCol = \
        #             -defaultBiasesCol[metricsSensOrder,0]/np.abs(normMetricValsCol[metricsSensOrder,0])
        #    normlzdResidBias = (-defaultBiasesApproxNonlin-defaultBiasesCol)[metricsSensOrder,0] \
        #                       / np.abs(normMetricValsCol[metricsSensOrder,0])
        #    metricsBarChart = createMetricsBarChart(metricsNames[metricsSensOrder],paramsNames,
        #                                            minusNormlzdDefaultBiasesCol, normlzdResidBias, minusNonlinMatrixDparamsOrdered,
        #                                            title='Removal of biases in each metric by each parameter')
        metricsBarChart = createMetricsBarChart(metricsNamesMaskedOrdered, paramsNames,
                                                defaultBiasesColMaskedOrdered, defaultBiasesApproxNonlinMaskedOrdered,
                                                normMetricValsColMaskedOrdered,
                                                minusNonlinMatrixDparamsOrderedMasked,
                                                title='Removal of biases in each metric by each parameter')

    if plot_paramsIncrsBarChart:

        print("Creating paramsIncrBarChart . . .")

        u, s, vh = \
            np.linalg.svd(normlzdLinplusSensMatrixPoly,
                          full_matrices=False)

        sumParamsIncrs = np.mean(minusNonlinMatrixDparamsOrdered, axis=0)
        calcBias = np.sum( sumParamsIncrs )
        sumParamsIncrsPlusBiases = np.append(sumParamsIncrs, calcBias)
        sumParamsIncrsPlusBiases = np.append(sumParamsIncrsPlusBiases, np.mean(normlzdDefaultBiasesCol))

        paramsAbbrvPlusBiases = np.append(paramsAbbrv, 'calc bias')
        paramsAbbrvPlusBiases = np.append(paramsAbbrvPlusBiases, 'obs bias')

        barColorsPlusObs = np.append( ["blue"]*len(paramsNames), np.array(["orange", "red"]) )
        colorDiscreteMap = dict(zip(paramsAbbrvPlusBiases,
                                    barColorsPlusObs.tolist() ))

        paramsIncrsBarChart = \
            createBarChart(sumParamsIncrsPlusBiases, index=paramsAbbrvPlusBiases,
                           columns=["Parameter<br>contribution"],
                           orientation='v',
                           colorDiscreteMap=colorDiscreteMap,
                           title="Mean parameter contributions to removal of biases <br>" \
                                   + "(= column means of sensMatrix*dParams+0.5*curvMatrix*dParams**2)",
                           xlabel="Parameter", ylabel="Contribution to bias removal",
                           width=800, height=500,
                           showLegend=False)

    if plot_paramsIncrsBarChart:

        absParamsIncrs = np.mean( np.abs(minusNonlinMatrixDparamsOrdered), axis=0 )
        absParamsIncrsPlusBias = np.append(absParamsIncrs, np.mean( np.abs(normlzdDefaultBiasesCol)))
        absParamsAbbrvPlusBias = np.append(paramsAbbrv, 'abs bias')
        barColors = np.append( ["blue"]*len(paramsNames), "red" )
        absColorDiscreteMap = dict(zip(absParamsAbbrvPlusBias.tolist(),
                                    barColors.tolist()))

        paramsAbsIncrsBarChart = \
            createBarChart(absParamsIncrsPlusBias, index=absParamsAbbrvPlusBias,
                           columns=["Parameter<br>contribution"],
                           orientation='v',
                           colorDiscreteMap=absColorDiscreteMap,
                           title="Size of abs parameter contributions <br>" \
                                   + "(= column means of abs sensMatrix*dParams+0.5*curvMatrix*dParams**2)",
                           xlabel="Parameter", ylabel="Contribution to bias removal",
                           width=800, height=500,
                           showLegend=False)




    if False:
        biasLinNlIndivContrbBarFig = \
            createBiasLinNlIndivContrbBarFig(normlzdSensParamsMatrixOrdered, curvParamsMatrixOrdered,
                                             metricsNamesOrdered, paramsNames)

    if False:
        biasVsBiasApproxScatterplot = \
            createBiasVsBiasApproxScatterplot(defaultBiasesApproxNonlin, defaultBiasesCol,
                                              normMetricValsCol,
                                              metricsNames)

    if plot_biasesVsDiagnosticScatterplot:
        diagnosticPrefix = ["U10", "U10"]
        biasVsDiagnosticScatterplot = \
            createBiasVsDiagnosticScatterplot(diagnosticPrefix, defaultBiasesCol,
                                              normMetricValsCol,
                                              defaultNcFilename)

    if False:
        biasSensMatrixScatterFig = \
            createBiasSensMatrixScatterFig(defaultBiasesCol, defaultBiasesApproxElastic,
                                           normMetricValsCol, metricsNames)

    if plot_dpMinMatrixScatterFig:
        dpMinMatrixScatterFig = \
            createDpMinMatrixScatterFig(defaultBiasesCol, normlzdSensMatrixPoly,
                                        normMetricValsCol, metricsNames)

    if False:
        maxSensMetricsFig = \
            createMaxSensMetricsFig(normlzdSensMatrixPoly, metricsNames)

    if plot_biasesVsSensMagScatterplot:
        print("Creating biasesVsSensMagScatterplot . . .")

        #sensCol = np.linalg.norm(normlzdLinplusSensMatrixPoly, axis=1)
        sensCol = np.linalg.norm(normlzdSensMatrixPoly, axis=1)
        # Find the index of the element with the largest magnitude
        maxSensIdx = np.argmax(sensCol)
        signSens = np.sign(normlzdLinplusSensMatrixPoly @ normlzdLinplusSensMatrixPoly[maxSensIdx, :].T)
        xCol = sensCol * signSens
        yCol = -defaultBiasesCol[:, 0] / np.abs(normMetricValsCol[:, 0])
        #yCol = (-defaultBiasesApproxNonlin - defaultBiasesCol)[:, 0] \
        #       / np.abs(normMetricValsCol[:, 0])

        biasesVsSensMagScatterplot = \
            createScatterplot(xCol=xCol, xColLabel='sens',
                              yCol=yCol, yColLabel='bias',
                              #colorCol=normlzdResid, colorColLabel='resid',
                              #colorScale='Spectral',
                              #colorCol=tunedLossChange[:, 0],
                              colorCol=1e3*np.sign(tunedLossChange[:,0])*np.sqrt(np.sqrt(np.abs(tunedLossChange[:,0]))),
                              colorColLabel='sqrtsqrttunedLossChange',
                              colorScale='Spectral_r',
                              plotBgColor='lightgrey',
                              pointLabels=metricsNamesNoprefix, pointLabelsHeader='Metric',
                              plotTitle="""Regional normalized biases vs. signed magnitude of sensitivity.""",
                              xaxisTitle="Signed magnitude of sensitivity of regional metrics to parameter changes",
                              yaxisTitle="Regional biases",
                              showLegend=False, hoverMode="closest",
                              plotWidth=700, plotHeight=500)

        #xCol = np.linalg.norm(normlzdLinplusSensMatrixPoly, axis=1)
        #yCol = normlzdResid

        biasVsSensMagResidScatterplot = \
            createScatterplot(xCol=xCol, xColLabel='sens',
                              #yCol=normlzdResid, yColLabel='normlzdResid',
                              yCol=yCol, yColLabel='bias',
                              colorCol=normlzdResid, colorColLabel='resid',
                              colorScale='Spectral',
                              #colorCol=normlzdResid, colorColLabel='normlzdResid',
                              #colorScale='Rainbow',
                              plotBgColor='lightgrey',
                              pointLabels=metricsNamesNoprefix, pointLabelsHeader='Metric',
                              plotTitle = """Regional normalized biases vs. signed magnitude of sensitivity.""",
                              xaxisTitle = "Signed magnitude of sensitivity of regional metrics to parameter changes",
                              yaxisTitle = "Regional biases",
                              #plotTitle="""Regional normalized residuals vs. signed magnitude of sensitivity.""",
                              #xaxisTitle="Signed magnitude of sensitivity of regional metrics to parameter changes",
                              #yaxisTitle="Regional normalized residuals",
                              showLegend=False, hoverMode="closest",
                              plotWidth=700, plotHeight=500)

        #biasesVsSensMagScatterplot = \
        #    createBiasesVsSensMagScatterplot(normlzdLinplusSensMatrixPoly,
        #                                     defaultBiasesCol, defaultBiasesApproxNonlin,
        #                                     normMetricValsCol, metricsNames)
        #biasesVsSensMagScatterplot = \
        #createBiasesVsSensMagScatterplot(normlzdLinplusSensMatrixPolyMetricsMasked, defaultBiasesColMasked,
        #                                 normMetricValsColMasked, metricsNamesMasked)

    if plot_biasesVsSvdScatterplot:
        print("Creating biasesVsSvdScatterplot . . .")

        # vh = V^T = transpose of right-singular vector matrix, V.
        #normlzdSensMatrixPolyCentered = normlzdSensMatrixPoly - np.mean(normlzdSensMatrixPoly,0)

        #normlzdLinplusSensMatrixPolyCentered = normlzdLinplusSensMatrixPoly - np.mean(normlzdLinplusSensMatrixPoly, 0)

        u, s, vh = \
            np.linalg.svd(normlzdLinplusSensMatrixPoly,
                          full_matrices=False)

        #normlzdLinplusSensMatrixPolyPlusBias = np.hstack((normlzdLinplusSensMatrixPoly, normlzdDefaultBiasesCol))
        #normlzdLinplusSensMatrixPolyPlusBiasCentered = \
        #    normlzdLinplusSensMatrixPolyPlusBias - np.mean(normlzdLinplusSensMatrixPolyPlusBias, 0)
        #u, s, vh = \
        #    np.linalg.svd(normlzdLinplusSensMatrixPolyPlusBiasCentered, full_matrices=False)

        RSquaredSvd = np.square(s) / np.sum(np.square(s))
        if False:
            #if beVerbose:
            print("Variance explained by each SVD component (R**2) = ", RSquaredSvd)
            print("paramsNames = ", paramsNames)
            print("vh[:,0] = ", vh[:, 0])
            print("vh[:,1] = ", vh[:, 1])
            print("u[:,0] = ", u[:, 0])
            print("u[:,1] = ", u[:, 1])
            print("s = ", s)
            print("eigenvalues = lambda = s**2/(n-1)", s * s / (u.shape[0] - 1))
            us = u @ np.diag(s)
            print("u@s[:,0] = ", us[:, 0])

        # Set xCol to first left singular vector
        xCol = u[:, 0] * normlzdDefaultBiasesCol[:,0]
        yCol = u[:, 1] * normlzdDefaultBiasesCol[:,0]
        #yCol = -defaultBiasesCol[:, 0] / np.abs(normMetricValsCol[:, 0])

        biasesVsSvdScatterplot = \
            createScatterplot(xCol=xCol, xColLabel='SV1*bias',
                              yCol=yCol, yColLabel='SV2*bias',
                              colorCol=tunedLossChange[:,0],
                              #colorCol=np.minimum(1, -normlzdDefaultBiasesCol[:, 0] ),
                              colorColLabel='loss change',
                              colorScale='Rainbow',
                              plotBgColor='lightgrey',
                              pointLabels=metricsNamesNoprefix, pointLabelsHeader='Region',
                              plotTitle=(
                                          "Biases (color) as a function of first and second left singular vector values<br>" \
                                          + "Variance explained by each SVD component (R**2) = <br>" \
                                          + np.array2string(RSquaredSvd)),
                              xaxisTitle="First left singular vector values",
                              yaxisTitle="Second left singular vector values",
                              showLegend=False, hoverMode="closest",
                              plotWidth=700, plotHeight=500)

        # Set xCol to first left singular vector
        #xCol = u[whitelistedMetricsMask, 0]
        #yCol = u[whitelistedMetricsMask, 1]

        residVsSvdScatterplot = \
            createScatterplot(xCol=xCol[whitelistedMetricsMask], xColLabel='SV1',
                              yCol=yCol[whitelistedMetricsMask], yColLabel='SV2',
                              colorCol=tunedLossChange[whitelistedMetricsMask, 0],
                              #colorCol=np.minimum(1, -defaultBiasesCol[whitelistedMetricsMask, 0] / np.abs(normMetricValsCol[whitelistedMetricsMask, 0])),
                              colorColLabel='loss change',
                              colorScale='Rainbow',
                              plotBgColor='lightgrey',
                              pointLabels=metricsNamesNoprefix[whitelistedMetricsMask],
                              pointLabelsHeader='Region',
                              plotTitle=(
                                      "Biases (color) as a function of first and second left singular vector values<br>" \
                                      + "Variance explained by each SVD component (R**2) = <br>" \
                                      + np.array2string(RSquaredSvd)),
                              xaxisTitle="First left singular vector values",
                              yaxisTitle="Second left singular vector values",
                              showLegend=False, hoverMode="closest",
                              plotWidth=700, plotHeight=500)

        #residVsSvdScatterplot = \
        #    createScatterplot(xCol=xCol, xColLabel='SV1',
        #                      yCol=yCol, yColLabel='SV2',
        #                      colorCol=np.minimum(1, normlzdResid),
        #                      colorColLabel='normlzdResid',
        #                      plotBgColor='lightgrey',
        #                      pointLabels=metricsNamesNoprefix, pointLabelsHeader='Region',
        #                      plotTitle="""Residuals (color) as a function of first and second left singular vector values""",
        #                      xaxisTitle="First left singular vector values",
        #                      yaxisTitle="Second left singular vector values",
        #                      showLegend=False, hoverMode="closest",
        #                      plotWidth=700, plotHeight=500)

    if False:
        biasesVsSensArrowFig = \
            createBiasesVsSensArrowFig(normlzdWeightedSensMatrixPoly, defaultBiasesCol,
                                       defaultBiasesApproxNonlin,
                                       normMetricValsCol, metricsNames)

    if False:
        normlzdSensMatrixColsFig = \
            createNormlzdSensMatrixColsFig(defaultBiasesCol, normlzdSensMatrixPoly,
                                           normMetricValsCol, metricsNames, paramsNames)

    if False:
        normlzdSensMatrixRowsFig = \
            createNormlzdSensMatrixRowsFig(normlzdSensMatrixPoly,
                                           metricsNames, paramsNames)

    if plot_sensMatrixAndBiasVecFig:
        print("Creating sensMatrixAndBiasVecFig . . .")
        # Create figure that shows the sensitivity matrix and bias column, both color coded.
        matrixDictKeyString = "normlzdLinplusSensMatrixPoly"
        matrixDict = {matrixDictKeyString: normlzdLinplusSensMatrixPoly}
        sensMatrixAndBiasVecFig = createMatrixPlusColFig(matrix=matrixDict[matrixDictKeyString],
                                                         matIndexLabel=metricsNames,
                                                         matColLabel=paramsAbbrv,
                                                         colVector=-np.around(
                                                             defaultBiasesCol / np.abs(normMetricValsCol), decimals=2),
                                                         colVectIndexLabel=metricsNames,
                                                         colVectColLabel=['-Normalized Biases'],
                                                         plotHeight=1400, plotWidth=1000,
                                                         #plotHeight=2800 displays all 20x20 regions
                                                         printCellText=False,
                                                         plotTitle='Color-coded normalized sensitivity matrix, ' + matrixDictKeyString,
                                                         reversedYAxis='reversed',
                                                         eqnAdd=True)

    # Needed for several plots:
    XT_dot_X_Linplus = normlzdLinplusSensMatrixPoly.T @ normlzdLinplusSensMatrixPoly
    #XT_dot_X_Linplus = normlzdSensMatrixPoly.T @ normlzdSensMatrixPoly
    #XT_dot_X_Linplus = normlzdWeightedLinplusSensMatrixPoly.T @ normlzdWeightedLinplusSensMatrixPoly

    if False:
        # Create figure that plots color-coded parameter correlation matrix plus parameter-bias correlation column.
        (XT_dot_X_Linplus_corr, stdMatrixInv) = covMatrix2corrMatrix(XT_dot_X_Linplus, returnStd=True)
        normlzdStdDefaultBiasesCol = stdMatrixInv @ normlzdLinplusSensMatrixPoly.T @ normlzdDefaultBiasesCol
        #normlzdStdDefaultBiasesCol = stdMatrixInv @ normlzdSensMatrixPoly.T @ normlzdDefaultBiasesCol
        #normlzdStdDefaultBiasesCol = stdMatrixInv @ normlzdWeightedLinplusSensMatrixPoly.T @ normlzdWeightedDefaultBiasesCol
        paramsCorrArrayBiasFig = createMatrixPlusColFig(matrix=XT_dot_X_Linplus_corr,
                                                        matIndexLabel=paramsNames,
                                                        matColLabel=paramsNames,
                                                        colVector=-np.around(normlzdStdDefaultBiasesCol, decimals=2),
                                                        colVectIndexLabel=paramsNames,
                                                        colVectColLabel=['Projection onto -biases'],
                                                        plotHeight=700, plotWidth=1000,
                                                        printCellText=True,
                                                        plotTitle='Cosines of angles between columns of sensitivity matrix',
                                                        reversedYAxis='reversed',
                                                        eqnAdd=False)

    if plot_projectionMatrixFigs:
        # Create figure that plots color-coded projection matrix plus bias column.
        XT_dot_X_Linplus_inv = np.linalg.inv(XT_dot_X_Linplus)
        fullProjectionMatrix = normlzdLinplusSensMatrixPoly @ XT_dot_X_Linplus_inv @ normlzdLinplusSensMatrixPoly.T
        allLeverages = np.diag(fullProjectionMatrix)
        if beVerbose:
            np.set_printoptions(linewidth=200)
            print("All leverages = ", allLeverages)
            print("Sum of leverages = ", fullProjectionMatrix.trace())
            print("Number of parameters =", len(paramsNames))
            print("A large leverage is > 3p/n, which = ",
                  3 * len(paramsNames) / len(metricsNames))
        #projectionMatrixFig = createMatrixPlusColFig( matrix = projectionMatrix,
        #                 matIndexLabel = metricsNames,
        #                 matColLabel = metricsNames,
        #                 colVector = -np.around(normlzdDefaultBiasesCol, decimals=2),
        #                 colVectIndexLabel = metricsNames,
        #                 colVectColLabel = ['-Normalized Biases'],
        #                 plotHeight=700, plotWidth=1000,
        #                 plotTitle='Projection matrix',
        #                 reversedYAxis = 'reversed',
        #                 eqnAdd = False)
        linplusProjectionMatrixMasked = normlzdLinplusSensMatrixPolyMetricsMasked \
                                        @ XT_dot_X_Linplus_inv \
                                        @ normlzdLinplusSensMatrixPolyMetricsMasked.T
        if beVerbose:
            print("Masked leverages = ", np.diag(linplusProjectionMatrixMasked))
            #print("projectionMatrix rows=", np.linalg.norm( projectionMatrix, axis=1))

        matrixDictKeyString = "linplusProjectionMatrixMasked"
        matrixDict = {matrixDictKeyString: linplusProjectionMatrixMasked}
        projectionMatrixFig = createMatrixPlusColFig(matrix=matrixDict[matrixDictKeyString],
                                                     matIndexLabel=metricsNamesMasked,
                                                     matColLabel=metricsNamesMasked,
                                                     colVector=-np.around(normlzdDefaultBiasesColMasked, decimals=2),
                                                     colVectIndexLabel=metricsNamesMasked,
                                                     colVectColLabel=['-Normalized Biases'],
                                                     plotHeight=700, plotWidth=1000,
                                                     printCellText=True,
                                                     plotTitle='Excerpt of projection matrix, ' + matrixDictKeyString,
                                                     reversedYAxis='reversed',
                                                     eqnAdd=False)

        xCol = allLeverages
        yCol = np.abs(-defaultBiasesCol[:, 0]) / np.abs(normMetricValsCol[:, 0])

        biasesVsLeveragesScatterplot = \
            createScatterplot(xCol=xCol, xColLabel='Lev',
                              yCol=yCol, yColLabel='bias',
                              colorCol=yCol, colorColLabel='bias',
                              colorScale='Rainbow',
                              plotBgColor='lightgrey',
                              pointLabels=metricsNames, pointLabelsHeader='Metric',
                              plotTitle="""Regional biases vs. leverages.""",
                              xaxisTitle="Leverages",
                              yaxisTitle="Regional biases",
                              showLegend=False, hoverMode="closest",
                              plotWidth=700, plotHeight=500)

    if plot_paramsCorrArrayFig:
        matrixDictKeyString = "normlzdSensMatrixPoly"
        matrixDict = {matrixDictKeyString: normlzdSensMatrixPoly}
        #matrixDictKeyString = "normlzdLinplusSensMatrixPoly"
        #matrixDict={matrixDictKeyString:normlzdLinplusSensMatrixPoly}
        paramsCorrArrayFig = \
            createParamsCorrArrayFig(matrix=matrixDict[matrixDictKeyString],
                                     biasesCol=normlzdDefaultBiasesCol,
                                     paramsNames=paramsNames,
                                     plotTitle='cos(angle) among parameters<br>\
                                               (i.e., X^T*X using debiased columns of sens matrix)<br>' \
                                               + matrixDictKeyString + '<br>' \
                                     )

    if plot_dpMin2PtFig:
        print("Creating dpMin2PtFig . . .")
        #    dpMin2PtFig = \
        #    createDpMin2PtFig( normlzdLinplusSensMatrixPoly, defaultBiasesCol,
        #                          normMetricValsCol, metricsNames )
        dpMin2PtFig = \
            createDpMin2PtFig(normlzdLinplusSensMatrixPolyMetricsMasked, defaultBiasesColMasked,
                              normMetricValsColMasked, metricsNamesMasked)

    # The following plot is redundant:
    #biasTotContrbBarFig = \
    #      createBarChart( minusNonlinMatrixDparamsOrdered, index=metricsNamesOrdered, columns=paramsNames,
    #       #               barBase = np.zeros(numMetrics),
    #                      #barBase = -defaultBiasesCol[metricsSensOrder]/np.abs(normMetricValsCol[metricsSensOrder]) @ np.ones((1,len(paramsNames))),
    #                      orientation = 'v',
    #                      title="""Linear + nonlinear contributions to removal of biases <br>
    #                               (dMetrics=sensMatrix*dParams)""",
    #                      xlabel="Regional metric", ylabel="Contribution to bias removal",
    #                      width=800, height=500 )

    if plot_PcaBiplot:
        print("Creating PcaBiplotFig . . .")

        paramsNamesAbbr = np.char.replace(paramsNames, 'clubb_', '')
        paramsNamesAbbr = np.char.replace(paramsNamesAbbr, 'c_invrs_tau_', '')
        paramsNamesAbbr = np.char.replace(paramsNamesAbbr, 'wpxp_n2', 'n2')
        paramsNamesAbbr = np.char.replace(paramsNamesAbbr, 'threshold', 'thresh')

        # Create scatterplot to look at outliers
        PcaBiplotFig = \
            createPcaBiplot(normlzdSensMatrixPoly, normlzdDefaultBiasesCol,
                            paramsNamesAbbr,
                            whitelistedMetricsMask,
                            xColLabel='SV1', yColLabel='SV2',
                            colorCol=np.minimum(1, -normlzdDefaultBiasesCol[:, 0]),
                            colorColLabel='bias',
                            pointLabels=metricsNames, pointLabelsHeader='Metric',
                            plotTitle=(
                                        "Biases (color) as a function of first and second left singular vector values<br>" \
                                        + "Variance explained by each SVD component (R**2) = <br>"),
                            xaxisTitle="First left singular vector values",
                            yaxisTitle="Second left singular vector values",
                            showLegend=False, hoverMode="closest",
                            plotWidth=700, plotHeight=500
                            )

    if plot_PcSensMap:

        print("Creating PcSensMap . . .")

        PcMapPanelBias = \
            createMapPanel(fieldToPlotCol=normlzdDefaultBiasesCol,
                           plotWidth=500,
                           plotTitle="normlzdDefaultBiasesCol",
                           boxSize=20)

        # Use same color range in residual plot as in bias plot
        minFieldBias = np.min(normlzdDefaultBiasesCol)
        maxFieldBias = np.max(normlzdDefaultBiasesCol)

        PcMapPanelResid = \
            createMapPanel(fieldToPlotCol=-normlzdResid,
                           plotWidth=500,
                           plotTitle="-normlzdResid",
                           boxSize=20,
                           minField=minFieldBias,
                           maxField=maxFieldBias)

        BiasParamsDashboardChildren = [html.Div(children=[
            dcc.Graph(id="PcMapPanelBias", figure=PcMapPanelBias,
                      style={'display': 'inline-block'}, config=downloadConfig),
            dcc.Graph(id="PcMapPanelResid", figure=PcMapPanelResid,
                      style={'display': 'inline-block'}, config=downloadConfig)
        ])]

        PcMapPanelDefaultLoss = \
            createMapPanel(fieldToPlotCol=1e3*np.sqrt(defaultLossCol),
                           plotWidth=500,
                           plotTitle="Sqrt Default Loss (x 1e3)",
                           boxSize=20)

        PcMapPanelTunedLossChange = \
            createMapPanel(fieldToPlotCol=1e3*np.sign(tunedLossChange)*np.sqrt(np.abs(tunedLossChange)),
                           plotWidth=500,
                           plotTitle="Sqrt Tuned Loss Change (x 1e3)",
                           boxSize=20)

        BiasParamsDashboardChildren.append(html.Div(children=[
                dcc.Graph(figure=PcMapPanelDefaultLoss, style={'display': 'inline-block'},
                          config=downloadConfig),
                dcc.Graph(figure=PcMapPanelTunedLossChange, style={'display': 'inline-block'},
                          config=downloadConfig)
            ]))


        paramsIdx = 0
        while paramsIdx < len(paramsNames):
            leftFig = \
                createMapPanel(fieldToPlotCol=normlzdLinplusSensMatrixPoly[:, paramsIdx],
                               plotWidth=500,
                               plotTitle=f"normlzdSensMatrixPoly[:,{paramsNames[paramsIdx]}]",
                               boxSize=20)
            if paramsIdx + 1 < len(paramsNames):
                rightFig = \
                    createMapPanel(fieldToPlotCol=normlzdLinplusSensMatrixPoly[:, paramsIdx + 1],
                                   plotWidth=500,
                                   plotTitle=f"normlzdSensMatrixPoly[:,{paramsNames[paramsIdx + 1]}]",
                                   boxSize=20)
            else:
                rightFig = go.Figure()
                plotWidth = 500
                plotHeight = np.rint(plotWidth * (490 / 700))
                rightFig.update_layout(width=plotWidth, height=plotHeight)

            BiasParamsDashboardChildren.append(html.Div(children=[
                dcc.Graph(figure=leftFig, style={'display': 'inline-block'},
                          config=downloadConfig),
                dcc.Graph(figure=rightFig, style={'display': 'inline-block'},
                          config=downloadConfig)
            ]))

            paramsIdx += 2

        #print("1e6*defaultLossCol = ", 1e6*np.sort(defaultLossCol[:,0]))

        u, s, vh = \
            np.linalg.svd(normlzdLinplusSensMatrixPoly, full_matrices=False)

        #u, s, vh = \
        #    np.linalg.svd(normlzdCurvMatrix, full_matrices=False)

        PcMapPanelU0 = \
            createMapPanel(fieldToPlotCol=u[:, 0],
                           plotWidth=500,
                           plotTitle="SVD 1",
                           boxSize=20)

        PcMapPanelU1 = \
            createMapPanel(fieldToPlotCol=u[:, 1],
                           plotWidth=500,
                           plotTitle="SVD 2",
                           boxSize=20)

        PcMapPanelU0bias = \
            createMapPanel(fieldToPlotCol=u[:, 0]*normlzdDefaultBiasesCol[:,0],
                           plotWidth=500,
                           plotTitle="SVD 1 * bias",
                           boxSize=20)

        PcMapPanelU1bias = \
            createMapPanel(fieldToPlotCol=u[:, 1]*normlzdDefaultBiasesCol[:,0],
                           plotWidth=500,
                           plotTitle="SVD 2 * bias",
                           boxSize=20)

        U0U3DashboardChildren = [
            html.Div(children=
                     [dcc.Graph(id="PcMapPanelU0", figure=PcMapPanelU0,
                                style={'display': 'inline-block'},
                                config=downloadConfig),
                      dcc.Graph(id="PcMapPanelU1", figure=PcMapPanelU1,
                                style={'display': 'inline-block'},
                                config=downloadConfig)
                      ]
                     ),
            html.Div(children=
                     [dcc.Graph(id="PcMapPanelU0bias", figure=PcMapPanelU0bias,
                                style={'display': 'inline-block'},
                                config=downloadConfig),
                      dcc.Graph(id="PcMapPanelU1bias", figure=PcMapPanelU1bias,
                                style={'display': 'inline-block'},
                                config=downloadConfig)
                      ]
                     )
        ]

        # Now combine the matrix and column sub-figures into one figure
        #PcMapFig = make_subplots(
        #    rows=1, cols=2,
        #    column_widths=[0.5, 0.5],
        #    horizontal_spacing=0.1,
        #    specs=[[{"type": "scattergeo"}, {"type": "scattergeo"}]]
        #)

        # subplot layouts
        #https://community.plotly.com/t/how-to-pass-figure-layout-info-into-subplot/72608/2
        #PcMapFig.add_trace(PcMapPanelU0.data[0], row=1, col=1)
        #PcMapFig.add_trace(PcMapPanelU1.data[0], row=1, col=2)

    if plot_vhMatrixFig:

        print("Creating SVD vh matrix figure . . .")

        u, s, vh = \
            np.linalg.svd(normlzdLinplusSensMatrixPoly, full_matrices=False)

        uCurv, sCurv, vhCurv = \
            np.linalg.svd(normlzdCurvMatrix, full_matrices=False)

        #normlzdLinplusSensMatrixPolyPlusBias = \
        #    np.hstack((normlzdLinplusSensMatrixPoly, normlzdDefaultBiasesCol))

        #u, s, vh = \
        #    np.linalg.svd(normlzdLinplusSensMatrixPolyPlusBias, full_matrices=False)

        svdLabel = (np.arange(vh.shape[0])+1).astype(str)
        svdLabel = np.char.add("SVD ", svdLabel)
        svdLabel = np.char.add(svdLabel, " ")

        vhMatrixFig = \
        createColoredMatrixFig(
            matrix=vh,
            matrixRowLabel=svdLabel,
            matrixColLabel=paramsAbbrv,
            plotTitle="SVD v^T Matrix",
            plotWidth=500,
            plotHeight=500,
            printCellText=True
        )


    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

    sensMatrixDashboard = dash.Dash(__name__, external_stylesheets=external_stylesheets)

    dashboardChildren = [
        html.H1(children='Sensitivity matrix diagnostics'),

        html.Div(children=''' ''')]

    if plot_PcSensMap:
        dashboardChildren.extend(BiasParamsDashboardChildren)
        dashboardChildren.extend(U0U3DashboardChildren)
        #dashboardChildren.append(dcc.Graph(id='PcMapFig', figure=PcMapFig))
    if plot_vhMatrixFig:
        dashboardChildren.append(dcc.Graph(id='vhMatrixFig', figure=vhMatrixFig,
                                           config=downloadConfig))
    if plot_PcaBiplot:
        dashboardChildren.append(dcc.Graph(id='PcaBiplotFig', figure=PcaBiplotFig,
                                           config=downloadConfig))
    if plot_paramsErrorBarsFig:
        dashboardChildren.append(dcc.Graph(id='paramsErrorBarsFig', figure=paramsErrorBarsFig,
                                           config=downloadConfig))
    if plot_biasesOrderedArrowFig:
        dashboardChildren.append(dcc.Graph(id='biasesOrderedArrowFig', figure=biasesOrderedArrowFig,
                                           config=downloadConfig))
    if plot_metricsBarChart:
        dashboardChildren.append(dcc.Graph(id='metricsBarChart', figure=metricsBarChart,
                                           config=downloadConfig))
    if plot_paramsTotContrbBarChart:
        dashboardChildren.append(dcc.Graph(id='paramsTotContrbBarChart', figure=paramsTotContrbBarChart,
                                           config=downloadConfig))
    if plot_paramsIncrsBarChart:
        dashboardChildren.append(dcc.Graph(id='paramsIncrsBarChart', figure=paramsIncrsBarChart,
                                           config=downloadConfig))
    if plot_paramsAbsIncrsBarChart:
        dashboardChildren.append(dcc.Graph(id='paramsAbsIncrsBarChart', figure=paramsAbsIncrsBarChart,
                                           config=downloadConfig))
    if False:
        dashboardChildren.append(dcc.Graph(id='linplusSensMatrixBarFig', figure=linplusSensMatrixBarFig,
                                           config=downloadConfig))
    if False:
        dashboardChildren.append(dcc.Graph(id='biasLinNlIndivContrbBarFig', figure=biasLinNlIndivContrbBarFig,
                                           config=downloadConfig))
    if plot_dpMin2PtFig:
        dashboardChildren.append(dcc.Graph(id='dpMin2PtFig', figure=dpMin2PtFig,
                                           config=downloadConfig))
    if False:
        dashboardChildren.append(dcc.Graph(id='biasVsBiasApproxScatterplot', figure=biasVsBiasApproxScatterplot,
                                           config=downloadConfig))
    #config= { 'toImageButtonOptions': { 'scale': 6 } }
    if plot_biasesVsDiagnosticScatterplot:
        dashboardChildren.append(dcc.Graph(id='biasVsDiagnosticScatterplot', figure=biasVsDiagnosticScatterplot,
                                           config=downloadConfig))
    if plot_sensMatrixAndBiasVecFig:
        dashboardChildren.append(dcc.Graph(id='sensMatrixAndBiasVecFig', figure=sensMatrixAndBiasVecFig,
                                           config=downloadConfig))
    if False:
        dashboardChildren.append(dcc.Graph(id='paramsCorrArrayBiasFig', figure=paramsCorrArrayBiasFig,
                                           config=downloadConfig))
    if plot_paramsCorrArrayFig:
        dashboardChildren.append(dcc.Graph(id='paramsCorrArrayFig', figure=paramsCorrArrayFig,
                                           config=downloadConfig))
    if False:
        dashboardChildren.append(dcc.Graph(id='metricsCorrArrayFig', figure=metricsCorrArrayFig,
                                           config=downloadConfig))
    if plot_projectionMatrixFigs:
        dashboardChildren.append(dcc.Graph(id='projectionMatrixFig', figure=projectionMatrixFig,
                                           config=downloadConfig))
        dashboardChildren.append(dcc.Graph(id='biasesVsLeveragesScatterplot', figure=biasesVsLeveragesScatterplot,
                                           config=downloadConfig))
    if plot_biasesVsSensMagScatterplot:
        dashboardChildren.append(dcc.Graph(id='biasesVsSensMagScatterplot',
                                           figure=biasesVsSensMagScatterplot,
                                           config=downloadConfig))
        dashboardChildren.append(dcc.Graph(id='biasVsSensMagResidScatterplot',
                                           figure=biasVsSensMagResidScatterplot,
                                           config=downloadConfig))
    if plot_biasesVsSvdScatterplot:
        dashboardChildren.append(dcc.Graph(id='biasesVsSvdScatterplot', figure=biasesVsSvdScatterplot,
                                           config=downloadConfig))
        dashboardChildren.append(dcc.Graph(id='residVsSvdScatterplot', figure=residVsSvdScatterplot,
                                           config=downloadConfig))
    if plot_threeDotFig:
        dashboardChildren.append(dcc.Graph(id='threeDotFig', figure=threeDotFig,
                                           config=downloadConfig))
    if False:
        dashboardChildren.append(dcc.Graph(id='biasSensScatterFig', figure=biasSensMatrixScatterFig,
                                           config=downloadConfig))
    if plot_dpMinMatrixScatterFig:
        dashboardChildren.append(dcc.Graph(id='dpMinMatrixScatterFig', figure=dpMinMatrixScatterFig,
                                           config=downloadConfig))
    if False:
        dashboardChildren.append(dcc.Graph(id='maxSensMetricsFig', figure=maxSensMetricsFig,
                                           config=downloadConfig))
    if False:
        dashboardChildren.append(dcc.Graph(id='normlzdSensMatrixColsFig', figure=normlzdSensMatrixColsFig,
                                           config=downloadConfig))
    if False:
        dashboardChildren.append(dcc.Graph(id='normlzdSensMatrixRowsFig', figure=normlzdSensMatrixRowsFig,
                                           config=downloadConfig))
    if False:
        dashboardChildren.append(dcc.Graph(id='biasesVsSensArrowFig', figure=biasesVsSensArrowFig,
                                           config=downloadConfig))

    print("At end of createFigs . . .")
    ##redundant dcc.Graph( id='biasTotContrbBarFig', figure=biasTotContrbBarFig ),

    sensMatrixDashboard.layout = html.Div(children=dashboardChildren)

    sensMatrixDashboard.run_server(debug=True, use_reloader=False)

    return


def createMapPanel(fieldToPlotCol,
                   plotWidth,
                   plotTitle,
                   boxSize,
                   minField=None, maxField=None):

    regionalMapPanel = go.Figure(go.Scattergeo())

    # Set boundaries for colored regional boxes
    regionalMapPanel.update_layout(xaxis=dict(range=[0, 360]),
                                   yaxis=dict(range=[-90, 90]))
    regionalMapPanel.update_xaxes(showgrid=False,
                                  tickmode="linear",
                                  dtick=60, tick0=0)
    regionalMapPanel.update_yaxes(showgrid=False,
                                  zeroline=False,
                                  tickmode="linear",
                                  dtick=30, tick0=-90)

    # Calculate number of regions in the east-west (X) and north-south (Y) directions
    numXBoxes = np.rint(360 / boxSize).astype(int)  # 18
    numYBoxes = np.rint(180 / boxSize).astype(int)  # 9

    # Draw longitudinal regional box number along top of plot
    xBoxNums = np.linspace(start=1, stop=numXBoxes, num=numXBoxes).astype(int)
    for xIdx, xBoxNum in np.ndenumerate(xBoxNums):
        xPos = int(0.5 * boxSize + xIdx[0] * boxSize)
        regionalMapPanel.add_annotation(x=xPos, yref="paper", y=1.1,
                                        text=xBoxNum.astype(str), showarrow=False)
    # Draw latitudinal box number along right-hand side of plot
    yBoxNums = np.linspace(start=1, stop=numYBoxes, num=numYBoxes).astype(int)
    for yIdx, yBoxNum in np.ndenumerate(yBoxNums):
        yPos = int(90 - 0.5 * boxSize - yIdx[0] * boxSize)
        regionalMapPanel.add_annotation(xref="paper", x=1.04,
                                        y=yPos,
                                        text=yBoxNum.astype(str), showarrow=False)

    plotHeight = np.rint(plotWidth * (490 / 700))
    regionalMapPanel.update_layout(width=plotWidth, height=plotHeight)
    #regionalMapPanel.update_layout(width=700, height=450)
    regionalMapPanel.update_layout(title=plotTitle, title_y=0.9, title_x=0.5)

    # Shift coastal boundaries to correct location on map
    regionalMapPanel.update_geos(lataxis_range=[-90, 90],
                                 lonaxis_range=[0, 360])

    fieldToPlotMatrix = fieldToPlotCol.reshape(numYBoxes, numXBoxes)
    # Shift from 0,360 to -180,180 degrees longitude
    #fieldToPlotMatrix = np.roll(fieldToPlotMatrix, -9, axis=1)

    # Set color scaling for colors of regional boxes
    latRange = range(90, -90, -boxSize)
    lonRange = range(0, 360, boxSize)
    normlzdColorMatrix = np.zeros_like(fieldToPlotMatrix)
    if minField == None or maxField == None:
        minField = np.min(fieldToPlotMatrix)
        maxField = np.max(fieldToPlotMatrix)
    rangeField = np.maximum( np.abs(maxField), np.abs(minField) )

    # If the data include 0, use a diverging colorscale.
    # Remap 0 to 0.5 and normalize the data so that
    #    the data lie within the range [0,1].
    #if ( False ):
    if ( maxField > 0) & (minField < 0):
        colorScale = 'RdBu_r'
        for latIdx, lat in np.ndenumerate(latRange):
            for lonIdx, lon in np.ndenumerate(lonRange):
                normlzdColorMatrix[latIdx, lonIdx] = \
                    0.5 * fieldToPlotMatrix[latIdx, lonIdx] / rangeField + 0.5
                #if fieldToPlotMatrix[latIdx][lonIdx] < 0:
                #    normlzdColorMatrix[latIdx, lonIdx] = \
                #        0.5 * (fieldToPlotMatrix[latIdx, lonIdx]-minField)/np.abs(minField)
                #else:
                #    normlzdColorMatrix[latIdx, lonIdx] = \
                #        0.5 * fieldToPlotMatrix[latIdx, lonIdx] / maxField + 0.5
    else:  # Don't use diverging colorscale
        colorScale = 'Bluered'
        normlzdColorMatrix = (fieldToPlotMatrix - minField) / \
                      (maxField - minField)

    # Draw a colored rectangle in each region in layer underneath
    #print("After setting colorScale, it = ", colorScale)
    #colorIdx = 0
    colorList = []
    for latIdx, lat in np.ndenumerate(latRange):
        for lonIdx, lon in np.ndenumerate(lonRange):
            # 'bluered'
            colorString = pc.sample_colorscale(colorscale=colorScale,
                                         samplepoints=normlzdColorMatrix[latIdx, lonIdx],
                                         )[0]
                                         #low = 0, high = 1)[0]

            colorList.append([normlzdColorMatrix[latIdx, lonIdx].item(),
                              colorString])

            regionalMapPanel.add_shape(
                type="rect",
                xref="x",
                yref="y",
                x0=lon,
                y0=lat,
                x1=lon + boxSize,
                y1=lat - boxSize,
                line=dict(color="black", width=1),
                #fillcolor=colorList[colorIdx],
                fillcolor=colorString,
                opacity=1.0,
                layer="below"
            )

            #colorIdx += 1

    # Draw map of land boundaries in layer on top
    regionalMapPanel.update_geos(showcoastlines=True,
                                 coastlinecolor='black',
                                 coastlinewidth=1,
                                 showlakes=False,
                                 showland=False,
                                 showocean=False,
                                 bgcolor='rgba(0,0,0,0)',
                                 lonaxis_showgrid=False,
                                 lataxis_showgrid=False
                                 )
    # Add colorbar by creating an invisible, fake scatterplot
    if (colorScale == 'RdBu_r'):
        tickVals = [np.min(normlzdColorMatrix),
                    0.5*np.min(normlzdColorMatrix)+0.5*np.max(normlzdColorMatrix),
                    np.max(normlzdColorMatrix)]
        tickText = [f"{-rangeField:.2f}",
                    '0.0',
                    f"{rangeField:.2f}"]
        #if (maxField > np.abs(minField)):
        #    #tickVals = [0.5 * minField / rangeField + 0.5,
        #    #            0.25 * minField / rangeField + 0.75,
        #    #            1.0]
        #    tickVals = [0.5 * np.min(fieldToPlotMatrix) / rangeField + 0.5,
        #                0.25 * np.min(fieldToPlotMatrix) / rangeField + 0.75,
        #                1.0]
        #    tickText = [f"{-maxField:.2f}",
        #                '0.0',
        #                f"{maxField:.2f}"]
        #else:
        #    #tickVals = [0.0,
        #    #            0.25 * maxField / rangeField + 0.25,
        #    #            0.5 * maxField / rangeField + 0.5]
        #    tickVals = [0.0,
        #                0.25 * np.max(fieldToPlotMatrix) / rangeField + 0.25,
        #                0.5 * np.max(fieldToPlotMatrix) / rangeField + 0.5]
        #    tickText = [f"{minField:.2f}",
        #                '0.0',
        #                f"{-minField:.2f}"]
    else:
        tickVals = [0.0, 1.0]
        tickText = ['0', f"{rangeField:.2f}"]
    colorbar_trace = go.Scatter(x=[None],
                                y=[None],
                                mode='markers',
                                marker=dict(
                                    #color=fieldToPlotCol,
                                    color=normlzdColorMatrix.reshape((numYBoxes*numXBoxes,)),
                                    colorscale=colorScale,
                                    #colorscale=colorList,
                                    showscale=True,
                                    colorbar=dict(thickness=15,
                                                  tickvals = tickVals,
                                                  ticktext = tickText
                                                  )
                                                  #outlinewidth=0)
                                )
                                )

    regionalMapPanel.add_trace(colorbar_trace)

    return regionalMapPanel


def covMatrix2corrMatrix(covMatrix, returnStd=False):
    # https://gist.github.com/wiso/ce2a9919ded228838703c1c7c7dad13b

    import numpy as np

    stdVector = np.sqrt(np.diag(covMatrix))
    stdMatrixInv = np.diag(1.0 / stdVector)
    corrMatrix = stdMatrixInv @ covMatrix @ stdMatrixInv
    if returnStd:
        return (corrMatrix, stdMatrixInv)
    else:
        return corrMatrix


def createMatrixPlusColFig(matrix, matIndexLabel, matColLabel,
                           colVector, colVectIndexLabel, colVectColLabel,
                           plotHeight, plotWidth,
                           printCellText,
                           plotTitle, reversedYAxis=None,
                           eqnAdd=False):
    '''Creates a figure that displays a color-coded matrix and an accompanying column vector.'''

    import numpy as np
    import pandas as pd
    import plotly.express as px

    # First create a sub-figure that displays color-coded matrix
    matSubfig = \
        createColoredMatrixFig(
            matrix=matrix,
            matrixRowLabel=matIndexLabel,
            matrixColLabel=matColLabel,
            plotTitle=plotTitle,
            plotWidth=None,
            plotHeight=None,
            printCellText=printCellText
        )

    roundedMatrix = np.around(matrix, decimals=2)
    matMaxAbs = np.max(np.abs(roundedMatrix))

    # Now create a sub-figure showing a color-coded column vector
    df_biasArray = pd.DataFrame(colVector,
                                index=colVectIndexLabel,
                                columns=colVectColLabel)
    colVectSubfig = px.imshow(
        df_biasArray.to_numpy(),
        x=df_biasArray.columns.tolist(),
        y=df_biasArray.index.tolist(),
        text_auto=printCellText
    )
    colVectSubfig.update_layout(
        title_text='',
        title_x=0.5,
        #width=10,
        #height=1400,
        xaxis_showgrid=False,
        yaxis_showgrid=False,
        xaxis_zeroline=False,
        yaxis_zeroline=False,
    )

    # Now combine the matrix and column sub-figures into one figure
    matrixPlusColFig = make_subplots(
        rows=1, cols=2,
        column_widths=[0.9, 0.1],
        horizontal_spacing=0.3,
    )
    matrixPlusColFig.add_trace(matSubfig.data[0], row=1, col=1)
    matrixPlusColFig.add_trace(colVectSubfig.data[0], row=1, col=2)
    matrixPlusColFig.update_layout(
        title_text=plotTitle,
        height=plotHeight,
        width=plotWidth,
        template='plotly_white')
    matrixPlusColFig.update_layout(coloraxis=dict(colorscale='RdBu_r', cmin=-matMaxAbs, cmax=matMaxAbs),
                                   showlegend=False)
    matrixPlusColFig.update_yaxes(autorange=reversedYAxis, row=1, col=2)
    matrixPlusColFig.update_yaxes(autorange=reversedYAxis, row=1, col=1)

    if eqnAdd == True:

        matrixPlusColFig.add_annotation(dict(font=dict(color="black", size=36),
                                             # x=x_loc,
                                             x=0.82,
                                             y=0.5,
                                             showarrow=False,
                                             text='<b>=</b>',
                                             textangle=0,
                                             xref="paper",
                                             yref="paper"
                                             ))

        #matrixPlusColFig.add_vline(x=-0.78,col=2)
        left = 0.67
        right = 0.76
        top = 0.67
        bottom = 0.35
        delta = 0.015
        matrixPlusColFig.add_shape(type='line', x0=right, y0=bottom, x1=right, y1=top,
                                   xref="paper", yref="paper")
        matrixPlusColFig.add_shape(type='line', x0=right, y0=top, x1=right - delta, y1=top,
                                   xref="paper", yref="paper")
        matrixPlusColFig.add_shape(type='line', x0=left, y0=top, x1=left + delta, y1=top,
                                   xref="paper", yref="paper")

        matrixPlusColFig.add_shape(type='line', x0=left, y0=bottom, x1=left, y1=top,
                                   xref="paper", yref="paper")
        matrixPlusColFig.add_shape(type='line', x0=right, y0=bottom, x1=right - delta, y1=bottom,
                                   xref="paper", yref="paper")
        matrixPlusColFig.add_shape(type='line', x0=left, y0=bottom, x1=left + delta, y1=bottom,
                                   xref="paper", yref="paper")

        # print column vector of parameter names from y=0.38 to 0.62
        labelSpacing = 0.3 / len(matColLabel)
        for counter, colLabel in enumerate(reversed(matColLabel)):
            matrixPlusColFig.add_annotation(dict(font=dict(color="black", size=12),
                                                 x=0.75,
                                                 y=0.38 + labelSpacing * counter,
                                                 showarrow=False,
                                                 text=colLabel,
                                                 textangle=0,
                                                 xref="paper",
                                                 yref="paper"
                                                 ))

    return (matrixPlusColFig)


def createColoredMatrixFig(
        matrix,
        matrixRowLabel, matrixColLabel,
        plotTitle, plotWidth, plotHeight,
        printCellText
                            ):


    # Create a figure that displays a color-coded matrix
    roundedMatrix = np.around(matrix, decimals=2)
    df_matrix = pd.DataFrame(roundedMatrix,
                             index=matrixRowLabel,
                             columns=matrixColLabel)
    # matMaxAbs = np.max(np.abs(roundedNormlzdSensMatrix))
    coloredMatrixFig = px.imshow(
        df_matrix.to_numpy(),
        x=df_matrix.columns.tolist(),
        y=df_matrix.index.tolist(),
        text_auto=printCellText,
        color_continuous_scale=px.colors.diverging.balance
    )
    coloredMatrixFig.update_xaxes(side="bottom")
    coloredMatrixFig.update_layout(
        title_text=plotTitle,
        title_x=0.5,
        width=plotWidth,
        height=plotHeight,
        xaxis_showgrid=False,
        yaxis_showgrid=False,
        xaxis_zeroline=False,
        yaxis_zeroline=False,
        # template='plotly_white',
        coloraxis_showscale=False
    )

    return coloredMatrixFig


def createScatterplot(xCol, xColLabel,
                      yCol, yColLabel,
                      colorCol, colorColLabel,
                      colorScale,
                      plotBgColor,
                      pointLabels, pointLabelsHeader,
                      plotTitle,
                      xaxisTitle,
                      yaxisTitle,
                      showLegend, hoverMode,
                      plotWidth, plotHeight):
    # Helper function that plots a column of length numMetrics (yCol) vs. xCol
    #df = pd.DataFrame({
    #                   xColLabel: xCol,
    #                   yColLabel: yCol,
    #                   colorColLabel: colorCol,
    #                  }, index=pointLabels )

    #colorScale='RdBu'
    #colorScale='Rainbow'

    df = pd.DataFrame({
        xColLabel: xCol,
        yColLabel: yCol,
        colorColLabel: colorCol,
        pointLabelsHeader: pointLabels
    }, index=pointLabels)
    scatterplot = px.scatter(df, x=xColLabel, y=yColLabel,
                             hover_data=pointLabelsHeader,
                             title=plotTitle,
                             color=colorColLabel,
                             color_continuous_scale=colorScale,
                             color_continuous_midpoint=0,
                             range_color=[np.min(colorCol), -np.min(colorCol)]
                             )
    scatterplot.update_traces(opacity=0.0)
    # Add annotations with color-scaled text
    maxField = np.max(colorCol)
    minField = np.min(colorCol)
    rangeField = np.maximum(np.abs(maxField), np.abs(minField))
    normlzdColorCol = 0.5*(1*np.sort(colorCol))/rangeField + 0.5
    #normlzdColorCol = (colorCol - np.min(colorCol)) / \
    #                  (np.max(colorCol) - np.min(colorCol))
    dfSorted = df.sort_values(by=[colorColLabel], ascending=True)
    j = 0
    for i, row in dfSorted.iterrows():
        scatterplot.add_annotation(
            x=row[xColLabel],
            y=row[yColLabel],
            text=i,
            font=dict(color=pc.sample_colorscale( colorScale, normlzdColorCol[j])[0]),
            showarrow=False
        )
        j = j + 1
    #scatterplot.update_traces(textfont_color=col3)
    #scatterplot.data[0].marker.color
    #pc.sample_colorscale('Rainbow', 0.5)[0]
    #scatterplot.update_traces(textfont = dict(color='#f0f921'))
    #scatterplot.update_traces(textfont=dict(color='red'))
    #scatterplot.for_each_trace(lambda t: \
    #    t.update(textfont_color=pc.sample_colorscale('Rainbow', normlzdColorCol[t])[0],
    #             textposition='top center'))
    scatterplot.update_xaxes(title=xaxisTitle)
    scatterplot.update_yaxes(title=yaxisTitle)
    scatterplot.update_layout(showlegend=showLegend)
    scatterplot.update_layout(hovermode=hoverMode)
    scatterplot.update_layout(width=plotWidth, height=plotHeight)
    scatterplot.update_layout(plot_bgcolor=plotBgColor)

    return scatterplot


#def createMetricsBarChart( metricsNames, paramsNames, biases, normlzdResidBias, sensMatrix, title ):
def createMetricsBarChart(metricsNames, paramsNames,
                          defaultBiasesCol, defaultBiasesApproxNonlin,
                          normMetricValsCol,
                          sensMatrix,
                          title):
    import plotly.graph_objects as go
    import numpy as np
    import pdb

    #metricsNames = metricsNames[metricsSensOrder]

    biases = \
        -defaultBiasesCol / np.abs(normMetricValsCol)
    normlzdResidBias = (-defaultBiasesApproxNonlin - defaultBiasesCol) \
                       / np.abs(normMetricValsCol)

    biases = np.reshape(biases, (-1, 1))
    barBase = np.copy(biases)  # np.copy prevents biases variable from changing
    rightEnd = np.copy(biases)
    leftEnd = np.copy(biases)
    barsData = []
    for col in range(len(paramsNames)):
        #print("paramsNames[col]=", paramsNames[col])
        sensCol = sensMatrix[:, [col]]
        #print("sensCol=", sensCol )
        #print("rightEnd=", rightEnd )
        for row in range(len(sensCol)):
            if (np.sign(sensCol[row]) > 0):
                barBase[row] = rightEnd[row]
            else:
                barBase[row] = leftEnd[row]

        #print("barBase=", barBase)
        #print("biases during=", biases)
        barsData.append(go.Bar(name=paramsNames[col], y=metricsNames, x=sensCol[:, 0],
                               base=barBase[:, 0], orientation="h"))
        rightEnd = rightEnd + np.maximum(np.zeros_like(sensCol), sensCol)
        leftEnd = leftEnd + np.minimum(np.zeros_like(sensCol), sensCol)

    # Insert a narrow black horizontal line in each bar to denote the improvement wrought by tuning
    normlzdResidBias = np.reshape(normlzdResidBias, (-1, 1))
    barsData.append(go.Bar(name='+ tuning correction',
                           y=metricsNames, x=-normlzdResidBias[:, 0] + biases[:, 0], base=normlzdResidBias[:, 0],
                           orientation="h",
                           width=0.2,
                           marker_line_color='black', marker_color='black', marker_line_width=2,
                           opacity=1.0
                           )
                    )

    # Insert a black vertical line in each bar to denote default biases that we want to remove
    barsData.append(go.Bar(name='default bias',
                           y=metricsNames, x=np.zeros(len(metricsNames)), base=biases[:, 0],
                           orientation="h",
                           marker_line_color='black', marker_color='black', marker_line_width=5
                           )
                    )

    metricsBarChart = go.Figure(data=barsData)

    # Change the bar mode
    metricsBarChart.update_layout(title=title)
    metricsBarChart.update_layout(barmode='stack')
    metricsBarChart.update_xaxes(visible=True, zeroline=True, zerolinewidth=4, zerolinecolor='gray')  # Plot y axis
    metricsBarChart.update_layout(width=800, height=200 + 50 * len(metricsNames))
    metricsBarChart.update_xaxes(title="-Normalized biases")

    #pdb.set_trace()

    return metricsBarChart


def createBarChart(matrix, index, columns,
                   #                   barBase,
                   colorDiscreteMap=None,
                   orientation="v",
                   title=None,
                   xlabel=None, ylabel=None,
                   width=800, height=500,
                   showLegend=True):
    import plotly.express as px
    import plotly.graph_objects as go
    import pandas as pd

    df = pd.DataFrame(matrix,
                      index=index,
                      columns=columns)
    barChart = px.bar(df, x=df.index, y=df.columns,
                      #base=barBase,
                      #offset=1,
                      color=index,
                      color_discrete_map=colorDiscreteMap,
                      orientation=orientation,
                      title=title)
    barChart.update_xaxes(title=xlabel)
    barChart.update_yaxes(title=ylabel)
    barChart.update_layout(hovermode="closest")
    barChart.update_layout(showlegend=showLegend)
    barChart.update_yaxes(visible=True, zeroline=True, zerolinewidth=1, zerolinecolor='gray')  # Plot x axis
    barChart.update_layout(width=width, height=height)
    #barChart.update_layout(barmode='relative')

    return barChart


def createThreeDotFig(metricsNames, paramsNames, transformedParamsNames,
                      paramsAbbrv,
                      paramsSolnNonlin,
                      nonlinMatrixDparamsMasked,
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

    #    if ( len(paramsNames) != len(sens1NcFilenames)   ):
    #        print("Number of parameters must equal number of netcdf files.")
    #        quit()

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
    defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1, numParams))

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

    threeDotFig = make_subplots(rows=numMetrics, cols=numParams,
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
                           marker=dict(color='black', size=16)),
                row=arrayRow + 1,
                col=arrayCol + 1
            )

            # Plot orange symbol at tuned value of (parameter, metric)
            threeDotFig.add_trace(
                go.Scatter(x=paramsSolnNonlin[arrayCol,:],
                           y=[nonlinMatrixDparamsMasked[arrayRow,arrayCol]],
                           mode='markers',
                           marker_symbol='x',
                           marker=dict(color='orange', size=16)),
                row=arrayRow + 1,
                col=arrayCol + 1
            )

            # Plot quadratic curve passing through 3 dots
            paramPts = np.linspace(np.min(paramVals), np.max(paramVals), num=60)
            #print("paramPts =", paramPts)
            dnormlzdParamPts = np.linspace(normlzdOrdDparamsMin[arrayRow, arrayCol],
                                           normlzdOrdDparamsMax[arrayRow, arrayCol], num=60)
            metricPts = np.zeros_like(dnormlzdParamPts)
            for idx, dnormlzdParamPt in enumerate(dnormlzdParamPts):
                metricPts[idx] = 0.5 * normlzdCurvMatrix[arrayRow, arrayCol] * dnormlzdParamPt ** 2 \
                                 + normlzdSensMatrixPoly[arrayRow, arrayCol] * dnormlzdParamPt \
                                 + normlzdConstMatrix[arrayRow, arrayCol]
                metricPts[idx] = metricPts[idx] * np.abs(normMetricValsCol[arrayRow])

            #print("metricPts =", metricPts)
            threeDotFig.add_trace(
                go.Scatter(x=paramPts, y=metricPts,
                           mode='lines',
                           line=dict(color='blue', width=4)),
                row=arrayRow + 1,
                col=arrayCol + 1
            )

            # Calculate horizontal line at observed value
            threeObsMetricVals = np.squeeze(obsMetricValsCol[arrayRow][0] * np.ones((3, 1)))
            #threeObsMetricValsList = threeObsMetricVals.tolist()
            #print("obsMetricVals=", threeObsMetricValsList)
            #print("paramVals=", paramVals)            
            threeDotFig.add_trace(
                go.Scatter(x=paramVals, y=threeObsMetricVals,
                           mode='lines',
                           line=dict(color='red', width=4)),
                row=arrayRow + 1,
                col=arrayCol + 1
            )

            # Label the metric and parameter for the subplots
            paramsFontSize = np.rint(144.0/numParams)
            if (arrayRow == numMetrics - 1):  # Put params labels along bottom of plot
                threeDotFig.update_xaxes(title_text=paramsAbbrv[arrayCol],
                                         title_font_size=paramsFontSize,
                                         #font=dict(size=50),
                                         tickangle=45,
                                         row=arrayRow + 1, col=arrayCol + 1
                                         )
            if (arrayRow == 0):  # Put params labels along top of plot
                threeDotFig.update_xaxes(title_text=paramsAbbrv[arrayCol],
                                         row=arrayRow + 1, col=arrayCol + 1,
                                         title_font_size=paramsFontSize,
                                         #side="top", title_standoff=100
                                         )
            metricsFontSize = np.rint(300.0 / numMetrics)
            if (arrayCol == 0):  # Insert metrics label only along left edge of plot
                threeDotFig.update_yaxes(title_text=metricsNames[arrayRow], row=arrayRow + 1, col=arrayCol + 1,
                                         title_font_size=metricsFontSize)
            threeDotFig.update_layout(showlegend=False,
                                      title_text="Simulated metric values vs. parameter values",
                                      height=2500)

    threeDotFig.update_xaxes(tickangle=45)  # Put params label at 45-degree angle

    return (threeDotFig)


def createCorrArrayFig(matrix, indexLabels, title):
    import numpy as np
    import pandas as pd
    import plotly.figure_factory as ff
    import plotly.express as px
    import pdb

    cosAnglesMatrix = calcMatrixAngles(matrix)
    #cosAnglesMatrix = np.copy( matrix )
    roundedCosAnglesMatrix = np.around(cosAnglesMatrix, decimals=2)

    df = pd.DataFrame(roundedCosAnglesMatrix,
                      index=indexLabels,
                      columns=indexLabels)
    # Display only the lower-triangular elements of the matrix
    upTriMask = np.logical_not(np.tril(np.ones_like(roundedCosAnglesMatrix, dtype=bool)))
    df_mask = df.mask(upTriMask)
    #maskedRoundedCosAnglesMatrix = np.ma.masked_array(roundedCosAnglesMatrix, mask=upTriMask)
    #maskedRoundedCosAnglesMatrix.filled(np.nan)
    #print("maskedAngles =", maskedRoundedCosAnglesMatrix)
    #print("cosAnglesMatrix =", cosAnglesMatrix)
    #print("upTriMask =", upTriMask)
    corrArrayFig = ff.create_annotated_heatmap(
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
    corrArrayFig.update_xaxes(side="bottom")
    corrArrayFig.update_layout(
        title_text=title,
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

    #pdb.set_trace()

    return (corrArrayFig)


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

    if (len(paramsNames) != len(sensNcFilenames)):
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

    defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1, numParams))

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

    normlzdWeightedSensMatrixDiff = normlzdWeightedSensMatrixExt - normlzdWeightedSensMatrix
    sensMatrixDiff = sensMatrixExt - sensMatrix

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
                       transformedParamsNames)

    paramsLowValsPCBound = defaultParamValsOrigRow.T - 0.5 * np.abs(dparamsSolnPCBound)
    paramsHiValsPCBound = defaultParamValsOrigRow.T + 0.5 * np.abs(dparamsSolnPCBound)

    #pdb.set_trace()

    return (paramsLowValsPCBound, paramsHiValsPCBound)


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
    defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1, numParams))

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

    return (defaultBiasesCol, sensMatrix, normlzdWeightedSensMatrix)


def abbreviateClubbParamsNames(paramsNames):
    paramsAbbrv = np.char.replace(paramsNames, 'clubb_', '')
    paramsAbbrv = np.char.replace(paramsAbbrv, 'c_invrs_tau_', '')
    paramsAbbrv = np.char.replace(paramsAbbrv, 'wpxp_n2', 'n2')
    paramsAbbrv = np.char.replace(paramsAbbrv, 'altitude', 'alt')
    paramsAbbrv = np.char.replace(paramsAbbrv, 'threshold', 'thres')
    paramsAbbrv = np.char.replace(paramsAbbrv, 'thresh', 'thres')

    return paramsAbbrv


def createParamsErrorBarsFig(paramsAbbrv, defaultParamValsOrigRow, paramsScales,
                             paramsLowValsPCBound, paramsHiValsPCBound,
                             paramsSolnLin, dnormlzdParamsSolnLin,
                             paramsSolnNonlin, dnormlzdParamsSolnNonlin,
                             paramsSolnElastic, dnormlzdParamsSolnElastic):
    # Plot box and whiskers plot of optimal parameter values.
    # Multiply in the user-designated scale factors before plotting.
    df = pd.DataFrame(np.hstack(defaultParamValsOrigRow[0, :] * paramsScales),
                      index=paramsAbbrv, columns=["Default plus error bars"])
    df["err_minus"] = (defaultParamValsOrigRow[0, :] - paramsLowValsPCBound[:, 0]) * paramsScales
    df["err_plus"] = (paramsHiValsPCBound[:, 0] - defaultParamValsOrigRow[0, :]) * paramsScales
    paramsErrorBarsFig = px.scatter(df, x=df.index, y=df.columns,
                                    error_y="err_plus", error_y_minus="err_minus",
                                    title="""Best-fit parameter values""")
    paramsErrorBarsFig.update_traces(go.Scatter(
        mode='markers',
        marker=dict(color='black', size=14),
        error_y=dict(color='black', thickness=2, width=10)
    ))
    #paramsErrorBarsFig.add_trace(go.Scatter(x=paramsNames, y=paramsLowValsPCBound[:,0]*paramsScales,
    #                               name=r'$paramsSolnPC - \sigma$',
    #                               line=dict(color='white', width=0), mode='lines', showlegend=False))
    #paramsErrorBarsFig.add_trace(go.Scatter(x=paramsNames, y=paramsHiValsPCBound[:,0]*paramsScales, fill='tonexty',
    #                           name='Default Parameter Values +- sigma', mode='none',
    #                               fillcolor='rgba(253,253,150,1.0)'))
    #paramsErrorBarsFig.add_trace(go.Scatter(x=paramsNames, y=defaultParamValsOrigRow[0,:]*paramsScales,
    #                               name='Default Parameter Values', line=dict(color='black', width=6) ))
    paramsErrorBarsFig.add_trace(go.Scatter(x=paramsAbbrv, y=paramsSolnLin[:, 0] * paramsScales,
                                            mode='markers',
                                            marker=dict(color='green', size=8),
                                            name='Linear regression, |dp|='
                                                 + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSolnLin))))
    paramsErrorBarsFig.add_trace(go.Scatter(x=paramsAbbrv, y=paramsSolnNonlin[:, 0] * paramsScales,
                                            mode='markers',
                                            marker_symbol='x',
                                            marker=dict(color='orange', size=12),
                                            name='paramsSolnNonlin, |dpPC|='
                                                 + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSolnNonlin))))
    paramsErrorBarsFig.add_trace(go.Scatter(x=paramsAbbrv, y=paramsSolnElastic[:, 0] * paramsScales,
                                            mode='markers',
                                            marker_symbol='square',
                                            marker=dict(color='cyan', size=8),
                                            name='Lasso regression, |dpLasso|='
                                                 + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSolnElastic)),
                                            line=dict(color='red', width=2)))
    #paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnPCBound[:,0]*paramsScales,
    #                                name='paramsSolnPCBound, |dpBound|='
    #                               + '{:.2e}'.format(0.0) ))
    paramsErrorBarsFig.update_yaxes(title="User-scaled parameter value")
    paramsErrorBarsFig.update_xaxes(title="Parameter Name")
    paramsErrorBarsFig.update_layout(hovermode="closest")
    paramsErrorBarsFig.update_layout(width=1000, height=500)

    return paramsErrorBarsFig


def createParamsCorrArrayFig(matrix,
                             biasesCol,
                             paramsNames,
                             plotTitle):
    # Create color-coded matrix that displays correlations among parameter vectors
    normlzdSensMatrixConcatBiases = np.hstack((matrix, biasesCol))
    normlzdSensMatrixConcatBiasesDebiased = normlzdSensMatrixConcatBiases \
                                            - np.mean(normlzdSensMatrixConcatBiases, axis=0)
    #cosAnglesMatrix = calcMatrixAngles( normlzdSensMatrixConcatBiases.T )
    cosAnglesMatrix = calcMatrixAngles(normlzdSensMatrixConcatBiasesDebiased.T)
    roundedCosAnglesMatrix = np.around(cosAnglesMatrix, decimals=2)
    df = pd.DataFrame(roundedCosAnglesMatrix,
                      index=np.concatenate((paramsNames, ['bias vector'])),
                      columns=np.concatenate((paramsNames, ['bias vector'])))
    upTriMask = np.logical_not(np.tril(np.ones_like(roundedCosAnglesMatrix, dtype=bool)))
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
        title_text=plotTitle,
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

    return paramsCorrArrayFig


def minimize2ptDp(metricsNames, normMetricValsCol,
                  normlzdSensMatrix, defaultBiasesCol):
    # normlzdCurvMatrix,
    # reglrCoef,
    # beVerbose):

    import numpy as np
    # import pdb
    from scipy.optimize import minimize

    # from scipy.optimize import Bounds

    numMetrics = len(metricsNames)

    # pdb.set_trace()

    # Don't let parameter values go negative
    # lowerBoundsCol =  -defaultParamValsOrigRow[0]/magParamValsRow[0]

    # Perform nonlinear optimization
    normlzdDefaultBiasesCol = defaultBiasesCol / np.abs(normMetricValsCol)

    # dnormlzdParamsSolnNonlin = minimize(objFnc,x0=np.ones_like(np.transpose(defaultParamValsOrigRow)), \
    # dnormlzdParamsSolnNonlin = minimize(objFnc,x0=np.zeros_like(np.transpose(defaultParamValsOrigRow[0])), \
    # dnormlzdParamsSolnNonlin = minimize(objFnc,dnormlzdParamsSoln, \
    #                           args=(normlzdSensMatrix, normlzdDefaultBiasesCol, metricsWeights,
    #                           normlzdCurvMatrix, reglrCoef, numMetrics),\
    #                           method='Powell', tol=1e-12,
    #                           bounds=Bounds(lb=lowerBoundsCol) )
    # dnormlzdParamsSolnNonlin = np.atleast_2d(dnormlzdParamsSolnNonlin.x).T

    # dparamsSolnNonlin = dnormlzdParamsSolnNonlin * np.transpose(magParamValsRow)
    # paramsSolnNonlin = np.transpose(defaultParamValsOrigRow) + dparamsSolnNonlin
    # if beVerbose:
    #    print("paramsSolnNonlin.T=", paramsSolnNonlin.T)
    #    print("normlzdSensMatrix@dnPS.x.T=", normlzdSensMatrix @ dnormlzdParamsSolnNonlin)
    #    print("normlzdDefaultBiasesCol.T=", normlzdDefaultBiasesCol.T)
    #    print("normlzdSensMatrix=", normlzdSensMatrix)

    # normlzdWeightedDefaultBiasesApproxNonlin = \
    #         fwdFnc(dnormlzdParamsSolnNonlin, normlzdSensMatrix, normlzdCurvMatrix, numMetrics) \
    #         * metricsWeights

    # defaultBiasesApprox = (forward model soln - default soln)
    # defaultBiasesApproxNonlin = normlzdWeightedDefaultBiasesApproxNonlin \
    #                            * np.reciprocal(metricsWeights) * np.abs(normMetricValsCol)

    # To provide error bars, calculate solution with no nonlinear term and double the nonlinear term
    # defaultBiasesApproxNonlinNoCurv = \
    #         fwdFnc(dnormlzdParamsSolnNonlin, normlzdSensMatrix, 0*normlzdCurvMatrix, numMetrics) \
    #         * np.abs(normMetricValsCol)

    def dparamsSqdFnc(dnormlzdParamsRow):
        return (np.dot(dnormlzdParamsRow, dnormlzdParamsRow.T))

    def dpSqdJac(dnormlzdParamsRow):
        return (2. * dnormlzdParamsRow)

    # Two = 2
    # RowIdxs = slice(0, Two) # RowIdxs means "[0:Two]"

    # idx1 = 0
    # idx2 = 1

    dpMagMin2Row = np.zeros((numMetrics, numMetrics))

    for idx1 in np.arange(numMetrics - 1):  # each row in correlation matrix
        for idx2 in np.arange(idx1 + 1, numMetrics):  # each col in correlation matrix
            rowIdxs = [idx1, idx2]
            normlzdSensMatrix2Rows = normlzdSensMatrix[rowIdxs, :]
            normlzdDefaultBiases2Rows = normlzdDefaultBiasesCol[rowIdxs].flatten()

            #    normlzdSensRowEqns = {'type': 'eq',
            #            'fun': lambda x: -normlzdDefaultBiasesCol[0:2] - np.dot(normlzdSensMatrix2Rows, x.T),
            #            'jac': lambda x: -normlzdSensMatrix2Rows}

            normlzdSensRowEqns = {'type': 'eq',
                                  'fun': lambda x: -normlzdDefaultBiases2Rows - np.dot(normlzdSensMatrix2Rows, x)
                                  }

            dnormlzdParamsMin2pt = minimize(dparamsSqdFnc,
                                            # x0=np.zeros_like(np.transpose(defaultParamValsOrigRow[0])),
                                            x0=np.zeros(normlzdSensMatrix.shape[1]),
                                            jac=dpSqdJac,
                                            constraints=normlzdSensRowEqns,
                                            method='SLSQP')

            dpMagMin2Row[idx2, idx1] = np.sqrt(dnormlzdParamsMin2pt.fun)
            # print("dpMagMin2Row=", idx2, idx1, dpMagMin2Row[idx2,idx1])

    return (dpMagMin2Row)


def createDpMin2PtFig(normlzdSensMatrixPoly, defaultBiasesCol,
                      normMetricValsCol, metricsNames):
    cosAnglesMatrix = calcMatrixAngles(normlzdSensMatrixPoly)
    invrsCosFactorMinusMatrix = np.power(np.maximum(np.finfo(float).eps, 2. * (1. - cosAnglesMatrix)), -0.5)
    invrsCosFactorPlusMatrix = np.power(np.maximum(np.finfo(float).eps, 2. * (1. + cosAnglesMatrix)), -0.5)
    # dbOnAbsSensVector = column vector = dbias / magnitude_of_each_row_of_sensitivity_matrix
    dbOnAbsSensVector = \
        -defaultBiasesCol / np.abs(normMetricValsCol) \
        / np.linalg.norm(normlzdSensMatrixPoly, axis=1).reshape(-1, 1)
    # Create matrix whose (i,l) element = cosine_factor_il * ( (dbias/mag_sens)_i - (dbias/mag_sens)_l )
    dbOnAbsSensMatrix1 = np.ones((len(metricsNames), 1)) @ dbOnAbsSensVector.T
    dbOnAbsSensMatrix2 = dbOnAbsSensVector @ np.ones((1, len(metricsNames)))
    dpMin2PtMinusMatrix = invrsCosFactorMinusMatrix * \
                          np.abs(dbOnAbsSensMatrix2 - dbOnAbsSensMatrix1)
    dpMin2PtPlusMatrix = invrsCosFactorPlusMatrix * \
                         np.abs(dbOnAbsSensMatrix2 + dbOnAbsSensMatrix1)
    dpMin2PtMatrix = np.maximum(dpMin2PtMinusMatrix, dpMin2PtPlusMatrix)

    # Should this be based on the nonlinear matrix??
    dpMin2PtMatrix = minimize2ptDp(metricsNames, normMetricValsCol,
                                   normlzdSensMatrixPoly, defaultBiasesCol)

    roundedDpMin2PtMatrix = np.around(dpMin2PtMatrix, decimals=2)
    #dpMin2PtMatrix = np.fill_diagonal(roundedDpMin2PtMatrix, np.nan)
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
        title_text='Minimum size of parameter perturbation between 2 metrics',
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

    return dpMin2PtFig


def createNormlzdSensMatrixColsFig(defaultBiasesCol, normlzdSensMatrixPoly,
                                   normMetricValsCol, metricsNames, paramsNames):
    # Plot each column of normalized sensitivity matrix as a separate line.
    # Each column tells us how all metrics vary with a single parameter.
    df = pd.DataFrame(np.hstack((-defaultBiasesCol / np.abs(normMetricValsCol), normlzdSensMatrixPoly)),
                      index=metricsNames,
                      columns=np.append('Norm bias', paramsNames))
    normlzdSensMatrixColsFig = px.line(df, x=df.index, y=df.columns,
                                       title="""Columns of normalized, unweighted sensitivity matrix (plus the bias vector).<br>
                       Each column (line) shows how sensitive the metrics are to a change in a single parameter value.<br>
                       (A positive value means that an increase in parameter value brings the default simulation closer to obs.)""")
    normlzdSensMatrixColsFig.update_yaxes(title="Norml sens, (|param|/|obsmetric|) * dmetric/dparam")
    normlzdSensMatrixColsFig.update_xaxes(title="Regional metric")
    normlzdSensMatrixColsFig.layout.legend.title = "Parameter"
    normlzdSensMatrixColsFig.update_layout(hovermode="closest")
    #pdb.set_trace()

    return normlzdSensMatrixColsFig


def createNormlzdSensMatrixRowsFig(normlzdSensMatrixPoly,
                                   metricsNames, paramsNames):
    # Plot each row of normalized sensitivity matrix as a separate line.
    # Each row tells us how a single metric varies with all parameters.
    df = pd.DataFrame(np.transpose(normlzdSensMatrixPoly),
                      index=paramsNames,
                      columns=metricsNames)
    normlzdSensMatrixRowsFig = px.line(df, x=df.index, y=df.columns,
                                       title="""Rows of normalized, unweighted sensitivity matrix.<br>
                       Each row (line) tells us how a single metric varies with all parameters.<br>
                       (A positive value means that an increase in parameter value brings the default simulation closer to obs.)""")
    normlzdSensMatrixRowsFig.update_yaxes(title="Norml sens, (|param|/|obsmetric|) * dmetric/dparam")
    normlzdSensMatrixRowsFig.update_xaxes(title="Parameter")
    normlzdSensMatrixRowsFig.layout.legend.title = "Metric"
    normlzdSensMatrixRowsFig.update_layout(hovermode="closest")

    return normlzdSensMatrixRowsFig


def createBiasesOrderedArrowFig(metricsSensOrder, metricsNamesOrdered,
                                defaultBiasesCol, normMetricValsCol,
                                defaultBiasesApproxNonlinNoCurv, defaultBiasesApproxNonlin2xCurv,
                                defaultBiasesApproxNonlin,
                                linSolnBiasesCol):
    # Plot a black dot for each default-run bias
    biasesOrderMatrix = np.dstack((-defaultBiasesCol[metricsSensOrder])).squeeze()
    fracBiasesOrderMatrix = np.diagflat(np.reciprocal(np.abs(normMetricValsCol[metricsSensOrder]))) \
                            @ np.expand_dims(np.atleast_1d(biasesOrderMatrix), axis=1)
    #                           @ biasesOrderMatrix
    df = pd.DataFrame(fracBiasesOrderMatrix,
                      index=metricsNamesOrdered,
                      columns=['fracDefBias'])
    biasesOrderedArrowFig = px.line(df, x=df.index, y=df.columns,
                                    title="""<span style='color:blue'>Predicted</span> and <span style='color:red'>actual</span> removal of regional biases""")
    biasesOrderedArrowFig.update_yaxes(title="-Normalized bias")
    biasesOrderedArrowFig.update_xaxes(title="Regional metric")
    biasesOrderedArrowFig.update_layout(hovermode="closest")
    biasesOrderedArrowFig.update_layout(showlegend=False)
    biasesOrderedArrowFig.update_traces(mode='markers', line_color='black')  # Plot default biases as black dots
    biasesOrderedArrowFig.update_yaxes(visible=True, zeroline=True, zerolinewidth=1,
                                       zerolinecolor='gray')  # Plot x axis
    biasesOrderedArrowFig.update_layout(width=700, height=500)
    # Now plot an arrow for each region that points from default-run bias to new bias after tuning
    xArrow = np.arange(len(metricsNamesOrdered))  # x-coordinate of arrows
    yArrow = -defaultBiasesCol[metricsSensOrder, 0] / np.abs(normMetricValsCol[metricsSensOrder, 0])
    gap = 0.2  # horizontal spacing between arrows
    # Plot error bar on prediction arrow.  Bar runs between 0- and 2x-curvature solns.
    for i, item in enumerate(metricsNamesOrdered):
        biasesOrderedArrowFig.add_annotation(
            x=xArrow[i] - gap,  # ith arrow's head
            # ith arrow's head:
            y=(-defaultBiasesApproxNonlinNoCurv - defaultBiasesCol)[metricsSensOrder[i], 0] / np.abs(
                normMetricValsCol[metricsSensOrder[i], 0]),
            ax=xArrow[i] - gap,  # ith arrow's head
            # ith arrow's head:
            ay=(-defaultBiasesApproxNonlin2xCurv - defaultBiasesCol)[metricsSensOrder[i], 0] / np.abs(
                normMetricValsCol[metricsSensOrder[i], 0]),
            font=dict(family="bold", color="blue", size=30),
            showarrow=True,
            xref='x',
            yref='y',
            axref='x',
            ayref='y',
            text='',  # blank because we want only the arrow
            arrowhead=0,
            arrowsize=1,
            arrowwidth=6,
            arrowcolor='lightskyblue'
            # https://stackoverflow.com/questions/72496150/user-friendly-names-for-plotly-css-colors
        )
    #biasesOrderedArrowFig.add_scatter(x=df.index, y=df.columns, line_color='pink')  # attempt to make black dot appear on top
    biasesOrderedArrowFig.update_traces(mode='markers', line_color='black')
    # Plot arrows showing the tuner's nonlinear predicted bias removal
    for i, item in enumerate(metricsNamesOrdered):
        biasesOrderedArrowFig.add_annotation(
            x=xArrow[i] - gap,  # ith arrow's head
            # ith arrow's head:
            y=(-defaultBiasesApproxNonlin - defaultBiasesCol)[metricsSensOrder[i], 0] / np.abs(
                normMetricValsCol[metricsSensOrder[i], 0]),
            #y= (-defaultBiasesApproxNonlinNoCurv-defaultBiasesCol)[metricsSensOrder[i],0]/np.abs(normMetricValsCol[metricsSensOrder[i],0]),
            ax=xArrow[i] - gap,  # ith arrow's tail
            ay=yArrow[i],  # ith arrow's tail
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
            x=xArrow[i] + gap,  # ith arrow's head
            # ith arrow's head:
            #y= (-linSolnBiasesCol-defaultBiasesCol)[metricsSensOrder[i],0]/np.abs(normMetricValsCol[metricsSensOrder[i],0]),
            y=(-linSolnBiasesCol)[metricsSensOrder[i], 0] / np.abs(normMetricValsCol[metricsSensOrder[i], 0]),
            ax=xArrow[i] + gap,  # ith arrow's tail
            ay=yArrow[i],  # ith arrow's tail
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
    #         y = (-defaultBiasesApproxNonlinNoCurv-defaultBiasesCol)[metricsSensOrder[i],0]/np.abs(normMetricValsCol[metricsSensOrder[i],0]),
    #         text ='-',  # plot horizontal line
    #         font = dict(family="bold", color="blue", size=30),
    #         showarrow=False
    #       )
    #     # Plot 2x-curvature error bars on prediction arrow
    #     for i, item in enumerate(metricsNamesOrdered):
    #         biasesOrderedArrowFig.add_annotation(
    #         x =  xArrow[i]-gap,  # ith arrow's head
    #         # ith arrow's head:
    #         y = (-defaultBiasesApproxNonlin2xCurv-defaultBiasesCol)[metricsSensOrder[i],0]/np.abs(normMetricValsCol[metricsSensOrder[i],0]),
    #         text ='-',  # plot horizontal line
    #         font = dict(family="bold", color="blue", size=30),
    #         showarrow=False
    #       )

    #biasesOrderedArrowFig.add_trace(go.Scatter(x=xArrow, y=yArrow,
    #                      name='Region of improvement', mode='markers',
    #                       marker=dict(color='green', size=14)))
    #pdb.set_trace()

    #biasesOrderedArrowFig.write_image('biasesOrderedArrowFig.png', scale=6)

    return biasesOrderedArrowFig


def createBiasLinNlIndivContrbBarFig(normlzdSensParamsMatrixOrdered, curvParamsMatrixOrdered,
                                     metricsNamesOrdered, paramsNames):
    #dnormlzdParamsSolnNonlinMatrix = np.ones((len(metricsNames),1)) @ dnormlzdParamsSolnNonlin.T
    #curvParamsMatrixOrdered = 0.5 * normlzdCurvMatrix[metricsSensOrder,:] * dnormlzdParamsSolnNonlinMatrix**2
    #print("Sum rows=", np.sum(-normlzdSensParamsMatrixOrdered-curvParamsMatrixOrdered, axis=1))
    dfLin = pd.DataFrame(-1 * normlzdSensParamsMatrixOrdered,
                         index=metricsNamesOrdered,
                         columns=paramsNames)
    # biasLinNlIndivContrbBarFig = go.Figure()
    # biasLinNlIndivContrbBarFig.add_trace( go.Bar(x=df.index, y=list(curvParamsMatrixOrdered[:,2:3])) )
    #df_long = pd.wide_to_long( df, i=df.index, j=df.columns, stubnames=[''] )
    dfLin = dfLin.reset_index()
    dfLin.rename(columns={'index': 'metricsNamesOrdered'}, inplace=True)
    #print("biasContrGrouped df=", df.to_string())
    #print("df.columns.values=", df.columns.values)
    #df.columns[1] = ['metricsNamesOrdered']

    dfLin_long = dfLin.melt(id_vars='metricsNamesOrdered',
                            var_name='paramsNames', value_name='Contribution to bias removal')
    dfLin_long.insert(0, 'isNonlin', ['linear'] * len(paramsNames) * len(metricsNamesOrdered))
    #print("df_long=", df_long.to_string())

    dfNonlin = pd.DataFrame(-1 * curvParamsMatrixOrdered,
                            index=metricsNamesOrdered,
                            columns=paramsNames)
    # biasLinNlIndivContrbBarFig = go.Figure()
    # biasLinNlIndivContrbBarFig.add_trace( go.Bar(x=df.index, y=list(curvParamsMatrixOrdered[:,2:3])) )
    #df_long = pd.wide_to_long( df, i=df.index, j=df.columns, stubnames=[''] )
    dfNonlin = dfNonlin.reset_index()
    dfNonlin.rename(columns={'index': 'metricsNamesOrdered'}, inplace=True)
    #print("biasContrGrouped df=", df.to_string())
    #print("df.columns.values=", df.columns.values)
    #df.columns[1] = ['metricsNamesOrdered']

    dfNonlin_long = dfNonlin.melt(id_vars='metricsNamesOrdered',
                                  var_name='paramsNames', value_name='Contribution to bias removal')
    dfNonlin_long.insert(0, 'isNonlin', ['nonlinear'] * len(paramsNames) * len(metricsNamesOrdered))
    #print("df_long=", df_long.to_string())    

    dfLinNonlin_long = pd.concat([dfLin_long, dfNonlin_long], ignore_index=True)

    biasLinNlIndivContrbBarFig = px.bar(dfLinNonlin_long,
                                        facet_col='metricsNamesOrdered', y='Contribution to bias removal',
                                        facet_col_spacing=0.005,  # default 0.03
                                        x='isNonlin', color='paramsNames')  #,
    #title = """Long: Linear ++ nonlinear contributions to actual removal of regional biases""")
    #biasLinNlIndivContrbBarFig.update_yaxes(title="-(Def-Sim) / abs(obs metric value)")
    #biasLinNlIndivContrbBarFig.update_xaxes(title="Regional metric")
    ##biasLinNlIndivContrbBarFig.update_xaxes(visible=False)
    ##biasLinNlIndivContrbBarFig.update_yaxes(visible=False)
    biasLinNlIndivContrbBarFig.update_layout(hovermode="closest")
    biasLinNlIndivContrbBarFig.update_layout(showlegend=True)
    ##biasLinNlIndivContrbBarFig.update_xaxes(showticklabels=False).update_yaxes(showticklabels=False)
    #biasLinNlIndivContrbBarFig.for_each_annotation(lambda a: a.update(text=''))
    biasLinNlIndivContrbBarFig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    biasLinNlIndivContrbBarFig.update_annotations(textangle=-90)
    biasLinNlIndivContrbBarFig.update_layout(margin=dict(t=160))
    for axis in biasLinNlIndivContrbBarFig.layout:
        #if type(biasLinNlIndivContrbBarFig.layout[axis]) == go.layout.YAxis:
        #    biasLinNlIndivContrbBarFig.layout[axis].title.text = 'Contribution to bias removal'
        if type(biasLinNlIndivContrbBarFig.layout[axis]) == go.layout.XAxis:
            biasLinNlIndivContrbBarFig.layout[axis].title.text = ''
    biasLinNlIndivContrbBarFig.layout.title.text = "Linear and nonlinear contributions to removal of mean bias <br> (dMetrics=sensMatrix*dParams) <br><br><br>"
    #biasLinNlIndivContrbBarFig.update_yaxes(visible=True,zeroline=True,zerolinewidth=1,zerolinecolor='gray') # Plot x axis
    biasLinNlIndivContrbBarFig.update_layout(title_y=0.95)
    biasLinNlIndivContrbBarFig.update_layout(width=1000, height=450)
    #print("curvParams =", -1*curvParamsMatrixOrdered)
    #print("normlzdSens =", -1*normlzdSensParamsMatrixOrdered)

    return biasLinNlIndivContrbBarFig


def createBiasVsBiasApproxScatterplot(defaultBiasesApproxNonlin, defaultBiasesCol,
                                      normMetricValsCol,
                                      metricsNames):
    # Plot a scatterplot of default-simulation bias and SVD approximation of that bias.
    # Each column tells us how all metrics vary with a single parameter.
    biasSensDirMatrix = np.concatenate((defaultBiasesApproxNonlin / np.abs(normMetricValsCol),
                                        (-defaultBiasesCol / np.abs(normMetricValsCol))), axis=1)
    biasAndParamsNames = ["biasApproxNonlin", "bias"]
    #biasAndParamsNames = np.append(["bias", "bias_approx_pc"], paramsNames)
    df = pd.DataFrame(biasSensDirMatrix,
                      index=metricsNames,
                      columns=biasAndParamsNames)
    biasSensDirMatrixScatter = px.scatter(df, x="biasApproxNonlin", y="bias",
                                          text=metricsNames, title="Bias approx vs bias")
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
    biasVsBiasApproxScatterplot.update_yaxes(visible=True, zeroline=True, zerolinewidth=2,
                                             zerolinecolor='lightblue')  # Plot x axis
    biasVsBiasApproxScatterplot.update_layout(width=800, height=500)
    biasVsBiasApproxScatterplot.update_layout(title="Bias approx vs bias")

    return biasVsBiasApproxScatterplot


def createBiasVsDiagnosticScatterplot(diagnosticPrefix, defaultBiasesCol,
                                      normMetricValsCol,
                                      defaultNcFilename):
    from set_up_dashboard_inputs import setUp_x_MetricsList, \
        setupDefaultMetricValsCol

    diagnosticNamesWeightsAndNorms, dfdiagnosticMetricGlobalVal = setUp_x_MetricsList(diagnosticPrefix,
                                                                                      defaultNcFilename)
    # Split up the list above into metric names and the corresponding weights.
    dfdiagnosticNamesWeightsAndNorms = \
        pd.DataFrame(diagnosticNamesWeightsAndNorms,
                     columns=['diagnosticNames', 'diagnosticWeights', 'diagnosticNorms'])
    diagnosticNames = dfdiagnosticNamesWeightsAndNorms[['diagnosticNames']].to_numpy().astype(str)[:, 0]

    # Set up a column vector of metric values from the default simulation
    diagnosticValsCol = \
        setupDefaultMetricValsCol(diagnosticNames, defaultNcFilename)

    # Plot a scatterplot of default-simulation bias and SVD approximation of that bias.
    # Each column tells us how all metrics vary with a single parameter.
    biasDiagnosticMatrix = np.concatenate((diagnosticValsCol,
                                           (-defaultBiasesCol / np.abs(normMetricValsCol))), axis=1)
    biasAndDiagnosticNames = [diagnosticPrefix, "bias"]
    #biasAndParamsNames = np.append(["bias", "bias_approx_pc"], paramsNames)
    df = pd.DataFrame(biasDiagnosticMatrix,
                      index=diagnosticNames,
                      columns=biasAndDiagnosticNames)
    biasVsDiagnosticScatterplot = px.scatter(df, x=diagnosticPrefix, y="bias",
                                             text=diagnosticNames, title="Bias vs diagnostic")
    #biasSensDirMatrixOneOneLine = px.line(df, x="bias", y="bias")
    #biasSensDirMatrixOneMOneLine = px.line(df, x="bias", y=-df.loc[:,"bias"])
    #biasVsDiagnosticScatterplot = go.Figure(data=biasSensDirMatrixScatter.data
    #                                          + biasSensDirMatrixOneOneLine.data)
    #biasRange = (max(df.loc[:,"bias"]), min(df.loc[:,"bias"]))
    #biasVsDiagnosticScatterplot.add_trace(go.Scatter(x=biasRange, y=biasRange, fill='tozeroy',
    #                           name='Region of improvement', mode='none',
    #                           fillcolor='rgba(253,253,150,0.7)'))
    biasVsDiagnosticScatterplot.update_xaxes(title=diagnosticPrefix)
    biasVsDiagnosticScatterplot.update_yaxes(title="-defaultBiasesCol/obs")
    biasVsDiagnosticScatterplot.update_traces(textposition='top center')
    biasVsDiagnosticScatterplot.update_yaxes(visible=True, zeroline=True, zerolinewidth=2,
                                             zerolinecolor='lightblue')  # Plot x axis
    biasVsDiagnosticScatterplot.update_layout(width=800, height=500)
    biasVsDiagnosticScatterplot.update_layout(title="Bias vs diagnostic")

    return biasVsDiagnosticScatterplot


def createBiasSensMatrixScatterFig(defaultBiasesCol, defaultBiasesApproxElastic,
                                   normMetricValsCol, metricsNames):
    # Plot a scatterplot of default-simulation bias and SVD approximation of that bias.
    # Each column tells us how all metrics vary with a single parameter.
    biasSensMatrix = np.concatenate((-defaultBiasesCol / np.abs(normMetricValsCol),
                                     (-defaultBiasesApproxElastic - defaultBiasesCol) / np.abs(normMetricValsCol)),
                                    axis=1)
    #defaultBiasesApproxElastic/np.abs(normMetricValsCol)), axis=1)
    biasAndParamsNames = ["bias", "bias_approx_pc"]
    #biasAndParamsNames = np.append(["bias", "bias_approx_pc"], paramsNames)
    df = pd.DataFrame(biasSensMatrix,
                      index=metricsNames,
                      columns=biasAndParamsNames)
    biasSensMatrixScatter = px.scatter(df, x="bias", y="bias_approx_pc", text=metricsNames,
                                       #biasSensMatrixScatter = px.scatter(df, x=np.append(["bias_approx_pc"], paramsNames), y="bias",
                                       title="""Columns of normalized sensitivity matrix.<br>
                       vs. bias vector.<br>
                       """)
    biasSensMatrixOneOneLine = px.line(df, x="bias", y="bias")
    biasSensMatrixOneMOneLine = px.line(df, x="bias", y=-df.loc[:, "bias"])
    biasSensMatrixScatterFig = go.Figure(data=biasSensMatrixScatter.data
                                              + biasSensMatrixOneOneLine.data
                                              + biasSensMatrixOneMOneLine.data)
    biasRange = (max(df.loc[:, "bias"]), min(df.loc[:, "bias"]))
    biasSensMatrixScatterFig.add_trace(go.Scatter(x=biasRange, y=biasRange, fill='tozeroy',
                                                  name='Region of improvement', mode='none',
                                                  fillcolor='rgba(253,253,150,0.7)'))
    biasSensMatrixScatterFig.update_yaxes(title="(-defaultBiasesApproxElastic-defaultBiasesCol)/obs")
    biasSensMatrixScatterFig.update_xaxes(title="-defaultBiasesCol/obs")
    biasSensMatrixScatterFig.update_traces(textposition='top center')

    return biasSensMatrixScatterFig


def createDpMinMatrixScatterFig(defaultBiasesCol, normlzdSensMatrixPoly,
                                normMetricValsCol, metricsNames):
    # Plot a scatterplot of minimum parameter perturbation vs. fractional default bias approximation
    # Calculate lower bound on normalized parameter perturbations
    #normlzdDefaultBiasesCol = ( metricsWeights * (-defaultBiasesCol) /
    normlzdDefaultBiasesCol = ((defaultBiasesCol) /
                               np.abs(normMetricValsCol))
    #sensMatrixRowMag = np.linalg.norm(normlzdWeightedSensMatrixPoly, axis=1)
    sensMatrixRowMag = np.linalg.norm(normlzdSensMatrixPoly, axis=1)
    #sensMatrixRowMag = np.amax(np.abs(normlzdSensMatrixPoly), axis=1)
    dpMin = np.abs(normlzdDefaultBiasesCol) / np.atleast_2d(sensMatrixRowMag).T
    #u_dot_b = np.atleast_2d(sensMatrixRowMag).T * -normlzdDefaultBiasesCol
    dpMinMatrix = np.dstack((np.reciprocal(dpMin),
                             #dpMinMatrix = np.dstack((np.abs(u_dot_b),
                             ##dpMinMatrix = np.dstack((np.atleast_2d(sensMatrixRowMag).T,
                             #                      np.abs(defaultBiasesApproxElastic)/np.abs(normMetricValsCol)
                             np.abs(defaultBiasesCol) / np.abs(normMetricValsCol)
                             )).squeeze()
    biasAndParamsNames = ["dpMinInvrs", "bias_approx"]
    df = pd.DataFrame(dpMinMatrix,
                      index=metricsNames,
                      columns=biasAndParamsNames)
    dpMinMatrixScatter = px.scatter(df, x="dpMinInvrs", y="bias_approx", text=metricsNames,
                                    title="""dpMinInvrs  vs. |approx bias vector|.<br>
                       """)
    dpMinMatrixScatterFig = go.Figure(data=dpMinMatrixScatter.data)
    dpMinMatrixScatterFig.update_yaxes(title="|defaultBiasesApproxElastic|")
    dpMinMatrixScatterFig.update_xaxes(title="dpMinInvrs")
    dpMinMatrixScatterFig.update_traces(textposition='top center')

    return dpMinMatrixScatterFig


def createBiasesVsSensMagScatterplot(normlzdSensMatrixPoly,
                                     defaultBiasesCol, defaultBiasesApproxNonlin,
                                     normMetricValsCol, metricsNames):
    ## Plot the biases versus sensitivity of each regional metric.
    ##    More specifically, plot the maximum magnitude value of each row of the sensitivity matrix.
    #df = pd.DataFrame({
    ##df = pd.DataFrame({'Max abs normlzd sensitivity': np.max(np.abs(normlzdSensMatrixPoly), axis=1), # max |row elements|
    ##df = pd.DataFrame({'Max abs normlzd sensitivity': np.sum(normlzdWeightedSensMatrixPoly, axis=1), # sum of row elements
    ##df = pd.DataFrame({'Max abs normlzd sensitivity': np.linalg.norm(normlzdWeightedSensMatrixPoly, axis=1), # rms of row elements
    ##df = pd.DataFrame({'Max abs normlzd sensitivity': np.linalg.norm(normlzdSensMatrixPoly, axis=1), # rms of row elements
    ##df = pd.DataFrame({'Max abs normlzd sensitivity':
    ##                    -defaultBiasesCol[:,0]/np.abs(normMetricValsCol[:,0])*np.linalg.norm(normlzdWeightedSensMatrixPoly, axis=1), # sum of row elements
    #                   'Max abs normlzd lin+ sensitivity': np.linalg.norm(normlzdLinplusSensMatrixPoly, axis=1), # sum of row elements
    #                   'Default biases': np.abs(-defaultBiasesCol[:,0])/np.abs(normMetricValsCol[:,0]),
    ##                   'revised tuning': (-defaultBiasesApproxElastic-defaultBiasesCol)[:,0]/np.abs(normMetricValsCol[:,0])
    #                  }, index=metricsNames )
    ##biasesVsSensMagScatterplot = px.scatter(df, x='Max abs normlzd sensitivity', y=df.columns[1:2],
    ##biasesVsSensMagScatterplot = px.scatter(df, x='Max abs normlzd sensitivity', y='default tuning',
    #biasesVsSensMagScatterplot = px.scatter(df, x=['Max abs normlzd lin+ sensitivity'], y='Default biases',
    #                             text=metricsNames,
    #                             title = """Regional biases vs. magnitude of sensitivity.""")
    #biasesVsSensMagScatterplot.update_yaxes(title="Regional biases")
    #biasesVsSensMagScatterplot.update_xaxes(title="Magnitude of sensitivity of regional metrics to parameter changes")
    #biasesVsSensMagScatterplot.update_layout(showlegend=False)
    #biasesVsSensMagScatterplot.update_layout(hovermode="closest")
    #biasesVsSensMagScatterplot.update_layout( width=700, height=500  )

    xCol = np.linalg.norm(normlzdSensMatrixPoly, axis=1)
    #yCol = -defaultBiasesCol[:,0]/np.abs(normMetricValsCol[:,0])
    yCol = (-defaultBiasesApproxNonlin - defaultBiasesCol)[:, 0] \
           / np.abs(normMetricValsCol[:, 0])

    biasesVsSensMagScatterplot = \
        createScatterplot(xCol=xCol, xColLabel='sens',
                          yCol=yCol, yColLabel='bias',
                          colorCol=yCol, colorColLabel='bias',
                          colorScale='Rainbow',
                          pointLabels=metricsNames, pointLabelsHeader='Metric',
                          plotTitle="""Regional biases vs. magnitude of sensitivity.""",
                          xaxisTitle="Magnitude of sensitivity of regional metrics to parameter changes",
                          yaxisTitle="Regional biases",
                          showLegend=False, hoverMode="closest",
                          plotWidth=700, plotHeight=500)

    return biasesVsSensMagScatterplot


def createBiasesVsSensArrowFig(normlzdWeightedSensMatrixPoly, defaultBiasesCol,
                               defaultBiasesApproxNonlin,
                               normMetricValsCol, metricsNames):
    # Compute length of arrows between default and tuned biases
    #metricsNamesPadded = ",,".join(metricsNames).split(",")
    #metricsNamesPadded = ",,".join(metricsNamesPadded).split(",")
    #metricsNamesPadded = np.append(metricsNamesPadded, ["", "", ""], axis=0)
    xArrow = np.linalg.norm(normlzdWeightedSensMatrixPoly, axis=1)
    yArrow = -defaultBiasesCol[:, 0] / np.abs(normMetricValsCol[:, 0])
    #uArrow = np.zeros_like(xArrow)
    #vArrow = (-defaultBiasesApproxElasticNonlin)[:,0]/np.abs(normMetricValsCol[:,0])
    #arrowFig = create_quiver(xArrow, yArrow, uArrow, vArrow,
    #                         scale=1,text=metricsNamesPadded)
    #arrowFig.update_yaxes(title="Regional biases")
    #arrowFig.update_traces(mode='lines+text')  # make labels appear in plot, not just hovermode

    # Plot biases vs. sensitivity, but with arrows indicating the degree of bias reduction
    df = pd.DataFrame(
        {'Max abs normlzd sensitivity': np.linalg.norm(normlzdWeightedSensMatrixPoly, axis=1),  # sum of row elements

         'default tuning': -defaultBiasesCol[:, 0] / np.abs(normMetricValsCol[:, 0]),
         'revised tuning': (-defaultBiasesApproxNonlin - defaultBiasesCol)[:, 0] / np.abs(normMetricValsCol[:, 0])
         }, index=metricsNames)
    biasesVsSensArrowFig = px.scatter(df, x='Max abs normlzd sensitivity', y=df.columns[1:2],
                                      text=metricsNames,
                                      title="""Regional biases with default and nonlin tuning versus sensitivity, with arrows.""")
    biasesVsSensArrowFig.update_traces(textposition="middle right")
    for i, item in enumerate(metricsNames):
        biasesVsSensArrowFig.add_annotation(
            x=xArrow[i],  # ith arrow's head
            y=(-defaultBiasesApproxNonlin - defaultBiasesCol)[i, 0] / np.abs(normMetricValsCol[i, 0]),
            # ith arrow's length
            #y= (-defaultBiasesApproxNonlin2x-defaultBiasesCol)[i,0]/np.abs(normMetricValsCol[i,0]),  # ith arrow's length
            ax=xArrow[i],  # ith arrow's tail
            ay=yArrow[i],  # ith arrow's tail
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
    biasesVsSensArrowFig.update_layout(hovermode="closest")
    biasesVsSensArrowFig.update_traces(cliponaxis=False)
    biasesVsSensArrowFig.update_yaxes(automargin=True)

    return biasesVsSensArrowFig


def createMaxSensMetricsFig(normlzdSensMatrixPoly, metricsNames):
    # Plot the sensitivity of each regional metric.
    #    More specifically, plot the maximum magnitude value of each row of the sensitivity matrix.
    df = pd.DataFrame(np.max(np.abs(normlzdSensMatrixPoly), axis=1),  # max of absolute val of each row
                      index=metricsNames,
                      columns=['Max abs normlzd sensitivity'])
    maxSensMetricsFig = px.line(df, x=df.index, y=df.columns,
                                title="""Maximum normalized sensitivity of each metric with respect to parameters.<br>
                       (Low sensitivity means that the metric is unbudgeable by these parameters.)""")
    maxSensMetricsFig.update_yaxes(title="Max |sens row|")
    maxSensMetricsFig.update_xaxes(title="Regional metric")
    maxSensMetricsFig.update_layout(hovermode="closest")
    maxSensMetricsFig.update_traces(mode='lines+markers')

    return maxSensMetricsFig


def calcMatrixVectorAngles(matrix, row):
    '''Calculate cos(angle) between one row of a matrix and all rows of the same matrix.
       Returns a column vector, with length equal to the number of rows in the matrix.'''

    import sklearn

    normed_matrix = sklearn.preprocessing.normalize(matrix, axis=1, norm='l2')

    cosAngles = normed_matrix @ normed_matrix[row, :].T

    return cosAngles


def calcMatrixAngles(matrix):
    '''Calculate cos(angle) among all rows of the same matrix.'''

    import sklearn

    normed_matrix = sklearn.preprocessing.normalize(matrix, axis=1, norm='l2')

    cosAnglesMatrix = normed_matrix @ normed_matrix.T

    return cosAnglesMatrix


def createPcaBiplot(normlzdSensMatrix, normlzdDefaultBiasesCol,
                    paramsNames,
                    whitelistedMetricsMask,
                    xColLabel, yColLabel,
                    colorCol,
                    colorColLabel,
                    pointLabels, pointLabelsHeader,
                    plotTitle,
                    xaxisTitle,
                    yaxisTitle,
                    showLegend, hoverMode,
                    plotWidth, plotHeight
                    ):
    import numpy as np

    from sklearn.decomposition import PCA

    # The code below is based on https://plotly.com/python/pca-visualization/

    df = pd.DataFrame(normlzdSensMatrix, columns=paramsNames)

    df['bias'] = normlzdDefaultBiasesCol

    # Remove the column-mean from each column,
    #   although the PCA code below seems to do this automatically.
    df = df - df.mean()

    #u, s, vh = np.linalg.svd(df.to_numpy(), full_matrices=False)
    #np.square(s) / np.sum(np.square(s))

    paramsNamesPlusBias = np.append(paramsNames, 'bias')

    #pca = PCA(n_components=2)
    pca = PCA()
    components = pca.fit_transform(df)

    # Let's relate the contents of pca to a standard SVD.
    #    See https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
    # If the SVD of the matrix is denoted matrix ~= u*diag(s)*vh, then
    # components = pca.fit_transform(df) = u@diag(s)
    # pca.components_ = vh (or first 2 rows of vh)
    # pca.singular_values_ = s
    # pca.explained_variance_ = eigenval = lambda = s**2/(numMetrics-1)
    # pca.explained_variance_ratio_ = RSquaredSvd = np.square(s)/np.sum(np.square(s))
    # loadings = vh @ s/sqrt(numMetrics-1)
    # loadings.T @ loadings is orthogonal
    # loadings @ loadings.T = covariance matrix
    # covMatrix2corrMatrix( loadings @ loadings.T ) = correlation matrix
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)

    PcaBiplotFig = \
        createScatterplot(xCol=components[:, 0], xColLabel=xColLabel,
                          yCol=components[:, 1], yColLabel=yColLabel,
                          colorCol=colorCol,
                          colorColLabel=colorColLabel,
                          colorScale='Rainbow',
                          plotBgColor='lightgrey',
                          pointLabels=pointLabels, pointLabelsHeader=pointLabelsHeader,
                          plotTitle=(plotTitle + np.array2string(pca.explained_variance_ratio_)),
                          xaxisTitle=xaxisTitle,
                          yaxisTitle=yaxisTitle,
                          showLegend=showLegend, hoverMode=hoverMode,
                          plotWidth=plotWidth, plotHeight=plotHeight)

    print("Whitelisted metrics = ", pointLabels[whitelistedMetricsMask])
    print("SV1, SV2 = ", components[whitelistedMetricsMask,0:2])

    ## Bold the labels of masked metrics
    #dfMasked = pd.DataFrame({
    #    xColLabel: components[whitelistedMetricsMask, 0],
    #    yColLabel: components[whitelistedMetricsMask, 1],
    #    colorColLabel: colorCol[whitelistedMetricsMask],
    #    pointLabelsHeader: pointLabels[whitelistedMetricsMask]
    #}, index=pointLabels[whitelistedMetricsMask])
    #normlzdColorCol = (colorCol[whitelistedMetricsMask] - np.min(colorCol)) / \
    #                  (np.max(colorCol) - np.min(colorCol))
    #j = 0
    #for i, row in dfMasked.iterrows():
    #    PcaBiplotFig.add_annotation(
    #        x=row[xColLabel],
    #        y=row[yColLabel],
    #        text=i,
    #        font=dict(color=pc.sample_colorscale('Rainbow', normlzdColorCol[j])[0],
    #                  weight="bold"),
    #        showarrow=False
    #    )
    #    j = j + 1

    # Display arrows to represent loading vectors
    scaleArrowLength = len(df.columns)
    for i, paramsName in enumerate(paramsNamesPlusBias):
        PcaBiplotFig.add_annotation(
            ax=0, ay=0,
            axref="x", ayref="y",
            x=scaleArrowLength * loadings[i, 0],
            y=scaleArrowLength * loadings[i, 1],
            showarrow=True,
            arrowsize=2,
            arrowhead=2,
            xanchor="right",
            yanchor="top"
        )
        # Label each loading vector with a parameter name
        PcaBiplotFig.add_annotation(
            x=scaleArrowLength * loadings[i, 0],
            y=scaleArrowLength * loadings[i, 1],
            ax=0, ay=0,
            xanchor="center",
            yanchor="bottom",
            text=paramsName,
            yshift=5,
            font=dict(size=14, weight="bold")
        )

    #from pca import pca
    #import pdb

    # https://plotly.com/python/pca-visualization/

    # reduce the data towards 2 PCs
    #model = pca(n_components=2, detect_outliers='ht2')

    # Augmented array with LHS and RHS
    #augMatrix = np.concatenate((normlzdSensMatrix, -defaultBiasesCol / np.abs(normMetricValsCol) ), axis=1)

    #paramsList = list(paramsNames)
    #paramsList.append('dbias')
    #augParamsNames = np.asarray(paramsList)

    # Fit transform
    #results = model.fit_transform(augMatrix, row_labels=metricsNames, col_labels=augParamsNames)
    #results = model.fit_transform(normlzdSensMatrix, row_labels=metricsNames, col_labels=paramsNames)

    #PC_test = model.transform(augMatrix)
    ##PC_test = model.transform(normlzdSensMatrix)
    #outliers, outliers_params = model.compute_outliers(PC=PC_test)
    #print("PCA outliers = ", outliers)

    # Plot explained variance
    #fig, ax = model.plot()

    # Scatter first 2 PCs
    #fig, ax = model.scatter()

    # Make biplot with the number of features
    #PcaBiplotFig, ax = model.biplot(n_feat=paramsNames.size+1)
    #fig, ax = model.biplot(n_feat=paramsNames.size)

    #pdb.set_trace()

    return PcaBiplotFig


def oldCode():
    # Calculate the fraction of the default-sim bias that remains after tuning.
    # This is unweighted and hence is not necessarily less than one.
    # defaultBiasesApprox = J*delta_p = ( fwd - def )
    # numerator = ( fwd - def ) + ( def - obs ) = ( fwd - obs )
    #    Bias = ( defaultBiasesApprox + defaultBiasesCol )
    # defaultBiasesCol = delta_b = ( default - obs ) = denominator
    #    BiasMagRatio = np.linalg.norm(Bias/np.abs(normMetricValsCol))**2 / \
    #                   np.linalg.norm(defaultBiasesCol/np.abs(normMetricValsCol))**2

    # Calculate the global-model bias relative to the default-sim bias.
    # This is unweighted and hence is not necessarily less than one.
    # defaultBiasesApprox = J*delta_p = ( fwd - def )
    # numerator = ( linSoln - def ) + ( def - obs ) = ( linSoln - obs )
    #linSolnBias = ( linSolnBiasesCol + defaultBiasesCol )
    # defaultBiasesCol = delta_b = ( default - obs ) = denominator
    #linSolnBiasMagRatio = np.linalg.norm(linSolnBias/np.abs(normMetricValsCol))**2 / \
    #                      np.linalg.norm(defaultBiasesCol/np.abs(normMetricValsCol))**2

    # weightedBiasLin = metricsWeights * ( lin - obs ) = numerator
    #weightedBiasLinSoln = metricsWeights * ( linSolnBiasesCol + defaultBiasesCol ) / np.abs(normMetricValsCol)
    #weightedBiasLinSolnMagRatio = np.linalg.norm(weightedBiasLinSoln)**2 / np.linalg.norm(normlzdMDeltaB)**2

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
    #    vhNormlzdColsFig.update_layout(hovermode="closest")
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
    #    uNormlzdColsFig.update_layout(hovermode="closest")
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
    #    uNormlzdWeightedColsFig.update_layout(hovermode="closest")
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
    #    uNormlzdBiasColsFig.update_layout(hovermode="closest")
    #    for idx, val in np.ndenumerate(sNormlzd):
    #        uNormlzdBiasColsFig.data[idx[0]].name = "{}".format(idx[0]+1) + ", " + "{:.2e}".format(val)
    #pdb.set_trace()

    # Plot the residual*sensitivity vs. bias*sensitivity
    # The goal is to separate out which regional are amenable to tuning
    normlzdBias = -defaultBiasesCol[:, 0] / np.abs(normMetricValsCol[:, 0])
    normlzdResid = (-defaultBiasesApproxElastic - defaultBiasesCol)[:, 0] / np.abs(normMetricValsCol[:, 0])
    df = pd.DataFrame(
        {'bias times sensitivity': normlzdBias * np.sum(normlzdSensMatrixPoly, axis=1),  # sum of row elements
         'residual times sensitivity': normlzdResid * np.sum(normlzdSensMatrixPoly, axis=1)  # sum of row elements
         }, index=metricsNames)
    residVsBiasScatter = px.scatter(df, x='bias times sensitivity', y=df.columns[1:],
                                    text=metricsNames,
                                    title="""Residual times sensitivity versus bias times sensitivity.""")
    residVsBiasOneOneLine = px.line(df, x="bias times sensitivity", y="bias times sensitivity")
    residVsBiasFig = go.Figure(data=residVsBiasScatter.data
                                    + residVsBiasOneOneLine.data)
    residVsBiasFig.update_yaxes(title="Residual times sensitivity")
    residVsBiasFig.update_xaxes(title="Bias times sensitivity")
    residVsBiasFig.update_layout(hovermode="closest")

    # Plot the relative biases versus sensitivity of each regional metric.
    #    More specifically, plot the maximum magnitude value of each row of the sensitivity matrix.
    #pdb.set_trace()
    relBiasNumerator = np.abs(-defaultBiasesApproxElastic - defaultBiasesCol)[:, 0] / np.abs(normMetricValsCol[:, 0])
    relBiasDenom = np.maximum(0.02, np.abs(-defaultBiasesCol[:, 0] / np.abs(normMetricValsCol[:, 0])))
    df = pd.DataFrame(
        {'Max abs normlzd sensitivity': np.max(np.abs(normlzdSensMatrixPoly), axis=1),  # max |row elements|
         'relative bias': relBiasNumerator / relBiasDenom
         }, index=metricsNames)
    relBiasesVsSensFig = px.scatter(df, x='Max abs normlzd sensitivity', y=df.columns[1:],
                                    text=metricsNames,
                                    title="""Ratio of abs(approx_bias)/abs(default_bias) versus sensitivity.""")
    relBiasesVsSensFig.update_yaxes(title="Relative biases")
    relBiasesVsSensFig.update_xaxes(title="Max abs normlzd sensitivity")
    relBiasesVsSensFig.update_layout(hovermode="closest")

    # Plot the relative biases versus sensitivity of each regional metric.
    #    More specifically, plot the maximum magnitude value of each row of the sensitivity matrix.
    absBiasTuned = np.abs(-defaultBiasesApproxElastic - defaultBiasesCol)[:, 0] / np.abs(normMetricValsCol[:, 0])
    absBiasDefault = np.abs(-defaultBiasesCol[:, 0] / np.abs(normMetricValsCol[:, 0]))
    df = pd.DataFrame(
        {'Max abs normlzd sensitivity': np.max(np.abs(normlzdSensMatrixPoly), axis=1),  # max |row elements|
         'bias difference': absBiasDefault - absBiasTuned
         }, index=metricsNames)
    diffBiasesVsSensFig = px.scatter(df, x='Max abs normlzd sensitivity', y=df.columns[1:],
                                     text=metricsNames,
                                     title="""abs(default_bias) - abs(approx_bias) versus sensitivity.""")
    diffBiasesVsSensFig.update_yaxes(title="Bias difference")
    diffBiasesVsSensFig.update_xaxes(title="Max abs normlzd sensitivity")
    diffBiasesVsSensFig.update_layout(hovermode="closest")

    # Plot the parameter values recommended by SVD.
    # Multiply in the user-designated scale factors before plotting.
    paramsFig = go.Figure()
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsLowValsPCBound[:, 0] * paramsScales,
                                   name=r'$paramsSolnPC - \sigma$',
                                   line=dict(color='white', width=0), mode='lines', showlegend=False))
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsHiValsPCBound[:, 0] * paramsScales, fill='tonexty',
                                   name='Default Parameter Values +- sigma', mode='none',
                                   fillcolor='rgba(253,253,150,1.0)'))
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=defaultParamValsOrigRow[0, :] * paramsScales,
                                   name='Default Parameter Values', line=dict(color='black', width=6)))
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnLin[:, 0] * paramsScales,
                                   name='Linear regression, |dp|=' + '{:.2e}'.format(
                                       np.linalg.norm(dnormlzdParamsSolnLin))))
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnNonlin[:, 0] * paramsScales,
                                   name='paramsSolnNonlin, |dpNonlin|='
                                        + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSolnNonlin))))
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnElastic[:, 0] * paramsScales,
                                   name='Lasso regression, |dpLasso|='
                                        + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSolnElastic)),
                                   line=dict(color='red', width=2)))
    #paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnPCBound[:,0]*paramsScales,
    #                                name='paramsSolnPCBound, |dpBound|='
    #                               + '{:.2e}'.format(0.0) ))
    paramsFig.update_yaxes(title="User-scaled parameter value")
    paramsFig.update_xaxes(title="Parameter Name")
    paramsFig.update_layout(hovermode="closest")
    paramsFig.update_layout(width=1000, height=500)

    #    # Plot the biases of the default simulation and the SVD approximation of that
    #    biasesMatrix = np.dstack((-defaultBiasesCol,
    #                          #defaultBiasesApprox,
    #                          #defaultBiasesApproxPC,
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
    #    biasesFig.update_layout(hovermode="closest")
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
    #

    df = pd.DataFrame(-1 * normlzdSensParamsMatrixOrdered,
                      index=metricsNamesOrdered,
                      columns=paramsNames)
    biasContrOrderFig = px.bar(df, x=df.index, y=df.columns,
                               title="""Linear contributions to actual removal of regional biases""")
    biasContrOrderFig.update_yaxes(title="-(Def-Sim) / abs(obs metric value)")
    biasContrOrderFig.update_xaxes(title="Regional metric")
    biasContrOrderFig.update_layout(hovermode="closest")
    biasContrOrderFig.update_layout(showlegend=True)
    #biasContrOrderFig.update_traces(mode='markers', line_color='black')  # Plot default biases as black dots
    biasContrOrderFig.update_yaxes(visible=True, zeroline=True, zerolinewidth=1, zerolinecolor='gray')  # Plot x axis
    biasContrOrderFig.update_layout(width=800, height=500)
    # Now plot an arrow for each region that points from default-run bias to new bias after tuning
    xArrow = np.arange(len(metricsNamesOrdered))  # x-coordinate of arrows
    yArrow = -defaultBiasesCol[metricsSensOrder, 0] / np.abs(normMetricValsCol[metricsSensOrder, 0])
    gap = 0.1  # horizontal spacing between arrows
    #Plot arrows showing the tuner's nonlinear predicted bias removal
    for i, item in enumerate(metricsNamesOrdered):
        biasContrOrderFig.add_annotation(
            x=xArrow[i] - gap,  # ith arrow's head
            # ith arrow's head:
            y=(-defaultBiasesApproxNonlin - defaultBiasesCol)[metricsSensOrder[i], 0] / np.abs(
                normMetricValsCol[metricsSensOrder[i], 0]),
            #y= (-defaultBiasesApproxNonlinNoCurv-defaultBiasesCol)[metricsSensOrder[i],0]/np.abs(normMetricValsCol[metricsSensOrder[i],0]),
            ax=xArrow[i] - gap,  # ith arrow's tail
            ay=yArrow[i],  # ith arrow's tail
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
    biasContrOrderFig.add_annotation(text='tuner prediction of bias removal',
                                     font=dict(color='blue'),
                                     align='left', xref='paper', yref='paper', x=0.05, y=0.9, showarrow=False)
    biasContrOrderFig.add_annotation(text='realized E3SM bias removal',
                                     font=dict(color='red'),  #'rgba(255,0,0,0.0)'),
                                     align='left', xref='paper', yref='paper', x=0.05, y=0.8, showarrow=False)

    return
