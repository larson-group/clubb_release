# -*- coding: utf-8 -*-

# Run this app with `python3 sens_matrix_dashboard.py` and
# view the plots at http://127.0.0.1:8050/ in your web browser.
# (To open a web browser on a larson-group computer,
# login to malan with `ssh -X` and then type `firefox &`.)

def main():

    import dash
    import dash_core_components as dcc
    import dash_html_components as html
    import plotly.express as px
    import plotly.graph_objects as go
    import pandas as pd

    import numpy as np
    import pdb
    import sklearn
    from plotly.figure_factory import create_quiver
    from itertools import chain

    from analyze_sensitivity_matrix import \
            analyzeSensMatrix, setupObsCol, setupDefaultMetricValsCol, \
            findOutliers, findParamsUsingElastic
    from test_analyzeSensMatrix import write_test_netcdf_files
    from set_up_dashboard_inputs import setUpInputs


    # The user should input all tuning data into file set_up_dashboard_inputs.py
    metricsNames, metricsWeights, \
    paramsNames, paramsScales, \
    transformedParamsNames, \
    sensNcFilenames, sensNcFilenamesExt, \
    defaultNcFilename, linSolnNcFilename, \
    obsMetricValsDict, obsMetricValsCol, \
    reglrCoef \
    = \
        setUpInputs()


    # Estimate non-linearity of the global model to perturbations in parameter values.
    # To do so, calculate radius of curvature of the three points from the default simulation
    #   and the two sensitivity simulations.
    #calcNormlzdRadiusCurv(metricsNames, paramsNames, transformedParamsNames, paramsScales,
    #                      metricsWeights, obsMetricValsCol,
    #                      sensNcFilenames, sensNcFilenamesExt, defaultNcFilename)

    # Calculate changes in parameter values needed to match metrics.
    defaultMetricValsCol, defaultBiasesCol, \
    defaultBiasesApprox, defaultBiasesApproxLowVals, defaultBiasesApproxHiVals, \
    defaultBiasesApproxPC, defaultBiasesApproxLowValsPC, defaultBiasesApproxHiValsPC, \
    normlzdWeightedDefaultBiasesApprox, normlzdWeightedDefaultBiasesApproxPC, \
    defaultBiasesOrigApprox, defaultBiasesOrigApproxPC, \
    sensMatrixOrig, sensMatrix, normlzdSensMatrix, \
    normlzdWeightedSensMatrix, biasNormlzdSensMatrix, svdInvrsNormlzdWeighted, \
    vhNormlzd, uNormlzd, sNormlzd, \
    vhNormlzdWeighted, uNormlzdWeighted, sNormlzdWeighted, \
    magParamValsRow, \
    defaultParamValsOrigRow, dparamsSoln, dnormlzdParamsSoln, \
    dparamsSolnPC, dnormlzdParamsSolnPC, \
    paramsSoln, paramsLowVals, paramsHiVals, \
    paramsSolnPC, paramsLowValsPC, paramsHiValsPC = \
          analyzeSensMatrix(metricsNames, paramsNames, transformedParamsNames,
                            metricsWeights,
                            sensNcFilenames, sensNcFilenamesExt, defaultNcFilename,
                             obsMetricValsDict)


    paramsLowValsPCBound, paramsHiValsPCBound = \
        calcParamsBounds(metricsNames, paramsNames, transformedParamsNames,
                     metricsWeights, obsMetricValsCol,
                     magParamValsRow,
                     sensNcFilenames, sensNcFilenamesExt, defaultNcFilename)


    #print("dnormlzdParamsSoln.T=", dnormlzdParamsSoln.T)
    #print("normlzdSensMatrix=", normlzdSensMatrix)
    #print("normlzdSensMatrix@dnormlzdParamsSoln=", normlzdSensMatrix @ dnormlzdParamsSoln)

    normlzdCurvMatrix = \
        constructNormlzdCurvMatrix(metricsNames, paramsNames, transformedParamsNames, \
                                   metricsWeights, obsMetricValsCol, magParamValsRow, \
                                   sensNcFilenames, sensNcFilenamesExt, defaultNcFilename)

    defaultBiasesApproxNonlin, \
    dnormlzdParamsSolnNonlin, paramsSolnNonlin, \
    defaultBiasesApproxNonlin2x = \
        solveUsingNonlin(metricsNames, paramsNames, transformedParamsNames, \
                         metricsWeights, obsMetricValsCol, magParamValsRow, \
                         sensNcFilenames, sensNcFilenamesExt, defaultNcFilename, \
                         defaultParamValsOrigRow, \
                         dnormlzdParamsSoln, normlzdSensMatrix, defaultBiasesCol, \
                         normlzdCurvMatrix, \
                         reglrCoef)

    # Create scatterplot to look at outliers
    #createPcaBiplot(normlzdSensMatrix, defaultBiasesCol, obsMetricValsCol, metricsNames, paramsNames)

    ## Find outliers by use of the ransac algorithm
    #outlier_mask, defaultBiasesApproxRansac, normlzdWeightedDefaultBiasesApproxRansac, \
    #dnormlzdParamsSolnRansac, paramsSolnRansac = \
    #    findOutliers(normlzdSensMatrix, normlzdWeightedSensMatrix, \
    #                 defaultBiasesCol, obsMetricValsCol, magParamValsRow, defaultParamValsOrigRow)
    #print( "ransac_outliers = ", metricsNames[outlier_mask] )
    #print( "ransac_inliers = ", metricsNames[~outlier_mask] )
    ##pdb.set_trace()

    # Find best-fit params by use of the Elastic Net algorithm
    defaultBiasesApproxElastic, defaultBiasesApproxElasticNonlin, \
    dnormlzdParamsSolnElastic, paramsSolnElastic = \
        findParamsUsingElastic(normlzdSensMatrix, normlzdWeightedSensMatrix, \
                     defaultBiasesCol, obsMetricValsCol, metricsWeights, \
                     magParamValsRow, defaultParamValsOrigRow, \
                     normlzdCurvMatrix)

    defaultBiasesApproxElasticCheck = ( normlzdWeightedSensMatrix @ dnormlzdParamsSolnElastic ) \
                            * np.reciprocal(metricsWeights) * np.abs(obsMetricValsCol)

    print("defaultBiasesApproxElastic = ", defaultBiasesApproxElastic)
    print("defaultBiasesApproxElasticCheck = ", defaultBiasesApproxElasticCheck)

    #pdb.set_trace()

    # Set up a column vector of metric values from the default simulation
    defaultMetricValsCol = setupDefaultMetricValsCol(metricsNames, defaultNcFilename)

    # Set up a column vector of metric values from the global simulation based on optimized
    #     parameter values.
    linSolnMetricValsCol = setupDefaultMetricValsCol(metricsNames, linSolnNcFilename)

    # Store biases in default simulation, ( global_model - default )
    linSolnBiasesCol = np.subtract(linSolnMetricValsCol, defaultMetricValsCol)

    # Calculate the fraction of the default-sim bias that remains after tuning.
    # This is unweighted and hence is not necessarily less than one.
    # defaultBiasesApprox = J*delta_p = ( fwd - def )
    # numerator = ( fwd - def ) + ( def - obs ) = ( fwd - obs )
    Bias = ( defaultBiasesApprox + defaultBiasesCol )
    # defaultBiasesCol = delta_b = ( default - obs ) = denominator
    BiasMagRatio = np.linalg.norm(Bias/np.abs(obsMetricValsCol))**2 / \
                   np.linalg.norm(defaultBiasesCol/np.abs(obsMetricValsCol))**2

    # Calculate the fraction of the default-sim bias that remains after tuning,
    #    but using a truncated PC observation.
    # This is unweighted and hence is not necessarily less than one.
    # defaultBiasesApproxPC = J*delta_p = ( fwd - def )
    # numerator = ( fwd - def ) + ( def - obs ) = ( fwd - obs )
    BiasPC = ( defaultBiasesApproxPC + defaultBiasesCol )
    # defaultBiasesCol = delta_b = ( default - obs ) = denominator
    BiasPCMagRatio = np.linalg.norm(BiasPC/np.abs(obsMetricValsCol))**2 / \
                     np.linalg.norm(defaultBiasesCol/np.abs(obsMetricValsCol))**2

    ## Calculate the fraction of the default-sim bias that remains after tuning,
    ##    but using a truncated PC observation.
    ## This is unweighted and hence is not necessarily less than one.
    ## defaultBiasesApproxRansac = J*delta_p = ( fwd - def )
    ## numerator = ( fwd - def ) + ( def - obs ) = ( fwd - obs )
    #BiasRansac = ( defaultBiasesApproxRansac + defaultBiasesCol )
    ## defaultBiasesCol = delta_b = ( default - obs ) = denominator
    #BiasRansacMagRatio = np.linalg.norm(BiasRansac/np.abs(obsMetricValsCol))**2 / \
    #                 np.linalg.norm(defaultBiasesCol/np.abs(obsMetricValsCol))**2

    # Calculate the fraction of the default-sim bias that remains after tuning,
    #    but using a truncated PC observation.
    # This is unweighted and hence is not necessarily less than one.
    # defaultBiasesApproxElastic = J*delta_p = ( fwd - def )
    # numerator = ( fwd - def ) + ( def - obs ) = ( fwd - obs )
    BiasElastic = ( defaultBiasesApproxElastic + defaultBiasesCol )
    # defaultBiasesCol = delta_b = ( default - obs ) = denominator
    BiasElasticMagRatio = np.linalg.norm(BiasElastic/np.abs(obsMetricValsCol))**2 / \
                     np.linalg.norm(defaultBiasesCol/np.abs(obsMetricValsCol))**2

    # Calculate the global-model bias relative to the default-sim bias.
    # This is unweighted and hence is not necessarily less than one.
    # defaultBiasesApprox = J*delta_p = ( fwd - def )
    # numerator = ( linSoln - def ) + ( def - obs ) = ( linSoln - obs )
    linSolnBias = ( linSolnBiasesCol + defaultBiasesCol )
    # defaultBiasesCol = delta_b = ( default - obs ) = denominator
    linSolnBiasMagRatio = np.linalg.norm(linSolnBias/np.abs(obsMetricValsCol))**2 / \
                          np.linalg.norm(defaultBiasesCol/np.abs(obsMetricValsCol))**2

    # Calculate the fraction of bias removed by the non-PC soln, but normalized and weighted,
    # like the equations that the SVD actually solves, so that according to theory,
    # the value should be < 1.
    # But I'm not sure if it will be < 1 if the parameters are transformed to log space.
    normlzdMDeltaB = metricsWeights * defaultBiasesCol / np.abs(obsMetricValsCol) # right-hand side
    weightedBiasNumer = normlzdWeightedDefaultBiasesApprox + normlzdMDeltaB
    weightedBiasDenom = normlzdMDeltaB
    weightedBiasMagRatio = np.linalg.norm(weightedBiasNumer)**2 / np.linalg.norm(weightedBiasDenom)**2

    # Calculate the fraction of bias removed by PC soln, but normalized and weighted,
    # like the equations that the SVD actually solves, so that according to theory,
    # the value should be < 1.
    # But I'm not sure if it will be < 1 if the parameters are transformed to log space.
    weightedBiasPCNumer = normlzdWeightedDefaultBiasesApproxPC + normlzdMDeltaB
    weightedBiasPCDenom = normlzdMDeltaB
    weightedBiasPCMagRatio = np.linalg.norm(weightedBiasPCNumer)**2 / np.linalg.norm(weightedBiasPCDenom)**2

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
    weightedBiasLinSoln = metricsWeights * ( linSolnBiasesCol + defaultBiasesCol ) / np.abs(obsMetricValsCol)
    weightedBiasLinSolnMagRatio = np.linalg.norm(weightedBiasLinSoln)**2 / np.linalg.norm(normlzdMDeltaB)**2


    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

    sensMatrixDashboard = dash.Dash(__name__, external_stylesheets=external_stylesheets)

    # Plot the biases of the default simulation and the SVD approximation of that
    biasesMatrix = np.dstack((-defaultBiasesCol,
                          defaultBiasesApprox,
                          defaultBiasesApproxPC,
                          #defaultBiasesApproxElastic,
                          defaultBiasesApproxElasticNonlin,
                          defaultBiasesApproxNonlin,
                          linSolnBiasesCol
                         )).squeeze()
    fracBiasesMatrix = np.diagflat(np.reciprocal(np.abs(obsMetricValsCol))) @ biasesMatrix
    df = pd.DataFrame(fracBiasesMatrix,
                    index=metricsNames,
                  columns= ['fracDefBias',
                            'fracDefBiasesApprox',
                            'fracDefBiasesApproxPC',
                            'fracDefBiasesApproxElasticNonlin',
                            'fracDefBiasesApproxNonlin',
                            'fracLinSolnBiasesCol'
                           ])
    biasesFig = px.line(df, x=df.index, y=df.columns,
              title = """Fractional biases of default simulation and approximations thereof.<br>
                    Plotted quantities have the structure -(def-obs), -(def-fwd), -(def-lin)""")
    biasesFig.update_yaxes(title="-(Def-Sim) / abs(obs metric value)")
    biasesFig.update_xaxes(title="Metric and region")
    biasesFig.layout.legend.title = "Default or which approximation"
    biasesFig.update_layout(hovermode="x")
    biasesFig.data[1].name = "fracDefBiasesApprox, " \
                         + "{:.2f}".format(weightedBiasMagRatio) \
                        + ", {:.2f}".format(BiasMagRatio)
    biasesFig.data[2].name = "fracDefBiasesApproxPC, " \
                         + "{:.2f}".format(weightedBiasPCMagRatio) \
                         + ", {:.2f}".format(BiasPCMagRatio)
    biasesFig.data[3].name = "fracDefBiasesApproxElastic, " \
                         + "{:.2f}".format(weightedBiasElasticMagRatio) \
                         + ", {:.2f}".format(BiasElasticMagRatio)
    biasesFig.data[4].name = "fracDefBiasesApproxNonlin, " \
                         + "{:.2f}".format(-99) \
                         + ", {:.2f}".format(-99)
    biasesFig.data[5].name = "fracLinSolnBiasesCol, " \
                          + "{:.2f}".format(weightedBiasLinSolnMagRatio) \
                          + ", {:.2f}".format(linSolnBiasMagRatio)


    # Create plot showing how well the regional biases are actually removed
    metricsSens = np.linalg.norm(normlzdWeightedSensMatrix, axis=1) # measure of sensitivity of each metric
    # metricsSensOrdered = (rankdata(metricsSens) - 1).astype(int)  # this ordering doesn't work as an index
    metricsSensOrdered = metricsSens.argsort()
    metricsNamesOrdered = metricsNames[metricsSensOrdered]  # list of metrics names, ordered from least to most sensitive
    # Plot a black dot for each default-run bias
    biasesOrderMatrix = np.dstack(( -defaultBiasesCol[metricsSensOrdered] )).squeeze()
    fracBiasesOrderMatrix = np.diagflat(np.reciprocal(np.abs(obsMetricValsCol[metricsSensOrdered]))) @ biasesOrderMatrix
    df = pd.DataFrame(fracBiasesOrderMatrix,
                      index=metricsNamesOrdered,
                      columns= ['fracDefBias'])
    biasesOrderFig = px.line(df, x=df.index, y=df.columns,
              title = """Predicted and actual removal of regional biases""")
    biasesOrderFig.update_yaxes(title="-(Def-Sim) / abs(obs metric value)")
    biasesOrderFig.update_xaxes(title="Metric and region")
    biasesOrderFig.update_layout(hovermode="x")
    biasesOrderFig.update_layout(showlegend=False)
    biasesOrderFig.update_traces(mode='markers', line_color='black')  # Plot default biases as black dots
    biasesOrderFig.update_yaxes(visible=True,zeroline=True,zerolinewidth=1,zerolinecolor='gray') # Plot x axis
    biasesOrderFig.update_layout( width=500, height=500  )
    # Now plot an arrow for each region that points from default-run bias to new bias after tuning
    xArrow = np.arange(len(metricsNamesOrdered)) # x-coordinate of arrows
    yArrow = -defaultBiasesCol[metricsSensOrdered,0]/np.abs(obsMetricValsCol[metricsSensOrdered,0])
    gap = 0.1  # horizontal spacing between arrows
    # Plot arrows showing the tuner's nonlinear predicted bias removal
    for i, item in enumerate(metricsNamesOrdered):
        biasesOrderFig.add_annotation(
        x=  xArrow[i] - gap,  # ith arrow's head
        # ith arrow's length:
        y= (-defaultBiasesApproxNonlin-defaultBiasesCol)[metricsSensOrdered[i],0]/np.abs(obsMetricValsCol[metricsSensOrdered[i],0]),
        #y= (-defaultBiasesApproxNonlin2x-defaultBiasesCol)[metricsSensOrdered[i],0]/np.abs(obsMetricValsCol[metricsSensOrdered[i],0]),
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
    biasesOrderFig.add_annotation(text='tuner prediction of bias removal',
                                  font=dict(color='blue'),
                                  align='left', xref='paper', yref='paper', x=0.05, y=0.9, showarrow=False)
    biasesOrderFig.add_annotation(text='realized E3SM bias removal',
                                  font=dict(color='red'), #'rgba(255,0,0,0.0)'),
                                  align='left', xref='paper', yref='paper', x=0.05, y=0.8, showarrow=False)
    # Plot arrows showing the bias removal of E3SM's solution
    for i, item in enumerate(metricsNamesOrdered):
        biasesOrderFig.add_annotation(
        x=  xArrow[i]+gap,  # ith arrow's head
        # ith arrow's length:
        y= (-linSolnBiasesCol-defaultBiasesCol)[metricsSensOrdered[i],0]/np.abs(obsMetricValsCol[metricsSensOrdered[i],0]),
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
        arrowcolor='red' #,
        #opacity=0.0
	    )
    #pdb.set_trace()

    # Plot a scatterplot of default-simulation bias and SVD approximation of that bias.
    # Each column tells us how all metrics vary with a single parameter.
    biasSensMatrix = np.concatenate((-defaultBiasesCol/np.abs(obsMetricValsCol),
                                    (-defaultBiasesApproxElastic-defaultBiasesCol)/np.abs(obsMetricValsCol)), axis=1)
                                     #defaultBiasesApproxElastic/np.abs(obsMetricValsCol)), axis=1)
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
    #normlzdSensMatrixColsFig.layout.legend.title = "Parameter"
    #pdb.set_trace()


    # Plot a scatterplot of minimum parameter perturbation vs. fractional default bias approximation
    # Calculate lower bound on normalized parameter perturbations
    #normlzdDefaultBiasesCol = ( metricsWeights * (-defaultBiasesCol) /
    normlzdDefaultBiasesCol = ( (-defaultBiasesCol) /
                                np.abs(obsMetricValsCol) )
    sensMatrixRowMag = np.linalg.norm(normlzdWeightedSensMatrix, axis=1)
    dpMin = np.abs(normlzdDefaultBiasesCol) / np.atleast_2d(sensMatrixRowMag).T
    u_dot_b = np.atleast_2d(sensMatrixRowMag).T * normlzdDefaultBiasesCol
    #dpMinMatrix = np.dstack((np.reciprocal(dpMin),
    dpMinMatrix = np.dstack((np.abs(u_dot_b),
    #dpMinMatrix = np.dstack((np.atleast_2d(sensMatrixRowMag).T,
                          np.abs(defaultBiasesApproxElastic)/np.abs(obsMetricValsCol)
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
    df = pd.DataFrame(np.max(np.abs(normlzdSensMatrix), axis=1), # max of absolute val of each row
                  index=metricsNames,
                  columns=['Max abs normlzd sensitivity'])
    maxSensMetricsFig = px.line(df, x=df.index, y=df.columns,
              title = """Maximum normalized sensitivity of each metric with respect to parameters.<br>
                       (Low sensitivity means that the metric is unbudgeable by these parameters.)""" )
    maxSensMetricsFig.update_yaxes(title="Max |sens row|")
    maxSensMetricsFig.update_xaxes(title="Metric and region")
    maxSensMetricsFig.update_layout(hovermode="x")
    maxSensMetricsFig.update_traces(mode='lines+markers')


    # Plot the biases versus sensitivity of each regional metric.
    #    More specifically, plot the maximum magnitude value of each row of the sensitivity matrix.
    #df = pd.DataFrame({'Max abs normlzd sensitivity': np.max(np.abs(normlzdSensMatrix), axis=1), # max |row elements|
    #df = pd.DataFrame({'Max abs normlzd sensitivity': np.sum(normlzdWeightedSensMatrix, axis=1), # sum of row elements
    df = pd.DataFrame({'Max abs normlzd sensitivity': np.linalg.norm(normlzdWeightedSensMatrix, axis=1), # sum of row elements
    #df = pd.DataFrame({'Max abs normlzd sensitivity':
    #                    -defaultBiasesCol[:,0]/np.abs(obsMetricValsCol[:,0])*np.linalg.norm(normlzdWeightedSensMatrix, axis=1), # sum of row elements
                       'default tuning': -defaultBiasesCol[:,0]/np.abs(obsMetricValsCol[:,0]),
                       'revised tuning': (-defaultBiasesApproxElastic-defaultBiasesCol)[:,0]/np.abs(obsMetricValsCol[:,0])
                      }, index=metricsNames )
    biasesVsSensFig = px.scatter(df, x='Max abs normlzd sensitivity', y=df.columns[1:],
                                 text=metricsNames, title = """Regional biases with default and revised tuning, versus sensitivity.""" )
    biasesVsSensFig.update_yaxes(title="Regional biases")
    biasesVsSensFig.update_xaxes(title="Sensitivity of regional metrics to parameter changes")
    biasesVsSensFig.update_layout(hovermode="x")

    # Compute length of arrows between default and tuned biases
    #metricsNamesPadded = ",,".join(metricsNames).split(",")
    #metricsNamesPadded = ",,".join(metricsNamesPadded).split(",")
    #metricsNamesPadded = np.append(metricsNamesPadded, ["", "", ""], axis=0)
    xArrow = np.linalg.norm(normlzdWeightedSensMatrix, axis=1)
    yArrow = -defaultBiasesCol[:,0]/np.abs(obsMetricValsCol[:,0])
    #uArrow = np.zeros_like(xArrow)
    #vArrow = (-defaultBiasesApproxElasticNonlin)[:,0]/np.abs(obsMetricValsCol[:,0])
    #arrowFig = create_quiver(xArrow, yArrow, uArrow, vArrow,
    #                         scale=1,text=metricsNamesPadded)
    #arrowFig.update_yaxes(title="Regional biases")
    #arrowFig.update_traces(mode='lines+text')  # make labels appear in plot, not just hovermode


    # Plot biases vs. sensitivity, but with arrows indicating the degree of bias reduction
    df = pd.DataFrame({'Max abs normlzd sensitivity': np.linalg.norm(normlzdWeightedSensMatrix, axis=1), # sum of row elements

                       'default tuning': -defaultBiasesCol[:,0]/np.abs(obsMetricValsCol[:,0]),
                       'revised tuning': (-defaultBiasesApproxNonlin-defaultBiasesCol)[:,0]/np.abs(obsMetricValsCol[:,0])
                      }, index=metricsNames )
    biasesVsSensArrowFig = px.scatter(df, x='Max abs normlzd sensitivity', y=df.columns[1:2],
                                 text=metricsNames,
                                 title = """Regional biases with default and nonlin tuning versus sensitivity, with arrows.""" )
    biasesVsSensArrowFig.update_traces(textposition="middle right")
    for i, item in enumerate(metricsNames):
        biasesVsSensArrowFig.add_annotation(
        x=  xArrow[i],  # ith arrow's head
        y= (-defaultBiasesApproxNonlin-defaultBiasesCol)[i,0]/np.abs(obsMetricValsCol[i,0]),  # ith arrow's length
        #y= (-defaultBiasesApproxNonlin2x-defaultBiasesCol)[i,0]/np.abs(obsMetricValsCol[i,0]),  # ith arrow's length
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
    normlzdBias = -defaultBiasesCol[:,0]/np.abs(obsMetricValsCol[:,0])
    normlzdResid = (-defaultBiasesApproxElastic-defaultBiasesCol)[:,0]/np.abs(obsMetricValsCol[:,0])
    df = pd.DataFrame({'bias times sensitivity': normlzdBias*np.sum(normlzdSensMatrix, axis=1), # sum of row elements
                       'residual times sensitivity': normlzdResid*np.sum(normlzdSensMatrix, axis=1) # sum of row elements
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
    relBiasNumerator = np.abs(-defaultBiasesApproxElastic-defaultBiasesCol)[:,0]/np.abs(obsMetricValsCol[:,0])
    relBiasDenom = np.maximum(0.02, np.abs(-defaultBiasesCol[:,0]/np.abs(obsMetricValsCol[:,0])) )
    df = pd.DataFrame({'Max abs normlzd sensitivity': np.max(np.abs(normlzdSensMatrix), axis=1), # max |row elements|
                       'relative bias': relBiasNumerator / relBiasDenom
                      }, index=metricsNames )
    relBiasesVsSensFig = px.scatter(df, x='Max abs normlzd sensitivity', y=df.columns[1:],
                                 text=metricsNames, title = """Ratio of abs(approx_bias)/abs(default_bias) versus sensitivity.""" )
    relBiasesVsSensFig.update_yaxes(title="Relative biases")
    relBiasesVsSensFig.update_xaxes(title="Max abs normlzd sensitivity")
    relBiasesVsSensFig.update_layout(hovermode="x")

    # Plot the relative biases versus sensitivity of each regional metric.
    #    More specifically, plot the maximum magnitude value of each row of the sensitivity matrix.
    absBiasTuned = np.abs(-defaultBiasesApproxElastic-defaultBiasesCol)[:,0]/np.abs(obsMetricValsCol[:,0])
    absBiasDefault = np.abs(-defaultBiasesCol[:,0]/np.abs(obsMetricValsCol[:,0]))
    df = pd.DataFrame({'Max abs normlzd sensitivity': np.max(np.abs(normlzdSensMatrix), axis=1), # max |row elements|
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
    df = pd.DataFrame( np.hstack( (-defaultBiasesCol/np.abs(obsMetricValsCol),normlzdSensMatrix) ),
                  index=metricsNames,
                  columns=np.append('Norm bias', paramsNames) )
    normlzdSensMatrixColsFig = px.line(df, x=df.index, y=df.columns,
              title =  """Columns of normalized, unweighted sensitivity matrix (plus the bias
                                       vector).<br>
                       Each column (line) shows how sensitive the metrics are to a change in a single parameter value.<br>
                       (A positive value means that an increase in parameter value brings the default simulation closer to obs.)""" )
    normlzdSensMatrixColsFig.update_yaxes(title="Norml sens, (|param|/|obsmetric|) * dmetric/dparam")
    normlzdSensMatrixColsFig.update_xaxes(title="Metric and region")
    normlzdSensMatrixColsFig.layout.legend.title = "Parameter"
    normlzdSensMatrixColsFig.update_layout(hovermode="x")
    #pdb.set_trace()

    # Plot each row of normalized sensitivity matrix as a separate line.
    # Each row tells us how a single metric varies with all parameters.
    df = pd.DataFrame(np.transpose(normlzdSensMatrix),
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
    rightSingVectorNums = (np.arange(paramsNames.shape[0])+1).astype(str)
    df = pd.DataFrame(np.transpose(vhNormlzd),
                  index=paramsNames,
                  columns=rightSingVectorNums)
    vhNormlzdColsFig = px.line(df, x=df.index, y=df.columns,
              title = """Columns of normalized, unweighted right-singular vector matrix, V.<br>
                        Each column (line) is a vector of parameter values associated with a singular value.<br>""" )
    vhNormlzdColsFig.update_yaxes(title="Right-singular vector")
    vhNormlzdColsFig.update_xaxes(title="Parameter")
    vhNormlzdColsFig.layout.legend.title = "Singular vector"
    vhNormlzdColsFig.update_layout(hovermode="x")
    for idx, val in np.ndenumerate(sNormlzd):
        vhNormlzdColsFig.data[idx[0]].name = "{}".format(idx[0]+1) + ", " + "{:.2e}".format(val)

    # Plot each column of normalized, unweighted left-singular vector matrix, U.
    df = pd.DataFrame(uNormlzd,
                  index=metricsNames,
                  columns=rightSingVectorNums)
    uNormlzdColsFig = px.line(df, x=df.index, y=df.columns,
              title = """Columns of normalized, unweighted left-singular vector matrix, U.<br>
                       Each column (line) is a vector of metric values associated with a singular value.<br>""" )
    uNormlzdColsFig.update_yaxes(title="Left-singular vector")
    uNormlzdColsFig.update_xaxes(title="Parameter")
    uNormlzdColsFig.layout.legend.title = "Singular vector"
    uNormlzdColsFig.update_layout(hovermode="x")
    for idx, val in np.ndenumerate(sNormlzd):
        uNormlzdColsFig.data[idx[0]].name = "{}".format(idx[0]+1) + ", " + "{:.2e}".format(val)

    # Plot each column of normalized, weighted left-singular vector matrix, U.
    df = pd.DataFrame(uNormlzdWeighted,
                  index=metricsNames,
                  columns=rightSingVectorNums)
    uNormlzdWeightedColsFig = px.line(df, x=df.index, y=df.columns,
              title = """Columns of normalized, weighted left-singular vector matrix, U.<br>
                       Each column (line) is a vector of metric values associated with a singular value.<br>""" )
    uNormlzdWeightedColsFig.update_yaxes(title="Left-singular vector")
    uNormlzdWeightedColsFig.update_xaxes(title="Parameter")
    uNormlzdWeightedColsFig.layout.legend.title = "Singular vector"
    uNormlzdWeightedColsFig.update_layout(hovermode="x")
    for idx, val in np.ndenumerate(sNormlzdWeighted):
        uNormlzdWeightedColsFig.data[idx[0]].name = "{}".format(idx[0]+1) + ", " + "{:.2e}".format(val)

    # Plot each column of left-singular vector matrix, U, multiplied by biases.
    # Plot each column of left-singular vector matrix, U, multiplied by biases.
    df = pd.DataFrame(uNormlzd*defaultBiasesCol/np.abs(obsMetricValsCol),

                   index=metricsNames,
                  columns=rightSingVectorNums)
    uNormlzdBiasColsFig = px.line(df, x=df.index, y=df.columns,
              title = """Columns of normalized, unweighted left-singular vector matrix, U, dotted with delta_b (defaultBiasesCol).<br>
                       Each column (line) is a vector of metric values associated with a singular value.<br>""" )
    uNormlzdBiasColsFig.update_yaxes(title="uNormlzd dot delta_b")
    uNormlzdBiasColsFig.update_xaxes(title="Parameter")
    uNormlzdBiasColsFig.layout.legend.title = "Singular vector"
    uNormlzdBiasColsFig.update_layout(hovermode="x")
    for idx, val in np.ndenumerate(sNormlzd):
        uNormlzdBiasColsFig.data[idx[0]].name = "{}".format(idx[0]+1) + ", " + "{:.2e}".format(val)
    #pdb.set_trace()

    # Plot box and whiskers plot of optimal parameter values.
    # Multiply in the user-designated scale factors before plotting.
    df = pd.DataFrame( np.hstack( defaultParamValsOrigRow[0,:]*paramsScales ),
                  index=paramsNames, columns=["Default plus error bars"] )
    df["err_minus"] = ( defaultParamValsOrigRow[0,:] -  paramsLowValsPCBound[:,0] ) * paramsScales
    df["err_plus"]  = ( paramsHiValsPCBound[:,0] - defaultParamValsOrigRow[0,:] ) * paramsScales
    paramsBar = px.scatter(df, x=df.index, y=df.columns,
              error_y="err_plus", error_y_minus="err_minus",
              title =  """Best-fit parameter values""" )
    paramsBar.update_traces( go.Scatter(
                            mode='markers',
                            marker=dict(color='black', size=14),
                            error_y=dict( color='black', thickness=2, width=10 )
                                    ))
    #paramsBar.add_trace(go.Scatter(x=paramsNames, y=paramsLowValsPCBound[:,0]*paramsScales,
    #                               name=r'$paramsSolnPC - \sigma$',
    #                               line=dict(color='white', width=0), mode='lines', showlegend=False))
    #paramsBar.add_trace(go.Scatter(x=paramsNames, y=paramsHiValsPCBound[:,0]*paramsScales, fill='tonexty',
    #                           name='Default Parameter Values +- sigma', mode='none',
    #                               fillcolor='rgba(253,253,150,1.0)'))
    #paramsBar.add_trace(go.Scatter(x=paramsNames, y=defaultParamValsOrigRow[0,:]*paramsScales,
    #                               name='Default Parameter Values', line=dict(color='black', width=6) ))
    paramsBar.add_trace(go.Scatter(x=paramsNames, y=paramsSoln[:,0]*paramsScales,
                                   mode='markers',
                                   marker=dict(color='green', size=8),
                                   name='Linear regression, |dp|=' 
                                       + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSoln)) ))
    paramsBar.add_trace(go.Scatter(x=paramsNames, y=paramsSolnNonlin[:,0]*paramsScales,
                                   mode='markers',
                                   marker_symbol='x',
                                   marker=dict(color='orange',  size=12),
                                   name='paramsSolnNonlin, |dpPC|='
                                       + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSolnNonlin)) ))
    paramsBar.add_trace(go.Scatter(x=paramsNames, y=paramsSolnElastic[:,0]*paramsScales,
                                   mode='markers',
                                   marker_symbol='square',
                                   marker=dict(color='cyan', size=8),
                                   name='Lasso regression, |dpLasso|='
                                        + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSolnElastic)) ,
                                   line=dict(color='red', width=2)  ))
    #paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnPCBound[:,0]*paramsScales,
    #                                name='paramsSolnPCBound, |dpBound|='
    #                               + '{:.2e}'.format(0.0) ))
    paramsBar.update_yaxes(title="User-scaled parameter value")
    paramsBar.update_xaxes(title="Parameter Name")
    paramsBar.update_layout(hovermode="x")
    paramsBar.update_layout( width=1000, height=500  )


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
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSoln[:,0]*paramsScales,
                                   name='Linear regression, |dp|=' + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSoln)) ))
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnNonlin[:,0]*paramsScales,
                                   name='paramsSolnNonlin, |dpPC|='
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

    sensMatrixDashboard.layout = html.Div(children=[
        html.H1(children='Sensitivity matrix diagnostics'),

        html.Div(children=''' '''),
        
        dcc.Graph( id='biasesOrderFig', figure=biasesOrderFig ),
        dcc.Graph( id='paramsBar', figure=paramsBar ),
        dcc.Graph( id='paramsFig', figure=paramsFig ),
        dcc.Graph( id='biasesFig', figure=biasesFig ),
        dcc.Graph( id='biasesSensScatterFig', figure=biasSensMatrixScatterFig ),
        dcc.Graph( id='dpMinScatterFig', figure=dpMinMatrixScatterFig ),
        dcc.Graph( id='maxSensMetricsFig', figure=maxSensMetricsFig ),
        dcc.Graph( id='biasesVsSensFig', figure=biasesVsSensFig ),
        dcc.Graph( id='biasesVsSensArrowFig', figure=biasesVsSensArrowFig ),
        dcc.Graph( id='residVsBiasFig', figure=residVsBiasFig ),
        dcc.Graph( id='diffBiasesVsSensFig', figure=diffBiasesVsSensFig ),
        dcc.Graph( id='relBiasesVsSensFig', figure=relBiasesVsSensFig ),
        dcc.Graph( id='normlzdSensMatrixColsFig', figure=normlzdSensMatrixColsFig ),
        dcc.Graph( id='normlzdSensMatrixRowsFig', figure=normlzdSensMatrixRowsFig ),
        dcc.Graph( id='vhNormlzdColsFig', figure=vhNormlzdColsFig ),
        dcc.Graph( id='uNormlzdWeightedColsFig', figure=uNormlzdWeightedColsFig ),
        dcc.Graph( id='uNormlzdColsFig', figure=uNormlzdColsFig ),
        dcc.Graph( id='uNormlzdBiasColsFig', figure=uNormlzdBiasColsFig )
    ])

    sensMatrixDashboard.run_server(debug=True)

    return

# Calculate forward nonlinear solution, normalized but not weighted
def fwdFnc(dnormlzdParams, normlzdSensMatrix, normlzdCurvMatrix):
    import numpy as np

    normlzdDefaultBiasesApproxNonlin = \
            normlzdSensMatrix @ dnormlzdParams \
            + 0.5 * normlzdCurvMatrix @ (dnormlzdParams * dnormlzdParams) 

    return normlzdDefaultBiasesApproxNonlin

def solveUsingNonlin(metricsNames, paramsNames, transformedParamsNames, \
                     metricsWeights, obsMetricValsCol, magParamValsRow, \
                     sensNcFilenames, sensNcFilenamesExt, defaultNcFilename, \
                     defaultParamValsOrigRow, \
                     dnormlzdParamsSoln, normlzdSensMatrix, defaultBiasesCol, \
                     normlzdCurvMatrix, \
                     reglrCoef):

    import numpy as np
    import pdb
    from scipy.optimize import minimize

    # Construct numMetrics x numParams matrix of second derivatives, d2metrics/dparams2.
    # The derivatives are normalized by observed metric values and max param values.
    #normlzdCurvMatrix = \
    #    constructNormlzdCurvMatrix(metricsNames, paramsNames, transformedParamsNames, \
    #                               metricsWeights, obsMetricValsCol, magParamValsRow, \
    #                               sensNcFilenames, sensNcFilenamesExt, defaultNcFilename)

    # Define objective function that is to be minimized.
    def objFnc(dnormlzdParams, normlzdSensMatrix, normlzdDefaultBiasesCol, metricsWeights,
               normlzdCurvMatrix, reglrCoef):
        import numpy as np
        import pdb

        dnormlzdParams = np.atleast_2d(dnormlzdParams).T # convert from 1d row array to 2d column array
        chisqd = np.linalg.norm( (-normlzdDefaultBiasesCol \
                                  - fwdFnc(dnormlzdParams, normlzdSensMatrix, normlzdCurvMatrix) \
                                    ) * metricsWeights \
                                , ord=2 \
                               )**1  \
                + reglrCoef * np.linalg.norm( dnormlzdParams, ord=1 )
        #chisqdOrig = np.linalg.norm( (-normlzdDefaultBiasesCol ) * np.reciprocal(metricsWeights) )**2

        return chisqd


    # Perform nonlinear optimization
    normlzdDefaultBiasesCol = defaultBiasesCol/np.abs(obsMetricValsCol)
    dnormlzdParamsSolnNonlin = minimize(objFnc,x0=np.zeros_like(dnormlzdParamsSoln), \
    #dnormlzdParamsSolnNonlin = minimize(objFnc,dnormlzdParamsSoln, \
                               args=(normlzdSensMatrix, normlzdDefaultBiasesCol, metricsWeights,
                               normlzdCurvMatrix, reglrCoef),\
                               method='Powell')
    dnormlzdParamsSolnNonlin = np.atleast_2d(dnormlzdParamsSolnNonlin.x).T

    # Check whether the minimizer actually minimizes chisqd
    # Initial value of chisqd, which assumes parameter perturbations are zero
    chisqdZero = objFnc(np.zeros_like(dnormlzdParamsSoln).T, normlzdSensMatrix, \
                        normlzdDefaultBiasesCol, metricsWeights, \
                        normlzdCurvMatrix, reglrCoef)
    # Optimized value of chisqd, which uses optimal values of parameter perturbations
    chisqdMin = objFnc(dnormlzdParamsSolnNonlin.T, normlzdSensMatrix, \
                        normlzdDefaultBiasesCol, metricsWeights, \
                        normlzdCurvMatrix, reglrCoef)
    print("chisqdZero =", chisqdZero)  
    print("chisqdMin =", chisqdMin)      
    

    dparamsSolnNonlin = dnormlzdParamsSolnNonlin * np.transpose(magParamValsRow)
    paramsSolnNonlin = np.transpose(defaultParamValsOrigRow) + dparamsSolnNonlin
    #print("paramsSoln.T=", paramsSoln.T)
    print("paramsSolnNonlin.T=", paramsSolnNonlin.T)
    print("normlzdSensMatrix@dnPS.x.T=", normlzdSensMatrix @ dnormlzdParamsSolnNonlin)
    print("normlzdSensMatrix@dnormlzdParamsSoln=", normlzdSensMatrix @ dnormlzdParamsSoln)
    print("normlzdDefaultBiasesCol.T=", normlzdDefaultBiasesCol.T)
    #print("normlzdSensMatrix=", normlzdSensMatrix)
    #pdb.set_trace()

    normlzdWeightedDefaultBiasesApproxNonlin = \
             fwdFnc(dnormlzdParamsSolnNonlin, normlzdSensMatrix, normlzdCurvMatrix) \
             * metricsWeights

    scale = 2
    normlzdWeightedDefaultBiasesApproxNonlin2x = \
             fwdFnc(scale*dnormlzdParamsSolnNonlin, normlzdSensMatrix, 1*normlzdCurvMatrix) \
             * metricsWeights

    # defaultBiasesApprox = (forward model soln - default soln)
    defaultBiasesApproxNonlin = normlzdWeightedDefaultBiasesApproxNonlin \
                                * np.reciprocal(metricsWeights) * np.abs(obsMetricValsCol)

    defaultBiasesApproxNonlin2x = normlzdWeightedDefaultBiasesApproxNonlin2x \
                                * np.reciprocal(metricsWeights) * np.abs(obsMetricValsCol)

    return (defaultBiasesApproxNonlin, \
            dnormlzdParamsSolnNonlin, paramsSolnNonlin, \
            defaultBiasesApproxNonlin2x \
           )

def constructNormlzdCurvMatrix(metricsNames, paramsNames, transformedParamsNames,
                               metricsWeights, obsMetricValsCol, magParamValsRow,
                               sens1NcFilenames, sens2NcFilenames, defaultNcFilename):

    """
    For nonlinear 2nd-order term of Taylor series: 0.5*dp^2*d2m/dp2+...,
    construct a numMetrics x numParams matrix of 2nd-order derivatives, d2m/dp2.
    Each row is a different metric.  Each column is a different parameter.
    The matrix is nondimensionalized by the observed values of metrics
    and maximum values of parameters.
    """
    import numpy as np
    import sys
    import netCDF4
    #import matplotlib.pyplot as plt
    import pdb

    from analyze_sensitivity_matrix import setupDefaultMetricValsCol, setupDefaultParamVectors, \
                                           setupSensArrays
    from scipy.interpolate import UnivariateSpline

    if ( len(paramsNames) != len(sens1NcFilenames)   ):
        print("Number of parameters must equal number of netcdf files.")
        quit()

    # Number of tunable parameters
    numParams = len(paramsNames)

    # Number of metrics
    numMetrics = len(metricsNames)

    # For use in normalizing metrics matrices
    invrsObsMatrix = np.reciprocal(np.abs(obsMetricValsCol)) @ np.ones((1,numParams))

    # Set up a column vector of metric values from the default simulation
    defaultMetricValsCol = \
        setupDefaultMetricValsCol(metricsNames, defaultNcFilename)

    # Based on the default simulation,
    #    set up a column vector of metrics and a row vector of parameter values.
    defaultParamValsRow, defaultParamValsOrigRow = \
        setupDefaultParamVectors(metricsNames, paramsNames, transformedParamsNames,
                                 numMetrics, numParams,
                                 defaultNcFilename)
    normlzdDefaultParamValsRow = defaultParamValsRow * np.reciprocal(magParamValsRow)
    defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1,numParams))
    normlzdDefaultMetricValsMatrix = defaultMetricValsMatrix * invrsObsMatrix

    # Based on the numParams sensitivity simulations,
    #    set up a row vector of modified parameter values.
    # Also set up numMetrics x numParams matrix,
    #    each column of which lists the metrics
    #    from one of the sensitivity simulations
    sens1MetricValsMatrix, sens1ParamValsRow, sens1ParamValsOrigRow = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sens1NcFilenames)
    normlzdSens1ParamValsRow = sens1ParamValsRow * np.reciprocal(magParamValsRow)
    normlzdSens1MetricValsMatrix = sens1MetricValsMatrix * invrsObsMatrix

    # Set up sensitivity-simulation matrices from the extended sensitivity simulation
    sens2MetricValsMatrix, sens2ParamValsRow, sens2ParamValsOrigRow = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sens2NcFilenames)
    normlzdSens2ParamValsRow = sens2ParamValsRow * np.reciprocal(magParamValsRow)
    normlzdSens2MetricValsMatrix = sens2MetricValsMatrix * invrsObsMatrix

    # Initialize matrix to store second derivatives of metrics w.r.t. parameters
    normlzdCurvMatrix = np.full_like(sens1MetricValsMatrix, 0.0)
    normlzdCurvMatrix2 = np.full_like(sens1MetricValsMatrix, 0.0)  # 2nd way of calculating derivs

    # Calculate differences in parameter values between default, sensitivity,
    #    and extended sensitivity runs.
#    delta_params_def_sens = sensParamValsRow - defaultParamValsRow
#    delta_params_def_sensExt = sensParamValsRowExt - defaultParamValsRow
#    delta_params_sens_sensExt = sensParamValsRowExt - sensParamValsRow

    # Calculate numMetrics x numParams matrix of metric values.
#    delta_metrics_def_sens = sensMetricValsMatrix - defaultMetricValsMatrix
#    delta_metrics_def_sensExt = sensMetricValsMatrixExt - defaultMetricValsMatrix
#    delta_metrics_sens_sensExt = sensMetricValsMatrixExt - sensMetricValsMatrix
    for col in np.arange(numParams):
        for row in np.arange(numMetrics):

            # Set up three (x,y) points whose 2nd-order derivative we wish to calculate.
            # For the spline code below, the x points need to be ordered from least to greatest.
            if normlzdSens1ParamValsRow[0,col] < normlzdSens2ParamValsRow[0,col]:
                params = [ normlzdSens1ParamValsRow[0,col],
                       normlzdDefaultParamValsRow[0,col],
                       normlzdSens2ParamValsRow[0,col] ]
                metrics = [ normlzdSens1MetricValsMatrix[row,col],
                        normlzdDefaultMetricValsMatrix[row,col],
                        normlzdSens2MetricValsMatrix[row,col] ]
            else:
                params = [ normlzdSens2ParamValsRow[0,col],
                       normlzdDefaultParamValsRow[0,col],
                       normlzdSens1ParamValsRow[0,col] ]
                metrics = [ normlzdSens2MetricValsMatrix[row,col],
                        normlzdDefaultMetricValsMatrix[row,col],
                        normlzdSens1MetricValsMatrix[row,col] ]

            # Calculate second-order spline based on three given (x,y) points.
            metricValsSpline = UnivariateSpline(params,metrics,s=0,k=2)
            # Based on spline, find 2nd derivative at arbitrary point (1).
            # I hope that the derivative has the same value at all points,
            #    since it is a parabola.
            normlzdCurvMatrix[row,col] = metricValsSpline.derivative(n=2)(1)

            # Check results using a different calculation
            coefs = np.polyfit(params, metrics, 2)
            normlzdCurvMatrix2[row,col] = 2*coefs[0]
            #pdb.set_trace()

    print( 'normlzdCurvMatrix=', normlzdCurvMatrix )
    print( 'normlzdCurvMatrix2=', normlzdCurvMatrix2 )

            # Distance between points in simulations = sqrt(dparam**2 + dmetric**2)
#            length_def_sens = np.linalg.norm([delta_params_def_sens[0][col],
#                                              delta_metrics_def_sens[row][col]])
#            length_def_sensExt = np.linalg.norm([delta_params_def_sensExt[0][col],
#                                                 delta_metrics_def_sensExt[row][col]])
#            length_sens_sensExt = np.linalg.norm([delta_params_sens_sensExt[0][col],
#                                                  delta_metrics_sens_sensExt[row][col]])
#            semi_perim = 0.5 * ( length_def_sens + length_def_sensExt + length_sens_sensExt )
#            # area of triangle formed by points.  Use Heron's formula.
#            area = np.sqrt( semi_perim *
#                           (semi_perim-length_def_sens) *
#                           (semi_perim-length_def_sensExt) *
#                           (semi_perim-length_sens_sensExt)
#                          )
#            if (area == 0.0):
#                print( '\nIn calcNormlzdRadiusCurv, area == 0.0 for param ', paramsNames[col],
#                        'and metric ', metricsNames[row] )

            # Greatest distance between parameter values in the 3 simulations:
#            max_params_width = \
#            np.max(np.abs([delta_params_def_sens[0][col],
#                        delta_params_def_sensExt[0][col],
#                        delta_params_sens_sensExt[0][col]]))
#            if (max_params_width == 0.0):
#                print( '\nIn calcNormlzdRadiusCurv, max_params_width == 0.0 for param ', paramsNames[col],
#                        'and metric ', metricsNames[row] )

            # Calculate Menger curvature from triangle area and distance between points:
#            normlzd_radius_of_curv[row][col] = 0.25 * length_def_sens*length_def_sensExt*length_sens_sensExt \
#                                                / area / max_params_width



#    fig, axs = plt.subplots(numMetrics, numParams)
#    for col in np.arange(numParams):
#        for row in np.arange(numMetrics):

#            paramVals = [defaultParamValsRow[0][col], sensParamValsRow[0][col], sensParamValsRowExt[0][col]]
#            metricVals = [defaultMetricValsMatrix[row][col], sensMetricValsMatrix[row][col],
#                  sensMetricValsMatrixExt[row][col]]

#            axs[row, col].scatter( paramVals, metricVals )
#            axs[row, col].set_xlabel(paramsNames[col])
#            axs[row, col].set_ylabel(metricsNames[row])
#            fig.show()

    #pdb.set_trace()

    return normlzdCurvMatrix


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

    from analyze_sensitivity_matrix import setupDefaultMetricValsCol, setupDefaultParamVectors, \
                                           setupSensArrays

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
        setupDefaultParamVectors(metricsNames, paramsNames, transformedParamsNames,
                                 numMetrics, numParams,
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
                        sensNcFilenames)

    # Set up sensitivity-simulation matrices from the extended sensitivity simulation
    sensMetricValsMatrixExt, sensParamValsRowExt, sensParamValsOrigRowExt = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sensNcFilenamesExt)

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
                          metricsWeights, obsMetricValsCol,
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

    from analyze_sensitivity_matrix import setupDefaultMetricValsCol, setupDefaultParamVectors, \
                                           setupSensArrays, calcSvdInvrs, calcParamsSoln

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
        setupDefaultParamVectors(metricsNames, paramsNames, transformedParamsNames,
                                 numMetrics, numParams,
                                 defaultNcFilename)

    defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1,numParams))

    defaultBiasesCol, sensMatrix, normlzdWeightedSensMatrix = \
        calcParamsBoundsHelper(metricsNames, paramsNames, transformedParamsNames,
                           metricsWeights, obsMetricValsCol,
                           numMetrics, numParams,
                           magParamValsRow,
                           defaultMetricValsCol, defaultParamValsRow, defaultParamValsOrigRow,
                           sensNcFilenames)

    defaultBiasesColExt, sensMatrixExt, normlzdWeightedSensMatrixExt = \
        calcParamsBoundsHelper(metricsNames, paramsNames, transformedParamsNames,
                           metricsWeights, obsMetricValsCol,
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
         calcSvdInvrs(normlzdWeightedSensMatrixDiff, sValsRatio)

    paramsSolnPC, paramsLowValsPC, paramsHiValsPC, dparamsSolnPCBound, dnormlzdParamsSolnPC, \
    defaultBiasesApproxPC, defaultBiasesApproxLowValsPC, \
    defaultBiasesApproxHiValsPC = \
             calcParamsSoln(svdInvrsNormlzdWeightedPC, metricsWeights, magParamValsRow, \
                            sensMatrixDiff, normlzdWeightedSensMatrixDiff, \
                            obsMetricValsCol, defaultBiasesCol,
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
                           metricsWeights, obsMetricValsCol,
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

    from analyze_sensitivity_matrix import setupDefaultMetricValsCol, setupSensArrays, \
                                           constructSensMatrix, calcSvdInvrs, calcParamsSoln


    # Based on the numParams sensitivity simulations,
    #    set up a row vector of modified parameter values.
    # Also set up numMetrics x numParams matrix,
    #    each column of which lists the metrics
    #    from one of the sensitivity simulations
    sensMetricValsMatrix, sensParamValsRow, sensParamValsOrigRow = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sensNcFilenames)

    # Set up sensitivity-simulation matrices from the extended sensitivity simulation
    #sensMetricValsMatrixExt, sensParamValsRowExt, sensParamValsOrigRowExt = \
    #    setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
    #                    numMetrics, numParams,
    #                    sensNcFilenamesExt)

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
                                 obsMetricValsCol,
                                 numMetrics, numParams,
                                 beVerbose=True)

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
    #                        obsMetricValsCol, defaultBiasesCol,
    #                        defaultParamValsOrigRow, \
    #                        sValsTruncInvNormlzdWeightedPC,
    #                        vhNormlzdWeighted, \
    #                        numParams, paramsNames,
    #                        transformedParamsNames )


    #pdb.set_trace()

    return ( defaultBiasesCol, sensMatrix, normlzdWeightedSensMatrix )

def createPcaBiplot(normlzdSensMatrix, defaultBiasesCol, obsMetricValsCol, metricsNames, paramsNames):

    import numpy as np
    from pca import pca
    import pdb

    # reduce the data towards 2 PCs
    model = pca(n_components=2, detect_outliers='ht2')

    # Augmented array with LHS and RHS
    augMatrix = np.concatenate((normlzdSensMatrix, -defaultBiasesCol / np.abs(obsMetricValsCol) ), axis=1)

    paramsList = list(paramsNames)
    paramsList.append('dbias')
    augParamsNames = np.asarray(paramsList)

def createPcaBiplot(normlzdSensMatrix, defaultBiasesCol, obsMetricValsCol, metricsNames, paramsNames):


    import numpy as np
    from pca import pca
    import pdb

    # reduce the data towards 2 PCs
    model = pca(n_components=2, detect_outliers='ht2')

    # Augmented array with LHS and RHS
    augMatrix = np.concatenate((normlzdSensMatrix, -defaultBiasesCol / np.abs(obsMetricValsCol) ), axis=1)

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
