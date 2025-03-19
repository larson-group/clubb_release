import numpy as np

def bootstrap_calculations(numSamples,
                           metricsWeights,
                           metricsNames,
                           paramsNames,
                           numMetrics,
                           numMetricsToTune,
                           normMetricValsCol,
                           magParamValsRow,
                           defaultParamValsOrigRow,
                           normlzdSensMatrixPoly,
                           normlzdDefaultBiasesCol,
                           normlzdCurvMatrix,
                           reglrCoef,
                           defaultBiasesCol):
    from quadtune_driver import solveUsingNonlin, lossFncMetrics, fwdFnc
    # In order to do sampling with replacement of the metrics we sample the corresponding indices with replacement
    metricsSampleIdxMatrix = np.random.randint(0, numMetricsToTune, size=(numSamples, numMetricsToTune))
    paramsSolnNonlin = np.zeros((numSamples, len(paramsNames), 1))
    defaultBiasesApproxNonlinMatrix = np.full((numSamples, numMetrics, 1), np.nan)
    lossesLeftOut = np.full((numSamples, numMetricsToTune, 1), np.nan)
    lossesInSample = np.full((numSamples, numMetricsToTune, 1), np.nan)
    # Loop through every bootstrap sample and calculate parameter values
    for sampleIdx, metricsSampleIdxRow in enumerate(metricsSampleIdxMatrix):
        if sampleIdx % 10 == 0:
            print(sampleIdx)
        defaultBiasesApproxNonlin, dnormlzdParamsSolnNonlin, paramsSolnNonlin[sampleIdx], \
            dnormlzdParamsSolnLin, paramsSolnLin, defaultBiasesApproxNonlin2x, \
            defaultBiasesApproxNonlinNoCurv, defaultBiasesApproxNonlin2xCurv = (
            solveUsingNonlin(metricsNames[metricsSampleIdxRow],
                             metricsWeights[metricsSampleIdxRow],
                             normMetricValsCol[metricsSampleIdxRow],
                             magParamValsRow,
                             defaultParamValsOrigRow,
                             normlzdSensMatrixPoly[metricsSampleIdxRow, :],
                             normlzdDefaultBiasesCol[metricsSampleIdxRow],
                             normlzdCurvMatrix[metricsSampleIdxRow, :],
                             reglrCoef,
                             beVerbose=False))

        # Find the loss on the left out regions
        metricsLeftOut = np.setdiff1d(range(numMetricsToTune), metricsSampleIdxRow)
        lossesLeftOut[sampleIdx, metricsLeftOut] = lossFncMetrics(dnormlzdParamsSolnNonlin,
                                                                  normlzdSensMatrixPoly[metricsLeftOut, :],
                                                                  normlzdDefaultBiasesCol[metricsLeftOut],
                                                                  metricsWeights[metricsLeftOut],
                                                                  normlzdCurvMatrix[metricsLeftOut, :],
                                                                  metricsLeftOut.size)
        metricsInSample = np.unique(metricsSampleIdxRow)
        lossesInSample[sampleIdx, metricsInSample] = lossFncMetrics(dnormlzdParamsSolnNonlin,
                                                                    normlzdSensMatrixPoly[metricsInSample, :],
                                                                    normlzdDefaultBiasesCol[metricsInSample],
                                                                    metricsWeights[metricsInSample],
                                                                    normlzdCurvMatrix[metricsInSample, :],
                                                                    metricsInSample.size)

        defaultBiasesApproxNonlinMatrix[sampleIdx, :] = fwdFnc(dnormlzdParamsSolnNonlin, normlzdSensMatrixPoly,
                                                               normlzdCurvMatrix, numMetrics) * np.abs(normMetricValsCol)

    # get solution of the full dataset
    tunedBiasesFullData, dnormlzdParamsSolnNonlinFullData, *_ = solveUsingNonlin(metricsNames,
                                                                                 metricsWeights,
                                                                                 normMetricValsCol,
                                                                                 magParamValsRow,
                                                                                 defaultParamValsOrigRow,
                                                                                 normlzdSensMatrixPoly,
                                                                                 normlzdDefaultBiasesCol,
                                                                                 normlzdCurvMatrix,
                                                                                 reglrCoef,
                                                                                 beVerbose=False)
    residualsFullDataCol = -tunedBiasesFullData - defaultBiasesCol
    residualsBootstrapMatrix = -defaultBiasesApproxNonlinMatrix[:, :, 0] - defaultBiasesCol.T

    # lower and upper bounds for error pars plot
    lower_bounds_boot = np.quantile(paramsSolnNonlin[:, :, 0], 0.025, axis=0)
    upper_bounds_boot = np.quantile(paramsSolnNonlin[:, :, 0], 1 - 0.025, axis=0)
    median_boot = np.quantile(paramsSolnNonlin[:, :, 0], 0.5, axis=0)
    param_bounds_boot = np.vstack((upper_bounds_boot, median_boot, lower_bounds_boot))

    # Process lossesLeftOut
    lossesLeftOut = np.nanmean(lossesLeftOut[:, :, 0], axis=0)
    metricsToNotDrop = np.argsort(lossesLeftOut)[:-2]
    numRowsLoss = int(np.sqrt(lossesLeftOut.size / 2))
    lossesLeftOut = lossesLeftOut.reshape(numRowsLoss, 2 * numRowsLoss)

    # Process lossesInSample
    lossesInSample = np.nanmean(lossesInSample[:, :, 0], axis=0)
    lossesInSample = lossesInSample.reshape(numRowsLoss, 2 * numRowsLoss)

    # Process losses after dropping
    defaultBiasesApproxNonlinDrop, dnormlzdParamsSolnNonlinDrop, *_ = solveUsingNonlin(metricsNames[metricsToNotDrop],
                                                                                       metricsWeights[metricsToNotDrop],
                                                                                       normMetricValsCol[
                                                                                           metricsToNotDrop],
                                                                                       magParamValsRow,
                                                                                       defaultParamValsOrigRow,
                                                                                       normlzdSensMatrixPoly[
                                                                                       metricsToNotDrop, :],
                                                                                       normlzdDefaultBiasesCol[
                                                                                           metricsToNotDrop],
                                                                                       normlzdCurvMatrix[
                                                                                       metricsToNotDrop, :],
                                                                                       reglrCoef,
                                                                                       beVerbose=False)

    lossesDrop = lossFncMetrics(dnormlzdParamsSolnNonlinDrop,
                                normlzdSensMatrixPoly[:numMetricsToTune, :],
                                normlzdDefaultBiasesCol[:numMetricsToTune],
                                metricsWeights[:numMetricsToTune],
                                normlzdCurvMatrix[:numMetricsToTune, :],
                                numMetricsToTune)
    lossesDrop = lossesDrop.reshape(numRowsLoss, 2 * numRowsLoss)

    # Process losses full data
    lossesFullData = lossFncMetrics(dnormlzdParamsSolnNonlinFullData,
                                    normlzdSensMatrixPoly[:numMetricsToTune, :],
                                    normlzdDefaultBiasesCol[:numMetricsToTune],
                                    metricsWeights[:numMetricsToTune],
                                    normlzdCurvMatrix[:numMetricsToTune, :],
                                    numMetricsToTune)
    lossesFullData = lossesFullData.reshape(numRowsLoss, 2 * numRowsLoss)
    return paramsSolnNonlin, residualsFullDataCol, residualsBootstrapMatrix, param_bounds_boot, lossesDrop, lossesFullData, lossesLeftOut, lossesInSample