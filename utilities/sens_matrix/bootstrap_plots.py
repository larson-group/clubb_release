import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
from sklearn.decomposition import PCA

def bootstrap_plots(numSamples,
                    numMetrics,
                    numMetricsToTune,
                    metricsNames,
                    residualsBootstrapMatrix,
                    residualsFullDataCol,
                    defaultBiasesCol,
                    lossesLeftOut,
                    lossesInSample,
                    paramsNames,
                    paramsSolnNonlin,
                    lossesDrop,
                    lossesFullData):
    # create folders to save in
    folderName = get_folder_name(numMetricsToTune)
    create_folder(folderName)
    # Plot biases
    for i in range(numMetrics):
        residualsSample = residualsBootstrapMatrix[:, i]
        sns.displot(data=residualsSample, kde=True)
        plt.axvline(-defaultBiasesCol[i, 0], color='red', linestyle='--', label="default residuals")
        plt.axvline(residualsFullDataCol[i, 0], color='green', linestyle='--', label="tuned residuals")
        plt.xlabel(metricsNames[i])
        plt.title(f"mean: {np.mean(residualsSample):.2g}, std: {np.std(residualsSample):.2g}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"Outputs/{folderName}/biases/{metricsNames[i]}_distribution.png", dpi=300)
        plt.show()

    # PCA to 1 component
    pca = PCA(n_components=1)
    pca_result = pca.fit_transform(residualsBootstrapMatrix)
    sns.displot(data=pca_result.flatten(), kde=True)
    plt.axvline(pca.transform(-defaultBiasesCol.T)[0, 0], color='red', linestyle='--', label="default residuals")
    plt.axvline(pca.transform(residualsFullDataCol.T)[0, 0], color='green', linestyle='--', label="tuned residuals")
    plt.title(f"mean: {np.mean(pca_result):.2g}, std: {np.std(pca_result):.2g}")
    plt.xlabel("First pca component of metric samples")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"Outputs/{folderName}/biases/PCA_distribution.png", dpi=300)
    plt.show()


    # remove untuned metrics
    residualsBootstrapMatrix = residualsBootstrapMatrix[:, :numMetricsToTune]
    residualsFullDataCol = residualsFullDataCol[:numMetricsToTune]

    # MSR plot
    msrBiases = np.mean(np.power(residualsBootstrapMatrix-residualsFullDataCol.T, 2), axis=0)
    numRowsMsr = int(np.sqrt(msrBiases.size/2))
    msrBiases = msrBiases.reshape(numRowsMsr, 2*numRowsMsr)
    plt.imshow(msrBiases, extent=(0, 360, -90, 90))
    plt.colorbar(label="MSR")
    plt.title(fr"$\frac{{1}}{{{numSamples}}}\sum_{{i=1}}^{{{numSamples}}}(y_{{boot,i}}-y_{{full}})^2$")
    plt.tight_layout()
    plt.savefig(f"Outputs/{folderName}/biases/MSR_bias.png", dpi=300)
    plt.show()

    # Losses plot
    # Create a 2x3 grid of plots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))  # 2 rows, 3 columns

    # Determine a common color scale for all losses
    vmin = min(lossesLeftOut.min(), lossesInSample.min(), lossesDrop.min(), lossesFullData.min())
    vmax = max(lossesLeftOut.max(), lossesInSample.max(), lossesDrop.max(), lossesFullData.max())

    # Plot lossesLeftOut
    im1 = axes[0, 0].imshow(lossesLeftOut, extent=(0, 360, -90, 90))
    axes[0, 0].set_title(f"Avg out of Sample Loss, avg increase {np.mean((lossesInSample - lossesLeftOut)):.2g}")
    axes[0, 0].set_xlabel("Longitude")
    axes[0, 0].set_ylabel("Latitude")
    plt.colorbar(im1, ax=axes[0, 0])

    # Plot lossesInSample
    im2 = axes[0, 1].imshow((lossesInSample - lossesLeftOut) / lossesInSample, extent=(0, 360, -90, 90))
    axes[0, 1].set_title(r"$(e_{in}-e_{out})/e_{in}$")
    axes[0, 1].set_xlabel("Longitude")
    axes[0, 1].set_ylabel("Latitude")
    plt.colorbar(im2, ax=axes[0, 1])

    # Plot lossesInSample
    im3 = axes[0, 2].imshow((lossesInSample - lossesLeftOut) / lossesInSample > 0, extent=(0, 360, -90, 90))
    axes[0, 2].set_title("1 = better out of sample")
    axes[0, 2].set_xlabel("Longitude")
    axes[0, 2].set_ylabel("Latitude")
    plt.colorbar(im3, ax=axes[0, 2])

    # Plot lossesDrop
    im4 = axes[1, 0].imshow(lossesDrop, extent=(0, 360, -90, 90))
    axes[1, 0].set_title(f"Losses after dropping, avg improvement: {np.mean((lossesFullData - lossesDrop)):.2g}")
    axes[1, 0].set_xlabel("Longitude")
    axes[1, 0].set_ylabel("Latitude")
    plt.colorbar(im4, ax=axes[1, 0])

    # Plot lossesFullData
    im5 = axes[1, 1].imshow((lossesFullData - lossesDrop) / lossesFullData, extent=(0, 360, -90, 90), vmin=-5)
    axes[1, 1].set_title(r"$(e_{full}-e_{drop})/e_{full}$")
    axes[1, 1].set_xlabel("Longitude")
    axes[1, 1].set_ylabel("Latitude")
    plt.colorbar(im5, ax=axes[1, 1])

    # Plot lossesFullData
    im6 = axes[1, 2].imshow((lossesFullData - lossesDrop) / lossesFullData > 0, extent=(0, 360, -90, 90))
    axes[1, 2].set_title("1 = better after drop")
    axes[1, 2].set_xlabel("Longitude")
    axes[1, 2].set_ylabel("Latitude")
    plt.colorbar(im6, ax=axes[1, 2])

    # Adjust layout
    plt.tight_layout()
    plt.savefig(f"Outputs/{folderName}/losses/LossesMap.png", dpi=300)
    plt.show()

    # Generate a histogram and kernel density estimation for every parameter in the samples
    for i in range(len(paramsNames)):
        paramSample = paramsSolnNonlin[:, i].flatten()
        sns.displot(data=paramSample, kde=True)
        plt.xlabel(paramsNames[i])
        plt.title(f"mean: {np.mean(paramSample):.2g}, std: {np.std(paramSample):.2g}")
        plt.tight_layout()
        plt.savefig(f"Outputs/{folderName}/params/{paramsNames[i]}_distribution.png", dpi=300)
        plt.show()

    # Generate a scatter plot for every pair of parameters
    paramSamplesDf = pd.DataFrame(data=paramsSolnNonlin[:, :, 0], columns=paramsNames)
    sns.pairplot(paramSamplesDf, diag_kind='hist')
    plt.savefig(f"Outputs/{folderName}/params/Scatter.png", dpi=300)
    plt.show()


def get_folder_name(numMetrics):
    # Determine the subfolder name based on numMetrics
    folder_name = ""
    if numMetrics == 162:
        folder_name = "20x20"
    elif numMetrics == 72:
        folder_name = "30x30"
    elif numMetrics == 32:
        folder_name = "45x45"
    else:
        folder_name = "unknownGrid"
    return folder_name

def create_folder(folder_name):
    # Create the 'Outputs' folder if it doesn't exist
    outputs_folder = "Outputs"
    os.makedirs(outputs_folder, exist_ok=True)

    # Create the subfolder
    subfolder_path = os.path.join(outputs_folder, folder_name)
    os.makedirs(subfolder_path, exist_ok=True)

    # Create inner folders for biases, losses, and params
    inner_folders = ["biases", "losses", "params"]
    for folder in inner_folders:
        inner_folder_path = os.path.join(subfolder_path, folder)
        os.makedirs(inner_folder_path, exist_ok=True)
