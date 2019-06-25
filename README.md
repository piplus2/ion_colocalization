R scripts used for the data analysis of the paper:

[Inglese, P., Correia, G., Pruski, P., Glen, R. C., & Takats, Z. (2019). Colocalization features for classification of tumors using desorption electrospray ionization mass spectrometry imaging. Analytical chemistry.](https://pubs.acs.org/doi/10.1021/acs.analchem.8b05598)

Data is available on Mendeley (DOI: 10.17632/wpy848vsfy.1)

https://data.mendeley.com/datasets/wpy848vsfy/1

## Data preparation

1_recalibrate_imzml.R
Performs the single point re-calibration of the negative ion mode data using Palmitic acid theoretical mass as reference.

2_preprocess_within_sample.R
Pre-processing of the spectra within each tissue section sample.

3_preprocess_between_sample.R
Pre-processing to match peaks between multiple tissue section samples.

4_generate_multi_data.R
Load and arrange the spectra intensities for the following calculation of correlations.

5_correlations.R
Calculate the Spearman's correlations using the `multi_data` matrices.

    6_correlations_offset.R
Calculate the Spearman's correlations from the simulated offset datasets.

## Classification

1_classification.R
Main script for PLS-DA modelling of cancer type using the Spearman's correlations.

    2_classification_offset.R
Perform classification using the correlation features from the simulated offset data.

    3_classification_mean.R
Perform classification using the mean peaks intensity.

    4_classification_mean_offset.R
Perform classification using the mean peaks intensity from the simulated offset data.

    _accuracy_ scripts
Save the accuracies of the models.

## Univariate analysis

    1_univariate_analysis.R
Univariate Kruskal-Wallis test for significant different correlations between the contrasts.

	2_graph_adjacencies.R
Write the adjacencies of the N most significant correlations.

    3_select_top_correlations.R
Select the N most significant correlations for model interpretation purposes.
