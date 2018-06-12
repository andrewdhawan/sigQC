30/07/2017/Version 0.1.16: Fixed minor issues + changed display names of boxplot metrics (in interim Alessandro added negative controls).
21/06/2017/Version 1.13: Fixed issue with radarplot matrix output - standardise number of columns even when such values are not computed.
21/06/2017/Version 0.1.12: Alessandro fixed error with long file names in windows, adding virtual drive solution, and other bug fixes.
20/06/2017/Version 1.12: Fixed a tiny bug causing the warning for the heatmaps, removed borders around the barplots, added more colours (dynamic based on how many datasets there are) for the radarplot, no longer outputting the threshold value for binarization as it is just the median.
20/06/2017/Version 1.11: Alessandro added bootstrap resampling + bug fixes.
10/06/2017/Version 1.10: Fixed outptutting of radarchart_table.txt for consistency. Rowname format now standardised to have '.' replacing spaces in the signature name or dataset name and have an '_' between the signature and dataset name.
31/05/2017/Version 1.09: Added error checking to ensure that every gene signature has at least 2 elements, and that every dataset contains at least 2 expressed elemetns of each signature, otherwise the signature/dataset is removed.
30/05/2017/Version 1.08: Numerous bug fixes, added underscores to headers of output tables, rank product now outputting for all genes, fix for complexheatmap issue with differing numbers of rows.
29/05/2017/Version 1.07: Numerous bug fixes, automated legend placement/margins resizing for autocorrelation heatmaps, added proportion of variance explained by first PCA to radarplot.
17/05/2017/Version 1.06: Updates to rankProduct calculation, origin parameter error checking, removal of logged parameter, and improvements to legend placement in radarplot.
16/05/2017/Version 1.05: Further bug fix re: Rankproduct function; added parameters to main function: logged - T/F if data was log transformed or not (required) and origin - character vector of different origins of datasets, default is assumption that all datasets are same origin (both used for rankproduct).
16/05/2017/Version 1.04: Small bug fix re: Rankproduct function
15/05/2017/Version 1.03: Small bug fix re: handling of matrices with high concentrations of NA values - able to bypass these values to prevent errors.
13/05/2017/Version 1.02: Small bug fix re: Coefficient of variance ratio in radar plot; added absolute values to maintain range in [0,1].
12/05/2017/Version 1.01: Small bug fix re: calculating proportion of expressed/NA genes, accounting for genes not present as rownnames in the dataset.
11/05/2017/Version 1.00: Base version.
