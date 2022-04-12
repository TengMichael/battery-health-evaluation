# Project of battery-health-evaluation
This project includes code and battery data related to Zhongwei Deng et al., "Battery health evaluation using a short random segment of constant current charging". 

# Datasets
In this study, four types of batteries are investigated, namely, LiNixCoyAl1−x−yO2 (NCA), LiNixMnyCo1−x−yO2 (NMC), LiFePO4 (LFP), and a blend of NMC and LiCoO2 (NMC-LCO). The dataset of NMC-LCO is contributed by Hawaii Natural Energy Institute (HNEI), and the other three are contributed by Sandia National Laboratories (SNL). The raw datasets of batteries come from [BatteryArchive](https://www.batteryarchive.org), which is a website to present battery data in a uniform format. People can download the battery datasets from the website or find it in the [Release](https://github.com/TengMichael/battery-health-evaluation/releases) of this project. 

# Release
In the [Release](https://github.com/TengMichael/battery-health-evaluation/releases) of this project, HNEI_raw_data.mat and SNL_raw_data.mat are the raw data extracted from .cvs files, HNEI_cell.mat, LFP_cell.mat, NCA_cell.mat and NMC_cell.mat include battery charging data and capacity of each cycle. All the above data can be extracted by running extract_data.m 

# Correlation analysis
Correlations between battery state of health (SOH) and features are analysed.
correlation_analysis_12seg.m is used to analyse correlation at a fixed number of segments, e.g. 12 in this study.
correlation_analysis_diff_seg.m is used to analyse correlation at different numbers of segments in this study.

# SOH estimation
soh_MLR.m, soh_SGPR.m, and soh_DCNN.m are used to estimate battery SOH by using multiple linear regression (MLR), sparse Gaussian process regression (SGPR), and deep convolutional neural network (DCNN).

# Requirement
Matlab>=2018a

Gaussian Process Regression toolbox in http://www.GaussianProcess.org/gpml/code
