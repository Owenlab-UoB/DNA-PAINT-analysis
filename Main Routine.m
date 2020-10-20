%% Reconstruction of proteins maps

% Step 1: Convert single molecule coordinates to x,y molecular coordinates

% Generates x,y molecular coordinates from the DNA-PAINT source data of x,y,sd single
% molecule (SM) localizations for each region of interest (ROI) using the
% information obtained after running Bayesian cluster analysis method for single-molecule 
% localization microscopy data. To convert SM localizations into molecular coordinate 
% the function uses k-means clustering and the average value of number of localizations 
% expected per binding site.
%  

data_file_name_in = 'data_channel2.txt'
data_file_name_out= 'ClusterCentroids_Ch2_v1.txt' % this file corresponds to molecular coordinates

summary_Bayesian_analysis = 'summary_channel2.txt'

nCalibration = 52; %num of localization per binding site 
nCalibrationLimit = 27;% cut-off value 

SMLocalizations_to_Molecular_Coordinates(data_file_name_in,summary_Bayesian_analysis,nCalibration,nCalibrationLimit,data_file_name_out)

%% Step 2: Plot SM localizations and molecular coordinates in the same graph

% Visualize the x,y SM localizations and x,y molecular coordinates for a specific ROI
% denoted by 1, 2, 3...
clear all

data_file_name_in = 'data_channel2.txt'
data_file_name_out= 'ClusterCentroids_Ch2_v1.txt' % this file corresponds to molecular coordinates

n_ROI = 1; %Change to the number of ROI that you would like to visualize

Plotting_SMLoc_and_Molecular_Coordinates(data_file_name_in,data_file_name_out,n_ROI)

%% Next steps correspond to Mixed protein cluster analysis

% Step 1: Generate merged data set after running Reconstruction of proteins
% maps for each of the different imaged proteins (here named as channel2,
% channel3, channel4).  

% Insert the Merged_Channel_ClusterCentroids.m function inside the folder that
% contains the ROIs folder (exactly like the case of the provided Example data folder).

data_file_name_in_1 = 'ClusterCentroids_Ch2_v1.txt'
data_file_name_in_2 = 'ClusterCentroids_Ch3_v1.txt'
data_file_name_in_3 = 'ClusterCentroids_Ch4_v1.txt'

data_file_name_out= 'ClusterCentroids_merged_v1.txt' % this file contains all the molecular coordinates of Ch2, Ch3 and Ch4.

Merged_Channel_ClusterCentroids(data_file_name_in_1,data_file_name_in_2,data_file_name_in_3,data_file_name_out);


%% Step 2: Select mixed protein cluster that contain at least two molecular units of each
% cluster combination type (i.e. RG, RB, GB, RGB) in the cluster. 
% Run this step, after running the after running Bayesian cluster analysis method for each of 
% the different sets of ROIs corresponding to each Cell of the same type and exposed to the same stimulii 

% Insert the Protein_Cluster_percentage.m function inside the folder that
% contains the ROIs folder (exactly like the case of the provided Example data folder).

clear all

data_file_name_in = 'list_Centroids_merged_v1.csv'
data_file_name_out= 'pProteins_Cell1.txt'

Protein_Cluster_percentage(data_file_name_in,data_file_name_out)


%% Step2: Combining data from Cell 1, Cell 2, etc...
 
clear all
load pProteins_Cell1.txt; dataset1 = pProteins_Cell1';
load pProteins_Cell2.txt; dataset2 = pProteins_Cell2';

Data_combined = [dataset1 dataset2]; Data_combined=Data_combined';

dlmwrite('pProteins_Total.txt',Data_combined,'delimiter',',','precision',10,'-append');


%% Step 3: Plotting ternary plot graphs

clear all

%load file that contains the percentual composition of each protein in the
%cluster (i.e. column 1: % of protein 1; column 2: % protein 2; column 3: %
%of protein 3) and remember to copy the functions TernaryPlot_color.m and 
% TernaryPlot_color_optional.m in the same folder that contains this data file.

load pProteins_Total.txt; 

data = pProteins_Total;

% Plot ternary plot of data
TernaryPlot_color(data);

% Enable TernaryPlot_color_optional function to plot most likely cluster 
% composition for each pair combination in the ternary plot. 
TernaryPlot_color_optional(data); 
