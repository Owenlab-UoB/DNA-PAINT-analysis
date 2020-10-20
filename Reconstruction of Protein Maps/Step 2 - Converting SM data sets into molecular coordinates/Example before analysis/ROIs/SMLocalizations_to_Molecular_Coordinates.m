function out = SMLocalizations_to_Molecular_Coordinates (data_file_name_in,summary_Bayesian_analysis,nCalibration,nCalibrationLimit,data_file_name_out)

% Generates x,y molecular coordinates from the DNA-PAINT source data of x,y single
% molecule (SM) localizations for each region of interest (ROI) labelled as Folder 1, 2,.... 
% The x,y SM localization data ('data_channelN.txt') should be structured as:
% x,y,sd with the headers printed in the first line of the file. 

% The molecular coordinates output, saved with the file name as given in data_file_name_out,
% provides the x,y molecular coordinates positions in the first and second
% column and the diameter of the cluster of SM localizations related to
% that molecular coordinate.
% 
% The function uses the information obtained in the 'labels' folder generated 
% when running the Bayesian cluster analysis method for single-molecule 
% localization microscopy data (Nature Protocols, 11, 2499?2514, 2016) on
% the the x,y SM localization data set. To access the 'labels' folder, input
% the 'summary_channelN.txt' file in the summary_Bayesian_analysis. See
% example provided with this function. 
%
% To convert SM localization data into molecular coordinates this function 
% uses the input value of nCalibration and nCalibrationLimit. nCalibration
% correspond to the average value of number of localizations expected per
% DNA-docking strand, while the nCalibrationLimit correspond to the cut-off
% value that is used to discard clusters of SM localizations that are likely
% to correspond to noise.

    files = dir;
    directoryNames = {files([files.isdir]).name};
    directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));

for k=1:length(directoryNames)
    subdirpath = directoryNames{k};
    
    % Open data files for Channel N, with N = 1, 2, etc. 
        
    data_channel1=reshape(textread(fullfile(subdirpath,data_file_name_in),'%f','delimiter',',','headerlines',1),3,[]).';
    data_channel1=data_channel1(2:end,:);
    lenght_c1 =length(data_channel1);
    
    % Read best cluster proposal string from the summary files. 

    strC1= fileread(fullfile(subdirpath,summary_Bayesian_analysis));
    newStrC1 =  strtrim(extractBetween(strC1,11,44));
    labelsC1 = load (fullfile(subdirpath,'labels', char(newStrC1)))';

    % Combine x,y positions with cluster labels in a single matrix.

    data_c1_labels = [data_channel1 labelsC1(2:end)];
    
    % Compute k-means clustering to each cluster of single molecule localizations
    % found in the 'data_channelN.txt' data to recover the molecular 
    % coordinates for each of the clusters, per ROI. 

    data_c1_labels_sort = sortrows(data_c1_labels,4);
    count_c1 = accumarray(data_c1_labels_sort(:,4) ,ones(size(data_c1_labels_sort(:,4))));
    filtered_count_c1 = count_c1(count_c1>1);
    filtered_count_c1 =  [0 filtered_count_c1']; filtered_count_c1 =filtered_count_c1';
    finalPositions = cumsum(filtered_count_c1);
    opts = statset('Display', 'final');
    temp_counts_positions = [filtered_count_c1 finalPositions];
    filtered_temp = temp_counts_positions(temp_counts_positions(:,1)>=nCalibrationLimit , :);
    filtered_count_c1 = filtered_temp(:,1); filtered_count_c1 =  [0 filtered_count_c1']; filtered_count_c1 =filtered_count_c1';
    finalPositions =  filtered_temp(:,2); finalPositions=  [0 finalPositions']; finalPositions =finalPositions';
    nClusters =  length(filtered_count_c1)-1; %calculate all the detected cluster  
    
    for i=1:1:nClusters;
        initial = 1 + finalPositions(i);
        final = finalPositions(i+1); 
        Cluster_xy = [data_c1_labels_sort(initial:final,1) data_c1_labels_sort(initial:final,2)];
        NumLocCluster = filtered_count_c1(i+1); 
        if NumLocCluster> nCalibration
        BindingSites =round(NumLocCluster./nCalibration);
        else
        BindingSites =1;
        end
        [idx,ClusterCentroids, sumd, D] = kmeans(Cluster_xy,BindingSites,'Distance','cityblock','Replicates',5,'Options',opts);
        D2 = min(D,[],2);
        D3 = [idx, D2];  Distances= sortrows(D3,1);
        ClusterDiameter = accumarray(Distances(:,1),Distances(:,2),[],@max); 
        ClusterCentroids = [ClusterCentroids ClusterDiameter];
        dlmwrite(fullfile(subdirpath,data_file_name_out), ClusterCentroids,'-append','delimiter',' ','precision', '%f');
    end
end

