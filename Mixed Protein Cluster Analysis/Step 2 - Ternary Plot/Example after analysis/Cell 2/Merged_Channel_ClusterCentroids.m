function Merged_Channel_ClusterCentroids(data_file_name_in_1,data_file_name_in_2,data_file_name_in_3,data_file_name_out)

files = dir;
directoryNames = {files([files.isdir]).name};
directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));

for i=1:1:length(directoryNames)
    subdirpath = directoryNames{i};
    
    % Open data files for Channel 1, Channel 2 and Merged Channel. 
    
    data_channel2=reshape(textread(fullfile(subdirpath,'ClusterCentroids_Ch2_v1.txt'),'%f','delimiter',',','headerlines',0),3,[]).';
    data_channel3=reshape(textread(fullfile(subdirpath,'ClusterCentroids_Ch3_v1.txt'),'%f','delimiter',',','headerlines',0),3,[]).';
    data_channel4=reshape(textread(fullfile(subdirpath,'ClusterCentroids_Ch4_v1.txt'),'%f','delimiter',',','headerlines',0),3,[]).';
  
   
    data_merged = [data_channel2' data_channel3' data_channel4'];
    data_merged = data_merged';     
    dlmwrite( fullfile(subdirpath,'ClusterCentroids_merged_v1.txt'),data_merged, 'delimiter',',','precision',10 );
    
end