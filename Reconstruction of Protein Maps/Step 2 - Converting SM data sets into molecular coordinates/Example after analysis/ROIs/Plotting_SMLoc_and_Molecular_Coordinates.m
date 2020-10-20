function out = Plotting_SMLoc_and_Molecular_Coordinates(data_file_name_in,data_file_name_out,n_ROI)
% Plot SM localizations and x,y molecular coordinates for a specific ROI
% denoted by 1, 2, 3...

    files = dir;
    directoryNames = {files([files.isdir]).name};
    directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));

    
    subdirpath = directoryNames{n_ROI};
    
    % Open simulated data files for Channel 1, Channel 2 and Merged Channel. 
    
    data_channel1=reshape(textread(fullfile(subdirpath,data_file_name_in),'%f','delimiter',',','headerlines',1),3,[]).';
    data_channel1=data_channel1(2:end,:);
    
    data_BS=load(fullfile(subdirpath,data_file_name_out));
    
    sz = 30; tz =2; 
    figure(1)
    scatter (data_channel1(:,1),data_channel1(:,2),tz,'b','filled');pbaspect([1 1 1])
    hold on
    scatter (data_BS(:,1),data_BS(:,2),sz,'r','filled')
end 
