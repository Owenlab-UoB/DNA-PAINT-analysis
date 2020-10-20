function Protein_Cluster_percentage(data_file_name_in,data_file_name_out)

files = dir;
directoryNames = {files([files.isdir]).name};
directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));

 for k=1:1:length(directoryNames)
    subdirpath = directoryNames{k};
    
    results_merged=csvread(fullfile(subdirpath,data_file_name_in),1,1);
    results_merged= results_merged(:,5:7); % these columns corresponds to the num of proteins per cluster;
    
    % When PAG & TRAF3 is 1 or 0, then its a pure Csk cluster -> remove
    % all of those rows
    TF1 = results_merged(:,2)==0 & results_merged(:,3)==0 ;
    results_merged(TF1,:) = [];
    TF2 = results_merged(:,2)==1 & results_merged(:,3)==0 ;
    results_merged(TF2,:) = [];
    TF3 = results_merged(:,2)==0 & results_merged(:,3)==1 ;
    results_merged(TF3,:) = [];
    TF4 = results_merged(:,2)==1 & results_merged(:,3)==1; 
    results_merged(TF4,:) = [];
    
    %When PAG & Csk is 1 or 0, then its a pure TRAF3 cluster -> remove
    % all of those rows
    
    TF5 = results_merged(:,1)==0 & results_merged(:,2)==0 ;
    results_merged(TF5,:) = [];
    TF6 = results_merged(:,1)==1 & results_merged(:,2)==0 ;
    results_merged(TF6,:) = [];
    TF7 = results_merged(:,1)==0 & results_merged(:,2)==1 ;
    results_merged(TF7,:) = [];
    TF8 = results_merged(:,1)==1 & results_merged(:,2)==1; 
    results_merged(TF8,:) = [];
    
    %When TRAF3 & Csk is 1 or 0, then its a pure PAG cluster -> remove
    % all of those rows
    
    TF9 = results_merged(:,1)==0 & results_merged(:,3)==0 ;
    results_merged(TF9,:) = [];
    TF10 = results_merged(:,1)==1 & results_merged(:,3)==0 ;
    results_merged(TF10,:) = [];
    TF11 = results_merged(:,1)==0 & results_merged(:,3)==1 ;
    results_merged(TF11,:) = [];
    TF12 = results_merged(:,1)==1 & results_merged(:,3)==1; 
    results_merged(TF12,:) = [];
    
    total_cluster = sum(results_merged,2);
    norm_results_merged =results_merged./total_cluster.*100;
     
    dlmwrite(data_file_name_out,norm_results_merged, 'delimiter',',','precision',10,'-append');
    end