function TernaryPlot_color(data)

%Protein contribution in each cluster

Protein_1 = data(:,1); 
Protein_2 = data(:,2); 
Protein_3 = data(:,3); 

% Converting protein contribution into scatter data set to plot ternary
% graph
x = Protein_1.*0.5+ Protein_2;
y = Protein_1;

sample_pts=[x y];

% Calculating most likely cluster composition

 PAG = data(:,2);  meanCSK_TRAF3comp = data;
 indices_1 = find(abs(PAG)<1); % nPAG proteins in clusters is 0 or 1 
 nCSK_TRAF3 = length(indices_1);pCSK_TRAF3 =nCSK_TRAF3 /length(data)*100;
 
 idx1 = setdiff(1:length(meanCSK_TRAF3comp(:,1)),  indices_1);
 meanCSK_TRAF3comp(idx1,:) = []; CSK_TRAF = mean(meanCSK_TRAF3comp);

 TRAF3 = data(:,3);   meanCSK_PAGcomp = data;
 indices_2 = find(abs(TRAF3)<1); % nTRAF3 proteins in clusters is 0 or 1 
 nCSK_PAG = length(indices_2);pCSK_PAG =nCSK_PAG/length(data)*100;
 
 idx2 = setdiff(1:length(meanCSK_PAGcomp(:,1)),  indices_2);
 meanCSK_PAGcomp(idx2,:) = []; CSK_PAG = mean(meanCSK_PAGcomp);
 
 CSK = data(:,1);  meanTRAF3_PAGcomp = data;
 indices_3 = find(abs(CSK)<1); % nPAG proteins in clusters is 0 or 1 
 nTRAF3_PAG = length(indices_3); pTRAF_PAG =  nTRAF3_PAG/length(data)*100;
 
 idx3 = setdiff(1:length( meanTRAF3_PAGcomp(:,1)),  indices_3);
 meanTRAF3_PAGcomp(idx3,:) = []; TRAF3_PAG = mean( meanTRAF3_PAGcomp);

 all =data; 
 indices_4 = [indices_1' indices_2' indices_3']
 idx4 = setdiff(1:length(meanTRAF3_PAGcomp(:,1)),  indices_4);
 all(idx4,:) = [];  all_comp = mean(all);
 
 nallproteins = length(data) - nCSK_TRAF3 - nCSK_PAG -  nTRAF3_PAG; pall =  nallproteins/length(data)*100;
 
 percentage_Proteins = [pCSK_PAG pCSK_TRAF3 pTRAF_PAG pall];
 
 data2 = [CSK_PAG' CSK_TRAF' TRAF3_PAG' all_comp']; data2 = data2';



% Defining color map. Here Protein 1 = Red (Csk), Protein 2 = Blue (PAG), Protein 3 = Green (TRAF3)

sample_color = [1 0.4 0.4].*Protein_1./100 + [0.302 0.745 0.933].*Protein_2./100 + [0 0.8 0.2].*Protein_3./100 ;
sample_color_2 = [1 0 0].*Protein_1./100 + [0 0 1].*Protein_2./100 + [0 1 0].*Protein_3./100 ;
% Defining triangle

P = [0 0; 50 100; 100 0];
% Vertices
scatter(P(:, 1), P(:, 2), 100, 'k');
hold on;
% Edges
for idx = 1:size(P, 1)-1
    plot([P(idx, 1) P(idx+1, 1)], [P(idx, 2) P(idx+1, 2)], 'k');
end
plot([P(end, 1) P(1, 1)], [P(end, 2) P(1, 2)], 'k');

% Colors for demo
C = ones(size(sample_pts, 1), 1).*sample_pts(:, 1);
s = 70;
% Scatter sample points
scatter(sample_pts(:, 1), sample_pts(:, 2),s, sample_color_2, 'filled');pbaspect([1 1 1])






end