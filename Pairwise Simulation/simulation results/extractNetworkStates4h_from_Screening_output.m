% Jeff 8/31/2018
% loading data from pairwise screen at 4 h for plotting in cytoscape
% pairwiseScreeningData = readtable('Screening_output.txt');

% 4 hours corresponds to column yEnd12
nodeNames = pairwiseScreeningData{1:137,2};
LPSIFNg = pairwiseScreeningData{contains(pairwiseScreeningData.Row,'LPS+IFNg'),'yEnd12'};
IL4 = pairwiseScreeningData{5618:5754,'yEnd12'}; % hard coded, regexp was taking too long
IFNg = pairwiseScreeningData{5070:5206,'yEnd12'}; % hard coded, regexp was taking too long
IFNgIL4 = pairwiseScreeningData{contains(pairwiseScreeningData.Row,'IFNg+IL4'),'yEnd12'};
T = table(nodeNames,LPSIFNg,IL4,IFNg,IFNgIL4);
writetable(T,'network_states_4h.csv')

