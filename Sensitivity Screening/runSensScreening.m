% Code orginally from Astor Liu runSensScreening.m (11/2/2017)
% Jingyuan Zhang edited version
% May 15, 2018

tic
% Run sensitivity analysis for M1 (LPS + IFNg) and M2 (IL-4) macrophages
clc;
clear all;
close all;

% read the validation sheet
    [~, txt, raw] = xlsread('Stimuli chart.xlsx');
    input1 = txt(2:end, 2); %second column, if 2 inputs used will need to alter this code
    input2 = txt(2:end, 3);
    inputCode = txt(2:end, 5); % this is the column with the code that will alter the input
    for ki=1:length(inputCode)
    code1 = strsplit(inputCode{ki}, ';');
%         for s=1:length(code1)
%             if length(strfind(code1{s},'w'))>0
%                 code1{s} = [code1{s}(1:end-1) '0.7']; % Set input levels
%             end
%         end
     inputCode{ki}=strjoin(code1,';');
    end

validationIDs = txt(2:end, 1);

% Loop through inputs
for i=1:length(inputCode)
    perturb=inputCode{i};
    dataset=strcat(input1{i},'+',input2{i});
    
    [sensM1,speciesM1] = sensAnalysisScr(perturb,dataset);

% % Input: IL-4
% [sensM2,~,~,speciesM2] = sensAnalysis('w(6)=1;',2);

% Convert the resultChart into table and save as a tab-delimited txt file
table = array2table(real(sensM1)); % Convert the validation outputs into a table
table.Properties.VariableNames = speciesM1; % Use the output file headings as the table labels
table.Properties.RowNames = speciesM1;

writetable(table,['./simulation results/macmodelSens_' dataset '_0.7.txt'],'Delimiter','\t','WriteRowNames',true); % Write the table variable into a txt file
end
toc