% ./Model Simulation/IFNg_IL4_Validation.m
% 
% This script generates model responses to IFNg and IL4 combined stimulation.
% 
% 	Input: 	./Input/
% 			IFNg validation GSE84520.xlsx
% 			IFNg+IL4 validation GSE84520.xlsx
%     			IL4 validation Illum-GSE84520.xlsx
%     			modelODE.m
% 			modelParams.m
% 			y0.mat
% 		./Model Simulation/QuantValidation_3inputs.m
% 
% 	Output: ./Model Simulation/simulation results/
% 			IFNg validation GSE84520_act.txt
% 			IFNg validation GSE84520_raw.txt
% 			IFNg validation GSE84520_validation.txt
% 			IFNg+IL4 validation GSE84520_act.txt
% 			IFNg+IL4 validation GSE84520_raw.txt
% 			IFNg+IL4 validation GSE84520_validation.txt
% 			IL4 validation Illum-GSE84520_act.txt
% 			IL4 validation Illum-GSE84520_raw.txt
% 			IL4 validation Illum-GSE84520_validation.txt
% Code orginally from Astor Liu runValidation.m (11/2/2017)
% Jingyuan Zhang reduced version to only generate IFNg IL4 combined simulation
% May 1, 2018
clear all;
close all;

% Set your file directory path here

cd('Model Simulation/')
addpath('../Input')

val_names={'IFNg validation GSE84520.xlsx' 'IFNg+IL4 validation GSE84520.xlsx' 'IL4 validation Illum-GSE84520.xlsx'};

inputlevel=0.7;

for val= 1:length(val_names)
    
    validationfname = sprintf('%s',val_names{val}); % Validation reference file name
    fname = strsplit(validationfname,'.');
    perturb = 'clc'; % Perturb value
    thresh = 2; % Threshold value
    
    [percentMatch(val),activityChange,resultChart,rawResult,matchL,yAll] = QuantValidation_3inputs(validationfname,perturb,thresh,inputlevel);
    
    % Convert the resultChart into table and save as a tab-delimited txt file
    table = cell2table(resultChart(2:end,:)); % Convert the validation outputs into a table
    table.Properties.VariableNames = resultChart(1,:); % Use the output file headings as the table labels
    writetable(table,['simulation results/' fname{1} '_validation.txt'],'Delimiter','\t','WriteRowNames',true); % Write the table variable into a txt file
    
    % Convert the rawResult into table and save as a tab-delimited txt file
    table2 = cell2table(rawResult(2:end,:));
    table2.Properties.VariableNames = rawResult(1,:);
    writetable(table2,['simulation results/' fname{1} '_raw.txt'],'Delimiter','\t','WriteRowNames',true);
    
    % Write simulation results
    table3 = cell2table(yAll); % Convert the validation outputs into a table
    %     table3.Properties.VariableNames = {'species' 'yStart' 'yEnd'}; % Use the output file headings as the table labels
    writetable(table3,['simulation results/' fname{1} '_act.txt'],'Delimiter','\t','WriteRowNames',true); % Write the table variable into a txt file
 end

percentMatch = percentMatch/100

cd('..')