% Code orginally from Astor Liu runValidation.m (11/2/2017)
% Jingyuan Zhang reduced version to only generate IFNg IL4 combined simulation
% May 1, 2018
clear all;
close all;

% Set file directory path

% cd('/Users/jingyuan/Documents/Academic/research/Macrophage Model/Astor macrophage/macmodel_test_original6.3p_JZ_edits_only_necessary_files/Model Simulation/')

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