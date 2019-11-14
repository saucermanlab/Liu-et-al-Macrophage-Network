% Code orginally from Astor Liu pairstimulation.m (9/21/2017)
% Jingyuan Zhang reduced version (no figures, no comparison to classic
% markers)
% May 8, 2018

% Running macrophage model validation
clc;
close all;

% Set file directory path
cd('/Users/jingyuan/Documents/Academic/research/Macrophage Model/Astor macrophage/macmodel_test_original6.3p_JZ_edits_only_necessary_files/Pairwise Simulation')

% Choose validation reference file; 0: mac validation spreadsheet model
% v5.1.xlsx; 1:
val_names={'M0 vs M1 validation exp automated_mRNA only.xlsx'...
    'M0 vs M2 validation exp automated_mRNA only.xlsx'};

% Read the stimuli combinations to be tested
[~, txt, raw] = xlsread('Stimuli chart.xlsx');
combo1 = txt(2:end, 2); %second column, if 2 inputs used will need to alter this code
combo2 = txt(2:end, 3);
combo3 = txt(2:end, 4);
comboCode = txt(2:end, 5); % this is the column with the code that will alter the input

for val=1:length(val_names)
    validationfname = sprintf('%s',val_names{val}); % Validation reference file name
    perturb = 'clc'; % Perturb value
    thresh = 0.1; % Threshold value
    
    % tableResult=[];
    %     tableRaw=[];
    for com=1:length(combo1)
        
        % Use these inputs as the inputs in validation instead of the inputs in
        % the validation spreedsheets
        sti1=sprintf('%s',combo1{com});
        sti2=sprintf('%s',combo2{com});
        sti3=sprintf('%s',combo3{com});
        stiCode=sprintf('%s',comboCode{com});
        
        [percentMatch(val,com),activityChange,resultChart,...
            rawResult,matchL,speciesNames,yEnd] = PairstimuliScreening_3inputs(validationfname,...
            perturb,thresh,sti1,sti2,sti3,stiCode);
        
        % Convert the resultChart into table and save as a tab-delimited txt file
        table = cell2table(resultChart(2:end,:)); % Convert the validation outputs into a table
        table.Properties.VariableNames = resultChart(1,:); % Use the output file headings as the table labels
        % writetable(table,[validationfname num2str(val) '.txt'],'Delimiter','\t','WriteRowNames',true); % Write the table variable into a txt file
        
        % Convert the rawResult into table and save as a tab-delimited txt file
        table2 = cell2table(rawResult(2:end,:));
        table2.Properties.VariableNames = rawResult(1,:);
        % writetable(table2,['macmodelvalidation' num2str(val) '_raw.txt'],'Delimiter','\t','WriteRowNames',true);
        
        
        % Convert the simulated outputs into table and save as a tab-delimited txt file
        % table3 = array2table(real(yEnd));
        table3 = cell2table(yEnd);
        comboname{com} = strcat(sti1,{'+'},sti2);
        outRNs = strcat(repmat(comboname{com},137,1),num2str([1:137]'));
        table3.Properties.RowNames = outRNs;
        
        
        if com==1
            tableResult=table;
            tableRaw=table2;
            tableOut=table3;
        else
            tableResult=[tableResult;table];
            tableRaw=[tableRaw;table2];
            tableOut=[tableOut;table3];
        end
        
    end
    writetable(tableResult,['simulation results/' validationfname '_Screening.txt'],'Delimiter','\t','WriteRowNames',true); % Write the table variable into a txt file
    writetable(tableRaw,['simulation results/' validationfname '_Screening_raw.txt'],'Delimiter','\t','WriteRowNames',true);
    
    if val==1
        writetable(tableOut,['simulation results/' 'Screening_output.txt'],'Delimiter','\t','WriteRowNames',true);
    else
    end  
end

percentMatch=[percentMatch(1:2,:)];

for l=1:length(combo1)
    comboname{l}=[sprintf('%s',combo1{l}) '+' sprintf('%s',combo2{l})];
end
tableP = array2table(percentMatch');
tableP.Properties.VariableNames = {'Match_to_M1' 'Match_to_M2'};
tableP.Properties.RowNames = comboname;
    writetable(tableP,['simulation results/' 'Screening_percentMatch.txt'],'Delimiter','\t','WriteRowNames',true);
