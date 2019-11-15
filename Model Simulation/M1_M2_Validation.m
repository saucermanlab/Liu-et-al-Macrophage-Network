% ./Model Simulation/M1_M2_Validation.m
% 
% This script generates model responses to LPS+IFNg and IL4 stimulation.
% The output includes Fig 2A and 2B.
% 
% 	Input: ./Input/
% 			mac validation sum M1.xlsx
% 			mac validation sum M2.xlsx
%     			LPS+IFNg validation RNASeq.xlsx
%     			IL4 validation RNASeq.xlsx
%     			mac validation sum M1 PM.xlsx
%     			mac validation sum M2 PM.xlsx
%     			modelODE.m
% 			modelParams.m
% 			y0.mat
% 		./Model Simulation/QuantValidation_3inputs.m
% 
% 	Output: ./Model Simulation/simulation results/
% 			LPS+IFNg validation RNASeq.txt
% 			LPS+IFNg validation RNASeq.txt
% 			IL4 validation RNASeq_raw.txt
% 			IL4 validation RNASeq.txt
% 			mac validation sum M1in0.7_raw.txt
% 			mac validation sum M1in0.7.txt
% 			mac validation sum M2in0.7_raw.txt
% 			mac validation sum M2in0.7.txt
% 			macmodelvalidation_M2in0.7.txt
% 			macmodelvalidation_M1in0.7.txt
% 
% 		./Model Simulation/plots/
% 			Arg1_mrna0.7.tif - Fig 2B
% 			CCL17_mrna0.7.tif
% 			IKBa_mrna0.7.tif
% 			IL1_mrna0.7.tif - Fig 2B
% 			IL6_mrna0.7.tif
% 			IL10_mrna0.7.tif
% 			IL12_mrna0.7.tif
% 			iNOS_mrna0.7.tif - Fig 2B
% 			IRF1_mrna0.7.tif
% 			macmodelvalidation M1 in0.7.tif - Fig 2A
% 			macmodelvalidation M2 in0.7.tif - Fig 2A
% 			MMP3_mrna0.7.tif
% 			MMP7_mrna0.7.tif
% 			MMP9_mrna0.7.tif
% 			MMP12_mrna0.7.tif
% 			PPARg_mrna0.7.tif
% 			SOCS1_mrna0.7.tif - Fig 2B
% 			SOCS3_mrna0.7.tif
% 			TNFa_mrna0.7.tif
% 			VEGF_mrna0.7.tif
% Code adapted from Astor Liu runValidation.m (11/2/2017)
% Jingyuan Zhang reduced version to only generate IFNg+LPS and IL4 simulation
% May 1, 2018

% Running macrophage model validation
clear all;
close all;

% Set your file directory path

cd('Model Simulation/')
addpath('../Input')

val_names={'LPS+IFNg validation RNASeq.xlsx' 'IL4 validation RNASeq.xlsx'};


inputlevel=0.7;
percentMatch=[];
activityChange=[];
resultChart={};
rawResult={};
matchL=[];
yAll=[];
yTran1=[];
yTran2=[];

for val=1:4 % 1:22
    validationfname = sprintf('%s',val_names{val}); % Validation reference file name
    fname = strsplit(validationfname,'.');
    perturb = 'clc'; % Perturb value
    thresh = 2; % Threshold value
    
    [percentMatch(val),activityChange,resultChart,rawResult,matchL,yAll] = QuantValidation_3inputs(validationfname,perturb,thresh,inputlevel);
    
    % Convert the resultChart into table and save as a tab-delimited txt file
    table = cell2table(resultChart(2:end,:)); % Convert the validation outputs into a table
    table.Properties.VariableNames = resultChart(1,:); % Use the output file headings as the table labels
    writetable(table,['./simulation results/' fname{1} 'in' num2str(inputlevel) '.txt'],'Delimiter','\t','WriteRowNames',true); % Write the table variable into a txt file
    
    % Convert the rawResult into table and save as a tab-delimited txt file
    table2 = cell2table(rawResult(2:end,:));
    table2.Properties.VariableNames = rawResult(1,:);
    writetable(table2,['./simulation results/' fname{1} 'in' num2str(inputlevel) '_raw.txt'],'Delimiter','\t','WriteRowNames',true);
    
    if val == 1
        table3 = cell2table(yAll); % Convert the validation outputs into a table
        % table3.Properties.VariableNames = {'species' 'yStart' 'yEnd'}; % Use the output file headings as the table labels
        writetable(table3,['./simulation results/macmodelvalidation_M1' 'in' num2str(inputlevel) '.txt'],'Delimiter','\t','WriteRowNames',true); % Write the table variable into a txt file
        
        speciesNames2 = strrep(yAll(:,1),'_','\_');
        yTran1 = yAll(:,2:end);
        yTran1 = cell2mat(yTran1);
        matches = strfind(speciesNames2,'mrna');
        out = find(~cellfun(@isempty,matches));
        out = out([1:2,4,10:11,14:16,18:19,22:23,25:26,31:32,34,17]);
        ptitle = ['./plots/macmodelvalidation M1 in' num2str(inputlevel) '.tif'];
        
        figure
        colormap(bone);
        imagesc(yTran1(out,1:400),[0,1]);
        ax=gca;
        axis([70,240,1,length(out)]);% colorbar;
        title('LPS+IFN-\gamma','fontsize',24,'FontWeight','normal');
        xlabel('Time (h)','FontSize',24);
        ylabel('Node','FontSize',24);
        set(gca,'YTick',1:length(out));
        set(gca,'YTickLabel',speciesNames2(out),'FontSize',12);
        ax.XTick = [102:40:3100];
        ax.XTickLabel = {'\fontsize{24} 0','\fontsize{24} 4','\fontsize{24} 8',...
            '\fontsize{24} 12',' ',' ',' ','','','\fontsize{24} 24'};
        %     ax.XTickLabel.FontSize = '\fontsize{24}';
        set(gcf,'PaperSize',[5*1.5,3*1.5]);
        set(gcf,'PaperPosition',[0 0 5*1.5 3*1.5]);
        set(ax,'Position',[0.25 0.18 0.65 0.72]);
        h=colorbar('location','EastOutside');
        % xf=get(gca,'position');
        set(h,'FontSize',24);
        saveas(gcf,ptitle,'tiffn');
        close
        
    elseif val == 2
        table4 = cell2table(yAll); % Convert the validation outputs into a table
        % table4.Properties.VariableNames = {'species' 'yStart' 'yEnd'}; % Use the output file headings as the table labels
        writetable(table4,['./simulation results/macmodelvalidation_M2' 'in' num2str(inputlevel) '.txt'],'Delimiter','\t','WriteRowNames',true); % Write the table variable into a txt file
        
        speciesNames2 = strrep(yAll(:,1),'_','\_');
        yTran2 = yAll(:,2:end);
        yTran2 = cell2mat(yTran2);
        %     matches = strfind(speciesNames2,'mrna');
        %     out = find(~cellfun(@isempty,matches));
        %     out = out([1:2,4,6,10:12,14:16,18:20,22:23,25:26,29:34,17]);
        ptitle = ['./plots/macmodelvalidation M2 in' num2str(inputlevel) '.tif'];
        figure
        colormap(bone);
        ax=gca;
        imagesc(yTran2(out,:),[0,1]);
        axis([70,240,1,length(out)]);% colorbar;
        title('IL4','fontsize',24,'FontWeight','normal');
        xlabel('Time (h)','FontSize',24);
        ylabel('Node','FontSize',24);
        set(gca,'YTick',1:length(out));
        set(gca,'YTickLabel',speciesNames2(out),'FontSize',12);
        ax.XTick = [102:40:3100];
        ax.XTickLabel = {'\fontsize{24} 0','\fontsize{24} 4','\fontsize{24} 8',...
            '\fontsize{24} 12',' ',' ',' ','','','\fontsize{24} 24'};
        ax=gca;
        set(gcf,'PaperSize',[5*1.5,3*1.5]);
        set(gcf,'PaperPosition',[0 0 5*1.5 3*1.5]);
        set(ax,'Position',[0.25 0.18 0.65 0.72]);
        h=colorbar('location','EastOutside');
        % xf=get(gca,'position');
        set(h,'FontSize',24);
        saveas(gcf,ptitle,'tiffn');
        close
        
        % Plot for example
        otitle1 = yAll(out,1);
        for ko=1:length(out)
            otitle = otitle1{ko};
            figure
            plot([1:3102],yTran1(out(ko),:),'.-','Color',[1,0.65,0],'LineWidth',4);
            hold on
            plot([1:3102],yTran2(out(ko),:),'.-','Color',[0.13,0.55,0.13],'LineWidth',4);
            hold off
            axis([70,240,0,1.1]);% colorbar;
            line([1,1000],[1.1*0.999,1.1*0.999],'Color','k','LineWidth',2);
            line([240,240],[0,1.2],'Color','k','LineWidth',2);
            ax=gca;
            set(ax, 'box','off')
            title(speciesNames2(out(ko),:),'fontsize',24,'FontWeight','normal');
            xlabel('Time (h)','FontSize',24);
            ylabel('Norm. Activity','FontSize',24);
            set(gca,'YTick',0:0.5:1);
            set(gca,'YTickLabel',[0:0.5:1],'FontSize',24);
            ax.XTick = [102:40:3100];
            ax.XTickLabel = {'\fontsize{24} 0','\fontsize{24} 4','\fontsize{24} 8','\fontsize{24} 12',...
                '\fontsize{24} 16',' ',' ',' ','','','\fontsize{24} 24'};
            ax=gca;
            set(gcf,'PaperSize',[4,3]);
            set(gcf,'PaperPosition',[0 0 4 3]);
            set(ax,'Position',[0.3 0.3 0.6 0.55]);
            if ko == 1
                legend({sprintf( '%s\n%s', 'LPS+', 'IFN\gamma' ),'IL4'},'box','off','FontSize',24,'Location','southeast');
            else
            end
            cd('./plots')
            saveas(gcf,[otitle,num2str(inputlevel),'.tif'],'tiffn');
            cd ..
            close
        end
        
    end
    percentMatch = percentMatch/100
end

cd('..')