tic
% Run sensitivity analysis for M1 (LPS + IFNg) and M2 (IL-4) macrophages
clc;
clear all;
close all;

% Set file directory path

path = matlab.desktop.editor.getActiveFilename;
cd(fileparts(path))
% cd('J:\private\Netflux-master\Netflux\MacrophageModel\macmodel_test_original2')
% addpath('J:\private\Netflux-master\Netflux\MacrophageModel')
% cd('J:\private\Netflux-master\Netflux')

run = 0; % 0: analysis only; 1: run M1-M2 sens analysis
scr = 0; % 0: M1-M2 sens analysis 1: screening sens analysis

if run==1
%% Run sensitivity analysis

% Input: LPS and IFNg
[sensM1,~,~,speciesM1] = sensAnalysis('w(4)=0.7;w(5)=0.7;',1);

% Input: IL-4
[sensM2,~,~,speciesM2] = sensAnalysis('w(6)=0.7;',2);

% Convert the resultChart into table and save as a tab-delimited txt file
table = array2table(real(sensM1)); % Convert the validation outputs into a table
table.Properties.VariableNames = speciesM1; % Use the output file headings as the table labels
table.Properties.RowNames = speciesM1;
writetable(table,['macmodelSensLPS+IFNg.txt'],'Delimiter','\t','WriteRowNames',true); % Write the table variable into a txt file

% Convert the resultChart into table and save as a tab-delimited txt file
table2 = array2table(real(sensM2)); % Convert the validation outputs into a table
table2.Properties.VariableNames = speciesM2; % Use the output file headings as the table labels
table2.Properties.RowNames = speciesM2;
writetable(table2,['macmodelSensIL4.txt'],'Delimiter','\t','WriteRowNames',true); % Write the table variable into a txt file

else
    %% Plot only
if scr == 0

[params,y0] = modelParams; % this accesses the parameters from model Params (part of the ODE code)

%pull parameters out in order to alter them in the loop without negating
%the previous alteration
[rpar,tau,ymax,speciesNames]=params{:};   %this line unpackages the structure of params

speciesNames = strrep(speciesNames,'_','\_'); % Add a backslash before the underscores in the texts for correct display
matches = strfind(speciesNames,'mrna');
out = find(~cellfun(@isempty,matches));
labels = speciesNames(out);
for i=1:size(labels,2)
    labels{i} = lower(labels{i});
    labels{i}(1)=upper(labels{i}(1));
    speciesNames{out(i)} = labels{i};
end
speciesNames = strrep(speciesNames,'Il1\_mrna','Il1a(mRNA)');
speciesNames = strrep(speciesNames,'Inos\_mrna','Nos2(mRNA)');
speciesNames = strrep(speciesNames,'\_mrna','(mRNA)');

        

dataset = {'macmodelSensLPS+IFNg.txt' 'macmodelSensIL4.txt'};
dataset2 = {'macmodelSensIFN.txt' 'macmodelSensMix.txt'};
%      base = {'macmodelvalidation_M1.txt' 'macmodelvalidation_M2.txt'};  
    figname = {'Sensitive Matrix LPS+IFNg' 'Sensitive Matrix IL4'};
    figname2 = {'Sensitive Matrix IFNg' 'Sensitive Matrix IFNg+IL4'};

    R2=[];
for k = 1:2
    sens = dlmread(dataset{k},'\t',1,1);
    %plot of sensitivity matrix
figure
cmaprange = [0.5:0.01:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
%the above few lines define the color scheme for color map

% if k==2
% % maxVal = max(max(sens));
% % maxVal = max(max(abs(sens)));
% caxis([-0.15, 0.15]);
% imagesc(sens,[-0.15, 0.15]);
% else
    caxis([-1, 1]);
imagesc(sens,[-1, 1]);

% end
% text(50,150,'Perturbation','FontSize',20); % x label
% text(-20,90,'Output','FontSize',20,'Rotation',90); % y label

xlabel('Perturbation','fontsize',20);
set(gca,'YTick',1:length(speciesNames));
set(gca,'YTickLabel',{''},'fontsize',4);
% set(gca,'YTickLabel',speciesNames,'fontsize',4);
ylabel('Output','fontsize',20);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(speciesNames));
set(gca,'XTickLabel',{''},'fontsize',4);
% set(gca,'XTickLabel',speciesNames,'fontsize',4);
% xticklabel_rotate;
% title(['Sensitivity Analysis']);

ax=gca;
set(gcf,'PaperSize',[8,6]);
set(gcf,'PaperPosition',[0 0 8 6]);
set(ax,'Position',[0.15 0.15 0.7 0.8]);
h=colorbar('location','EastOutside');
% xf=get(gca,'position');
set(h,'FontSize',16);
saveas(gcf,[figname{k}],'pdf');
close


%  
%%
%Diminished Sensitivity Matrix
%  %finds the 10 most changed columns in the sensitivity matrix
S = [];
M = [];
Col = [];

R = [];
N = [];
Row = [];

Dimsens = [];

S = sum(abs(sens)); % sum gives a row vector suming each colume
[M,I]= sort(S); % sort vector elements in aceding order
Col = find(S>M(1,(end-10))); % indices in S that give rise to a larger change than the 11th change -> top 10 changed columns

if k==1
    I1=fliplr(I);
else
end
R2 = [R2,S(I1)'];
snames = speciesNames(I1);

%finds the 10 most changed rows in the sensitivity matrix
% R = sum(abs(sens(:,Col)),2);
R = sum(abs(sens),2);
N = sort(R);
Row = find(R>=N((end-9),1));
lenR = length(Row);
if lenR>10
    for lenk=[1:lenR-10]
        Row(find(min(R(Row)))) = []; % Make sure only 10 values are picked
    end
else
end

% %uses the most changed columns to generate a new diminished matrix
% Col=[4,80,81,87,88,19,6,105,106,107];
% %Col=[88, 93, 86, 89, 10, 66, 7, 30, 77, 78];
% Row=[126,116,2,59,29,33,97];


%generates a diminished sensitivity matrix based on the 10 most changed
%rows found in line 55
for i = 1:length(Row);
    for j = 1:length(Col);
        DimSens(i,j) = sens(Row(i),Col(j));
    end
end

for i = 1:length(Row)
    Resp(1,i) = speciesNames(Row(i));
end

for i = 1:length(Col)
    KO(1,i) = speciesNames(Col(i));
end


% makes a figure of the abbreviated snesitivity matrix
figure
cmaprange = [0.5:0.01:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
DimSens = real(DimSens);
% DimSens = log2(DimSens);
% maxVal = max(max(DimSens));
% maxVal = 0.25;
caxis([-1, 1]);
imagesc(DimSens,[-1,1]);
text(0,14.5,['Knockdown of Most Influential Nodes'],'FontSize',28); % x label
text(-4,11,['Most Sensitive Nodes'],'FontSize',28,'Rotation',90); % y label
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(KO));
set(gca,'XTickLabel',KO,'fontsize',16);
xticklabel_rotate;
% xlabel('Perturbation','FontSize',16);
set(gca,'YTick',1:length(Resp));
set(gca,'YTickLabel',Resp,'fontsize',16);
% ylabel('Change in Activity','FontSize',16);
% title(['Sensitivity Matrix (M' num2str(dataset) ' Stimulation)'],'fontsize',16);
ax=gca;
set(gcf,'PaperSize',[8,7]);
set(gcf,'PaperPosition',[0 0 8 7]);
set(ax,'Position',[0.3 0.3 0.6 0.6]);
h=colorbar('location','EastOutside');
% xf=get(gca,'position');
set(h,'FontSize',28);
saveas(gcf,[figname{k} ' Dim'],'pdf');
close

close
end

% Plot M1-M2 sens bar charts
    R2=R2/max(max(R2)); % Normalize to the maximum
    
    figure
    bar(R2,1,'histc');
    colormap([1,0.6,0;0,0.6,0]);
    axis([0,length(snames)+1,0,1.1]);
    line([0,200],[1.1*0.99,1.1*0.99],'Color','k','LineWidth',2);
    line([length(snames)+1,length(snames)+1],[0,1.2],'Color','k','LineWidth',2);
    ax=gca;
    set(ax,'box','off','LineWidth',2);
%     title('Accumulated Activity Change Post-knockdown','fontsize',24,'FontWeight','normal');
    xlabel('Knockdown Nodes','FontSize',24);
%     ylabel('Norm. Accum. Activity Change','FontSize',24);
    set(gca,'XTick',1:10:length(snames));
    set(gca,'XTickLabel',' ','fontsize',1);
    ax.YTick = [0:0.5:1];
    ax.YTickLabel = {'\fontsize{24} 0','\fontsize{24} 0.5','\fontsize{24} 1','\fontsize{24} 6','\fontsize{24} 8'};
%     xticklabel_rotate;
%     text(40,-0.8,'Knockdown Node','FontSize',24); % x label
    text(-10,0,sprintf('Activity'),'FontSize',24,'Rotation',90); % y label
    legend({'LPS+IFN','IL4'},'Position',[0.25,0.5,0.3,0.2],'box','off','FontSize',24);
    set(gcf,'PaperSize',[12,3]);
    set(gcf,'PaperPosition',[0 0 12 3]);
    set(ax,'Position',[0.15 0.15 0.8 0.7]);
    saveas(gcf,'Sens_Accum_M1-M2','pdf');
    close

for k = 1:2
    sens = dlmread(dataset2{k},'\t',1,1);
    %plot of sensitivity matrix
figure
cmaprange = [0.5:0.01:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
%the above few lines define the color scheme for color map

caxis([-1, 1]);
imagesc(sens,[-1, 1]);

% end
% text(50,150,'Perturbation','FontSize',20); % x label
% text(-20,90,'Output','FontSize',20,'Rotation',90); % y label

xlabel('Perturbation','fontsize',20);
set(gca,'YTick',1:length(speciesNames));
set(gca,'YTickLabel',{''},'fontsize',4);
% set(gca,'YTickLabel',speciesNames,'fontsize',4);
ylabel('Output','fontsize',20);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(speciesNames));
set(gca,'XTickLabel',{''},'fontsize',4);
% set(gca,'XTickLabel',speciesNames,'fontsize',4);
% xticklabel_rotate;
% title(['Sensitivity Analysis']);

ax=gca;
set(gcf,'PaperSize',[8,6]);
set(gcf,'PaperPosition',[0 0 8 6]);
set(ax,'Position',[0.15 0.15 0.7 0.8]);
h=colorbar('location','EastOutside');
% xf=get(gca,'position');
set(h,'FontSize',16);
saveas(gcf,[figname2{k}],'pdf');
close


%  
%%
%Diminished Sensitivity Matrix
%  %finds the 10 most changed columns in the sensitivity matrix
S = [];
M = [];
Col = [];

R = [];
N = [];
Row = [];

Dimsens = [];

S = sum(abs(sens)); % sum gives a row vector suming each colume
[M,I]= sort(S); % sort vector elements in aceding order
Col = find(S>M(1,(end-10))); % indices in S that give rise to a larger change than the 11th change -> top 10 changed columns

if k==1
    I1=fliplr(I);
else
end
R2 = [R2,S(I1)'];
snames = speciesNames(I1);

%finds the 10 most changed rows in the sensitivity matrix
% R = sum(abs(sens(:,Col)),2);
R = sum(abs(sens),2);
N = sort(R);
Row = find(R>=N((end-9),1));
lenR = length(Row);
if lenR>10
    for lenk=[1:lenR-10]
        Row(find(min(R(Row)))) = []; % Make sure only 10 values are picked
    end
else
end

% %uses the most changed columns to generate a new diminished matrix
% Col=[4,80,81,87,88,19,6,105,106,107];
% %Col=[88, 93, 86, 89, 10, 66, 7, 30, 77, 78];
% Row=[126,116,2,59,29,33,97];


%generates a diminished sensitivity matrix based on the 10 most changed
%rows found in line 55
for i = 1:length(Row);
    for j = 1:length(Col);
        DimSens(i,j) = sens(Row(i),Col(j));
    end
end

for i = 1:length(Row)
    Resp(1,i) = speciesNames(Row(i));
end

for i = 1:length(Col)
    KO(1,i) = speciesNames(Col(i));
end


% makes a figure of the abbreviated snesitivity matrix
figure
cmaprange = [0.5:0.01:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
DimSens = real(DimSens);
% DimSens = log2(DimSens);
% maxVal = max(max(DimSens));
% maxVal = 0.25;
caxis([-1, 1]);
imagesc(DimSens,[-1,1]);
text(0,14.5,['Knockdown of Most Influential Nodes'],'FontSize',28); % x label
text(-4,11,['Most Sensitive Nodes'],'FontSize',28,'Rotation',90); % y label
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(KO));
set(gca,'XTickLabel',KO,'fontsize',16);
xticklabel_rotate;
% xlabel('Perturbation','FontSize',16);
set(gca,'YTick',1:length(Resp));
set(gca,'YTickLabel',Resp,'fontsize',16);
% ylabel('Change in Activity','FontSize',16);
% title(['Sensitivity Matrix (M' num2str(dataset) ' Stimulation)'],'fontsize',16);
ax=gca;
set(gcf,'PaperSize',[8,7]);
set(gcf,'PaperPosition',[0 0 8 7]);
set(ax,'Position',[0.3 0.3 0.6 0.6]);
h=colorbar('location','EastOutside');
% xf=get(gca,'position');
set(h,'FontSize',28);
saveas(gcf,[figname2{k} ' Dim'],'pdf');
close

close
end

%%
else % Screening sens analysis

    cd('../Sensitivity Screening');
%parameters and initial conditions
[params,y0] = modelParams; % this accesses the parameters from model Params (part of the ODE code)

%pull parameters out in order to alter them in the loop without negating
%the previous alteration
[rpar,tau,ymax,speciesNames]=params{:};   %this line unpackages the structure of params

speciesNames = strrep(speciesNames,'_','\_'); % Add a backslash before the underscores in the texts for correct display
matches = strfind(speciesNames,'mrna');
out = find(~cellfun(@isempty,matches));
labels = speciesNames(out);
for i=1:size(labels,2)
    labels{i} = lower(labels{i});
    labels{i}(1)=upper(labels{i}(1));
    speciesNames{out(i)} = labels{i};
end
speciesNames = strrep(speciesNames,'Il1\_mrna','Il1a(mRNA)');
speciesNames = strrep(speciesNames,'Inos\_mrna','Nos2(mRNA)');
speciesNames = strrep(speciesNames,'\_mrna','(mRNA)');

    files = dir('simulation results/macmodel*.txt'); % List all .txt files
    fname=struct2cell(files); % Extract file names
    fnames=fname(1,1:end)';
    tend=24;
%     if tend==4
%         fnames=fname(1,2:5:length(files))';
%     elseif tend==5
%         fnames=fname(1,5:5:length(files))';
%     elseif tend==24
%         fnames=fname(1,1:5:length(files))';
%     elseif tend==48
%         fnames=fname(1,4:5:length(files))';
%     end

    S2=[];
    R2=[];
    cnames={};
for dataset=1:length(fnames)

%     sens = dlmread(fnames{dataset},'\t',1,1);
    sens = real(tblread(['simulation results/',fnames{dataset}],'\t')); % Read sens matrix files

figname = strsplit(fnames{dataset},'.'); % Extract figure names
ftitle = strsplit(figname{1},'_');
cnames{dataset} = [ftitle{2}];
% ftitle = strrep(ftitle,'_','\_');

    %plot of sensitivity matrix
figure
cmaprange = [0.5:0.01:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
%the above few lines define the color scheme for color map


% maxVal = max(max(sens));
% maxVal = max(max(abs(sens)));
caxis([-1, 1]);
imagesc(sens,[-1, 1]);

% text(50,150,'Perturbation','FontSize',20); % x label
% text(-20,90,'Output','FontSize',20,'Rotation',90); % y label

xlabel('Perturbed Nodes','fontsize',20);
set(gca,'YTick',1:length(speciesNames));
set(gca,'YTickLabel',{''},'fontsize',4);
% set(gca,'YTickLabel',speciesNames,'fontsize',4);
ylabel('Affected Nodes','fontsize',20);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(speciesNames));
set(gca,'XTickLabel',{''},'fontsize',4);
% set(gca,'XTickLabel',speciesNames,'fontsize',4);
% xticklabel_rotate;
title([ftitle{2} '.7'],'fontsize',20,'FontWeight','normal');

ax=gca;
set(gcf,'PaperSize',[4,3]);
set(gcf,'PaperPosition',[0 0 4 3]);
set(ax,'Position',[0.2 0.2 0.7 0.6]);
h=colorbar('location','EastOutside');
% xf=get(gca,'position');
set(h,'FontSize',16);
saveas(gcf,['plots/',strcat(figname{1},'.',figname{2}) '.pdf'],'pdf');
close
%  
%%
%Diminished Sensitivity Matrix
S = [];
M = [];
Col = [];

R = [];
N = [];
Row = [];

Dimsens = [];
%  %finds the 10 most changed columns in the sensitivity matrix
S = sum(abs(sens)); % sum gives a row vector suming each colume
M = sort(S); % sort vector elements in aceding order
Col = find(S>M(1,(end-10))); % indices in S that give rise to a larger change than the 11th change -> top 10 changed columns


%finds the 10 most changed rows in the sensitivity matrix
R = sum(abs(sens),2);
N = sort(R);
Row = find(R>=N((end-9),1));
lenR = length(Row);
if lenR>10
    for lenk=[1:lenR-10]
        Row(find(min(R(Row)))) = []; % Make sure only 10 values are picked
    end
else
end

% Accumulated influence
S2 = [S2,S'];
% Accumulated sensitivity
R2 = [R2,R];


% %uses the most changed columns to generate a new diminished matrix
% Col=[4,80,81,87,88,19,6,105,106,107];
% %Col=[88, 93, 86, 89, 10, 66, 7, 30, 77, 78];
% Row=[126,116,2,59,29,33,97];


%generates a diminished sensitivity matrix based on the 10 most changed
%rows found in line 55
for i = 1:length(Row);
    for j = 1:length(Col);
        DimSens(i,j) = sens(Row(i),Col(j));
    end
end

for i = 1:length(Row)
    Resp(1,i) = speciesNames(Row(i));
end

for i = 1:length(Col)
    KO(1,i) = speciesNames(Col(i));
end


% makes a figure of the abbreviated snesitivity matrix
figure
cmaprange = [0.5:0.01:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
DimSens = real(DimSens);
% DimSens = log2(DimSens);
% maxVal = max(max(DimSens));
% maxVal = 0.25;
caxis([-1, 1]);
imagesc(DimSens,[-1,1]);
text(2,14.5,'Most Influential Nodes','FontSize',28); % x label
text(-4,10,'Most Sensitive Nodes','FontSize',28,'Rotation',90); % y label
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(KO));
set(gca,'XTickLabel',KO,'fontsize',16);
xticklabel_rotate;
% xlabel('Perturbation','FontSize',16);
set(gca,'YTick',1:length(Resp));
set(gca,'YTickLabel',Resp,'fontsize',16);
% ylabel('Change in Activity','FontSize',16);
title([ftitle{2} '.7'],'fontsize',28,'FontWeight','normal');
ax=gca;
set(gcf,'PaperSize',[8,7]);
set(gcf,'PaperPosition',[0 0 8 7]);
set(ax,'Position',[0.3 0.3 0.6 0.6]);
h=colorbar('location','EastOutside');
% xf=get(gca,'position');
set(h,'FontSize',28);
saveas(gcf,['plots/',strcat(figname{1},'.',figname{2}) ' (Dim).pdf'],'pdf');
close
end
% Convert the total infl. and sens. into table and save as a tab-delimited txt file
tableS = array2table(real(S2)); % Convert the validation outputs into a table
% tableS.Properties.VariableNames = cnames; % Use the output file headings as the table labels
tableS.Properties.RowNames = speciesNames;
writetable(tableS,['plots/','screeningSens_influence' num2str(tend) '.txt'],'Delimiter','\t','WriteRowNames',true); % Write the table variable into a txt file


cgoS = clustergram(S2,'RowLabels',speciesNames,'ColumnLabels',cnames,'ColorMap','redbluecmap',...
    'Symmetric',0,'Cluster','row','DisplayRange',max(max(S2)),'DisplayRatio',0.1);

figure
plot(cgoS);
caxis([0, max(max(S2))]);
ax=gca;
set(gcf,'PaperSize',[12,18]);
set(gcf,'PaperPosition',[0 0 12 18]);
set(ax,'Position',[0.3 0.1 0.4 0.85]);
h=colorbar('location','NorthOutside');
% xf=get(gca,'position');
set(h,'FontSize',20);
saveas(gcf,['plots/','screeningSens_influence' num2str(tend) '.pdf'],'pdf');
close

tableR = array2table(real(R2)); % Convert the validation outputs into a table
% tableR.Properties.VariableNames = cnames; % Use the output file headings as the table labels
tableR.Properties.RowNames = speciesNames;
writetable(tableR,['plots/','screeningSens_sensitivity' num2str(tend) '.txt'],'Delimiter','\t','WriteRowNames',true); % Write the table variable into a txt file


cgo = clustergram(R2,'RowLabels',speciesNames,'ColumnLabels',cnames,'ColorMap','redbluecmap',...
    'Symmetric',0,'Cluster','row','DisplayRange',max(max(R2)),'DisplayRatio',0.1);
figure
plot(cgo);
caxis([0, max(max(R2))]);
ax=gca;
set(gcf,'PaperSize',[12,18]);
set(gcf,'PaperPosition',[0 0 12 18]);
set(ax,'Position',[0.3 0.1 0.4 0.85]);
h=colorbar('location','NorthOutside');
% xf=get(gca,'position');
set(h,'FontSize',20);
saveas(gcf,['plots/','screeningSens_sensitivity' num2str(tend) '.pdf'],'pdf');
close
end
end
toc