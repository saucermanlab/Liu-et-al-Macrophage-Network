% Code orginally from Astor Liu sensAnalysisScr.m (11/2/2017)
% Jingyuan Zhang edited version (got rid of Col and Row as outputs since
% the code ran into errors with them
% May 15, 2018

function [sens,speciesNames] = sensAnalysisScr(perturb,dataset)

% Sensitivity analysis
% ACZ 6.2015
%
% This function will perform a sensitivity analysis given a specific
% perturbation. A sensitivity analysis is a comprehensive and systematic
% knock down of each node (full knock down which can be altered on
% line 50) and a measurement of the change in activity of all nodes with
% that knock down.

% Input:
% perturb = the context within which the cell is modeled, ie the input
% state of the model. Should be expressed as 'w(i) = x' where i is the
% index of the input to be altered and x is the value of the input. If no
% perturbation is required perturb = 'clc' should be used. Note that
% perturb must be a string.

% Outputs:
% sens = full matrix of sensitivity analysis
% Col = most influential (most changed) columns of the matrix
% Row = most sensitive (most changed) rows of the matrix

%parameters and initial conditions
[params,y0] = modelParams; % this accesses the parameters from model Params (part of the ODE code)
tspan = [0 700]; %this is how long you want the simulations to be run
options = [];  % these are options for the ode solver

%pull parameters out in order to alter them in the loop without negating
%the previous alteration
[rpar,tau,ymax,speciesNames]=params{:};   %this line unpackages the structure of params
w = rpar(1,:);  % this is the weight parameter now as an array
n = rpar(2,:);  % this is the hill coefficient parameter now as an array
EC50 = rpar(3,:);   % this is the EC50 now as an array

%perturbations above the control conditions
eval(perturb)    % this evaluates the perturbation to change the input to the model
rpar = [w;n;EC50];   % this and the line below repackage the parameters so they can be used by the ODE function
params = {rpar,tau,ymax,speciesNames};


%default simulation creating a vector showing steady state activation of
%all species
[t,y] = ode15s(@modelODE,tspan,y0,options,params);  % this is an initial simulation using the modelODE, t = time, y = matrix with activity of each species at that time
yEnd0(:,1) = y(end,:);   % this is an array that takes the final activity state for all species and stores it to compare it to the activity with each knock down


%preturbed simulations
yEnd1 = zeros(length(y0));  % this and the line below it form blank arrays to fill later
sens = zeros(length(y0));
deltaP = -1; % this value is used in line 59, it's how much the ymax value is changed



%full sensitivity analysis of full KO
for i = 1:length(ymax);  %this is a for loop going through as many iterations as there are species
    disp(['Simulation # ',num2str(i),' of ', num2str(length(ymax))])  %This displays the iteration number on the command line
    
    ymaxNew = ymax;   % this creates a shadow variable so that you don't corrupt the original parameter values
    ymaxNew(i) = (1+deltaP)*ymax(i);  % this alters the value of the parameter for the species at index i 
    params = {rpar,tau,ymaxNew,speciesNames}; % repackaging the new parameters for the ODE file to read
    tspan = [0 700]; 
    options = []; 
    [t2,y2] = ode15s(@modelODE,tspan,yEnd0,options,params); % a new simulation with the new parameter value
    yEnd1(:,i) = y2(end,:)';   % pulling out the final activation state
    sens(:,i) = (yEnd0 - yEnd1(:,i))/(ymaxNew(i) - ymax(i))*ymax(i)./ymax'; %expressing sensitivity as the difference in activity = (knock down - control activity) * normalization factor
    % normalization factor is 1 if ymax is 1 to start with
end


% %plot of sensitivity matrix
% figure
% cmaprange = [0.5:0.01:1];
% blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
% myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
% colormap(myrgbcmap);
% %the above few lines define the color scheme for color map
% sens = real(sens);
% % maxVal = max(max(sens));
% maxVal = 1;
% caxis([-maxVal, maxVal]);
% imagesc(sens,[-maxVal,maxVal]);
% 
% xlabel('Perturbation','fontsize',16);
% set(gca,'YTick',1:length(speciesNames));
% set(gca,'YTickLabel',speciesNames,'fontsize',4);
% ylabel('Sensitivity of Outputs','fontsize',16);
% set(gca,'XAxisLocation','bottom');
% set(gca,'XTick',1:length(speciesNames));
% set(gca,'XTickLabel',speciesNames,'fontsize',4);
% xticklabel_rotate;
% % title(['Sensitivity Analysis']);
% 
% ax=gca;
% set(gcf,'PaperSize',[8,6]);
% set(gcf,'PaperPosition',[0 0 8 6]);
% set(ax,'Position',[0.15 0.15 0.7 0.8]);
% h=colorbar('location','EastOutside');
% % xf=get(gca,'position');
% set(h,'FontSize',16);
% saveas(gcf,['Sensitivity Matrix ' dataset],'tiffn');
% close
% %  
% %%
%Diminished Sensitivity Matrix
%  %finds the 10 most changed columns in the sensitivity matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S = sum(abs(sens)); % sum gives a row vector suming each colume
% M = sort(S); % sort vector elements in aceding order
% Col = find(S>M(1,(end-10))); % indices in S that give rise to a larger change than the 11th change -> top 10 changed columns
% 
% 
% %finds the 10 most changed rows in the sensitivity matrix
% R = sum(abs(sens),2);
% N = sort(R);
% Row = find(R>N((end-11),1));
% Row = Row(1:10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % %uses the most changed columns to generate a new diminished matrix
% % Col=[4,80,81,87,88,19,6,105,106,107];
% % %Col=[88, 93, 86, 89, 10, 66, 7, 30, 77, 78];
% % Row=[126,116,2,59,29,33,97];
% 
% 
% %generates a diminished sensitivity matrix based on the 10 most changed
% %rows found in line 55
% for i = 1:length(Row);
%     for j = 1:length(Col);
%         DimSens(i,j) = sens(Row(i),Col(j));
%     end
% end
% 
% for i = 1:length(Row)
%     Resp(1,i) = speciesNames(Row(i));
% end
% 
% for i = 1:length(Col)
%     KO(1,i) = speciesNames(Col(i));
% end
% 

% % makes a figure of the abbreviated snesitivity matrix
% figure
% cmaprange = [0.5:0.01:1];
% blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
% myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
% colormap(myrgbcmap);
% DimSens = real(DimSens);
% % DimSens = log2(DimSens);
% % maxVal = max(max(DimSens));
% maxVal = 1;
% caxis([-maxVal, maxVal]);
% imagesc(DimSens,[-maxVal,maxVal]);
% set(gca,'YTick',1:length(Resp));
% set(gca,'YTickLabel',Resp,'fontsize',16);
% set(gca,'XAxisLocation','bottom');
% 
% xlabel('Perturbation','fontsize',16);
% 
% % ylabel('Change in Activity');
% % title(['Sensitivity Matrix (M' num2str(dataset) ' Stimulation)'],'fontsize',16);
% ax=gca;
% set(gcf,'PaperSize',[8,6]);
% set(gcf,'PaperPosition',[0 0 8 6]);
% set(ax,'Position',[0.3 0.3 0.6 0.6]);
% h=colorbar('location','EastOutside');
% % xf=get(gca,'position');
% set(h,'FontSize',16);
% set(gca,'XTick',1:length(KO));
% set(gca,'XTickLabel',KO,'fontsize',16);
% xticklabel_rotate;
% saveas(gcf,['Sensitivity Matrix ' dataset 'Dim'],'tiffn');
% close
% % 
end