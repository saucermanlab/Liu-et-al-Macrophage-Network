function [percentMatch,activityChange,resultChart,rawResult,matchL,speciesNames,yTran] = PairstimuliScreening_3inputs(validationfname,perturb,thresh,sti1,sti2,sti3,stiCode)
    rawResult = {};    
    yTran = {};

    %parameters and initial conditions
    [params,y0] = modelParams; % this accesses the parameters from model Params (part of the ODE code)
    tspan0 = [0 5];
    tspan = [0 300]; %this is how long you want the simulations to be run
    options = [];  % these are options for the ode solver

    %pull parameters out in order to alter them in the loop without negating
    %the previous alteration
    [rpar,tau,ymax,speciesNames]=params{:};   %this line unpackages the structure of params
    w = rpar(1,:);  % this is the weight parameter now as an array
    n = rpar(2,:);  % this is the hill coefficient parameter now as an array
    EC50 = rpar(3,:);
    
    % read the validation sheet
    [~, txt, raw] = xlsread(validationfname);
    % remove rows without data
%     noData = cellfun(@(x)isequal(x,'No Data'), txt(1:end, 7));
%     txt(noData, :) = []; 
%     noData = cellfun(@isempty, txt(1:end, 7));
%     txt(noData, :) = [];
%     assignin('base', 'txt', txt);
    measurement = txt(2:end,8); % this is the column with the qualitative measurement
    outputSpec = txt(2:end,6);
    control = txt(2:end,7);
    validationIDs = txt(2:end, 1);
    input1 = cell(length(measurement),1); input1(:)={sti1}; % 1 stimulation
    input2 = cell(length(measurement),1); input2(:)={sti2}; % 2 stimulation
    input3 = cell(length(measurement),1); input3(:)={sti3}; % 3 stimulation
    inputCode = cell(length(measurement),1);
    inputCode(:)={stiCode};% this is the column with the code that will alter the input
  
    % set validation threshold change
 % threshold, Ryall et al., 2012 set to 0.001 for sensitivity analysis
    inc = {'Increase'};
    dec = {'Decrease'};
    noc = {'No Change'};
    
        % determination of number of predictions matching references
    numMatching = 0; % number of predictions consistent with the qualitative literature species behavior
    
    % find indices of output species
    outputSpeciesIndex = zeros(1, length(measurement));
    for k = 1:length(outputSpec)
        [~,outputSpeciesIndex(k)] = ismember(outputSpec{k},speciesNames);
    end
    
        [params,y0] = modelParams; %reset the model parameters
        [rpar,tau,ymax,speciesNames]=params{:}; 
        w = rpar(1,:);  
        n = rpar(2,:); 
        EC50 = rpar(3,:);
        
%         % evaluate alternate baseline inputs
%         eval(perturb)
%         if isempty(control{i,1}) == 0 %find out if there is an alternative control condition
%             TF = strncmp(control{i,1},validationIDs,10);
%             eval(inputCode{TF});
%         end
            t1=[];y1=[];
            yInt1 = [];
            yStart = [];

            rpar = [w;n;EC50];
            params = {rpar,tau,ymax,speciesNames};
            [t1,y1] = ode15s(@modelODE, tspan0, y0, options, params);
            tInt1 = [0:1:max(tspan0)];
            yInt1 = interp1(t1,y1,tInt1);
            yStart = yInt1(end,:)';
        
        eval(stiCode);
        
            t2=[];y2=[];
            yInt2 = [];

        rpar = [w;n;EC50];
        params = {rpar,tau,ymax,speciesNames};
        y0 = yInt1(end,:);
        [t2,y2] = ode15s(@modelODE, tspan, y0, options, params);
        tInt2 = [0:1:max(tspan)];
        yInt2 = interp1(t2,y2,tInt2);
        yEnd = yInt2(end,:)';
        
        yT = real([yInt1' yInt2']);
        yTran = horzcat(num2cell(yT));
        sP=speciesNames';
        yTran = horzcat(sP,yTran);
        
    % loop over all validation simulations read in from the excel sheet
    for i = 1:length(inputCode)
        disp(['Validation # ', num2str(i), ' of ',num2str(length(inputCode))]) % write the simulation number to the command line to track loop progress
        i
        activityChange(i,1) = real(yEnd(outputSpeciesIndex(i)))/real(yStart(outputSpeciesIndex(i)));
        temptr = vertcat({num2str(real(yStart(outputSpeciesIndex(i))))},{num2str(real(yEnd(outputSpeciesIndex(i))))});
        rawResult = horzcat(rawResult,temptr);
        
        if activityChange(i) > 1 * thresh % increase
            prediction{i} = 'Increase';
            predChange{i} = num2str(activityChange(i));
            if isequal(inc,measurement(i))
                numMatching = numMatching + 1;
                match{i} = 'yes'; %if the simulation matches the experimental validation put a 1 in the vector
            else
                match{i} = 'no'; %if the simulation does not match put a 0 in the matrix
            end

        elseif activityChange(i) < 1 / thresh % decrease
            prediction{i} = 'Decrease';
            predChange{i} = num2str(activityChange(i));% why num2str
            if isequal(dec,measurement(i))
                numMatching = numMatching + 1;
                match{i} = 'yes';
            else
                match{i} = 'no';
            end
        else % no change
            prediction{i} = 'No Change';
            predChange{i} = num2str(activityChange(i));
            if isequal(noc,measurement(i))
                numMatching = numMatching + 1;
                match{i} = 'yes';
            else
                match{i} = 'no';
            end
        end
    end
    
    
    
    percentMatch = numMatching/length(measurement)*100;
%     for i = 1:length(control)
%          if isempty(control{i,1})==0
%             TF = strncmp(control{i,1},validationIDs,10);
%             controlText{i,1} = strcat(input1(TF),',',input2(TF));
%          end
%     end
    resultChart = {input1,input2,input3,outputSpec,measurement,prediction',predChange',match'};
    resultChart = horzcat(resultChart{:});
    header = { 'input1' ,'input2','input3','output','measurement', 'prediction','predictedChange', 'match'};
    resultChart = vertcat(header, resultChart);
    matchL = strncmp('yes',match,3);
    matchL = double(matchL');
    %l = length(matchL)+1;
    %xlswrite(validationfname,num2str(matchL),['J2:J',num2str(l)]);
    rawResult = horzcat(validationIDs,input1,input2,input3,outputSpec,control,rawResult');
    header2 = {'ID','input1','input2','input3','output','control','yStart','yEnd'};
    rawResult = vertcat(header2,rawResult);
end
        
        
        
        