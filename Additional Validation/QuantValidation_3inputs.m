function [percentMatch,activityChange,resultChart,rawResult,matchL,yAll] = QuantValidation_3inputs(validationfname,perturb,thresh,inputlevel)
    inscr = 1; % Run input level screening: 1
    rawResult = {};    

    %parameters and initial conditions
    [params,y0] = modelParams; % this accesses the parameters from model Params (part of the ODE code)
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
    input1 = txt(2:end, 2); %second column, if 2 inputs used will need to alter this code
    input2 = txt(2:end, 3);
    input3 = txt(2:end, 4);
    inputCode = txt(2:end, 5); % this is the column with the code that will alter the input
    if inscr==1
        code1 = strsplit(inputCode{1}, ';');
        for s=1:length(code1)
            if length(strfind(code1{s},'w'))>0
                code1{s} = [code1{s}(1:end-1) num2str(inputlevel)]; % Set input levels
            end
        end
        inputCode{1}=strjoin(code1,';');
    else
    end
    measurement = txt(2:end,8); % this is the column with the qualitative measurement
    expression = raw(2:length(measurement)+1,10); % RNA expression fold-changes
    outputSpec = txt(2:end,6);
    control = txt(2:end,7);
    validationIDs = txt(2:end, 1);
  
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
        
        % evaluate alternate baseline inputs
        eval(perturb)
        
%         if isempty(control{i,1}) == 0 %find out if there is an alternative control condition
%             TF = strncmp(control{i,1},validationIDs,10);
%             eval(inputCode{TF});
%         end
            
            t1=[];y1=[];
            yInt1 = [];
            yStart = [];

            rpar = [w;n;EC50];
            params = {rpar,tau,ymax,speciesNames};
            [t1,y1] = ode15s(@modelODE, tspan, y0, options, params);
            tInt1 = [0:1:10];
            yInt1 = interp1(t1,y1,tInt1);
            yStart = yInt1(end,:)';
        
        eval(inputCode{1});
        
            t2=[];y2=[];
            yInt2 = [];
            yEnd = [];

            rpar = [w;n;EC50];
            params = {rpar,tau,ymax,speciesNames};
            [t2,y2] = ode15s(@modelODE, tspan, y0, options, params);
            tInt2 = [0:1:max(tspan)];
            yInt2 = interp1(t2,y2,tInt2);
            yEnd = yInt2(end,:)';
        
        % End and start activities of all species
        yall = real([yInt1' yInt2']);
%         size(yall)
        yAll = horzcat(num2cell(yall));
        sP = speciesNames';
        yAll = horzcat(sP,yAll);
        
% loop over all validation simulations read in from the excel sheet
    for i = 1:length(inputCode)
        disp(['Validation # ', num2str(i), ' of ',num2str(length(inputCode))]) % write the simulation number to the command line to track loop progress
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
    resultChart = {input1,input2,input3,outputSpec,measurement,prediction',predChange',match',expression};
    resultChart = horzcat(resultChart{:});
    header = { 'input1' ,'input2','input3','output','measurement', 'prediction','predictedChange', 'match', 'expression'};
    resultChart = vertcat(header, resultChart);
    matchL = strncmp('yes',match,3);
    matchL = double(matchL');
    %l = length(matchL)+1;
    %xlswrite(validationfname,num2str(matchL),['J2:J',num2str(l)]);
    rawResult = horzcat(validationIDs,input1,input2,input3,outputSpec,control,rawResult',expression);
    header2 = {'ID','input1','input2','input3','output','control','yStart','yEnd', 'expression'};
    rawResult = vertcat(header2,rawResult);
end
        
        
        
        