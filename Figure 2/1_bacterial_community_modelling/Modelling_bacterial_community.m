%% Initialize the COBRA Toolbox
initCobraToolbox();

%% load a community model
%% gut personal model path
modPath = '.../Models/Low_community_exclude'
%% model list name
% gut
[~,infoFile,~]=xlsread('...Tables/Average_Low_SampleID_exclude.xlsx');
PathToModels.name = 'PmergedModel';

%% Load the GEM models to be joined.
for i=2:size(infoFile,1)
    % model=readCbModel([modPath filesep infoFile{i,1} '.mat']);
    % inputModels{i-1,1}=model;
    model=load([modPath filesep infoFile{i,1}],PathToModels.name);
    inputModels{i-1,1}=model.PmergedModel;
end
%% save path
% saving directory path
% gut
saveDir = '.../simu_out';
cd(saveDir)

%% load diet table

% use model derived diet
dietpath = '...Tables/modelDerivedDiet_community.xlsx';
UKaverage = readtable(dietpath);

UKaverage=table2cell(UKaverage);
UKaverage(:,2)=cellstr(num2str(cell2mat(UKaverage(:,2))));
UKaverage(:,3)=cellstr(num2str(cell2mat(UKaverage(:,3))));

% Add 'FeEx', and lb * 24
dietTable = UKaverage;

% make the value in dietTable numeric
% dietTable{:, 2} = str2double(dietTable{:, 2});
for i = 1:size(dietTable, 1)
    % Try to convert the string to a double
    num = str2double(dietTable{i, 2});
    
    % Check if the conversion was successful
    if ~isnan(num)
        dietTable{i, 2} = num;
    else
        % Handle the case where the conversion failed (e.g., display a message or set a default value)
        fprintf('Error converting entry at row %d, column 2\n', i);
        dietTable{i, 2} = 0; % Set a default value, modify as needed
    end
end

dietTable{1:3,2}

% numric for column 3
for i = 1:size(dietTable, 1)
    % Try to convert the string to a double
    num = str2double(dietTable{i, 3});
    
    % Check if the conversion was successful
    if ~isnan(num)
        dietTable{i, 3} = num;
    else
        % Handle the case where the conversion failed (e.g., display a message or set a default value)
        fprintf('Error converting entry at row %d, column 2\n', i);
        dietTable{i, 3} = 0; % Set a default value, modify as needed
    end
end

dietTable{1:3,3}

%% Join the first two columns of each sheet by the first column (rxns)

% join rxns
model1=inputModels{1,1};
model2=inputModels{2,1};

rxns_12 = outerjoin(table(model1.rxns,'VariableNames',{'rxns'}),table(model2.rxns,'VariableNames',{'rxns'}),'MergeKeys',true);

for i=3:size(inputModels,1)
    model=inputModels{i,1};
    rxns_12_full = outerjoin(rxns_12,table(model.rxns,'VariableNames',{'rxns'}),'MergeKeys',true);
end

%% loop start

newInputModels = cell(size(inputModels));

% for j=68:size(inputModels,1) % 169
% for j=1:35
for j=1:size(inputModels,1)
    model=inputModels{j,1};
    
    %% set MSP_xxxx_Ex_rxns (+ to 1000, - to -1000, 0 keep 0)

    % Get the number of reactions
    num_rxns = numel(model.rxns);

    % Iterate through each reaction
    for i = 1:num_rxns
        % Check if the reaction starts with 'msp_' and contains 'Ex'
        if startsWith(model.rxns{i}, 'msp_') % && contains(model.rxns{i}, 'Ex')
            % Check the sign of the upper bound (ub)
            if model.ub(i) > 0
                % If positive, set it to 1000
                model.ub(i) = 1000;
                % If 0, keep it as 0
            end

            % Check the sign of the lower bound (lb)
            if model.lb(i) < 0
                % If negative, set it to -1000
                model.lb(i) = -1000;
                % If 0, keep it as 0
            end
        end
    end
    

    %% set FeEx_ only positive (0 ~ +10000), only out of the systems

    % Get the number of reactions
    num_rxns = numel(model.rxns);

    % Iterate through each reaction
    for i = 1:num_rxns
        % Check if the reaction starts with 'FeEx_'
        if strncmp(model.rxns{i}, 'FeEx_', 5)
            % Check if the reaction is 'FeEX_O2'
            if strcmp(model.rxns{i}, 'FeEx_O2')
                % Set the upper bound (ub) to 0
                model.ub(i) = 0;
            else
                % Update the upper bound (ub) based on the specified conditions
                if model.ub(i) >= 0
                    model.ub(i) = 1000;
                end
            end

            % Update the lower bound (lb) based on the specified conditions
            if model.lb(i) < 0
                % If negative, keep it 0
                model.lb(i) = 0;
                % If 0, keep it as 0
            end
        end
    end

    %% set FoEx all 0
    % Get the number of reactions
    num_rxns = numel(model.rxns);

    % Iterate through each reaction
    for i = 1:num_rxns
        % Check if the reaction starts with 'FoEx_'
        if strncmp(model.rxns{i}, 'FoEx_', 5)
            % Update the upper bound (ub) based on the specified conditions
            if model.ub(i) >= 0
                % If positive, set it to 0, no food taken
                model.ub(i) = 0;
                % If 0, keep it as 0
            end

            % Update the lower bound (lb) based on the specified conditions
            if model.lb(i) <= 0
                % If negative, set it to 0, irreversible
                model.lb(i) = 0;
                % If 0, keep it as 0
            end
        end
    end

    %% --- set the Diet LOWER bounds on FeEx (microbial metabolite pool overall uptake)

    for i = 1:size(dietTable, 1)
        % Find the indices of the reaction in model.rxns that match the first column of dietTable
        rxnIndices = find(strcmp(model.rxns, dietTable{i, 1}));

        % If matches are found, update the lower bound (model.lb) with the UKaverage value from dietTable
        if ~isempty(rxnIndices)
            model.lb(rxnIndices) = (dietTable{i, 3});
        else
            disp(['Reaction not found in model: ', dietTable{i, 1}]);
        end
    end


    % %% oral set oxygen
    % model=changeRxnBounds(model,'FoEx_O2',-240,'l');
    % rxnIndices = find(strcmp(model.rxns, 'FoEx_O2'));
    % % model.lb(rxnIndices)

    %% set 'BiomassAll' ub as 1000 * 24
    desiredValue = 'BiomassAll';
    BiomassIndex = find(strcmp(model.rxns, desiredValue));

    % disp(['Index of ' desiredValue ' in model.rxns: ' num2str(BiomassIndex) ', ub is:' num2str(model.ub(BiomassIndex))]);

    model.ub(BiomassIndex) = 24000;

    % disp(['Index of ' desiredValue ' in model.rxns: ' num2str(BiomassIndex) ', ub is: ' num2str(model.ub(BiomassIndex))]);
    
    % 'FeEx_Biomass'
    desiredValue2 = 'FeEx_Biomass';
    BiomassIndex2 = find(strcmp(model.rxns, desiredValue2));

    % disp(['Index of ' desiredValue ' in model.rxns: ' num2str(BiomassIndex) ', ub is:' num2str(model.ub(BiomassIndex))]);

    model.ub(BiomassIndex2) = 24000;

    %% end

    
    newInputModels{j, 1} = model;
end

%% simulation
join_withoutO2 = rxns_12_full;

FinalModels = cell(size(inputModels));

for j=1:size(newInputModels,1)
    model=newInputModels{j,1};

    % for simulation

    % set objective function
    model=changeObjective(model,model.rxns(find(strncmp(model.rxns,'BiomassAll',10))));

    % FBA
    FBA=optimizeCbModel(model,'max');
    % FBA.f
    % FBA.v

    model.flux = FBA.v;
    FinalModels{j, 1} = model;

    % model sampFBAle name
    model_name_1 = strcat(char(infoFile(j+1,1)), '');

    if isempty(FBA.v)
        disp(['Model not grow: ', model_name_1]);
    else
        %% save the predicted flux
        % Create a table with the specified columns
        resultTable = table(model.rxns, FBA.v, 'VariableNames', {'rxns', 'flux'});
        resultTable.Properties.VariableNames{2} = model_name_1;

        % join table
        join_withoutO2 = outerjoin(join_withoutO2, resultTable, 'Keys', 'rxns', 'MergeKeys', true);

        % Specify the filename
        % model_name_1 = strcat(char(infoFile(j+1,1)), '');
        filename = [model_name_1 '_simu.xlsx'];

        % Write the table to an Excel file
        % writetable(resultTable, [saveDir filesep filename]);
        disp(['save done: ', model_name_1]);

    end
end
%%
% Specify the filename
model_name_1 = strcat(char(infoFile(i+1,1)), '');
filename = [model_name_1 '_gut_simu2.mat'];
variableToSave = newInputModels{1,1}; % This is the data you want to save
save(filename, 'variableToSave');

% Save the variable FinalModels to a .mat file named 'FinalModels.mat'
save('FinalModels.mat', 'FinalModels');

save('FinalModels_LC_High_OGD_exc_add16_avg.mat', 'FinalModels', '-v7.3');


% Write the table to an Excel file

writetable(join_withoutO2, [saveDir filesep 'Top15_exc_add16.xlsx']);

%% save FeEx
% Define the prefix to match
prefix = 'FeEx_';

% Apply the startsWith function to the first column of the table
rowsToKeep = startsWith(join_withoutO2{:,1}, prefix);

% Filter the table to keep only the rows that match the prefix
filteredTable = join_withoutO2(rowsToKeep, :);

% Define the directory and filename where you want to save the filtered table

% Write the filtered table to a new file
writetable(filteredTable, [saveDir filesep 'FeEx_Top15_exc_add16.xlsx']);

%% End

