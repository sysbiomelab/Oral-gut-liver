% Initialize the COBRA Toolbox
%%
initCobraToolbox()
%% Prepare input data and models
% change directory to where the tutorial is located
%%
MIGRENEindivPath = '/Users/.../Single_modeling/MATLAB_out_HFD';
cd(MIGRENEindivPath);

%% gut model path
modPath = '/Users/.../GutModels';
%%
[~,infoFile,~]=xlsread('/Users/.../table_out/speceis_model_list.xlsx');
PathToModels.name = 'model';


%% Load the GEM models to be joined.
for i=2:size(infoFile,1)
    % model=readCbModel([modPath filesep infoFile{i,1} '.mat']);
    % inputModels{i-1,1}=model;
    model=load([modPath filesep infoFile{i,1}],PathToModels.name);
    inputModels{i-1,1}=model.model;
end

%% Join the first two columns of each sheet by the first column (rxns)


% join rxns tables
model1=inputModels{1,1};
model2=inputModels{2,1};

rxns_tmp = outerjoin(table(model1.rxns,'VariableNames',{'ReactionIDs'}),table(model2.rxns,'VariableNames',{'ReactionIDs'}),'MergeKeys',true);

for i=3:size(inputModels,1)
    model=inputModels{i,1};
    rxns_tmp_full = outerjoin(rxns_tmp,table(model.rxns,'VariableNames',{'ReactionIDs'}),'MergeKeys',true);
end
  
%%
%% run the FBA
% join_withoutO2 = rxns_tmp_full;
join_withO2 = rxns_tmp_full;
% join_noO2_O2 = rxns_tmp_full;

for i=1:size(inputModels,1)
    model=inputModels{i,1};

    %% set obj
    model=changeObjective(model,model.rxns(find(strncmp(model.rxns,'Biomass_Bacteria',7))));

    %% set Ex_rxns to 0

    % Get the number of reactions
    num_rxns = numel(model.rxns);

    % % Iterate through each reaction
    % for j = 1:num_rxns
    %     % Check if the reaction starts with 'msp_' and contains 'Ex'
    %     if startsWith(model.rxns{j}, 'Ex_') % && contains(model.rxns{i}, 'Ex')
    %         % Check the sign of the upper bound (ub)
    % 
    %         model.lb(j) = 0;
    % 
    %     end
    % end

    %% set diet
    dietNumber=4; % 1:high Fibre Plant Based, 2:high Fibre omnivore, 3:high Protein Plant based
              % 4:high protein omnivore, 5:UK average.

    [model]=DietConstrain(model,dietNumber);


    %% end of N constrains

    % % Enable uptake of oxygen (saliva)
    % model=changeRxnBounds(model,'Ex_O2',-10,'l');

    % % Enable no uptake/generate of oxygen (gut)
    % model=changeRxnBounds(model,'Ex_O2',0,'u');

    % FBA
    FBA_2=optimizeCbModel(model,'max');
    data_1(i,2)=FBA_2.f;

    % % FBA_1.v - FBA_2.v
    % gut_oral = FBA_1.v - FBA_2.v;

    % Convert the reaction IDs and flux values into a table
    % model_name_1 = strcat(char(infoFile(i+1,1)), '');
    model_name_2 = strcat(char(infoFile(i+1,1)), '');
    % model_name_3 = strcat(char(infoFile(i+1,1)), '_noO2_O2');

    % result_table_withoutO2 = table(model.rxns, FBA_1.v, 'VariableNames', {'ReactionIDs', model_name_1});
    result_table_withO2 = table(model.rxns, FBA_2.v, 'VariableNames', {'ReactionIDs', model_name_2});

    join_withO2 = outerjoin(join_withO2,result_table_withO2,'MergeKeys',true);


end

%% save predicted flux
cd(MIGRENEindivPath);
writetable(join_withO2, [MIGRENEindivPath filesep 'H_HPD_O_out.xlsx']);



%% run the FVA

% Initialize the COBRA Toolbox
%%
initCobraToolbox()
%% Prepare input data and models
% change directory to where the tutorial is located
%%
MIGRENEindivPath = '/Users/.../Single_modeling/MATLAB_out/FVA_analysis';
cd(MIGRENEindivPath);

%% gut model path
modPath = '/Users/.../GutModels';
%%
[~,infoFile,~] = xlsread('/Users/.../table_out/speceis_model_list.xlsx');
PathToModels.name = 'model';

%% Load the GEM models to be joined.
for i=2:size(infoFile,1)
    model=load([modPath filesep infoFile{i,1}], PathToModels.name);
    inputModels{i-1,1} = model.model;
end

%% Join the first two columns of each sheet by the first column (rxns)
model1 = inputModels{1,1};
model2 = inputModels{2,1};

rxns_tmp = outerjoin(table(model1.rxns, 'VariableNames', {'ReactionIDs'}), table(model2.rxns, 'VariableNames', {'ReactionIDs'}), 'MergeKeys', true);

for i=3:size(inputModels,1)
    model = inputModels{i,1};
    rxns_tmp_full = outerjoin(rxns_tmp, table(model.rxns, 'VariableNames', {'ReactionIDs'}), 'MergeKeys', true);
end

%% FVA
join_withO2 = rxns_tmp_full; % Prepare to join FVA results

for i=1:size(inputModels,1)
    model = inputModels{i,1};

    %% Set objective
    model = changeObjective(model, model.rxns(find(strncmp(model.rxns, 'Biomass_Bacteria', 7))));

    %% Set diet constraints
    dietNumber = 5; % UK average diet for example
    [model] = DietConstrain(model, dietNumber);

    %% Perform FVA
    [minFlux, maxFlux] = fluxVariability(model, 100, 'max', model.rxns);  % FVA over all reactions

    %% Save FVA results
    model_name_min = strcat(char(infoFile(i+1,1)), '_minFlux');
    model_name_max = strcat(char(infoFile(i+1,1)), '_maxFlux');

    % Convert the reaction IDs and min/max flux values into tables
    result_table_min = table(model.rxns, minFlux, 'VariableNames', {'ReactionIDs', model_name_min});
    result_table_max = table(model.rxns, maxFlux, 'VariableNames', {'ReactionIDs', model_name_max});

    % Join FVA results (min and max) to the main table
    join_withO2 = outerjoin(join_withO2, result_table_min, 'MergeKeys', true);
    join_withO2 = outerjoin(join_withO2, result_table_max, 'MergeKeys', true);
end

%% Save the predicted flux variability analysis results
cd(MIGRENEindivPath);
writetable(join_withO2, [MIGRENEindivPath filesep 'gut_in_oral_overlap_FVA_out.xlsx']);

