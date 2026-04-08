% Initialize the COBRA Toolbox
%%
initCobraToolbox()
%% Prepare input data and models
% change directory to where the tutorial is located
%%
MIGRENEindivPath = '/Users/jinyi/Documents/.../Single_modeling/MATLAB_out_HFD';
cd(MIGRENEindivPath);

%% gut model path
modPath = '/Users/jinyi/Documents/.../GutModels';
%%
[~,infoFile,~]=xlsread('/Users/jinyi/Documents/.../OG_coenriched_no_H.xlsx');
PathToModels.name = 'model';


%% Load the GEM models to be joined.
for i=2:size(infoFile,1)
    % model=readCbModel([modPath filesep infoFile{i,1} '.mat']);
    % inputModels{i-1,1}=model;
    model=load([modPath filesep infoFile{i,1}],PathToModels.name);
    inputModels{i-1,1}=model.model;
end

%% Join the first two columns of each sheet by the first column (rxns)


% join rxns
model1=inputModels{1,1};
model2=inputModels{2,1};

rxns_tmp = outerjoin(table(model1.rxns,'VariableNames',{'ReactionIDs'}),table(model2.rxns,'VariableNames',{'ReactionIDs'}),'MergeKeys',true);

for i=3:size(inputModels,1)
    model=inputModels{i,1};
    rxns_tmp_full = outerjoin(rxns_tmp,table(model.rxns,'VariableNames',{'ReactionIDs'}),'MergeKeys',true);
end
  
%%
%% FBA

join_withO2 = rxns_tmp_full;

FinalModels = cell(size(inputModels));

for i=1:size(inputModels,1)
    model=inputModels{i,1};

    %% set obj
    model=changeObjective(model,model.rxns(find(strncmp(model.rxns,'Biomass_Bacteria',7))));

    %% set diet
    dietNumber=5; % 1:high Fibre Plant Based, 2:high Fibre omnivore, 3:high Protein Plant based
              % 4:high protein omnivore, 5:UK average.

    [model]=DietConstrain(model,dietNumber);

    %% FBA
    FBA_2=optimizeCbModel(model,'max');
    data_1(i,2)=FBA_2.f;

    model.flux = FBA_2.v;
    FinalModels{i, 1} = model;

    % Convert the reaction IDs and flux values into a table
    % model_name_1 = strcat(char(infoFile(i+1,1)), '');
    model_name_2 = strcat(char(infoFile(i+1,1)), '');
    % model_name_3 = strcat(char(infoFile(i+1,1)), '_noO2_O2');

    % result_table_withoutO2 = table(model.rxns, FBA_1.v, 'VariableNames', {'ReactionIDs', model_name_1});
    result_table_withO2 = table(model.rxns, FBA_2.v, 'VariableNames', {'ReactionIDs', model_name_2});

    join_withO2 = outerjoin(join_withO2,result_table_withO2,'MergeKeys',true);

end

%% save path
saveDir = '/Users/jinyi/Documents/.../MATLAB_equation';
cd(saveDir)

%% extract equation from each model
equationStrings_name = constructEquations(model, model.rxns, false, false, false, false, false, false);
equationStrings_formular = constructEquations(model, model.rxns, false, false, false, false, true, false);

%% test for one model: table reaction + FBA_flux
Final_model = FinalModels{1,1};
combinedTable = table(equationStrings_name, equationStrings_formular, model.rxns, model.flux, ...
    'VariableNames', {'EquationName', 'EquationFormular', 'rxn', 'Flux'});

% Display the table
head(combinedTable);

% Optionally, save the table to an Excel file
filename = 'HighAvgInc_combinedTable.xlsx';
% writetable(combinedTable, filename);
writetable(combinedTable, [saveDir filesep filename]);

% Confirm saving
disp(['Table saved as ', filename]);

%% loop extract equation
% Assuming FinalModels is a 16x1 cell array of structs
for i = 1:length(FinalModels)
    % Extract the current model
    Final_model = FinalModels{i,1};
    
    % Extract the model name
    modelName = Final_model.modelName;
    
    % Replace spaces or special characters in the model name (if necessary) to make valid variable names
    modelName = matlab.lang.makeValidName(modelName);
    
    % extract equations in each model
    equationStrings_name = constructEquations(Final_model, Final_model.rxns, false, false, false, false, false, false);
    equationStrings_formular = constructEquations(Final_model, Final_model.rxns, false, false, false, false, true, false);

    % Generate the combined table for this model
    combinedTable = table(equationStrings_name, equationStrings_formular, Final_model.rxns, Final_model.flux, ...
        'VariableNames', {'EquationName', 'EquationFormular', 'rxn', 'Flux'});
    
    % Dynamically create a variable name using the model name
    eval([ 'combinedTable_' modelName ' = combinedTable;' ]);

    % Generate the Excel filename dynamically based on the model name
    filename = [saveDir filesep 'combinedTable_' modelName '.xlsx'];
    
    % Save the combined table to an Excel file
    writetable(combinedTable, filename);
    
    % Confirm saving
    disp(['Table for ' modelName ' saved as ', filename]);
end


%% end




