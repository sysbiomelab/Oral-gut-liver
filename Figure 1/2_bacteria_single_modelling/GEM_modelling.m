%% run_FBA_tMSPs.m
% Flux balance analysis (FBA) for selected tMSP GEM models
%
% This script:
% 1. Initializes the COBRA Toolbox
% 2. Loads GEM models listed in an input spreadsheet
% 3. Merges reaction IDs across all models
% 4. Applies diet constraints using MIGRENE
% 5. Runs FBA for each model
% 6. Exports reaction fluxes and model-level FBA summary tables
%
% Requirements:
% - COBRA Toolbox
% - MIGRENE toolbox/functions, including DietConstrain
% - GEM model .mat files
%
% Author: Yi Jin
% Date: 2026-04-08
%
% Notes:
% - Models with infeasible solutions are retained in the output.
% - Their flux vectors are filled with NaN values.
% - Output variable names are sanitized for valid MATLAB table headers.

clear;
clc;

%% ----------------------------- User settings -----------------------------
% Root working directory for this project
projectRoot = 'out';

% Directory containing GEM model .mat files
modelDir = 'GutModels';

% Excel file listing model filenames
modelListFile = 'Tables/tMSPs.xlsx';

% Name of the model variable inside each .mat file
modelVarName = 'model';

% Output files
fluxOutputFile   = fullfile(projectRoot, 'MSP_16_fluxes.xlsx');
summaryOutputFile = fullfile(projectRoot, 'MSP_16_FBA_summary.xlsx');

% Diet selection
% 1 = high-fibre plant-based
% 2 = high-fibre omnivore
% 3 = high-protein plant-based
% 4 = high-protein omnivore
% 5 = UK average
dietNumber = 4;

%% ------------------------- Initialize COBRA ------------------------------
fprintf('Initializing COBRA Toolbox...\n');
initCobraToolbox(false);

%% ------------------------- Read model list -------------------------------
fprintf('Reading model list from: %s\n', modelListFile);
[~, infoFile, ~] = xlsread(modelListFile);

if size(infoFile, 1) < 2
    error('The model list file does not contain any model entries.');
end

modelNames = infoFile(2:end, 1);
nModels = numel(modelNames);

fprintf('Number of models to process: %d\n', nModels);

%% --------------------------- Load models ---------------------------------
fprintf('Loading models...\n');

inputModels = cell(nModels, 1);

for i = 1:nModels
    currentModelFile = modelNames{i};

    modelPath = fullfile(modelDir, currentModelFile);

    if ~isfile(modelPath)
        error('Model file not found: %s', modelPath);
    end

    tmp = load(modelPath, modelVarName);

    if ~isfield(tmp, modelVarName)
        error('Variable "%s" not found in file: %s', modelVarName, modelPath);
    end

    inputModels{i} = tmp.(modelVarName);
end

fprintf('All models loaded successfully.\n');

%% -------------------- Merge reaction IDs across models -------------------
fprintf('Building union of reaction IDs across all models...\n');

if nModels < 2
    error('At least two models are required for reaction merging.');
end

allRxnsTable = table(inputModels{1}.rxns, 'VariableNames', {'ReactionIDs'});

for i = 2:nModels
    currentRxnTable = table(inputModels{i}.rxns, 'VariableNames', {'ReactionIDs'});
    allRxnsTable = outerjoin(allRxnsTable, currentRxnTable, ...
        'Keys', 'ReactionIDs', 'MergeKeys', true);
end

fluxTable = allRxnsTable;

fprintf('Total unique reactions across all models: %d\n', height(fluxTable));

%% ---------------------- Prepare summary output ---------------------------
summaryTable = table( ...
    strings(nModels,1), ...   % ModelName
    zeros(nModels,1), ...     % NumReactions
    zeros(nModels,1), ...     % NumMetabolites
    zeros(nModels,1), ...     % NumGenes
    strings(nModels,1), ...   % BiomassReaction
    zeros(nModels,1), ...     % SolverStatus
    nan(nModels,1), ...       % ObjectiveValue
    strings(nModels,1), ...   % Notes
    'VariableNames', { ...
        'ModelName', ...
        'NumReactions', ...
        'NumMetabolites', ...
        'NumGenes', ...
        'BiomassReaction', ...
        'SolverStatus', ...
        'ObjectiveValue', ...
        'Notes'});

%% ------------------------------ Run FBA ----------------------------------
fprintf('Running FBA under diet %d...\n', dietNumber);

for i = 1:nModels
    model = inputModels{i};
    rawModelName = char(modelNames{i});
    safeModelName = matlab.lang.makeValidName(rawModelName);

    fprintf('\n[%d/%d] Processing model: %s\n', i, nModels, rawModelName);

    % Fill basic summary information
    summaryTable.ModelName(i) = string(rawModelName);
    summaryTable.NumReactions(i) = numel(model.rxns);

    if isfield(model, 'mets')
        summaryTable.NumMetabolites(i) = numel(model.mets);
    else
        summaryTable.NumMetabolites(i) = NaN;
    end

    if isfield(model, 'genes')
        summaryTable.NumGenes(i) = numel(model.genes);
    else
        summaryTable.NumGenes(i) = NaN;
    end

    % ---------------------- Set biomass objective -------------------------
    biomassIdx = find(startsWith(model.rxns, 'Biomass_Bacteria'));

    if isempty(biomassIdx)
        warning('No biomass reaction found for model %s. Filling fluxes with NaN.', rawModelName);

        summaryTable.BiomassReaction(i) = "";
        summaryTable.SolverStatus(i) = -1;
        summaryTable.ObjectiveValue(i) = NaN;
        summaryTable.Notes(i) = "No biomass reaction found";

        fluxVec = nan(numel(model.rxns), 1);
        resultTable = table(model.rxns, fluxVec, ...
            'VariableNames', {'ReactionIDs', safeModelName});
        fluxTable = outerjoin(fluxTable, resultTable, ...
            'Keys', 'ReactionIDs', 'MergeKeys', true);

        continue;
    end

    if numel(biomassIdx) > 1
        warning('Multiple biomass reactions found for model %s. Using the first one.', rawModelName);
        biomassIdx = biomassIdx(1);
        summaryNote = "Multiple biomass reactions found; first used";
    else
        summaryNote = "";
    end

    biomassRxn = model.rxns{biomassIdx};
    summaryTable.BiomassReaction(i) = string(biomassRxn);

    model = changeObjective(model, biomassRxn);

    % ----------------------- Apply diet constraint ------------------------
    try
        model = DietConstrain(model, dietNumber);
    catch ME
        warning('DietConstrain failed for model %s: %s', rawModelName, ME.message);

        summaryTable.SolverStatus(i) = -2;
        summaryTable.ObjectiveValue(i) = NaN;
        summaryTable.Notes(i) = "DietConstrain failed: " + string(ME.message);

        fluxVec = nan(numel(model.rxns), 1);
        resultTable = table(model.rxns, fluxVec, ...
            'VariableNames', {'ReactionIDs', safeModelName});
        fluxTable = outerjoin(fluxTable, resultTable, ...
            'Keys', 'ReactionIDs', 'MergeKeys', true);

        continue;
    end

    % ----------------------------- Run FBA -------------------------------
    try
        sol = optimizeCbModel(model, 'max');
    catch ME
        warning('optimizeCbModel failed for model %s: %s', rawModelName, ME.message);

        summaryTable.SolverStatus(i) = -3;
        summaryTable.ObjectiveValue(i) = NaN;
        summaryTable.Notes(i) = "optimizeCbModel failed: " + string(ME.message);

        fluxVec = nan(numel(model.rxns), 1);
        resultTable = table(model.rxns, fluxVec, ...
            'VariableNames', {'ReactionIDs', safeModelName});
        fluxTable = outerjoin(fluxTable, resultTable, ...
            'Keys', 'ReactionIDs', 'MergeKeys', true);

        continue;
    end

    % ----------------------- Check solution validity ---------------------
    isValidSolution = isfield(sol, 'stat') && ...
                      isfield(sol, 'v') && ...
                      ~isempty(sol.v) && ...
                      numel(sol.v) == numel(model.rxns) && ...
                      sol.stat == 1;

    if isValidSolution
        fluxVec = sol.v;
        objVal = sol.f;

        summaryTable.SolverStatus(i) = sol.stat;
        summaryTable.ObjectiveValue(i) = objVal;

        if strlength(summaryNote) == 0
            summaryTable.Notes(i) = "Optimal solution";
        else
            summaryTable.Notes(i) = summaryNote;
        end

        fprintf('  Status: optimal\n');
        fprintf('  Objective value: %.6f\n', objVal);
    else
        fluxVec = nan(numel(model.rxns), 1);

        if isfield(sol, 'stat') && ~isempty(sol.stat)
            summaryTable.SolverStatus(i) = sol.stat;
            fprintf('  Status: infeasible or non-optimal (stat = %d)\n', sol.stat);
        else
            summaryTable.SolverStatus(i) = NaN;
            fprintf('  Status: unavailable\n');
        end

        if isfield(sol, 'f') && ~isempty(sol.f)
            summaryTable.ObjectiveValue(i) = sol.f;
        else
            summaryTable.ObjectiveValue(i) = NaN;
        end

        if strlength(summaryNote) == 0
            summaryTable.Notes(i) = "No valid flux solution returned";
        else
            summaryTable.Notes(i) = summaryNote + "; no valid flux solution returned";
        end

        warning('FBA infeasible or invalid for model %s under diet %d.', rawModelName, dietNumber);
    end

    % -------------------------- Merge flux table -------------------------
    resultTable = table(model.rxns, fluxVec, ...
        'VariableNames', {'ReactionIDs', safeModelName});

    fluxTable = outerjoin(fluxTable, resultTable, ...
        'Keys', 'ReactionIDs', 'MergeKeys', true);
end

%% ----------------------------- Save output -------------------------------
fprintf('\nSaving flux table to: %s\n', fluxOutputFile);
writetable(fluxTable, fluxOutputFile);

fprintf('Saving summary table to: %s\n', summaryOutputFile);
writetable(summaryTable, summaryOutputFile);

fprintf('\nFBA analysis completed successfully.\n');