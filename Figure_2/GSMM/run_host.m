
%% LOAD MODELS
load('data.mat')


%% SELECT MODEL AND RUN SIMULATION 

model = liver;
targetComp = 'x'; 
metIndicesToRemove = find(endsWith(model.mets, targetComp));
model = removeMets(model, model.mets(metIndicesToRemove)); % remove all mets in x compartment to "open" model

% close all exchanges
exch = model.rxns(findExcRxns(model));
model = setParam(model,'lb',exch,0);

% constrain with media
constraints = readtable('RPMI_humanGEM.csv'); % constraint with RPMI media
model = setParam(model,'lb',constraints.rxn,constraints.LB);
model = setParam(model,'ub',constraints.rxn,constraints.UB);
model = setParam(model,'lb',{'HMR_9404','HMR_9269','HMR_9153','HMR_9151'},-10); % additional sink reactions

% constrain with transcriptomics
% only set this for the liver model (cirrhosis or normal csv file)
% comment out for other tissues
constraints = readtable('data/eflux_cirrhosis.csv'); % constraint with expression
model = setParam(model,'lb',constraints.rxn,constraints.LB);
model = setParam(model,'ub',constraints.rxn,constraints.UB);

T = table(); % will contain all fluxes
count = 1;
fluxes = {0,0.1,0.5,1,2,5,10,20,30,40,50,100}; 
for i=1:length(fluxes)
    model2 = setParam(model,'eq',{'HMR_9073','HMR_9086'},-fluxes{i}); 
    fba = solveLP(model2);
    if fba.f > 0.01
        T.([num2str(i)]) = fba.x;
    end
    count = count + 1;
end