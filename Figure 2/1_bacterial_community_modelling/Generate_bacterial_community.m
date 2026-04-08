%% for community GEM
%% gut model path
modPath = 'GutModels'; % download from https://www.microbiomeatlas.org/
%% abundance table
ABUNDANCE='...Tables/Average_Low_SampleID.xlsx'; %

% define a directory to save microbiomeGEM
SAVEDIR='.../Average_model/Average_model_generate'; % 

cd(SAVEDIR)
% integer or a range specified as a 2-element vector of integers
numWorkers=4;
% for some functions such as FBA simulation, you need to install
% cobra toolbox
initCobraToolbox(false)
%%
%% 
%load microbiome (MSP) abundance profile. it could be metagenomics or 16s based 
[abundance,infoFile,~]=xlsread(ABUNDANCE);

% if abundance is cell array，need to run this
if iscell(abundance)
    abundance = cellfun(@(x) str2double(x), abundance, 'UniformOutput', false);
    
    abundance = cellfun(@(x) ifelse(isnan(x), NaN, x), abundance);
    
    abundance = cell2mat(abundance);
end

%name of models
modelList = infoFile(2:end,1);
%name of samples
sampleName = infoFile(1,2:end);
%check the samples,
%remove the MSP name if the abundance of bacteria in all samples are zero
% abundance=abundance(sum(abundance,2)~=0,:);
% modelList=modelList(sum(abundance,2)~=0,:);
abundance_sample0 = abundance(sum(abundance,2)~=0,:);
modelList_sample0 = modelList(sum(abundance,2)~=0,:);

abundance = abundance_sample0;
modelList = modelList_sample0;
% %remove the samples if the there is no bacterial abundance
% sampleName=sampleName(sum(abundance,1)~=0);
% abundance=abundance(:,sum(abundance,1)~=0);
% get the number of bacteria (bacterial richness) in each sample 
temp1=abundance;
temp1(find(temp1>0))=1;
BactrialRichness=table(sampleName',sum(temp1,1)');
BactrialRichness.Properties.VariableNames = {'sampleName','BactrialRichness'};

% give the path where the models are available and the name of model assgined in the .mat files
PathToModels.path=modPath;
PathToModels.name='model';
% generate gut microbiome reaction composition (reaction richness) of all individuals 
richness= RxnRichnessGenerator(modelList,PathToModels,abundance,sampleName);

% generate reaction abundance for all individuals; the function generates both reaction abundance
% and relative reaction abundance
[reactionRelativeAbun, rxnAbunPerSample]= ReactionAbundanceGenerator(modelList,PathToModels,abundance,sampleName);
% generate reactobiome for all individuals
countPerFiveBacteria= ReactobiomeGenerator(modelList,PathToModels,abundance,sampleName);

%%
%community modeling
% give the path where the models are available and the name of model assgined in the .mat files
PathToModels.path=modPath;
PathToModels.name='model';

% define the number of top abundant bacteria for community modeling 
%here we generate communities for top 10 bacteria 
top=16;
thre=[];
for i=1:size(abundance,2)
t1=sort(abundance(:,i),1,'descend');
thre(i,1)=t1(top,1);
abundance(find(abundance(:,i) < thre(i,1)),i)=0;
end
boxplot(thre)
median(thre)

%specify the metabolite ID and exchange reaction for biomass (optional)
biomass.EXrxn='Ex_Biomass';
biomass.mets='cpd11416ee[lu]';
% make a directory to save generated community models
if ~exist([SAVEDIR filesep 'community'],'dir')
mkdir([SAVEDIR filesep 'community']);
end
PathToSave=[SAVEDIR filesep 'community'];
% in report, "one" next to sample name shows that community model has been
% generated for the individuals in PathToSave directory
[report]= MakeCommunity2(modelList,PathToModels,abundance,sampleName,PathToSave,biomass);

%% end

