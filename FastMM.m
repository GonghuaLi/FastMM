%% FastMM
% step 1: read paramaters
addpath('./bin');
pars = parse_parsfile('./pars.txt');
load(pars.consModel);
consmodel = model;
%% Reconstruction by fastcore
fastcoreModel = reconstruction_by_fastcore(consmodel,pars.expressionFile,pars.cutoff);
%%  write model for FVA  and single gene deletion
modelMatrix = file2cell('./out/reconstructed_model.txt','\t');
%%
modelnames = modelMatrix(1,2:end);
matrix = logical(cell2float(modelMatrix(2:end,2:end)));
[N_rxns,N_samples] = size(matrix);
FVAmin = zeros(N_rxns,N_samples);
FVAmax = zeros(N_rxns,N_samples);
GeneKOout = zeros(length(consmodel.genes),N_samples);
outmodels = {};
for i = 1:N_samples
    rmrxns = consmodel.rxns(~matrix(:,i));
    outmodel = removeRxns(consmodel,rmrxns);
    outmodel.c = zeros(length(outmodel.rxns),1);
    outmodel.c(findRxnIDs(outmodel,'biomass_reaction')) = 1;
    outmodels{i} = outmodel;
end
%% multiple threading FVA and single gene knockout
allrxns = repmat(consmodel.rxns,1,N_samples);
CoreNum = pars.numCPU;
if matlabpool('size')<=0 
    matlabpool('open','local',CoreNum);
else  
    disp('matlab pool already started::reopen matlab pool');
    matlabpool close;
    matlabpool('open','local',CoreNum); 
end

parfor i = 1:N_samples
    % FVA
    [IA,IB] = ismember(allrxns(:,i),outmodels{i}.rxns);
    IB(IB==0) = length(outmodels{i}.rxns)+1;
    FVAout_tmp = FastMM_FVA(outmodels{i});
    FVAout_tmp = [FVAout_tmp;0,0];
    FVAout_d = FVAout_tmp(IB,:);
    FVAmin(:,i) = FVAout_d(:,1);
    FVAmax(:,i) = FVAout_d(:,2);
    % knockout
    GeneKOout_tmp = FastMM_singleGeneKO(outmodels{i});
    GeneKOout(:,i) = GeneKOout_tmp(2:end,3);
end
matlabpool close;

writetxt(['rxns',modelnames;[consmodel.rxns,num2cellstr(FVAmin)]],'./out/FVAmin.txt','\t');
writetxt(['rxns',modelnames;[consmodel.rxns,num2cellstr(FVAmax)]],'./out/FVAmax.txt','\t');
writetxt(['genes',modelnames;[consmodel.genes,num2cellstr(GeneKOout)]],'./out/singleGeneKO.txt','\t');
%% Fast mcmc using multithreading saved in ./out/mcmc
if strcmp(pars.ismcmc,'ON')
    if ~isdir('./out/mcmc/')
        mkdir('./out/mcmc/');
    end       
    for i  = 1:N_samples
        outfile = ['./out/mcmc/',modelnames{i},'_mcmc.mat'];
        outmodel = FastMCMC(outmodel,2000,1000,outfile,'needNullS');
    end
end
%%  
