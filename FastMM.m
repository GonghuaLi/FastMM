%% FastMM for modeling lung cancer
% step 1: read paramaters and setenv
addpath('./bin');
addpath(genpath('./bin/extern'))
binpath = which('FastMM_FVA.m');
binpath = binpath(1:end-13);
if isempty(strfind(getenv('PATH'),binpath))
    if ispc
        setenv('PATH', [getenv('PATH'), ';',binpath]);
    else
        setenv('PATH', [getenv('PATH'), ':',binpath]);
    end
end
pars = parse_parsfile('./pars.txt');
load(pars.consModel);
consmodel = model;
%% using cplex here for fastcore
changeCobraSolver('cplex')
%% Reconstruction by mCADRE the output saved in 'reconstructed_model_mCADRE', need cplex solver
%  mCADREModel = reconstruction_by_mCADRE(consmodel,pars.expressionFile,pars.cutoff);
%% Reconstruction by fastcore the output saved in './out/reconstructed_model.txt',need cplex solver
fastcoreModel = reconstruction_by_fastcore(consmodel,pars.expressionFile,pars.cutoff);
%%  write model for FVA  and gene deletion
modelMatrix = file2cell('./out/reconstructed_model.txt','\t');
modelnames = modelMatrix(1,2:end);
%% single gene KO
GeneKOout = FastMM_singleGeneKO_multi(consmodel,fastcoreModel,pars.numCPU);


%%  double gene KO: this will spend much time
%GeneKOout = FastMM_doubleGeneKO_multi(consmodel,fastcoreModel,pars.numCPU);


%% FVA mulit
[FVAmin,FVAmax] = FastMM_FVA_multi(consmodel,fastcoreModel,pars.numCPU);


%% single metebolite knockout analysis
MetKOout = FastMM_singleMetKO_multi(consmodel,fastcoreModel,pars.numCPU);

%%  double metabolic knockout: this will spend much time
%GeneKOout = FastMM_doubleMetKO_multi(consmodel,fastcoreModel,pars.numCPU);


%% write result
writetxt(['rxns',modelnames;[consmodel.rxns,num2cellstr(FVAmin)]],'./out/FVAmin.txt','\t');
writetxt(['rxns',modelnames;[consmodel.rxns,num2cellstr(FVAmax)]],'./out/FVAmax.txt','\t');
%write geneKO
KOout = zeros(length(consmodel.genes),size(fastcoreModel,2));
for i =1:length(GeneKOout)
    KOout(:,i) = GeneKOout{i}(2:end,3);
end
writetxt(['genes',modelnames;[consmodel.genes,num2cellstr(KOout)]],'./out/singleGeneKO.txt','\t');

%write MetKO
writetxt(['Mets',modelnames;[consmodel.mets,num2cellstr(MetKOout)]],'./out/singleMetKO.txt','\t');



%% Fast mcmc using multithreading saved in ./out/mcmc 
% this function required in linux and should compile locally
% please contact the author if you have some problem to compile mcmc
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
