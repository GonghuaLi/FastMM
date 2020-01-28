function [FVAmin,FVAmax] = FastMM_FVA_multi(consmodel,rxnsmatrix,numcpu)

matrix = logical(rxnsmatrix);
[N_rxns,N_samples] = size(matrix);
outmodels = {};
FVAmin = zeros(N_rxns,N_samples);
FVAmax = zeros(N_rxns,N_samples);
for i = 1:N_samples
    rmrxns = consmodel.rxns(~matrix(:,i));
    outmodel = removeRxns(consmodel,rmrxns);
    %outmodel.c = zeros(length(outmodel.rxns),1);
    %outmodel.c(findRxnIDs(outmodel,'biomass_reaction')) = 1;
    outmodels{i} = outmodel;
end

%% multiple threading FVA
%CoreNum = numcpu;
%if matlabpool('size')<=0 
%    matlabpool('open','local',CoreNum);
%else  
%    disp('matlab pool already started::reopen matlab pool');
%    matlabpool close;
%    matlabpool('open','local',CoreNum); 
%end

%multiple threading FVA
% using parpool, you should start matlab : matlab -nodisplay  or matlab
%                do not use -nojvm option
CoreNum = numcpu;
allrxns = repmat(consmodel.rxns,1,N_samples);
if isempty(gcp('nocreate')) 
    parpool('local',CoreNum);
else  
    disp('matlab pool already started::reopen matlab pool');
    delete(gcp('nocreate'));
    parpool('local',CoreNum); 
end
global CBTLPSOLVER
thesolver = CBTLPSOLVER;
parfor i = 1:N_samples
    % FVA
    [IA,IB] = ismember(allrxns(:,i),outmodels{i}.rxns);
    IB(IB==0) = length(outmodels{i}.rxns)+1;
    FVAout_tmp = FastMM_FVA_par(outmodels{i},thesolver);
    FVAout_tmp = [FVAout_tmp;0,0];
    FVAout_d = FVAout_tmp(IB,:);
    FVAmin(:,i) = FVAout_d(:,1);
    FVAmax(:,i) = FVAout_d(:,2);
end
%matlabpool close;
delete(gcp('nocreate'));
