function GeneKOout = FastMM_doubleGeneKO_multi(consmodel,rxnsmatrix,numcpu,varargin)

matrix = logical(rxnsmatrix);
N_samples = size(matrix,2);
N_gene = length(consmodel.genes);
GeneKOout = zeros(N_gene*(N_gene-1)/2+N_gene+1,N_samples+2);
outmodels = {};
for i = 1:N_samples
    rmrxns = consmodel.rxns(~matrix(:,i));
    outmodel = removeRxns(consmodel,rmrxns);
    %outmodel.c = zeros(length(outmodel.rxns),1);
    %outmodel.c(findRxnIDs(outmodel,'biomass_reaction')) = 1;
    outmodels{i} = outmodel;
end

%% multiple threading double gene knockout
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
    GeneKOout_tmp = FastMM_doubleGeneKO_par(outmodels{i},thesolver,varargin);
    GeneKOout(:,i+2) = GeneKOout_tmp(:,3);
end
%matlabpool close;
delete(gcp('nocreate'));
%% indx
ind1 = GeneKOout(:,1);
ind2 = GeneKOout(:,2);
k = 1;
for i  = 0:N_gene-1
    ind1(k+1:k+N_gene-i) = repmat(i,N_gene-i,1);
    ind2(k+1:k+N_gene-i) = i+1:N_gene;
    k = k+N_gene-i;
end
tmp = ind1(2:N_gene+1);
ind1(2:N_gene+1) = ind2(2:N_gene+1);
ind2(2:N_gene+1) = tmp;
GeneKOout(:,1) = ind1;
GeneKOout(:,2) = ind2;

