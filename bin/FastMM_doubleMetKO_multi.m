function MetKOout = FastMM_doubleMetKO_multi(consmodel,rxnsmatrix,numcpu,varargin)

matrix = logical(rxnsmatrix);
N_samples = size(matrix,2);
N_Mets = length(consmodel.mets);
MetKOout = zeros(N_Mets*(N_Mets-1)/2+N_Mets+1,N_samples+2);
outmodels = {};
for i = 1:N_samples
    rmrxns = consmodel.rxns(~matrix(:,i));
    %outmodela = removeRxns(consmodel,rmrxns);
    %outmodel.c = zeros(length(outmodel.rxns),1);
    %outmodel.c(findRxnIDs(outmodel,'biomass_reaction')) = 1;
    outmodela = consmodel;
    outmodela.lb(~matrix(:,i)) = 0;
    outmodela.ub(~matrix(:,i)) = 0;
    outmodels{i} = outmodela;
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
    MetKOout_tmp = FastMM_doubleMetKO_par(outmodels{i},thesolver,varargin);
    MetKOout(:,i+2) = MetKOout_tmp(:,3);
end
%matlabpool close;
delete(gcp('nocreate'));
%% indx
ind1 = MetKOout(:,1);
ind2 = MetKOout(:,2);
k = 1;
for i  = 0:N_Mets-1
    ind1(k+1:k+N_Mets-i) = repmat(i,N_Mets-i,1);
    ind2(k+1:k+N_Mets-i) = i+1:N_Mets;
    k = k+N_Mets-i;
end
tmp = ind1(2:N_Mets+1);
ind1(2:N_Mets+1) = ind2(2:N_Mets+1);
ind2(2:N_Mets+1) = tmp;
MetKOout(:,1) = ind1;
MetKOout(:,2) = ind2;

