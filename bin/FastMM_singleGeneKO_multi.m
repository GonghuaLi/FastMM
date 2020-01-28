function GeneKOout = FastMM_singleGeneKO_multi(inmodel,rxnsmatrix,numcpu,varargin)
consmodel = inmodel;
matrix = logical(rxnsmatrix);
N_samples = size(matrix,2);
N_gene = length(consmodel.genes);
GeneKOout = {};
outmodels = {};

isobj = 0;
if isempty(varargin)
    isobj = 0;
else
    obj = regexp(varargin{1},'-f\s+([^\s]+)','tokens');
    if ~isempty(obj)
        obj = obj{1};
        obj = obj{1};
        ot = file2cell(obj,'\t');
        ids = findRxnIDs(consmodel,ot(:,1));
        isobj = 1;
    end
end
if isobj == 1
    consmodel.c = zeros(size(consmodel.c));
    consmodel.c(ids) = 1;
    a =sort(ids);
    [ia,ib] = ismember(ids,a);
end

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
    parpool(CoreNum);
else  
    disp('matlab pool already started::reopen matlab pool');
    delete(gcp('nocreate'));
    parpool(CoreNum); 
end

global CBTLPSOLVER
thesolver = CBTLPSOLVER;
if isobj == 0
    parfor i = 1:N_samples
        GeneKOout_tmp = FastMM_singleGeneKO_par(outmodels{i},thesolver);
        %GeneKOout(:,i+2) = GeneKOout_tmp(:,3);
        GeneKOout{i} = GeneKOout_tmp;
    end
else
    parfor i = 1:N_samples
        GeneKOout_tmp = FastMM_singleGeneKO_par(outmodels{i},thesolver);
        %GeneKOout(:,i+2) = GeneKOout_tmp(:,3);
        GeneKOout{i} = [GeneKOout_tmp(:,1:2),GeneKOout_tmp(:,ib+2)];
    end
end
    
%matlabpool close;
delete(gcp('nocreate'));