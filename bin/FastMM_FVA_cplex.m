function flux = FastMM_FVA_cplex(model,varargin)
t = clock;

lb = model.lb;
ub = model.ub;
b = model.b;
c = zeros(length(model.rxns),1);
[bed,ind,val] = sparse_to_csr(model.S);
bed = bed-1;
ind = ind-1;

inmatname = ['fvaIn_cplex_',num2str(ceil(rand(1,1)*1000)),num2str(ceil(t(6)*100)),'.mat'];
outmatname = ['fvaOut_cplex_',num2str(ceil(rand(1,1)*1000)),num2str(ceil(t(6)*100)),'.mat'];

save(inmatname,'lb','ub','b','c','bed','ind','val');

programm_name = 'FVA_cplex_mat';
system([programm_name,' ',inmatname,' ',outmatname]);
load(outmatname);
flux = [minFlux,maxFlux];

delete(inmatname);
delete(outmatname);