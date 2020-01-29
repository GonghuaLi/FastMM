% test
addpath('./bin');
addpath(genpath('./bin/extern/cobratoolbox-3.0.4_base'))
binpath = which('FastMM_FVA.m');
binpath = binpath(1:end-13);
if ~contains(getenv('PATH'),binpath)
    if ispc
        setenv('PATH', [getenv('PATH'), ';',binpath]);
    else
        setenv('PATH', [getenv('PATH'), ':',binpath]);
    end
end
load('./data/consistRecon2_v3.mat');
%% test cobra solver
if changeCobraSolver('gurobi5')
    display('Using gurobi5 solver');
elseif changeCobraSolver('cplex')
    display('Using cplex solver');
end


%%
display('Test FVA...');
try
   flux = FastMM_FVA(model);
    isfva = 1;
catch
    isfva = 0;
end

%%
display('Test knockout...');
try
    flux = FastMM_singleGeneKO(model);
    isko = 1;
catch
    isko = 0;
end

%%
display('Test mcmc...');
try
    outmodel = FastMCMC(model,1000,500,'testouttmp');
    ismcmc  = 1;
    delete('testouttmp.mat');
catch
    ismcmc = 0;
end
%%
if isfva == 1
    display('Fast FVA: passed');
else
    display('Fast FVA: not pass');
    display('Please recompile the Fast FVA and knockout in your computer');
end

if isko == 1
    display('Fast Knockout: passed');
else
    display('Fast Knockout: not pass');
    display('Please recompile the Fast FVA and knockout in your computer');
end

if ismcmc == 1
    display('Fast MCMC: passed');
else
    display('Fast MCMC: not pass');
    display('In version 1.0: only linux Fast MCMC was supported, please wait for version 2.0');
end





