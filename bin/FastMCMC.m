function outmodel = FastMCMC(inmodel,numPoints,numSteps,outfile,varargin)
%
%   written by Gonghua Li. Ph.D
%   ligonghua@mail.kiz.ac.cn
%   Kunming Institute of Zoology, CAS, PR China.
%   20161018

% if use null(S) to checke the output Points, default is not
% nullS
tic;
t = clock;
cmps = ['mcmc',num2str(ceil(rand(1,1)*1000)),num2str(ceil(t(6)*100)),'.mps'];
ud =0;
if isempty(varargin)
    is_nullS = 0;
else
    is_nullS = 1;
    nullSname = varargin{1};
    if exist(nullSname,'file')
        is_del = 0;
    else
        nullS_varname = ['nullS',num2str(ceil(rand(1,1)*1000)),num2str(ceil(t(6)*100))];
        tmp_nullS = null(full(inmodel.S));
        eval([nullS_varname,' =  tmp_nullS ; clear tmp_nullS'])
        save(nullS_varname,nullS_varname);
        nullSname = [nullS_varname,'.mat'];     
        is_del = 1;
    end
    if length(varargin)>2
        ud =1;
        uTol = varargin{2};
        dTol = varargin{3};
    end
end
% cobra2mps
if length(outfile)>4 && strcmpi(outfile(end-3:end),'.mat')
    cout = outfile;
else
    cout = [cout,'.mat'];
end
cobra2mps(inmodel,cmps);

if is_nullS == 0
    ['mcmc_gurobi_nomod ',cmps,' ',num2str(numPoints),...
        ' ',num2str(numSteps),' ',cout]
    system(['mcmc_gurobi_nomod ',cmps,' ',num2str(numPoints),...
        ' ',num2str(numSteps),' ',cout]);
else
    fprintf('The null sapce of model were used to check points: this will spend much time!\n')
    if ud ==0
        ['mcmc_gurobi_need_nullS ',cmps,' ',nullSname,' ',num2str(numPoints),...
        ' ',num2str(numSteps),' ',cout]
        system(['mcmc_gurobi_need_nullS ',cmps,' ',nullSname,' ',num2str(numPoints),...
        ' ',num2str(numSteps),' ',cout]);
    else
        ['mcmc_gurobi_need_nullS ',cmps,' ',nullSname,' ',num2str(numPoints),...
        ' ',num2str(numSteps),' ',cout,' ',uTol,' ',dTol]
        system(['mcmc_gurobi_need_nullS ',cmps,' ',nullSname,' ',num2str(numPoints),...
        ' ',num2str(numSteps),' ',cout,' ',uTol,' ',dTol]);
    end
        
end
load(cout);
outmodel = inmodel;
outmodel.Points = Points;
if is_nullS == 1 && is_del == 1
    delete(nullSname);
end
delete(cmps);
toc




