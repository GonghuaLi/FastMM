function initCobraToolbox_simplified(updateToolbox)
%      _____   _____   _____   _____     _____     |
%     /  ___| /  _  \ |  _  \ |  _  \   / ___ \    |   COnstraint-Based Reconstruction and Analysis
%     | |     | | | | | |_| | | |_| |  | |___| |   |   The COBRA Toolbox - 2017
%     | |     | | | | |  _  { |  _  /  |  ___  |   |
%     | |___  | |_| | | |_| | | | \ \  | |   | |   |   Documentation:
%     \_____| \_____/ |_____/ |_|  \_\ |_|   |_|   |   http://opencobra.github.io/cobratoolbox
%                                                  |
%
%     initCobraToolbox Initialize COnstraint-Based Reconstruction and Analysis Toolbox
%
%     Defines default solvers and paths, tests SBML io functionality.
%     Function only needs to be called once per installation. Saves paths afer script terminates.
%
%     In addition add either of the following into startup.m (generally in MATLAB_DIRECTORY/toolbox/local/startup.m)
%
%     initCobraToolbox
%           -or-
%     changeCobraSolver('gurobi');
%     changeCobraSolver('gurobi', 'MILP');
%     changeCobraSolver('tomlab_cplex', 'QP');
%     changeCobraSolver('tomlab_cplex', 'MIQP');
%     changeCbMapOutput('svg');
%
%     Maintained by Ronan M.T. Fleming, Sylvain Arreckx, Laurent Heirendt

    % define GLOBAL variables
    global CBTDIR;
    global SOLVERS;
    global OPT_PROB_TYPES;
    global CBT_LP_SOLVER;
    global CBT_MILP_SOLVER;
    global CBT_QP_SOLVER;
    global CBT_MIQP_SOLVER;
    global CBT_NLP_SOLVER;
    global GUROBI_PATH;
    global ILOG_CPLEX_PATH;
    global TOMLAB_PATH;
    global MOSEK_PATH;
    global WAITBAR_TYPE;
    global ENV_VARS;
    global gitBashVersion;
    global CBT_MISSING_REQUIREMENTS_ERROR_ID;

    %if ~exist('updateToolbox','var')
    %    updateToolbox = true;
    % end
    % define a base version of gitBash that is tested
    gitBashVersion = '2.13.3';

     % Set the Requirements Error Message ID:
    CBT_MISSING_REQUIREMENTS_ERROR_ID = 'COBRA:RequirementsNotMet';

    % default waitbar is of type text
    if isempty(WAITBAR_TYPE)
        WAITBAR_TYPE = 1;
    end

    % Linux only - default save path location
    defaultSavePathLocation = '~/pathdef.m';

    % initialize the cell of solvers
    SOLVERS = {};

    % declare that the environment variables have not yet been configured
    ENV_VARS.STATUS = 0;

    % initialize the paths
    if exist('GUROBI_PATH', 'var') ~= 1
        GUROBI_PATH = '';
    end
    if exist('ILOG_CPLEX_PATH', 'var') ~= 1
        ILOG_CPLEX_PATH = '';
    end
    if exist('TOMLAB_PATH', 'var') ~= 1
        TOMLAB_PATH = '';
    end
    if exist('MOSEK_PATH', 'var') ~= 1
        MOSEK_PATH = '';
    end

    % print header
    ENV_VARS.printLevel = true;


    % retrieve the current directory
    currentDir = pwd;

    % define the root path of The COBRA Toolbox and change to it.
    CBTDIR = fileparts(which('initCobraToolbox'));
    cd(CBTDIR);

    % add the external install folder
    addpath(genpath([CBTDIR filesep 'external' filesep 'base' filesep 'install']));

    %And the rdir directory
    addpath(genpath([CBTDIR filesep 'external' filesep 'base' filesep 'utilities' filesep 'rdir']));

    % add the install folder
    addpath(genpath([CBTDIR filesep 'src' filesep 'base' filesep 'install']));



    % set the depth flag if the version of git is higher than 2.10.0
    depthFlag = '--depth=1';


    % change to the root of The COBRA Tooolbox
    cd(CBTDIR);

 

    % temporary disable ssl verification
   
    % check curl
   



    % add the folders of The COBRA Toolbox
    folders = {'tutorials', 'papers', 'binary', 'deprecated', 'src', 'test'};



    % add the root folder
    addpath(CBTDIR);

    % add the external folder
    addpath(genpath([CBTDIR filesep 'external']));

    % add specific subfolders
    for k = 1:length(folders)
        tmpDir = [CBTDIR, filesep, folders{k}];
        if exist(tmpDir, 'dir') == 7
            addpath(genpath(tmpDir));
        end
    end

    %Adapt the mac path depending on the mac version.
    if ismac
        adaptMacPath()
    end

    % add the docs/source/notes folder
    addpath(genpath([CBTDIR filesep 'docs' filesep 'source' filesep 'notes']));

    % print a success message


    % Define default CB map output

    for CbMapOutput = {'svg', 'matlab'}
        CbMapOutputOK = changeCbMapOutput(char(CbMapOutput));
        if CbMapOutputOK
            break;
        end
    end
    if CbMapOutputOK
        if ENV_VARS.printLevel
            fprintf(' set to %s.\n', char(CbMapOutput));
        end
    else
        if ENV_VARS.printLevel
            fprintf('FAILED.\n');
        end
    end

    % define xml test file
    xmlTestModel = 'ecoli_core_model.xml';
    xmlTestFile = 'ecoli_core_model.xml';
    % save the userpath
    originalUserPath = path;

    % Set global LP solution accuracy tolerance
    changeCobraSolverParams('LP', 'optTol', 1e-6);

    % Check that SBML toolbox is installed and accessible
    if ~exist('TranslateSBML', 'file')
        if ENV_VARS.printLevel
            warning('libsbml binaries not on Matlab path: COBRA Toolbox will be unable to read/write SBML files');
        end
    else
        % Test the installation with:
        try
            TranslateSBML(xmlTestFile);
            if ENV_VARS.printLevel
                fprintf(' > TranslateSBML is installed and working properly.\n');
            end
        catch
            if ENV_VARS.printLevel
                warning(' > TranslateSBML is installed but is not working properly.');
                fprintf([' > Try running\n   >> TranslateSBML(\''', strrep(xmlTestFile, '\', '\\'), '\'');\n   in order to debug.\n']);
            end
        end
    end

    if ENV_VARS.printLevel
        fprintf(' > Configuring solver environment variables ...\n');
        configEnvVars(1);
        fprintf('   Done.\n');
        fprintf(' > Checking available solvers and solver interfaces ...');
    end

    % define categories of solvers: LP, MILP, QP, MIQP, NLP
    OPT_PROB_TYPES = {'LP', 'MILP', 'QP', 'MIQP', 'NLP'};

    %Define a set of "use first" solvers, other supported solvers will also be added to the struct.
    %This allows to assign them in any order but keep the most commonly used ones on top of the struct.
    SOLVERS = struct('gurobi',struct(),...
                     'ibm_cplex',struct(),...
                     'tomlab_cplex',struct(),...
                     'glpk',struct(),...
                     'mosek',struct(),...
                     'matlab',struct(),...
                     'cplex',struct());
    % active support - supported solvers
    SOLVERS.cplex_direct.type = {'LP', 'MILP', 'QP', 'MIQP'};
    SOLVERS.dqqMinos.type = {'LP'};
    SOLVERS.glpk.type = {'LP', 'MILP'};
    SOLVERS.gurobi.type = {'LP', 'MILP', 'QP', 'MIQP'};
    SOLVERS.ibm_cplex.type = {'LP', 'MILP', 'QP'};
    SOLVERS.cplex.type = {'LP', 'MILP', 'QP'};
    SOLVERS.matlab.type = {'LP', 'NLP'};
    SOLVERS.mosek.type = {'LP', 'QP'};
    SOLVERS.pdco.type = {'LP', 'QP'};
    SOLVERS.quadMinos.type = {'LP'};
    SOLVERS.tomlab_cplex.type = {'LP', 'MILP', 'QP', 'MIQP'};

    % passive support - solver interfaces
    SOLVERS.qpng.type = {'QP'};
    SOLVERS.tomlab_snopt.type = {'NLP'};

    % legacy solvers
    %SOLVERS.gurobi_mex.type = {'LP', 'MILP', 'QP', 'MIQP'};
    %SOLVERS.lindo_old.type = {'LP'};
    %SOLVERS.lindo_legacy.type = {'LP'};
    SOLVERS.lp_solve.type = {'LP'};
    %SOLVERS.opti.type = {'LP', 'MILP', 'QP', 'MIQP', 'NLP'};

    % definition of category of solvers with active support
    SOLVERS.cplex_direct.categ = 'active';
    SOLVERS.dqqMinos.categ = 'active';
    SOLVERS.glpk.categ = 'active';
    SOLVERS.gurobi.categ = 'active';
    SOLVERS.ibm_cplex.categ = 'active';
    SOLVERS.cplex.categ = 'active';
    SOLVERS.matlab.categ = 'active';
    SOLVERS.mosek.categ = 'active';
    SOLVERS.pdco.categ = 'active';
    SOLVERS.quadMinos.categ = 'active';
    SOLVERS.tomlab_cplex.categ = 'active';

    % definition of category of solvers with passive support
    SOLVERS.qpng.categ = 'passive';
    SOLVERS.tomlab_snopt.categ = 'passive';

    % definition of category of solvers with legacy support
    %SOLVERS.gurobi_mex.categ = 'legacy';
    %SOLVERS.lindo_old.categ = 'legacy';
    %SOLVERS.lindo_legacy.categ = 'legacy';
    SOLVERS.lp_solve.categ = 'legacy';
    %SOLVERS.opti.categ = 'legacy';

    % definition of categories of solvers
    supportedSolversNames = fieldnames(SOLVERS);
    catSolverNames.LP = {}; catSolverNames.MILP = {}; catSolverNames.QP = {}; catSolverNames.MIQP = {}; catSolverNames.NLP = {};
    for i = 1:length(supportedSolversNames)
        SOLVERS.(supportedSolversNames{i}).installed = false;
        SOLVERS.(supportedSolversNames{i}).working = false;
        types = SOLVERS.(supportedSolversNames{i}).type;
        for j = 1:length(types)
            catSolverNames.(types{j}){end + 1} = supportedSolversNames{i};
        end
    end
   
end

function [installed, versionGit] = checkGit()
% Checks if git is installed on the system and throws an error if not
%
% USAGE:
%     versionGit = checkGit();
%
% OUTPUT:
%     installed:      boolean to determine whether git is installed or not
%     versionGit:     version of git installed
%

    global ENV_VARS

    % set the boolean as false (not installed)
    installed = false;

    if ENV_VARS.printLevel
        fprintf(' > Checking if git is installed ... ')
    end

   
   
    if status_gitVersion == 0 && ~isempty(index)

        % determine the version of git
        versionGitStr = result_gitVersion(length(searchStr)+1:end);

        % replace line breaks and white spaces
        versionGitStr = regexprep(versionGitStr(1:7),'\s+','');

        % replace the dots in the version number
        tmp = strrep(versionGitStr, '.', '');

        % convert the string of the version number to a number
        versionGit = str2num(tmp);

        % set the boolean to true
        installed = true;

        if ENV_VARS.printLevel
            fprintf([' Done (version: ' versionGitStr ').\n']);
        end
    else
        if ispc
            fprintf('(not installed).\n');
            installGitBash();
        else
            fprintf(result_gitVersion);
            fprintf(' > Please follow the guidelines on how to install git: https://opencobra.github.io/cobratoolbox/docs/requirements.html.\n');
            error(' > git is not installed.');
        end
    end
end

function [status_curl, result_curl] = checkCurlAndRemote(throwError)
% Checks if curl is installed on the system, can connect to the opencobra URL, and throws an error if not
%
% USAGE:
%     status_curl = checkCurlAndRemote(throwError)
%
% INPUT:
%     throwError:   boolean variable that specifies if an error is thrown or a message is displayed
%

    global ENV_VARS

    if nargin < 1
        throwError = true;
    end

    if ENV_VARS.printLevel
        fprintf(' > Checking if curl is installed ... ')
    end

    origLD = getenv('LD_LIBRARY_PATH');
    newLD = regexprep(getenv('LD_LIBRARY_PATH'), [matlabroot '/bin/' computer('arch') ':'], '');

    % check if curl is properly installed
    [status_curl, result_curl] = system('curl --version');

    if status_curl == 0 && ~isempty(strfind(result_curl, 'curl')) && ~isempty(strfind(result_curl, 'http'))
        if ENV_VARS.printLevel
            fprintf(' Done.\n');
        end
    elseif ((status_curl == 127 || status_curl == 48) && isunix)
        % status_curl of 48 is "An unknown option was passed in to libcurl"
        % status_curl of 127 is a bash/shell error for "command not found"
        % You can get either if there is mismatch between the library
        % libcurl and the curl program, which happens with matlab's
        % distributed libcurl. In order to avoid library mismatch we
        % temporarily fchange LD_LIBRARY_PATH
        setenv('LD_LIBRARY_PATH', newLD);
        [status_curl, result_curl] = system('curl --version');
        setenv('LD_LIBRARY_PATH', origLD);
        if status_curl == 0 && ~isempty(strfind(result_curl, 'curl')) && ~isempty(strfind(result_curl, 'http'))
            if ENV_VARS.printLevel
                fprintf(' Done.\n');
            end
        else
            if throwError
                fprintf(result_curl);
                fprintf(' > Please follow the guidelines on how to install curl: https://opencobra.github.io/cobratoolbox/docs/requirements.html.\n');
                error(' > curl is not installed.');
            else
                if ENV_VARS.printLevel
                    fprintf(' (not installed).\n');
                end
            end
        end
    else
        if throwError
            if ispc
                fprintf('(not installed).\n');
                installGitBash();
            else
                fprintf(result_curl);
                fprintf(' > Please follow the guidelines on how to install curl: https://opencobra.github.io/cobratoolbox/docs/requirements.html.\n');
                error(' > curl is not installed.');
            end
        else
            if ENV_VARS.printLevel
                fprintf(' (not installed).\n');
            end
        end
    end

    if ENV_VARS.printLevel
        fprintf(' > Checking if remote can be reached ... ')
    end

    % check if the remote repository can be reached
    [status_curl, result_curl] = system('curl -s -k --head https://github.com/opencobra/cobratoolbox');

    % check if the URL exists
    if status_curl == 0 && ~isempty(strfind(result_curl, ' 200'))
        if ENV_VARS.printLevel
            fprintf(' Done.\n');
        end
    elseif ((status_curl == 127 || status_curl == 48) && isunix)
        setenv('LD_LIBRARY_PATH', newLD);
        [status_curl, result_curl] = system('curl -s -k --head https://github.com/opencobra/cobratoolbox');
        setenv('LD_LIBRARY_PATH', origLD);
        if status_curl == 0 && ~isempty(strfind(result_curl, ' 200'))
            if ENV_VARS.printLevel
                fprintf(' Done.\n');
            end
        else
            if throwError
                fprintf(result_curl);
                error('The remote repository cannot be reached. Please check your internet connection.');
            else
                if ENV_VARS.printLevel
                    fprintf(' (unsuccessful - no internet connection).\n');
                end
            end
        end
    else
        if throwError
            fprintf(result_curl);
            error('The remote repository cannot be reached. Please check your internet connection.');
        else
            if ENV_VARS.printLevel
                fprintf(' (unsuccessful - no internet connection).\n');
            end
        end
    end
end
