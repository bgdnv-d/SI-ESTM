function systemData = setupProjectTransition(setup)

% SETUPPROJECT(PROJNAME,ROOTDIR) prepares input parameters. Data gets
% collected and stored in appropriate GMPL format. First input argument PROJNAME
% is the identifier for the project; it has to be the same specified in function
% NEWPROJECT. Second argument ROOTDIR specifies the parent directory of
% functions, projects, etc. (your git working directory).
%
%
%Dmitrii Bogdanov
%last change 24.07.2025


rootDir = setup.rootDir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Paths, Files and Names %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modFiles.paramsAndSets = 'paramsAndSets.mod';
modFiles.variables = 'varsLossless.mod';
modFiles.objAndCons = 'objAndCons.mod';
modFiles.specCons = 'specificCons.mod';
modFiles.grid = 'grid.mod';
modFiles.end = 'end.mod';

datFiles.cost = 'cost.dat';
datFiles.demand = 'demand.dat';
datFiles.feedIn = 'feedin.dat';
datFiles.hydro = 'hydro.dat';
datFiles.sets = 'sets.dat';
datFiles.physical = 'physical.dat';
datFiles.TLlengths = 'TL_lengths.dat';
datFiles.limits = 'limits.dat';
datFiles.misc = 'misc.dat';
datFiles.ramping = 'rampingCost';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input data from Excel files %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Modelparams
try
    files.financial = setup.files.financial;
catch
    warning('Modelparams financial file is not defined')
    files.financial = ['Modelparams financial_Base.xlsx'];
end

try
    files.limits=setup.files.limits;
catch
    warning('Modelparams limits file is not defined')
    files.limits=['Modelparams limits_Base.xlsx'];
end

try
    files.mobility=setup.files.mobility;
catch
    warning('Modelparams mobility file is not defined')
    files.mobility=['Modelparams transport_Base.xlsx'];
end

try
    files.physical=setup.files.physical;
catch
    warning('Modelparams physical file is not defined')
    files.physical=['Modelparams physical_Base.xlsx'];
end

try
    files.geographical=setup.files.geographical;
catch
    warning('Modelparams geographical file is not defined')
    files.geographical=['Modelparams geographical_Base.xlsx'];
end
try
    files.instalation=setup.files.instalation;
catch
    warning('Modelparams instalation file is not defined')
    files.instalation=['Modelparams instalation_Base.xlsx'];
end
try
    files.ramping=setup.files.ramping;
catch
    warning('Modelparams ramping file is not defined')
    files.ramping=['Modelparams rampingCost_Base.xlsx'];
end
try
    files.grid=setup.files.grid;
catch
    warning('Grid distances file is not defined')
    files.grid=['grid_distances_Base.xlsx'];
end

% Turn off warnings on missing COM-server
warning('off','MATLAB:xlsread:ActiveX')


systemParams = GetSystemFromExcelTransition(setup,files);

% % Get ramping cost data


% Grid distances
[~,systemParams.TLlength,gridMat,systemParams.TL_DC_LowLimits,systemParams.TL_DC_UpLimits,systemParams.TL_AC_LowLimits,systemParams.TL_AC_UpLimits] = GetGridFromExcel([rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep files.grid]);
systemParams.gridMat=gridMat;
% obtain structure of active system elements
activeElements = ActiveComponents(systemParams);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Matlab vars %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep 'simulation-input_' 'Base' '.mat'],'systemParams','activeElements', 'modFiles','datFiles')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

systemData.systemParams = systemParams;
systemData.activeElements = activeElements;

return

clear all setupProject
