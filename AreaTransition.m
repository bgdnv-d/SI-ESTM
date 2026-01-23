function [systemParams,results] = AreaTransition(setup,projName,name,systemParams,costYear,activeElements,fossilLimit,fossilGasLimit,SupportFossilsFlag)
%
% FUNCTION main_Area_TransitionNew(setup, projName, name, systemParams, costYear, activeElements, fossilLimit, fossilGasLimit, SupportFossilsFlag)
%
% Main function to run area-level transition calculations for the energy system.
%
%
% INPUT:
%            setup:              Structure that contains all necessary settings and data for processing.
%            projName:           Add!!!!!!!!!!!!!!!!!
%            name:               Add!!!!!!!!!!!!!!!!!
%            systemParams:       Structure with parameters for the energy system.
%            costYear:           Year to which all costs are adjusted.
%            activeElements:     Structure with active elements in the system
%            fossilLimit:        Upper limit on total fossil fuel use.
%            fossilGasLimit:     Upper limit on fossil gas use.
%            SupportFossilsFlag: Flag to allow or disallow support measures for fossil fuels.
%
% OUTPUT:
%            systemParams:       Updated structure with results from the transition calculation.
%
%Dmitrii Bogdanov
%last change 23.07.2025


costYearCheck = 2050;

modFiles.paramsAndSets = 'paramsAndSets.mod';
modFiles.variables = 'varsLossless.mod';
modFiles.objAndCons = 'objAndConsNew.mod';
modFiles.specCons = 'specificCons.mod';
modFiles.grid = ['grid_' name '.mod'];
modFiles.end = 'end.mod';

datFiles.cost = ['cost_' name '.dat'];
datFiles.demand = ['demand_' name '.dat'];
datFiles.feedIn = ['feedin_' name '.dat'];
datFiles.hydro = ['hydro_' name '.dat'];
datFiles.sets = ['sets_' name '.dat'];
datFiles.physical = ['physical_' name '.dat'];
datFiles.TLlengths = ['TL_lengths_' name '.dat'];
datFiles.limits = ['limits_' name '.dat'];
datFiles.misc = ['misc_' name '.dat'];
datFiles.ramping = ['rampingCost_' name '.dat'];


if isempty(name)
    save([setup.rootDir filesep 'projects' filesep projName filesep 'input-data' filesep 'simulation-input_' projName '.mat'],'systemParams','activeElements', 'modFiles','datFiles');
else
    save([setup.rootDir filesep 'projects' filesep projName filesep 'input-data' filesep 'simulation-input_' projName '_' name '.mat'],'systemParams','activeElements', 'modFiles','datFiles');
end

if isempty(name)

    tempFilename=[setup.rootDir filesep 'projects' filesep projName filesep 'output' filesep 'results.mat'];

else
    tempFilename=[setup.rootDir filesep 'projects' filesep projName filesep 'output' filesep 'results_' name '.mat'];

end
if ~isfile(tempFilename)
    disp('Prepearing .dat files for the model')
    WriteCostFile([setup.rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.cost],setup,systemParams,activeElements,costYear);
    WriteSets([setup.rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.sets],setup,systemParams,activeElements);
    WritePhysicalData([setup.rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.physical],setup,systemParams,activeElements,costYear);
    WriteHydro([setup.rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.hydro],setup,systemParams,activeElements,setup.endHour);
    WriteLengths([setup.rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.TLlengths],setup,systemParams);
    WriteMisc([setup.rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.misc],setup,systemParams,activeElements,'',fossilLimit,fossilGasLimit,setup.storePeriod,costYear);
    WriteLimits([setup.rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.limits],setup,systemParams,activeElements);
    WriteGridModwithLosses([setup.rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep modFiles.grid],setup,systemParams.gridMat,systemParams);
    WriteRampingCostData([setup.rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.ramping],setup,systemParams,activeElements);
    WriteDemand([setup.rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.demand],setup,systemParams,activeElements,costYear,setup.endHour);
    WriteFeedIn([setup.rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.feedIn],setup,systemParams,activeElements,setup.endHour);
    fprintf('All .dat files are written \n')
else
    fprintf('Data inputs are not redefined, existing .dat files are used \n')
end


if costYear == 2050;
    qq=1;
end

try
    setup.SolverSolveTool;
catch
    setup.SolverSolveTool = 'mtlb';
end

try
    setup.SolverNumThreads;
catch
    setup.SolverNumThreads = 5;
end

results = RunSOLVER(setup.Solver,projName,setup,setup.rootDir,costYear,SupportFossilsFlag,name,setup.endHour,setup.SolverSolveTool,setup.SolverNumThreads);
if ~(SupportFossilsFlag & (costYear== setup.startYear))
    for i = 1:length(systemParams.IndexID)

        if systemParams.Active(i)==1
            try
                eval(['instCap = results.OPT_SIZE_' systemParams.IndexID{i} ';']);

                systemParams.Instalations(:,find(systemParams.IndexYears==costYear),i)=systemParams.Instalations(:,find(systemParams.IndexYears==costYear),i)+max((instCap-shiftdim(systemParams.SizeLimits(i,1,:),1))',0);

            end
        end
    end
    try
        systemParams.TL_DC_LowLimits = results.OPT_SIZE_TRTL';
        systemParams.TL_AC_LowLimits = results.OPT_SIZE_THAO';
    catch
    end
end
