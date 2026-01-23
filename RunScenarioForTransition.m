%% RunScenario_ForTransition_fixNew
%
% specifies the sets with active elements of the centralise energy system
%
% last change Dmitrii Bogdanov 07.11.2023


% Create folder

projName = [setup.projType];

try
    setup.MacroMc;
catch
    setup.MacroMc = 0;
end

try
    setup.MacroReg;
catch
    setup.MacroReg = 0;
end

if setup.MacroMc
    projName = [projName '_Mc'];
end

if setup.MacroReg
    projName = [projName '_R'];
end

if setup.OvernightFlag
    projName = [projName '_Overnight'];
end

if setup.SC.Flag
    projName = [projName '_EL_SC'];
end

if setup.Heat.Flag
    projName = [projName '_HE'];
end

if setup.Mobility
    projName = [projName '_TR'];
end

if setup.IndustryFlag
    projName = [projName '_IND'];
end

if setup.GasFlag
    projName = [projName '_GAS'];
end

if setup.DesalinationFlag
    projName = [projName '_DES'];
end

if setup.OnlyFlag
    projName = [projName '_Only'];
end

projName = [projName '_' num2str(costYear)];

if exist([setup.rootDir filesep 'projects' filesep projName],'dir') == 0
    mkdir([setup.rootDir filesep 'projects' filesep projName filesep 'input-data']); % map, modelparam*.xls
    mkdir([setup.rootDir filesep 'projects' filesep projName filesep 'dat-files']);
    mkdir([setup.rootDir filesep 'projects' filesep projName filesep 'output']);
    copyfile([setup.rootDir filesep 'projects' filesep 'Base' filesep 'models'],[setup.rootDir filesep 'projects' filesep projName filesep 'models']);
    fprintf(['\n\nThe directory structure of your new project \n\t >> ' projName ' << \n has successfully been generated.\n\n'])
else
    fprintf(['\n\nThe directory structure of your new project \n\t >> ' projName ' << \n already exists.\n\n'])
end

if isempty(systemParams.gridMat)
    systemParams.GridMax = 0;
end

% Write data files
%     switch setup.projType
%         case 'Area'
switch setup.projType
    case 'Area'
        name = '';
    case {'Regions', 'Countries'}
        name = systemParams.name;
end

[systemParams,results] = AreaTransition(setup,projName,name,systemParams,costYear,activeElements,fossilLimit,fossilGasLimit,SupportFossilsFlag);

systemParams.Instalations(systemParams.Instalations<1)=0;
systemData.systemParams.Instalations = systemParams.Instalations;
systemData.systemParams.TL_DC_LowLimits = round(systemParams.TL_DC_LowLimits);
systemData.systemParams.TL_AC_LowLimits = round(systemParams.TL_AC_LowLimits);

% %% read results data
% if isempty(name)
%     load([setup.rootDir filesep 'projects' filesep projName filesep 'output' filesep 'results']);
% else
%     load([setup.rootDir filesep 'projects' filesep projName filesep 'output' filesep 'results_'  name]);
% end

%% calculate RE share in capcities
capacitiesREN = 0;

if costYear == setup.startYear
    for techI = 1:length(setup.REN_Tech)
        capacitiesREN = capacitiesREN + shiftdim(systemParams.SizeLimits(ismember(systemParams.IndexID,setup.REN_Tech(techI)),1,:),1);
    end
else
    for techI = 1:length(setup.REN_Tech)
        capacitiesREN = capacitiesREN + results.(['OPT_SIZE_' setup.REN_Tech{techI}]);
    end

end

capacitiesFOS = 0;


if costYear == setup.startYear
    for techI = 1:length(setup.FOS_Tech)
        capacitiesFOS = capacitiesFOS + shiftdim(systemParams.SizeLimits(ismember(systemParams.IndexID,setup.FOS_Tech(techI)),1,:),1);
    end
else
    for techI = 1:length(setup.FOS_Tech)
        capacitiesFOS = capacitiesFOS + results.(['OPT_SIZE_' setup.FOS_Tech{techI}]);
    end

end


REN_Share_Lim = capacitiesREN./(capacitiesREN + capacitiesFOS) + setup.REGrowthLim(1+(costYear-setup.startYear)/setup.stepYear); % value from 0 to 1, not in percents!
REN_Share_Lim_Area = sum(capacitiesREN)./sum(capacitiesREN + capacitiesFOS) + setup.REGrowthLim(1+(costYear-setup.startYear)/setup.stepYear);


%% save instalation data
InstData.IndexNodes = systemData.systemParams.IndexNodes;
InstData.IndexYears = systemData.systemParams.IndexYears;
InstData.IndexID = systemData.systemParams.IndexID;
InstData.Instalations = systemData.systemParams.Instalations;

if isempty(name)
    save([setup.rootDir filesep 'projects' filesep projName filesep 'output' filesep 'results_inst_'  projName],'InstData');
else
    save([setup.rootDir filesep 'projects' filesep projName filesep 'output' filesep 'results_inst_'  name],'InstData');
end



%     end

try setup.fixTech.Flag;
catch
    setup.fixTech.Flag = 0;
end

if setup.fixTech.Flag

    for tt = 1:length(setup.fixTech.Tech)

        systemData.systemParams.prevStepCap.(setup.fixTech.Tech{tt}) = results.(['OPT_SIZE_' setup.fixTech.Tech{tt}]);
        if strcmp(setup.fixTech.Tech{tt},'RWIN')
            systemData.systemParams.prevStepCap.RWIO = results.OPT_SIZE_RWIO;
        end
    end
end

clear('results')