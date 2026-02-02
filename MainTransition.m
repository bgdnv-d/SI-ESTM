function MainTransition(setup)
% preapres data for the centralised system simulations
%
% last change Dmitrii Bogdanov 07.11.2023


setup.endHour = 8760;

% try
%     setup.AvF.Coal(setup.AvF.Coal>0.95)=0.95;
% catch
%     setup.AvF.Coal=repmat([0.95],1,setup.totNumReg);
% end
% try
%     setup.AvF.Gas(setup.AvF.Gas>0.95)=0.95;
% catch
%     setup.AvF.Gas=repmat([0.95],1,setup.totNumReg);
% end
% try
%     setup.AvF.Nuc(setup.AvF.Nuc>0.95)=0.95;
% catch
%     setup.AvF.Nuc=repmat([0.95],1,setup.totNumReg);
% end
% try
%     setup.AvF.Oil(setup.AvF.Oil>0.95)=0.95;
% catch
%     setup.AvF.Oil=repmat([0.95],1,setup.totNumReg);
% end



%% read Base data folder

if setup.readData
    systemData = SetupProjectTransition(setup); % basic data only
    %systemParams = systemData.systemParams;
else
    try
        load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep 'simulation-input_' 'Base' '.mat'])
        systemData.activeElements=activeElements;
        systemData.systemParams=systemParams;
        clear('systemParams');
    catch
        error('Data cannot be found, please read data from files (set ''setup.readData=true'')');
    end
end
tempPath = strsplit(pwd,'\');
systemData.systemParams.modelVersion = tempPath{end};
%systemParams.modelVersion = tempPath{end};



%% run SC for sectors
if ~setup.OvernightFlag
    if setup.SC.Flag
        if setup.Heat.Flag
            projNameSC = ['SC_HE'];
            if exist([setup.rootDir filesep 'projects' filesep projNameSC],'dir') == 0
                setupSC = setup;
                setupSC.projName = 'SC_HE';
                SelfConsTransition(setupSC,systemData);
            end
        else
            projNameSC = ['SC'];
            if exist([setup.rootDir filesep 'projects' filesep projNameSC],'dir') == 0
                setupSC = setup;
                setupSC.projName = 'SC';
                SelfConsTransition(setupSC,systemData);
            end
        end
    end
else
    if setup.SC.Flag
        if setup.Heat.Flag
            projNameSC = ['SC_HE_Overnight_' num2str(setup.startYear)];
            if exist([setup.rootDir filesep 'projects' filesep projNameSC],'dir') == 0
                setupSC = setup;
                setupSC.projName =['SC_HE_Overnight_' num2str(setup.startYear)];
                SelfConsTransition(setupSC,systemData);
            end
        else
            projNameSC = ['SC_Overnight_' num2str(setup.startYear)];
            if exist([setup.rootDir filesep 'projects' filesep projNameSC],'dir') == 0
                setupSC = setup;
                setupSC.projName = ['SC_Overnight_' num2str(setup.startYear)];
                SelfConsTransition(setupSC,systemData);
            end
        end
    end
end


%% run centralised system
%systemParams = systemData.systemParams;
switch setup.projType
    case 'Regions'
        setup.Countries = [[1:length(systemData.systemParams.IndexNodes)]',[1:length(systemData.systemParams.IndexNodes)]'];

        PrepareCountriesScenario(systemData,setup);

    case 'Countries'
        try setup.Countries;
        catch
            error('%s scenario countries list is not presented, please provide list of countries in ''setup.Countries''',setup.projType)
        end
        PrepareCountriesScenario(systemData,setup);

    otherwise
        PrepareScenarioForTransition(systemData,setup);
end

end








