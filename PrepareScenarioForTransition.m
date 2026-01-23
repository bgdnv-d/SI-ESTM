function PrepareScenarioForTransition(systemData,setup)
%
% FUNCTION PrepareScenario_ForTransition(systemData, setup)
%
% Prepares system data for a scenario with transition settings.
%
%
% INPUT:
%            systemData:   Structure with input data for the energy system.
%            setup:        Structure that contains all necessary settings and data for processing.
%
% OUTPUT:
%            -
%
%Dmitrii Bogdanov
%last change 23.07.2025


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

try systemData.Reg1;
catch
    systemData.Reg1 = 1;
end

if isscalar(setup.REGrowthLim)
    setup.REGrowthLim = setup.REGrowthLim * ones(1+(setup.endYear-setup.startYear)/setup.stepYear,1);
end


origInstalations = systemData.systemParams.Instalations;

if ~setup.Heat.Flag
    systemData.systemParams.Instalations(:,:,systemData.activeElements.activeHeatTransformer) = 0;
    systemData.systemParams.Instalations(:,:,ismember(systemData.systemParams.IndexID,'RCSP')) = origInstalations(:,:,ismember(systemData.systemParams.IndexID,'RCSP'));
end

if ~setup.DesalinationFlag
    systemData.systemParams.Instalations(:,:,systemData.activeElements.activeDesalination) = 0;
end

if setup.Mobility
    systemData.systemParams = MobilityInstalation(setup,systemData.systemParams);
else
    systemData.systemParams = MobilityEfficiency(systemData.systemParams);
end


%% Hourly heat pumps COP
systemData.systemParams = HourlyCOP(setup, systemData.systemParams);

systemParams = systemData.systemParams;
%% fix instalation data for fossil fuels
% coal

try setup.ExtendStartCapacityLifetime
catch
    setup.ExtendStartCapacityLifetime = 5;
end


overuseMask = (repmat(systemParams.IndexYears,size(systemParams.IndexID,1),1)<=(setup.startYear-systemParams.Lifetime)); % before fix 23.09.24

OverUsedCapacity = zeros(size(overuseMask'));
for regD = 1:size(systemParams.IndexNodes,1)
    OverUsedCapacity(regD,:) = sum(shiftdim(systemParams.Instalations(regD,:,:),1).*overuseMask');
end




% Create annual scenario
for costYear = setup.startYear:setup.stepYear:setup.endYear

    if costYear <= setup.startYear
        SupportFossilsFlag = setup.SupportFossilsFlag;

    else
        SupportFossilsFlag = false;
    end
    %% no SHHS in main
    systemData.systemParams.Active(ismember(systemData.systemParams.IndexID,'SHHS'))=0;
    systemData.systemParams.Active(ismember(systemData.systemParams.IndexID,'IHHS'))=0;

    systemParams = systemData.systemParams;
    %% no SHHS in main



    activeElements = ActiveComponentsMain(systemParams);

    % temp desalination height and distances - no more used




    systemParams.shareOfSNG(:,systemParams.IndexYears==costYear) = setup.shareOfSNG(:,(costYear-setup.startYear)/setup.stepYear+1);
    %topYear = find(systemParams.IndexYears == costYear);

    instalMask = (repmat(systemParams.IndexYears,size(systemParams.IndexID,1),1)<=costYear).*(repmat(systemParams.IndexYears,size(systemParams.IndexID,1),1)>=(costYear+setup.stepYear-systemParams.Lifetime));


    %% Calculate parameters

    % Electricity demand
    if setup.SC.Flag
        if setup.MacroMc
            if ~setup.OvernightFlag
                if ~setup.Heat.Flag
                    load([setup.rootDir '\projects\Base\input-data\Macro\SelfConsumptionData' '_' num2str(costYear) '_MjR_' num2str(setup.Macro_MajorRegNum) '_SC.mat']);
                else
                    load([setup.rootDir '\projects\Base\input-data\Macro\SelfConsumptionData' '_' num2str(costYear) '_MjR_' num2str(setup.Macro_MajorRegNum) '_SC_HE.mat']);
                end
            else
                load([setup.rootDir filesep 'projects' filesep 'SC_Overnight_' num2str(costYear) filesep 'output' filesep 'SelfConsumptionData' '_' num2str(costYear) '.mat']);
            end
        else
            if ~setup.OvernightFlag
                if ~setup.Heat.Flag
                    load([setup.rootDir filesep 'projects' filesep 'SC' filesep 'output' filesep 'SelfConsumptionData' '_' num2str(costYear) '.mat']);
                else
                    load([setup.rootDir filesep 'projects' filesep 'SC_HE' filesep 'output' filesep 'SelfConsumptionData' '_' num2str(costYear) '.mat']);
                end
            else
                if ~setup.Heat.Flag
                    load([setup.rootDir filesep 'projects' filesep 'SC_Overnight_' num2str(costYear) filesep 'output' filesep 'SelfConsumptionData' '_' num2str(costYear) '.mat']);
                else
                    load([setup.rootDir filesep 'projects' filesep 'SC_HE_Overnight_' num2str(costYear) filesep 'output' filesep 'SelfConsumptionData' '_' num2str(costYear) '.mat']);
                end
                
            end
        end
        if strcmp(setup.projType,'Countries')|strcmp(setup.projType,'Regions')
            if strcmp(setup.CountriesType,'byCountries')
                SCDataFields = fieldnames(SCData);

                for i=1:length(SCDataFields)
                    temp = getfield(SCData,SCDataFields{i});
                    if ~setup.MacroReg
                        SCData = setfield(SCData,SCDataFields{i},temp(:,systemData.systemParams.Regions));
                    else
                        SCData = setfield(SCData,SCDataFields{i},temp(:,systemData.systemParams.IndexNumNodes));
                    end
                end

            end
        end
        systemParams.SCData=SCData;
        for numReg=1:size(systemParams.IndexNodes,1)
            try
                demandProf = load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LELE'),numReg} '_' num2str(costYear)]);
            catch
                demandProf = load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LELE'),numReg}]);
            end
            fn1 = fields(demandProf);
            demandProf.('loadProf')=(demandProf.(fn1{1}))';

            demandProf.loadProf(isnan(demandProf.loadProf))= 0;
            newLoad = demandProf.loadProf/(sum(demandProf.loadProf,2)).*(systemParams.Instalations(numReg,find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'LELE'))));
            if sum(isnan(newLoad))
                warning('Zero power sector demand profile')
            end
            newLoad(isnan(newLoad)) = 0;

            systemParams.ValueLoad(1,:,1,numReg) = round(((newLoad-SCData.SC_demand(:,numReg)'-SCData.ToGrid(:,numReg)'+SCData.FromGrid(:,numReg)').*((newLoad-SCData.SC_demand(:,numReg)'-SCData.ToGrid(:,numReg)'+SCData.FromGrid(:,numReg)')>0))./(1-systemParams.AC_losses(numReg,find(systemParams.IndexYears==costYear))/100));
            systemParams.ValueLoad(1,:,2,numReg) = round(((newLoad-SCData.SC_demand(:,numReg)'-SCData.ToGrid(:,numReg)'+SCData.FromGrid(:,numReg)').*((newLoad-SCData.SC_demand(:,numReg)'-SCData.ToGrid(:,numReg)'+SCData.FromGrid(:,numReg)')>0))./(1-systemParams.AC_losses(numReg,find(systemParams.IndexYears==costYear))/100));


            systemParams.ProsDemand_Gas(numReg,:) = SCData.RNGA_Cons(:,numReg)'.*((SCData.RNGA_Cons(:,numReg)')>0);
            systemParams.ProsDemand_Oil(numReg,:) = SCData.RPET_Cons(:,numReg)'.*((SCData.RPET_Cons(:,numReg)')>0);
            systemParams.ProsDemand_Coal(numReg,:) = SCData.RHAR_Cons(:,numReg)'.*((SCData.RHAR_Cons(:,numReg)')>0);

            systemParams.ProsDemand_BGA(numReg,:) = SCData.RBGA_Cons(:,numReg)'.*((SCData.RBGA_Cons(:,numReg)')>0);
            systemParams.ProsDemand_BMW(numReg,:) = SCData.RBMW_Cons(:,numReg)'.*((SCData.RBMW_Cons(:,numReg)')>0);
            systemParams.ProsDemand_WOO(numReg,:) = SCData.RWOO_Cons(:,numReg)'.*((SCData.RWOO_Cons(:,numReg)')>0);
            systemParams.ProsDemand_WWO(numReg,:) = SCData.RWWO_Cons(:,numReg)'.*((SCData.RWWO_Cons(:,numReg)')>0);

        end
        shareIndHeatElPros = SCData.shareIndHeatElPros(:,1:numReg)'.*(systemParams.SectorsCons(1:numReg,:))./repmat(sum(systemParams.SectorsCons(1:numReg,:),2),1,3);
        setup.Heat.shareIndHeatElPros = shareIndHeatElPros(:,1) + shareIndHeatElPros(:,2);
    else

        SCData.generation = zeros(8760,size(systemParams.IndexNodes,1));
        SCData.RPVO_EL = zeros(8760,size(systemParams.IndexNodes,1));

        SCData.ToGrid = zeros(8760,size(systemParams.IndexNodes,1));
        SCData.FromGrid = zeros(8760,size(systemParams.IndexNodes,1));

        SCData.RBGA_Cons = zeros(8760,size(systemParams.IndexNodes,1));
        SCData.RBMW_Cons = zeros(8760,size(systemParams.IndexNodes,1));
        SCData.RWWO_Cons = zeros(8760,size(systemParams.IndexNodes,1));
        SCData.RWOO_Cons = zeros(8760,size(systemParams.IndexNodes,1));

        SCData.RHAR_Cons = zeros(8760,size(systemParams.IndexNodes,1));
        SCData.RPET_Cons = zeros(8760,size(systemParams.IndexNodes,1));
        SCData.RNGA_Cons = zeros(8760,size(systemParams.IndexNodes,1));

        SCData.SC_demand = zeros(8760,size(systemParams.IndexNodes,1));

        SCData.SC_share = zeros(3,size(systemParams.IndexNodes,1));
        SCData.shareIndHeatElPros = zeros(3,size(systemParams.IndexNodes,1));

        systemParams.SCData=SCData;

        for numReg=1:size(systemParams.IndexNodes,1)

            try
                demandProf = load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LELE'),numReg} '_' num2str(costYear)]);
                fn = fieldnames(demandProf);
                demandProf.loadProf = reshape((demandProf.(fn{1})),1,8760);
            catch
                demandProf = load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LELE'),numReg}]);
                fn = fieldnames(demandProf);
                demandProf.loadProf = reshape((demandProf.(fn{1})),1,8760);
            end

            newLoad = demandProf.loadProf./(sum(demandProf.loadProf,2)).*(systemParams.Instalations(numReg,find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'LELE'))));
            newLoad(isnan(newLoad)) = 0;
            systemParams.ValueLoad(1,:,1,numReg) = round(newLoad./(1-systemParams.AC_losses(numReg,find(systemParams.IndexYears==costYear))/100));
            systemParams.ValueLoad(1,:,2,numReg) = round(newLoad./(1-systemParams.AC_losses(numReg,find(systemParams.IndexYears==costYear))/100));

            systemParams.ProsDemand_Gas(numReg,:) = SCData.RNGA_Cons(:,numReg)'.*((SCData.RNGA_Cons(:,numReg)')>0);
            systemParams.ProsDemand_BGA(numReg,:) = SCData.RBGA_Cons(:,numReg)'.*((SCData.RBGA_Cons(:,numReg)')>0);
            systemParams.ProsDemand_WOO(numReg,:) = SCData.RWOO_Cons(:,numReg)'.*((SCData.RWOO_Cons(:,numReg)')>0);
            systemParams.ProsDemand_WWO(numReg,:) = SCData.RWWO_Cons(:,numReg)'.*((SCData.RWWO_Cons(:,numReg)')>0);

        end
        setup.Heat.shareIndHeatElPros = zeros(0,length((systemParams.IndexNodes)));
    end


    order = [1:length(systemParams.ValueLoad(1,:,1,1))];
    systemParams.ValuesCheck = 0;

    systemParams.Heat.Flag = setup.Heat.Flag;
    try
        systemParams.Heat.shareOfDistrHeat=setup.Heat.shareOfDistrHeat;
        systemParams.Heat.shareOfLowHeatInd=setup.Heat.shareOfLowHeatInd;
        systemParams.Heat.shareOfHighHeatInd=setup.Heat.shareOfHighHeatInd;
        systemParams.Heat.shareIndHeatElPros=setup.Heat.shareIndHeatElPros;
    end
    if setup.Heat.Flag

        systemParams.Heat.shareOfDistrHeat=setup.Heat.shareOfDistrHeat;
        systemParams.Heat.shareOfLowHeatInd=setup.Heat.shareOfLowHeatInd;
        systemParams.Heat.shareOfHighHeatInd=setup.Heat.shareOfHighHeatInd;
        systemParams.Heat.shareIndHeatElPros=setup.Heat.shareIndHeatElPros;

        for numReg=1:size(systemParams.IndexNodes,1)
            try
                tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LHIN'),numReg} '_' num2str(costYear)]);
                tt = fieldnames(tempDem);
                tempDem.IHDprof = tempDem.(tt{1});
                newProfH=round(tempDem.IHDprof(order)/sum(tempDem.IHDprof)*systemParams.Instalations(numReg,find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'LHIN'))));
                newProfH(isnan(newProfH))=0;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHIN')),:,1,numReg) = newProfH;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHIN')),:,2,numReg) = newProfH;
            catch
                tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LHIN'),numReg}]);
                tt = fieldnames(tempDem);
                tempDem.IHDprof = tempDem.(tt{1});
                newProfH=round(tempDem.IHDprof(order)/sum(tempDem.IHDprof)*systemParams.Instalations(numReg,find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'LHIN'))));
                newProfH(isnan(newProfH))=0;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHIN')),:,1,numReg) = newProfH;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHIN')),:,2,numReg) = newProfH;
            end
            try
                tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LHSP'),numReg} '_' num2str(costYear)]);
                tt = fieldnames(tempDem);
                tempDem.SHDprof = tempDem.(tt{1});
                newProfH=round(tempDem.SHDprof(order)/sum(tempDem.SHDprof)*systemParams.Instalations(numReg,find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'LHSP'))));
                newProfH(isnan(newProfH))=0;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHSP')),:,1,numReg) = newProfH;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHSP')),:,2,numReg) = newProfH;
            catch
                tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LHSP'),numReg}]);
                tt = fieldnames(tempDem);
                tempDem.SHDprof = tempDem.(tt{1});
                newProfH=round(tempDem.SHDprof(order)/sum(tempDem.SHDprof)*systemParams.Instalations(numReg,find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'LHSP'))));
                newProfH(isnan(newProfH))=0;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHSP')),:,1,numReg) = newProfH;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHSP')),:,2,numReg) = newProfH;
            end
            try
                tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LHDW'),numReg} '_' num2str(costYear)]);
                tt = fieldnames(tempDem);
                tempDem.DHWprof = tempDem.(tt{1});
                newProfH=round(tempDem.DHWprof(order)/sum(tempDem.DHWprof)*systemParams.Instalations(numReg,find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'LHDW'))));
                newProfH(isnan(newProfH))=0;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHDW')),:,1,numReg) = newProfH;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHDW')),:,2,numReg) = newProfH;
            catch
                tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LHDW'),numReg}]);
                tt = fieldnames(tempDem);
                tempDem.DHWprof = tempDem.(tt{1});
                newProfH=round(tempDem.DHWprof(order)/sum(tempDem.DHWprof)*systemParams.Instalations(numReg,find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'LHDW'))));
                newProfH(isnan(newProfH))=0;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHDW')),:,1,numReg) = newProfH;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHDW')),:,2,numReg) = newProfH;
            end

            try
                tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LHBC'),numReg} '_' num2str(costYear)]);
                tt = fieldnames(tempDem);
                tempDem.BCHprof = tempDem.(tt{1});
                newProfH=round(tempDem.BCHprof(order)/sum(tempDem.BCHprof)*systemParams.Instalations(numReg,find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'LHBC'))));
                newProfH(isnan(newProfH))=0;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHBC')),:,1,numReg) = newProfH;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHBC')),:,2,numReg) = newProfH;
            catch
                tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LHBC'),numReg}]);
                tt = fieldnames(tempDem);
                tempDem.BCHprof = tempDem.(tt{1});
                newProfH=round(tempDem.BCHprof(order)/sum(tempDem.BCHprof)*systemParams.Instalations(numReg,find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'LHBC'))));
                newProfH(isnan(newProfH))=0;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHBC')),:,1,numReg) = newProfH;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHBC')),:,2,numReg) = newProfH;
            end









        end
    else
        for numReg=1:size(systemParams.IndexNodes,1)
            newLoad = systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHIN')),:,1,numReg)*0;
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHIN')),:,1,numReg) = round(newLoad);
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHIN')),:,2,numReg) = round(newLoad);

            newLoad = systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHSP')),:,1,numReg)*0;
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHSP')),:,1,numReg) = round(newLoad);
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHSP')),:,2,numReg) = round(newLoad);

            newLoad = systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHDW')),:,1,numReg)*0;
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHDW')),:,1,numReg) = round(newLoad);
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHDW')),:,2,numReg) = round(newLoad);

            newLoad = systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHBC')),:,1,numReg)*0;
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHBC')),:,1,numReg) = round(newLoad);
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHBC')),:,2,numReg) = round(newLoad);
        end
    end

    for numReg=1:size(systemParams.IndexNodes,1)
        systemParams.ValueLoad(1,:,1,numReg) = systemParams.ValueLoad(1,order,1,numReg);
    end

    if setup.GasFlag
        for numReg=1:size(systemParams.IndexNodes,1)
            newLoad = systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,1,numReg)./(sum(systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,1,numReg),2)./systemParams.Instalations(numReg,find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'LIGA'))));
            newLoad(isnan(newLoad)) = 0;
            if setup.SC.Flag
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,1,numReg) = round(newLoad);
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,2,numReg) = round(newLoad);
            else
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,1,numReg) = round(newLoad);
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,2,numReg) = round(newLoad);
            end
        end
    else
        for numReg=1:size(systemParams.IndexNodes,1)
            newLoad = systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,1,numReg)*0;
            newLoad(isnan(newLoad)) = 0;
            if setup.SC.Flag
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,1,numReg) = round(newLoad);
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,2,numReg) = round(newLoad);
            else
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,1,numReg) = round(newLoad);
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,2,numReg) = round(newLoad);
            end
        end
    end

    % Water demand (with conversion to m3/h)
    if setup.DesalinationFlag
        for numReg=1:size(systemParams.IndexNodes,1)
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LPOT')),:,1,numReg)=round(repmat(systemParams.Instalations(numReg,find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'LPOT')))/24,1,8760));
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LPOT')),:,2,numReg)=round(repmat(systemParams.Instalations(numReg,find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'LPOT')))/24,1,8760));
        end
    else
        for numReg=1:size(systemParams.IndexNodes,1)
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LPOT')),:,1,numReg)=zeros(1,8760);
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LPOT')),:,2,numReg)=zeros(1,8760);
        end
    end

    if setup.Mobility
        MobilityLoadNames = {'LMLL','LMLW','LMLB','LMLM','LMLH','LMRP','LMRF','LMMP','LMMF','LMAP','LMAF'};
        for numTech = 1:length(MobilityLoadNames)

            tech = MobilityLoadNames(numTech);

            for numReg=1:size(systemParams.IndexNodes,1)
                temp=round(systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,tech)),:,1,numReg)./...
                    sum(systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,tech)),:,1,numReg)).*...
                    systemParams.Instalations(numReg,find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,tech))));
                temp(isnan(temp)) = 0;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,tech)),:,1,numReg)= temp;

                temp=round(systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,tech)),:,2,numReg)./...
                    sum(systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,tech)),:,2,numReg)).*...
                    systemParams.Instalations(numReg,find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,tech))));
                temp(isnan(temp)) = 0;
                systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,tech)),:,2,numReg)=temp;
            end
        end
    end





    if ~setup.OvernightFlag





        % Capacities mix efficiency

        efficiency = sum(shiftdim(sum(systemParams.Instalations,1),1).*instalMask'.*systemParams.Efficiency')./sum(shiftdim(sum(systemParams.Instalations,1),1).*instalMask');
        efficiency(isnan(efficiency)) = systemParams.Efficiency(isnan(efficiency),find(systemParams.IndexYears==costYear));



        systemParams.Efficiency = efficiency;
        systemParams.EtaTrans(isnan(systemParams.EtaTrans))=0;
        for i = 1:length(systemParams.IndexIDS)
            if systemParams.Efficiency(ismember(systemParams.IndexID,systemParams.IndexIDS(i)))
                systemParams.EtaStorage(i,:) = systemParams.EtaStorage(i,:).*systemParams.Efficiency(ismember(systemParams.IndexID,systemParams.IndexIDS(i)));
            end
        end

        for i =  1:length(systemParams.IndexIDT)
            if systemParams.Efficiency(ismember(systemParams.IndexID,systemParams.IndexIDT(i)))
                systemParams.EtaTrans(i,1:end-1) = systemParams.EtaTrans(i,1:end-1).*systemParams.Efficiency(ismember(systemParams.IndexID,systemParams.IndexIDT(i)));
            end
        end

        % Installed capacities Online, resource efficiency
        for regD = 1:size(systemParams.IndexNodes,1)

            for tt = 1:length(systemParams.IndexIDR)
                if ~isempty(systemParams.ValueResourceTag{tt,regD})
                    temp=load([setup.path filesep 'projects\Base\input-data' filesep systemParams.ValueResourceTag{tt,regD} '_' num2str(costYear) '.mat']);
                    ss=fieldnames(temp);
                    temp=temp.(ss{1});
                    systemParams.ValueResource(tt,:,1,regD) = temp;
                    systemParams.ValueResource(tt,:,2,regD) = temp;
                    systemParams.ValueResourceTotal(tt,regD) = sum(temp);
                end
            end

            for tt = 1:length(systemParams.IndexIDH)
                if ~isempty(systemParams.ValueHydroDamTag{tt,regD})
                    temp=load([setup.path filesep 'projects\Base\input-data' filesep systemParams.ValueResourceTag{tt,regD} '_' num2str(costYear) '.mat']);
                    ss=fieldnames(temp);
                    temp=temp.(ss{1});
                    systemParams.ValueHydroDam(tt,:,1,regD) = temp;
                    systemParams.ValueHydroDam(tt,:,2,regD) = temp;
                end
            end

            for i = 1:length(systemParams.IndexIDR)

                systemParams.ValueResource(i,:,1,regD) = systemParams.ValueResource(i,:,1,regD) * efficiency(ismember(systemParams.IndexID,systemParams.IndexIDR(i)));
                systemParams.ValueResource(i,:,2,regD) = systemParams.ValueResource(i,:,2,regD) * efficiency(ismember(systemParams.IndexID,systemParams.IndexIDR(i)));

            end

            LowLimit = sum(shiftdim(systemParams.Instalations(regD,:,:),1).*instalMask');

            if((costYear-setup.startYear)<=(setup.ExtendStartCapacityLifetime))
                LowLimit=(LowLimit+round(OverUsedCapacity(regD,:)*(1+(setup.ExtendStartCapacityLifetime - costYear+setup.startYear)/setup.stepYear)/(1+(setup.ExtendStartCapacityLifetime)/setup.stepYear)));
            end


            % To keep all hydro in the system
            LowLimit(ismember(systemParams.IndexID,'RRRI')) = sum(shiftdim(systemParams.Instalations(regD,1:find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'RRRI'))),1));
            LowLimit(ismember(systemParams.IndexID,'HDAM')) = sum(shiftdim(systemParams.Instalations(regD,1:find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'HDAM'))),1));
            LowLimit(ismember(systemParams.IndexID,'SPHS')) = sum(shiftdim(systemParams.Instalations(regD,1:find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'SPHS'))),1));
            LowLimit(ismember(systemParams.IndexID,'IPHS')) = sum(shiftdim(systemParams.Instalations(regD,1:find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,'IPHS'))),1));


            LowLimit = floor(LowLimit);
            UpLimit = systemParams.SizeLimits(:,2,regD)';

            LowLimit(LowLimit>UpLimit) = UpLimit(LowLimit>UpLimit);


            systemParams.SizeLimits(:,1,regD)=floor(LowLimit); % for Morocco instalations are given in MWh

            if ~setup.OvernightFlag   %SupportFossilsFlag
                switch costYear
                    case setup.startYear

                        newValueResourceTotal_RBWM = LowLimit(ismember(systemParams.IndexID,'TMSW'))/...
                            systemParams.EtaTrans(find(ismember(systemParams.IndexIDT,'TMSW')),find(ismember(systemParams.IndexFuels,'FSOL')))*8760*0.9;
                        setup.biomassUseStep.RBMW2015(regD) = min(1,(newValueResourceTotal_RBWM+sum(SCData.RBMW_Cons(:,regD)))/systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBMW'),regD));


                        newValueResourceTotal_RBGA = LowLimit(ismember(systemParams.IndexID,'TCHP'))/...
                            systemParams.EtaTrans(find(ismember(systemParams.IndexIDT,'TCHP')),find(ismember(systemParams.IndexFuels,'FGAS')))*8760*0.8;
                        setup.biomassUseStep.RBGA2015(regD) = min(1,(newValueResourceTotal_RBGA+sum(SCData.RBGA_Cons(:,regD)))/systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBGA'),regD));

                        newValueResourceTotal_RWWO = LowLimit(ismember(systemParams.IndexID,'TBPP'))/...
                            systemParams.EtaTrans(find(ismember(systemParams.IndexIDT,'TBPP')),find(ismember(systemParams.IndexFuels,'FSOL')))*8760*0.8;
                        if setup.Mobility
                            setup.biomassUseStep.RWWO2015(regD) = min(1,(newValueResourceTotal_RWWO+0.5*systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBDS'),regD)+sum(SCData.RWWO_Cons(:,regD)))/systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWWO'),regD));
                        else
                            setup.biomassUseStep.RWWO2015(regD) = min(1,(newValueResourceTotal_RWWO+sum(SCData.RWWO_Cons(:,regD)))/systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWWO'),regD));
                        end
                        newValueResourceTotal_RWOO = LowLimit(ismember(systemParams.IndexID,'TBPP'))/...
                            systemParams.EtaTrans(find(ismember(systemParams.IndexIDT,'TBPP')),find(ismember(systemParams.IndexFuels,'FSOL')))*8760*0.8;
                        if setup.Mobility
                            setup.biomassUseStep.RWOO2015(regD) = min(1,1.5*(newValueResourceTotal_RWOO+0.5*systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBDS'),regD)+sum(SCData.RWOO_Cons(:,regD))+sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHBC'),:,1,regD)))/systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWOO'),regD));
                        else
                            setup.biomassUseStep.RWOO2015(regD) = min(1,1.5*(newValueResourceTotal_RWOO+sum(SCData.RWOO_Cons(:,regD))+sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHBC'),:,1,regD)))/systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWOO'),regD));
                        end
                        setup.biomassUseStep.RWOO2015(regD) = max(setup.biomassUseStep.RWOO2015(regD),0.1);

                        setup.biomassUseStep.RBMW(regD) = setup.biomassUseStep.RBMW2015(regD);
                        setup.biomassUseStep.RBGA(regD) = setup.biomassUseStep.RBGA2015(regD);
                        setup.biomassUseStep.RBGA(regD) = setup.biomassUseStep.RBGA2015(regD);
                        setup.biomassUseStep.RWWO(regD) = setup.biomassUseStep.RWWO2015(regD);
                        setup.biomassUseStep.RWOO(regD) = setup.biomassUseStep.RWOO2015(regD);

                    case setup.startYear+setup.stepYear
                        try setup.biomassUseStep.RBMW2020(regD);
                            setup.biomassUseStep.RBMW(regD) = setup.biomassUseStep.RBMW2020(regD);
                        catch
                            setup.biomassUseStep.RBMW2020(regD) = setup.biomassUseStep.RBMW2015(regD) + 0.33*(1-setup.biomassUseStep.RBMW2015(regD));
                            setup.biomassUseStep.RBMW(regD) = setup.biomassUseStep.RBMW2020(regD);
                        end

                        try setup.biomassUseStep.RBGA2020(regD);
                            setup.biomassUseStep.RBGA(regD) = setup.biomassUseStep.RBGA2020(regD);
                        catch
                            setup.biomassUseStep.RBGA2020(regD) = setup.biomassUseStep.RBGA2015(regD) + 0.33*(1-setup.biomassUseStep.RBGA2015(regD));
                            setup.biomassUseStep.RBGA(regD) = setup.biomassUseStep.RBGA2020(regD);
                        end

                        try setup.biomassUseStep.RWWO2020(regD);
                            setup.biomassUseStep.RWWO(regD) = setup.biomassUseStep.RWWO2020(regD);
                        catch
                            setup.biomassUseStep.RWWO2020(regD) = setup.biomassUseStep.RWWO2015(regD) + 0.33*(1-setup.biomassUseStep.RWWO2015(regD));
                            setup.biomassUseStep.RWWO(regD) = setup.biomassUseStep.RWWO2020(regD);
                        end

                        try setup.biomassUseStep.RWOO2020(regD);
                            setup.biomassUseStep.RWOO(regD) = setup.biomassUseStep.RWOO2020(regD);
                        catch
                            setup.biomassUseStep.RWOO2020(regD) = setup.biomassUseStep.RWOO2015(regD) + 0.33*(1-setup.biomassUseStep.RWOO2015(regD));
                            setup.biomassUseStep.RWOO(regD) = setup.biomassUseStep.RWOO2020(regD);
                        end


                    case setup.startYear+setup.stepYear*2
                        try setup.biomassUseStep.RBMW2025(regD);
                            setup.biomassUseStep.RBMW(regD) = setup.biomassUseStep.RBMW2025(regD);
                        catch
                            setup.biomassUseStep.RBMW2025(regD) = setup.biomassUseStep.RBMW2015(regD) + 0.66*(1-setup.biomassUseStep.RBMW2015(regD));
                            setup.biomassUseStep.RBMW(regD) = setup.biomassUseStep.RBMW2025(regD);
                        end

                        try setup.biomassUseStep.RBGA2025(regD);
                            setup.biomassUseStep.RBGA(regD) = setup.biomassUseStep.RBGA2025(regD);
                        catch
                            setup.biomassUseStep.RBGA2025(regD) = setup.biomassUseStep.RBGA2015(regD) + 0.66*(1-setup.biomassUseStep.RBGA2015(regD));
                            setup.biomassUseStep.RBGA(regD) = setup.biomassUseStep.RBGA2025(regD);
                        end

                        try setup.biomassUseStep.RWWO2025(regD);
                            setup.biomassUseStep.RWWO(regD) = setup.biomassUseStep.RWWO2025(regD);
                        catch
                            setup.biomassUseStep.RWWO2025(regD) = setup.biomassUseStep.RWWO2015(regD) + 0.66*(1-setup.biomassUseStep.RWWO2015(regD));
                            setup.biomassUseStep.RWWO(regD) = setup.biomassUseStep.RWWO2025(regD);
                        end

                        try setup.biomassUseStep.RWOO2025(regD);
                            setup.biomassUseStep.RWOO(regD) = setup.biomassUseStep.RWOO2025(regD);
                        catch
                            setup.biomassUseStep.RWOO2025(regD) = setup.biomassUseStep.RWOO2015(regD) + 0.66*(1-setup.biomassUseStep.RWOO2015(regD));
                            setup.biomassUseStep.RWOO(regD) = setup.biomassUseStep.RWOO2025(regD);
                        end

                    case setup.startYear+setup.stepYear*3
                        try setup.biomassUseStep.RBMW2030(regD);
                            setup.biomassUseStep.RBMW(regD) = setup.biomassUseStep.RBMW2030(regD);
                        catch
                            setup.biomassUseStep.RBMW2030(regD) = setup.biomassUseStep.RBMW2015(regD) + 1*(1-setup.biomassUseStep.RBMW2015(regD));
                            setup.biomassUseStep.RBMW(regD) = setup.biomassUseStep.RBMW2030(regD);
                        end

                        try setup.biomassUseStep.RBGA2030(regD);
                            setup.biomassUseStep.RBGA(regD) = setup.biomassUseStep.RBGA2030(regD);
                        catch
                            setup.biomassUseStep.RBGA2030(regD) = setup.biomassUseStep.RBGA2015(regD) + 1*(1-setup.biomassUseStep.RBGA2015(regD));
                            setup.biomassUseStep.RBGA(regD) = setup.biomassUseStep.RBGA2030(regD);
                        end

                        try setup.biomassUseStep.RWWO2030(regD);
                            setup.biomassUseStep.RWWO(regD) = setup.biomassUseStep.RWWO2030(regD);
                        catch
                            setup.biomassUseStep.RWWO2030(regD) = setup.biomassUseStep.RWWO2015(regD) + 1*(1-setup.biomassUseStep.RWWO2015(regD));
                            setup.biomassUseStep.RWWO(regD) = setup.biomassUseStep.RWWO2030(regD);
                        end

                        try setup.biomassUseStep.RWOO2030(regD);
                            setup.biomassUseStep.RWOO(regD) = setup.biomassUseStep.RWOO2030(regD);
                        catch
                            setup.biomassUseStep.RWOO2030(regD) = setup.biomassUseStep.RWOO2015(regD) + 1*(1-setup.biomassUseStep.RWOO2015(regD));
                            setup.biomassUseStep.RWOO(regD) = setup.biomassUseStep.RWOO2030(regD);
                        end

                end




                systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBMW'),regD) = setup.biomassUseStep.RBMW(regD).*systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBMW'),regD);% - sum(SCData.RBMW_Cons(:,regD));

                systemParams.ValueResource(find(ismember(systemParams.IndexIDR,'RBGA')),:,1,regD) = floor(max(max(SCData.RBGA_Cons(:,regD)),setup.biomassUseStep.RBGA(regD).*systemData.systemParams.ValueResource(find(ismember(systemParams.IndexIDR,'RBGA')),:,1,regD)));% - (SCData.RBGA_Cons(:,regD)')));
                systemParams.ValueResource(find(ismember(systemParams.IndexIDR,'RBGA')),:,2,regD) = floor(max(max(SCData.RBGA_Cons(:,regD)),setup.biomassUseStep.RBGA(regD).*systemData.systemParams.ValueResource(find(ismember(systemParams.IndexIDR,'RBGA')),:,2,regD)));% - (SCData.RBGA_Cons(:,regD)')));
                systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBGA'),regD) = sum(round(systemParams.ValueResource(find(ismember(systemParams.IndexIDR,'RBGA')),:,1,regD)))+1;

                systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWWO'),regD) = setup.biomassUseStep.RWWO(regD).*systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWWO'),regD);% - sum(SCData.RWWO_Cons(:,regD));

                systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWOO'),regD) = max(setup.biomassUseStep.RWOO(regD).*systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWOO'),regD),1.5*sum(SCData.RWOO_Cons(:,regD)));% - sum(SCData.RWOO_Cons(:,regD));


            end

            fossilGasLimit(regD) = (1+SupportFossilsFlag*99)/100*8760*0.95*systemParams.SizeLimits(ismember(systemParams.IndexID,'TOCG'),1,regD)/0.43;

            % RE capacities limit
            capacitiesREN = 0;
            for iii = 1:length(setup.REN_Tech)
                capacitiesREN = capacitiesREN + systemParams.SizeLimits(ismember(systemParams.IndexID,setup.REN_Tech{iii}),1,regD);
            end

            capacitiesFOS = 0;
            for iii = 1:length(setup.FOS_Tech)
                capacitiesFOS = capacitiesFOS + systemParams.SizeLimits(ismember(systemParams.IndexID,setup.FOS_Tech{iii}),1,regD);
            end

            if costYear == setup.startYear
                systemParams.REN_Share_Lim(regD,find(systemParams.IndexYears==costYear)) = 1;capacitiesREN/(capacitiesREN + capacitiesFOS) + setup.REGrowthLim(1+(costYear-setup.startYear)/setup.stepYear); % value from 0 to 1, not in percents!
                systemParams.REN_Share_Lim_Area(find(systemParams.IndexYears==costYear)) = 1;sum(capacitiesREN)/sum(capacitiesREN + capacitiesFOS) + setup.REGrowthLim(1+(costYear-setup.startYear)/setup.stepYear); % value from 0 to 1, not in percents!

            else
                systemParams.REN_Share_Lim(regD,find(systemParams.IndexYears==costYear)) = REN_Share_Lim(regD);
                systemParams.REN_Share_Lim_Area(find(systemParams.IndexYears==costYear)) = REN_Share_Lim_Area;% value from 0 to 1, not in percents!

            end


        end

        fossilGasLimit = 10^13*ones(size(systemParams.IndexNodes))';
        fossilLimit = 10^15;

        if costYear >= setup.fossilBanYear
            fossilLimit = 0;
            fossilGasLimit = zeros(size(systemParams.IndexNodes))';
        end

        if costYear == setup.startYear
            systemParams.SizeLimits(ismember(systemParams.IndexID,systemParams.IndexIDH),2,:)=systemParams.SizeLimits(ismember(systemParams.IndexID,systemParams.IndexIDH),1,:)+1;
            systemParams.SizeLimits(ismember(systemParams.IndexID,systemParams.IndexIDS),2,:)=systemParams.SizeLimits(ismember(systemParams.IndexID,systemParams.IndexIDS),1,:)+1;
            systemParams.SizeLimits(ismember(systemParams.IndexID,systemParams.IndexIDT),2,:)=systemParams.SizeLimits(ismember(systemParams.IndexID,systemParams.IndexIDT),1,:)+1;
            systemParams.SizeLimits(ismember(systemParams.IndexID,'SBAT'),2,systemParams.SizeLimits(ismember(systemParams.IndexID,'SBAT'),2,:) < systemParams.SizeLimits(ismember(systemParams.IndexID,'IBAT'),2,:)) = systemParams.SizeLimits(ismember(systemParams.IndexID,'IBAT'),2,systemParams.SizeLimits(ismember(systemParams.IndexID,'SBAT'),2,:) < systemParams.SizeLimits(ismember(systemParams.IndexID,'IBAT'),2,:));

            systemParams.SizeLimits(ismember(systemParams.IndexID,'TBGD'),2,:) = Inf;
            systemParams.SizeLimits(ismember(systemParams.IndexID,'SBGA'),2,:) = Inf;
            systemParams.SizeLimits(ismember(systemParams.IndexID,'IBGA'),2,:) = Inf;

            systemParams.SizeLimits(ismember(systemParams.IndexID,'TICG'),2,:) = Inf;
            systemParams.SizeLimits(ismember(systemParams.IndexID,'TOCG'),2,:) = Inf;
        end

    else


        if setup.OvernightFlag
            instalMask(ismember(systemParams.IndexID,systemParams.IndexIDR),:) = 0;
            instalMask(ismember(systemParams.IndexID,systemParams.IndexIDT),:) = 0;
            OverUsedCapacity(:,:) = 0;
        end
		
        for regD = 1:size(systemParams.IndexNodes,1)
		    LowLimit = sum(shiftdim(systemParams.Instalations(regD,:,:),1).*instalMask');
    
            if((costYear-setup.startYear)<=(setup.ExtendStartCapacityLifetime))
                LowLimit=(LowLimit+round(OverUsedCapacity(regD,:)*(1+(setup.ExtendStartCapacityLifetime - costYear+setup.startYear)/setup.stepYear)/(1+(setup.ExtendStartCapacityLifetime)/setup.stepYear)));
            end
    
            LowLimit = floor(LowLimit);
    
            systemParams.SizeLimits(ismember(systemParams.IndexID,systemParams.IndexIDS),1,regD)=max(systemParams.SizeLimits(ismember(systemParams.IndexID,systemParams.IndexIDS),1,regD),LowLimit(ismember(systemParams.IndexID,systemParams.IndexIDS))'); 
        end
        
        systemParams.REN_Share_Lim=100*ones(length(systemParams.IndexNodes),length(systemParams.IndexYears));
        systemParams.REN_Share_Lim_Area=100*ones(length(systemParams.IndexYears));
        systemParams.AvF_Coal = 0.4*ones(size(systemParams.IndexNodes))';
        systemParams.AvF_Oil = 0.4*ones(size(systemParams.IndexNodes))';
        systemParams.AvF_Gas = 0.4*ones(size(systemParams.IndexNodes))';
        systemParams.AvF_Nuc = 0.85*ones(size(systemParams.IndexNodes))';

        systemParams.EtaTrans(isnan(systemParams.EtaTrans))=0;
        for i = 1:length(systemParams.IndexIDS)
            if systemParams.Efficiency(ismember(systemParams.IndexID,systemParams.IndexIDS(i)))
                systemParams.EtaStorage(i,:) = systemParams.EtaStorage(i,:).*systemParams.Efficiency(ismember(systemParams.IndexID,systemParams.IndexIDS(i)),systemParams.IndexYears==costYear);
            end
        end

        for i =  1:length(systemParams.IndexIDT)
            if systemParams.Efficiency(ismember(systemParams.IndexID,systemParams.IndexIDT(i)))
                systemParams.EtaTrans(i,1:end-1) = systemParams.EtaTrans(i,1:end-1).*systemParams.Efficiency(ismember(systemParams.IndexID,systemParams.IndexIDT(i)),systemParams.IndexYears==costYear);
            end
        end

        %systemParams.SizeLimits(ismember(systemParams.IndexID,'TNUC'),:,:) = zeros(1,2,numReg);
        %systemParams.SizeLimits(ismember(systemParams.IndexID,'THPP'),:,:) = zeros(1,2,numReg);
        %systemParams.SizeLimits(ismember(systemParams.IndexID,'TICG'),:,:) = zeros(1,2,numReg);

        fossilGasLimit = zeros(size(systemParams.IndexNodes))';
        fossilLimit = 0;

    end


    setup.OnlyFlag = setup.Only.Flag;

    if setup.OnlyFlag
        % NOT ACTIVE!
        % individual sectors simulations and synergy effect evaluation.

        systemParamsOrig = systemParams;

        %Power only

        shareOfSector = setup.Only.Power./(setup.Only.Power+setup.Only.Gas+setup.Only.Desalination);
        for regD=1:size(systemParams.IndexNodes,1)
            systemParams.ValueResource(find(ismember(systemParams.IndexIDR,'RBGA')),:,:,regD) = round(systemParamsOrig.ValueResource(find(ismember(systemParams.IndexIDR,'RBGA')),:,:,regD) * shareOfSector,3);
            systemParams.ValueResource(find(ismember(systemParams.IndexIDR,'RMSW')),:,:,regD) = round(systemParamsOrig.ValueResource(find(ismember(systemParams.IndexIDR,'RMSW')),:,:,regD) * shareOfSector,3);
            systemParams.ValueResource(find(ismember(systemParams.IndexIDR,'RWOO')),:,:,regD) = round(systemParamsOrig.ValueResource(find(ismember(systemParams.IndexIDR,'RWOO')),:,:,regD) * shareOfSector,3);
            systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBGA'),regD) = round(systemParamsOrig.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBGA'),regD) * shareOfSector,3);
            systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBMW'),regD) = round(systemParamsOrig.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBMW'),regD) * shareOfSector,3);
            systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWOO'),regD) = round(systemParamsOrig.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWOO'),regD) * shareOfSector,3);
            systemParams.SizeLimits(:,:,regD) = systemParamsOrig.SizeLimits(:,:,regD) * shareOfSector;
        end
        systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,:,:) = zeros(size(systemParamsOrig.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,:,:)));
        systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LPOT')),:,:,:) = zeros(size(systemParamsOrig.ValueLoad(find(ismember(systemParams.IndexIDL,'LPOT')),:,:,:)));
        systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LELE')),:,:,:) = systemParamsOrig.ValueLoad(find(ismember(systemParams.IndexIDL,'LELE')),:,:,:);

        setup.SC.Flag = true;
        setup.GasFlag = false;
        setup.DesalinationFlag = false;

        RunScenarioForTransition

        %Gas only

        shareOfSector = setup.Only.Gas./(setup.Only.Power+setup.Only.Gas+setup.Only.Desalination);
        for regD=1:size(systemParams.IndexNodes,1)
            systemParams.ValueResource(find(ismember(systemParams.IndexIDR,'RBGA')),:,:,regD) = round(systemParamsOrig.ValueResource(find(ismember(systemParams.IndexIDR,'RBGA')),:,:,regD) * shareOfSector,3);
            systemParams.ValueResource(find(ismember(systemParams.IndexIDR,'RMSW')),:,:,regD) = round(systemParamsOrig.ValueResource(find(ismember(systemParams.IndexIDR,'RMSW')),:,:,regD) * shareOfSector,3);
            systemParams.ValueResource(find(ismember(systemParams.IndexIDR,'RWOO')),:,:,regD) = round(systemParamsOrig.ValueResource(find(ismember(systemParams.IndexIDR,'RWOO')),:,:,regD) * shareOfSector,3);
            systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBGA'),regD) = round(systemParamsOrig.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBGA'),regD) * shareOfSector,3);
            systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBMW'),regD) = round(systemParamsOrig.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBMW'),regD) * shareOfSector,3);
            systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWOO'),regD) = round(systemParamsOrig.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWOO'),regD) * shareOfSector,3);
            systemParams.SizeLimits(:,:,regD) = systemParamsOrig.SizeLimits(:,:,regD) * shareOfSector;
        end
        systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,:,:) = systemParamsOrig.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,:,:);
        systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LPOT')),:,:,:) = zeros(size(systemParamsOrig.ValueLoad(find(ismember(systemParams.IndexIDL,'LPOT')),:,:,:)));
        systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LELE')),:,:,:) = zeros(size(systemParamsOrig.ValueLoad(find(ismember(systemParams.IndexIDL,'LELE')),:,:,:)));

        setup.SC.Flag = false;
        setup.GasFlag = true;
        setup.DesalinationFlag = false;

        RunScenarioForTransition

        %Desalination only

        shareOfSector = setup.Only.Desalination./(setup.Only.Power+setup.Only.Gas+setup.Only.Desalination);
        for regD=1:size(systemParams.IndexNodes,1)
            systemParams.ValueResource(find(ismember(systemParams.IndexIDR,'RBGA')),:,:,regD) = round(systemParamsOrig.ValueResource(find(ismember(systemParams.IndexIDR,'RBGA')),:,:,regD) * shareOfSector,3);
            systemParams.ValueResource(find(ismember(systemParams.IndexIDR,'RMSW')),:,:,regD) = round(systemParamsOrig.ValueResource(find(ismember(systemParams.IndexIDR,'RMSW')),:,:,regD) * shareOfSector,3);
            systemParams.ValueResource(find(ismember(systemParams.IndexIDR,'RWOO')),:,:,regD) = round(systemParamsOrig.ValueResource(find(ismember(systemParams.IndexIDR,'RWOO')),:,:,regD) * shareOfSector,3);
            systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBGA'),regD) = round(systemParamsOrig.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBGA'),regD) * shareOfSector,3);
            systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBMW'),regD) = round(systemParamsOrig.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBMW'),regD) * shareOfSector,3);
            systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWOO'),regD) = round(systemParamsOrig.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWOO'),regD) * shareOfSector,3);
            systemParams.SizeLimits(:,:,regD) = systemParamsOrig.SizeLimits(:,:,regD) * shareOfSector;
        end
        systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,:,:) = zeros(size(systemParamsOrig.ValueLoad(find(ismember(systemParams.IndexIDL,'LIGA')),:,:,:)));
        systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LPOT')),:,:,:) = systemParamsOrig.ValueLoad(find(ismember(systemParams.IndexIDL,'LPOT')),:,:,:);
        systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LELE')),:,:,:) = zeros(size(systemParamsOrig.ValueLoad(find(ismember(systemParams.IndexIDL,'LELE')),:,:,:)));

        setup.SC.Flag = false;
        setup.GasFlag = false;
        setup.DesalinationFlag = true;

        clear('SCData')
        RunScenarioForTransition

    else
        clear('SCData')
        RunScenarioForTransition

    end



end