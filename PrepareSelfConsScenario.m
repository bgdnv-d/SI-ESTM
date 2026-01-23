function PrepareSelfConsScenario(setup,Reg,systemParamsReg,activeElements,produsersTypes)
%
% FUNCTION PrepareSelfConsScenario(setup, Reg, systemParamsReg, activeElements, produsersTypes)
%
% Prepares self-consumption scenario data for a specific region.
%
%
% INPUT:
%            setup:            Structure that contains all necessary settings and data for processing.
%            Reg:              Name or index of the region to prepare.
%            systemParamsReg:  Structure with system parameters specific to the region.
%            activeElements:   Structure with active elements in the system
%            produsersTypes:   Types of producers.
%
% OUTPUT:
%            -
%
%Dmitrii Bogdanov
%last change 23.07.2025


setup.Mobility = 0;


BAT_lim = zeros(length(produsersTypes));

produsersNames = {'RPVR','RPVC','RPVI'};
batteryNames = {'SBAR','SBAC','SBAI'};
interfaceNames = {'IBAR','IBAC','IBAI'};

Bat_max = 10^12*ones(3,1);
WACC = systemParamsReg.WACC_SC;

modFiles.paramsAndSets = 'paramsAndSets.mod';
modFiles.variables = 'varsLossless.mod';
modFiles.objAndCons = 'objAndConsProsumers.mod';
modFiles.specCons = 'specificCons.mod';
modFiles.end = 'end.mod';

fossilShare = 0;
storePeriod = 72;

setup.ModelType = 'SC';

existPV = sum(systemParamsReg.Instalations(1,:,ismember(systemParamsReg.IndexID,'RPVO')))+...
    sum(systemParamsReg.Instalations(1,:,ismember(systemParamsReg.IndexID,'RPVA'))) +...
    sum(systemParamsReg.Instalations(1,:,ismember(systemParamsReg.IndexID,'RPVR')))+...
    sum(systemParamsReg.Instalations(1,:,ismember(systemParamsReg.IndexID,'RPVC')))+...
    sum(systemParamsReg.Instalations(1,:,ismember(systemParamsReg.IndexID,'RPVI')));

systemParamsReg = HourlyCOP(setup, systemParamsReg);

produsersTypesIndeces = [3 1 2];
for i=produsersTypesIndeces % all produsers types

    prosumersShareStep = 1;
    systemParams = systemParamsReg;
    systemParams.EtaTransOrig = systemParams.EtaTrans;
    systemParams.EtaStorageOrig = systemParams.EtaStorage;
    systemParams.Active = ActiveProsumers(systemParams);
    switch i
        case 1
            setup.SC.share_of_SC = setup.SC.share_of_SC_RES;
            existPV_tech = sum(systemParamsReg.Instalations(1,:,ismember(systemParamsReg.IndexID,'RPVR')));
        case 2
            setup.SC.share_of_SC = setup.SC.share_of_SC_COM;
            existPV_tech = sum(systemParamsReg.Instalations(1,:,ismember(systemParamsReg.IndexID,'RPVC')));
        case 3
            setup.SC.share_of_SC = setup.SC.share_of_SC_IND;
            existPV_tech = sum(systemParamsReg.Instalations(1,:,ismember(systemParamsReg.IndexID,'RPVI')));
    end

    for ii = 1:length(systemParams.IndexID)
        year1 = find(ismember(systemParamsReg.IndexYears,setup.startYear));
        shareIndHeatPros = max(0.8,setup.SC.share_of_SC(prosumersShareStep)/100*(1-setup.Heat.shareOfDistrHeat(year1)))*min(1,setup.SC.share_of_SC(prosumersShareStep)/100/(1-setup.Heat.shareOfDistrHeat(year1)));

        if systemParams.Active(ii)==1 & ~sum(ismember({'RPVO','SBAT','IBAT'},systemParams.IndexID{ii}))

            systemParams.Instalations(1,:,ii) = systemParams.Instalations(1,:,ii)*setup.Heat.shareOfIndivHeatSectors(i)/sum(setup.Heat.shareOfIndivHeatSectors);%*shareIndHeatPros;
        end

    end

    for costYear = setup.startYear:setup.stepYear:setup.endYear

        if exist([setup.rootDir filesep 'projects' filesep setup.projName],'dir') == 0
            mkdir([setup.rootDir filesep 'projects' filesep setup.projName filesep 'input-data']); % map, modelparam*.xls
            mkdir([setup.rootDir filesep 'projects' filesep setup.projName filesep 'dat-files']);
            mkdir([setup.rootDir filesep 'projects' filesep setup.projName filesep 'output']);
            copyfile([setup.rootDir filesep 'projects' filesep 'Base' filesep 'models'],[setup.rootDir filesep 'projects' filesep setup.projName filesep 'models']);
            fprintf(['\n\nThe directory structure of your new project \n\n\t >> ' setup.projName ' << \n\n has successfully been generated.\n\n'])
        else
            fprintf(['\n\nThe directory structure of your new project \n\n\t >> ' setup.projName ' << \n\n already exist.\n\n'])
        end


        year = find(ismember(systemParamsReg.IndexYears,costYear));


        systemParams.WACC = WACC(i);

        %% Demand
        try
            fn = [setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParamsReg.ValueLoadTag{ismember(systemParams.IndexIDL,'LELE'),systemParams.IndexNumNodes} '_' num2str(costYear) '.mat'];

            if exist(fn)>1
                demandProf = load(fn);
            else
                if costYear<2050
                    fn = [setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParamsReg.ValueLoadTag{ismember(systemParams.IndexIDL,'LELE'),systemParams.IndexNumNodes} '_' num2str(2020) '.mat'];
                    if exist(fn)>1
                        demandProf = load(fn);
                        warning('No annual profile data, 2020 data used')
                    else
                        fn = [setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParamsReg.ValueLoadTag{ismember(systemParams.IndexIDL,'LELE'),systemParams.IndexNumNodes} '.mat'];
                        demandProf = load(fn);
                        warning('No annual profile data')
                    end
                else
                    fn = [setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParamsReg.ValueLoadTag{ismember(systemParams.IndexIDL,'LELE'),systemParams.IndexNumNodes} '_' num2str(2050) '.mat'];
                    if exist(fn)>1

                        demandProf = load(fn);
                        warning('No annual profile data, 2050 data used')
                    else
                        fn = [setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParamsReg.ValueLoadTag{ismember(systemParams.IndexIDL,'LELE'),systemParams.IndexNumNodes} '.mat'];
                        demandProf = load(fn);
                        warning('No annual profile data')
                    end
                end
            end

        catch
            warning('No annual profile data')
            demandProf = load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParamsReg.ValueLoadTag{ismember(systemParams.IndexIDL,'LELE'),systemParams.IndexNumNodes}]);
        end
        fn1 = fields(demandProf);
        demandProf.('loadProf')=(demandProf.(fn1{1}))';

        systemParams.shareOfSectorEl = systemParamsReg.SectorsCons(1,i)/sum(systemParamsReg.SectorsCons(1,:));
        systemParams.ValueLoad(find(ismember(systemParamsReg.IndexIDL,'LELE')),:,1,1) = demandProf.loadProf/(sum(demandProf.loadProf,2)./systemParamsReg.Instalations(1,find(systemParamsReg.IndexYears==costYear),find(ismember(systemParamsReg.IndexID,'LELE'))))*...
            systemParamsReg.SectorsCons(1,i)/sum(systemParamsReg.SectorsCons(1,:))*setup.SC.share_of_SC(prosumersShareStep)/100;
        systemParams.ValueLoad(find(ismember(systemParamsReg.IndexIDL,'LELE')),:,2,1) = demandProf.loadProf./(sum(demandProf.loadProf,2)./systemParamsReg.Instalations(1,find(systemParamsReg.IndexYears==costYear),find(ismember(systemParamsReg.IndexID,'LELE'))))*...
            systemParamsReg.SectorsCons(1,i)/sum(systemParamsReg.SectorsCons(1,:))*setup.SC.share_of_SC(prosumersShareStep)/100;
        systemParams.ValueLoad(find(ismember(systemParamsReg.IndexIDL,'LIGA')),:,1,1)=systemParams.ValueLoad(find(ismember(systemParamsReg.IndexIDL,'LELE')),:,1,1)*0;
        systemParams.ValueLoad(find(ismember(systemParamsReg.IndexIDL,'LIGA')),:,2,1)=systemParams.ValueLoad(find(ismember(systemParamsReg.IndexIDL,'LELE')),:,2,1)*0;

        if setup.Heat.Flag
            try
                tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LHIN'),systemParams.IndexNumNodes} '_' num2str(costYear)]);
            catch
                tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LHIN'),systemParams.IndexNumNodes}]);
            end
            tt = fieldnames(tempDem);
            tempDem.IHDprof = tempDem.(tt{1});
            tempDem = round(tempDem.IHDprof/sum(tempDem.IHDprof)*systemParamsReg.Instalations(1,find(systemParamsReg.IndexYears==costYear),find(ismember(systemParamsReg.IndexID,'LHIN'))));
            tempDem(isnan(tempDem))=0;
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHIN')),:,1,1) = tempDem;
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHIN')),:,2,1) = tempDem;%tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep 'SHDprof_' num2str(systemParams.IndexNumNodes) '_' num2str(costYear)]);

            try
                tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LHSP'),systemParams.IndexNumNodes} '_' num2str(costYear)]);
            catch
                tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LHSP'),systemParams.IndexNumNodes}]);
            end
            tt = fieldnames(tempDem);
            tempDem.SHDprof = tempDem.(tt{1});
            tempDem = round(tempDem.SHDprof/sum(tempDem.SHDprof)*systemParamsReg.Instalations(1,find(systemParamsReg.IndexYears==costYear),find(ismember(systemParamsReg.IndexID,'LHSP'))));
            tempDem(isnan(tempDem))=0;
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHSP')),:,1,1) = tempDem;
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHSP')),:,2,1) = tempDem;

            try
                tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LHDW'),systemParams.IndexNumNodes} '_' num2str(costYear)]);
            catch
                tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LHDW'),systemParams.IndexNumNodes}]);
            end
            tt = fieldnames(tempDem);
            tempDem.DHWprof = tempDem.(tt{1});
            tempDem = round(tempDem.DHWprof/sum(tempDem.DHWprof)*systemParamsReg.Instalations(1,find(systemParamsReg.IndexYears==costYear),find(ismember(systemParamsReg.IndexID,'LHDW'))));
            tempDem(isnan(tempDem))=0;
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHDW')),:,1,1) = tempDem;
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHDW')),:,2,1) = tempDem;

            try
                tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LHBC'),systemParams.IndexNumNodes} '_' num2str(costYear)]);
            catch
                tempDem=load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep systemParams.ValueLoadTag{ismember(systemParams.IndexIDL,'LHBC'),systemParams.IndexNumNodes}]);
            end
            tt = fieldnames(tempDem);
            tempDem.BCHprof = tempDem.(tt{1});
            tempDem = round(tempDem.BCHprof/sum(tempDem.BCHprof)*systemParamsReg.Instalations(1,find(systemParamsReg.IndexYears==costYear),find(ismember(systemParamsReg.IndexID,'LHBC'))));
            tempDem(isnan(tempDem))=0;
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHBC')),:,1,1) = tempDem;
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHBC')),:,2,1) = tempDem;

        else

            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHIN')),:,1,1) = 0*systemParams.ValueLoad(find(ismember(systemParamsReg.IndexIDL,'LELE')),:,1,1);
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHIN')),:,2,1) = 0*systemParams.ValueLoad(find(ismember(systemParamsReg.IndexIDL,'LELE')),:,1,1);

            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHSP')),:,1,1) = 0*systemParams.ValueLoad(find(ismember(systemParamsReg.IndexIDL,'LELE')),:,1,1);
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHSP')),:,2,1) = 0*systemParams.ValueLoad(find(ismember(systemParamsReg.IndexIDL,'LELE')),:,1,1);

            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHDW')),:,1,1) = 0*systemParams.ValueLoad(find(ismember(systemParamsReg.IndexIDL,'LELE')),:,1,1);
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHDW')),:,2,1) = 0*systemParams.ValueLoad(find(ismember(systemParamsReg.IndexIDL,'LELE')),:,1,1);

            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHBC')),:,1,1) = 0*systemParams.ValueLoad(find(ismember(systemParamsReg.IndexIDL,'LELE')),:,1,1);
            systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LHBC')),:,2,1) = 0*systemParams.ValueLoad(find(ismember(systemParamsReg.IndexIDL,'LELE')),:,1,1);
        end
        %% El and tech costs



        switch i
            case 1
                el_cost = systemParamsReg.ElCostRES(Reg,year);

                systemParams.FeedInEfficiencies(ismember(systemParamsReg.IndexIDR,'RPVO')) = systemParamsReg.FeedInEfficiencies(ismember(systemParamsReg.IndexIDR,'RPVR'));

                systemParams.Capex(ismember(systemParamsReg.IndexID,'RPVO'),year) = systemParamsReg.Capex(ismember(systemParamsReg.IndexID,'RPVR'),year);
                systemParams.Opex_fix(ismember(systemParamsReg.IndexID,'RPVO'),year) = systemParamsReg.Opex_fix(ismember(systemParamsReg.IndexID,'RPVR'),year);
                systemParams.Opex_var(ismember(systemParamsReg.IndexID,'RPVO'),year) = systemParamsReg.Opex_var(ismember(systemParamsReg.IndexID,'RPVR'),year);

                systemParams.Capex(ismember(systemParamsReg.IndexID,'SBAT'),year) = systemParamsReg.Capex(ismember(systemParamsReg.IndexID,'SBAR'),year);
                systemParams.Opex_fix(ismember(systemParamsReg.IndexID,'SBAT'),year) = systemParamsReg.Opex_fix(ismember(systemParamsReg.IndexID,'SBAR'),year);
                systemParams.Opex_var(ismember(systemParamsReg.IndexID,'SBAT'),year) = systemParamsReg.Opex_var(ismember(systemParamsReg.IndexID,'SBAR'),year);
                systemParams.Capex(ismember(systemParamsReg.IndexID,'IBAT'),year) = systemParamsReg.Capex(ismember(systemParamsReg.IndexID,'IBAR'),year);
                systemParams.Opex_fix(ismember(systemParamsReg.IndexID,'IBAT'),year) = systemParamsReg.Opex_fix(ismember(systemParamsReg.IndexID,'IBAR'),year);
                systemParams.Opex_var(ismember(systemParamsReg.IndexID,'IBAT'),year) = systemParamsReg.Opex_var(ismember(systemParamsReg.IndexID,'IBAR'),year);

            case 2
                el_cost = systemParamsReg.ElCostCOM(Reg,year);

                systemParams.FeedInEfficiencies(ismember(systemParamsReg.IndexIDR,'RPVO')) = systemParamsReg.FeedInEfficiencies(ismember(systemParamsReg.IndexIDR,'RPVC'));

                systemParams.Capex(ismember(systemParamsReg.IndexID,'RPVO'),year) = systemParamsReg.Capex(ismember(systemParamsReg.IndexID,'RPVC'),year);
                systemParams.Opex_fix(ismember(systemParamsReg.IndexID,'RPVO'),year) = systemParamsReg.Opex_fix(ismember(systemParamsReg.IndexID,'RPVC'),year);
                systemParams.Opex_var(ismember(systemParamsReg.IndexID,'RPVO'),year) = systemParamsReg.Opex_var(ismember(systemParamsReg.IndexID,'RPVC'),year);

                systemParams.Capex(ismember(systemParamsReg.IndexID,'SBAT'),year) = systemParamsReg.Capex(ismember(systemParamsReg.IndexID,'SBAC'),year);
                systemParams.Opex_fix(ismember(systemParamsReg.IndexID,'SBAT'),year) = systemParamsReg.Opex_fix(ismember(systemParamsReg.IndexID,'SBAC'),year);
                systemParams.Opex_var(ismember(systemParamsReg.IndexID,'SBAT'),year) = systemParamsReg.Opex_var(ismember(systemParamsReg.IndexID,'SBAC'),year);
                systemParams.Capex(ismember(systemParamsReg.IndexID,'IBAT'),year) = systemParamsReg.Capex(ismember(systemParamsReg.IndexID,'IBAC'),year);
                systemParams.Opex_fix(ismember(systemParamsReg.IndexID,'IBAT'),year) = systemParamsReg.Opex_fix(ismember(systemParamsReg.IndexID,'IBAC'),year);
                systemParams.Opex_var(ismember(systemParamsReg.IndexID,'IBAT'),year) = systemParamsReg.Opex_var(ismember(systemParamsReg.IndexID,'IBAC'),year);

            case 3
                el_cost = systemParamsReg.ElCostIND(Reg,year);


                systemParams.FeedInEfficiencies(ismember(systemParamsReg.IndexIDR,'RPVO')) = systemParamsReg.FeedInEfficiencies(ismember(systemParamsReg.IndexIDR,'RPVI'));

                systemParams.Capex(ismember(systemParamsReg.IndexID,'RPVO'),year) = systemParamsReg.Capex(ismember(systemParamsReg.IndexID,'RPVI'),year);
                systemParams.Opex_fix(ismember(systemParamsReg.IndexID,'RPVO'),year) = systemParamsReg.Opex_fix(ismember(systemParamsReg.IndexID,'RPVI'),year);
                systemParams.Opex_var(ismember(systemParamsReg.IndexID,'RPVO'),year) = systemParamsReg.Opex_var(ismember(systemParamsReg.IndexID,'RPVI'),year);

                systemParams.Capex(ismember(systemParamsReg.IndexID,'SBAT'),year) = systemParamsReg.Capex(ismember(systemParamsReg.IndexID,'SBAI'),year);
                systemParams.Opex_fix(ismember(systemParamsReg.IndexID,'SBAT'),year) = systemParamsReg.Opex_fix(ismember(systemParamsReg.IndexID,'SBAI'),year);
                systemParams.Opex_var(ismember(systemParamsReg.IndexID,'SBAT'),year) = systemParamsReg.Opex_var(ismember(systemParamsReg.IndexID,'SBAI'),year);
                systemParams.Capex(ismember(systemParamsReg.IndexID,'IBAT'),year) = systemParamsReg.Capex(ismember(systemParamsReg.IndexID,'IBAI'),year);
                systemParams.Opex_fix(ismember(systemParamsReg.IndexID,'IBAT'),year) = systemParamsReg.Opex_fix(ismember(systemParamsReg.IndexID,'IBAI'),year);
                systemParams.Opex_var(ismember(systemParamsReg.IndexID,'IBAT'),year) = systemParamsReg.Opex_var(ismember(systemParamsReg.IndexID,'IBAI'),year);

        end

        systemParams.Lifetime(ismember(systemParamsReg.IndexID,'SBAT')) = systemParamsReg.Lifetime(ismember(systemParamsReg.IndexID,batteryNames{i}));
        systemParams.Lifetime(ismember(systemParamsReg.IndexID,'IBAT')) = systemParamsReg.Lifetime(ismember(systemParamsReg.IndexID,interfaceNames{i}));
        systemParams.Lifetime(ismember(systemParamsReg.IndexID,'RPVO')) = systemParamsReg.Lifetime(ismember(systemParamsReg.IndexID,produsersNames{i}));

        instalMask = (repmat(systemParamsReg.IndexYears,size(systemParamsReg.IndexID,1),1)<=costYear).*(repmat(systemParamsReg.IndexYears,size(systemParamsReg.IndexID,1),1)>=(costYear+setup.stepYear-systemParams.Lifetime));


        for ii = 1:length(systemParams.IndexID)

            if systemParams.Active(ii)==1 & ~sum(ismember({'RPVO','SBAT','IBAT'},systemParams.IndexID{ii}))
                if setup.Heat.Flag
                    if ~(setup.OvernightFlag)
                        if (costYear <= setup.startYear) & ~setup.OvernightFlag
                            systemParams.SizeLimits(ii,1) = round(sum(shiftdim(systemParams.Instalations(1,:,ii),1)));
                            systemParams.SizeLimits(ii,2) = round(sum(shiftdim(systemParams.Instalations(1,:,ii),1))+1);
                        else
                            systemParams.SizeLimits(ii,1) = round(sum(shiftdim(systemParams.Instalations(1,:,ii),1).*instalMask(ii,:)'));
                            systemParams.SizeLimits(ii,2) = systemParamsReg.SizeLimits(ii,2);
                        end
                    end
                else
                    systemParams.SizeLimits(ii,1) = 0;
                    systemParams.SizeLimits(ii,2) = 0;
                end
            end
        end


        LowLimitSBAT = round(sum(shiftdim(systemParams.Instalations(1,:,find(ismember(systemParams.IndexID,batteryNames{i}))),1).*instalMask(ismember(systemParams.IndexID,batteryNames{i}),:)'));
        LowLimitIBAT = round(sum(shiftdim(systemParams.Instalations(1,:,find(ismember(systemParams.IndexID,interfaceNames{i}))),1).*instalMask(ismember(systemParams.IndexID,interfaceNames{i}),:)'));
        LowLimitPV = round(sum(shiftdim(systemParams.Instalations(1,:,find(ismember(systemParams.IndexID,produsersNames{i}))),1).*instalMask(ismember(systemParams.IndexID,produsersNames{i}),:)'));



        systemParams.SizeLimits(ismember(systemParams.IndexID,'SBAT'),1) = LowLimitSBAT;
        systemParams.SizeLimits(ismember(systemParams.IndexID,'IBAT'),1) = LowLimitIBAT;
        systemParams.SizeLimits(ismember(systemParams.IndexID,'RPVO'),1) = LowLimitPV;

        if (costYear <= setup.startYear) & ~setup.OvernightFlag
            systemParams.SizeLimits(ismember(systemParams.IndexID,'SBAT'),2) = LowLimitSBAT+1;
            systemParams.SizeLimits(ismember(systemParams.IndexID,'IBAT'),2) = LowLimitIBAT+1;
            systemParams.SizeLimits(ismember(systemParams.IndexID,'RPVO'),2) = LowLimitPV+1;
        else
            systemParams.SizeLimits(ismember(systemParams.IndexID,'SBAT'),2) = systemParams.SizeLimits(ismember(systemParamsReg.IndexID,batteryNames{i}),2);
            systemParams.SizeLimits(ismember(systemParams.IndexID,'IBAT'),2) = systemParams.SizeLimits(ismember(systemParamsReg.IndexID,interfaceNames{i}),2);
            if existPV_tech+(systemParams.SizeLimits(ismember(systemParamsReg.IndexID,produsersNames{i}),2)-existPV)*systemParams.shareOfSectorEl*0.6 > LowLimitPV
                systemParams.SizeLimits(ismember(systemParams.IndexID,'RPVO'),2) = existPV_tech+(systemParams.SizeLimits(ismember(systemParamsReg.IndexID,produsersNames{i}),2)-existPV)*systemParams.shareOfSectorEl*0.6;
            else
                systemParams.SizeLimits(ismember(systemParams.IndexID,'RPVO'),2) = LowLimitPV*1.05;
            end

        end

        if setup.Heat.Flag
            setup.shareOfSector = setup.Heat.shareOfIndivHeatSectors(i)/sum(setup.Heat.shareOfIndivHeatSectors);
            setup.shareIndHeatPros = max(0.8,setup.SC.share_of_SC(prosumersShareStep)/100*(1-setup.Heat.shareOfDistrHeat(year)))*min(1,setup.SC.share_of_SC(prosumersShareStep)/100/(1-setup.Heat.shareOfDistrHeat(year)));
        else
            setup.shareOfSector = 0;
            setup.shareIndHeatPros = 0;
        end
        %%efficiency
        efficiency = sum(shiftdim(sum(systemParams.Instalations,1),1).*instalMask'.*systemParams.Efficiency')./sum(shiftdim(sum(systemParams.Instalations,1),1).*instalMask');
        efficiency(isnan(efficiency)) = systemParams.Efficiency(isnan(efficiency),find(systemParams.IndexYears==costYear));
        efficiency(isinf(efficiency)) = systemParams.Efficiency(isinf(efficiency),find(systemParams.IndexYears==costYear));

        for jj = 1:length(systemParams.IndexIDS)
            if systemParams.Efficiency(ismember(systemParams.IndexID,systemParams.IndexIDS(jj)))
                systemParams.EtaStorage(jj,:) = systemParams.EtaStorageOrig(jj,:).*efficiency(ismember(systemParams.IndexID,systemParams.IndexIDS(jj)));
            end
        end

        for jj =  1:length(systemParams.IndexIDT)
            if systemParams.Efficiency(ismember(systemParams.IndexID,systemParams.IndexIDT(jj)))
                systemParams.EtaTrans(jj,1:end-1) = systemParams.EtaTransOrig(jj,1:end-1).*efficiency(ismember(systemParams.IndexID,systemParams.IndexIDT(jj)));
            end
        end


        %% biomass


        systemParams.ValueResourceTotal(ismember(systemParamsReg.IndexIDR,'RBMW'),1) = systemParamsReg.ValueResourceTotal(ismember(systemParamsReg.IndexIDR,'RBMW'),1)*setup.SC.Biomass_WastesShare;%.*setup.shareIndHeatPros;
        systemParams.Opex_var_reg(ismember(systemParamsReg.IndexID,'RBMW'),1,year) = systemParamsReg.Opex_var_reg(ismember(systemParamsReg.IndexID,'RBMW'),1,year).*setup.SC.Biomass_WastesAddCost(year);

        systemParams.ValueResourceTotal(ismember(systemParamsReg.IndexIDR,'RWOO'),1) = systemParamsReg.ValueResourceTotal(ismember(systemParamsReg.IndexIDR,'RWOO'),1)*setup.SC.Biomass_BiomassShare;%.*setup.shareIndHeatPros;
        systemParams.Opex_var_reg(ismember(systemParamsReg.IndexID,'RWOO'),1,year) = systemParamsReg.Opex_var_reg(ismember(systemParamsReg.IndexID,'RWOO'),1,year).*setup.SC.Biomass_BiomassAddCost(year);
        systemParams.ValueResourceTotal(ismember(systemParamsReg.IndexIDR,'RWWO'),1) = 0;
        if costYear <=setup.startYear
            systemParams.ValueResourceTotal(ismember(systemParamsReg.IndexIDR,'RWOO'),1) = 5*10^10;%systemParamsReg.ValueResourceTotal(ismember(systemParamsReg.IndexIDR,'RWOO'),1)*10;%.*setup.shareIndHeatPros;
        end

        systemParams.ValueResource(find(ismember(systemParamsReg.IndexIDR,'RBGA')),:,1,1) = systemParamsReg.ValueResource(find(ismember(systemParamsReg.IndexIDR,'RBGA')),:,2,1)*setup.SC.Biomass_BiogasShare;%.*setup.shareIndHeatPros;
        systemParams.ValueResourceTotal(ismember(systemParamsReg.IndexIDR,'RBGA'),1)=sum(systemParamsReg.ValueResource(find(ismember(systemParamsReg.IndexIDR,'RBGA')),:,1,1))*setup.SC.Biomass_BiogasShare;%.*setup.shareIndHeatPros;
        systemParams.Opex_var_reg(ismember(systemParamsReg.IndexID,'RBGA'),1,year) = systemParamsReg.Opex_var_reg(ismember(systemParamsReg.IndexID,'RBGA'),1,year).*setup.SC.Biomass_BiogasAddCost(year);

        systemParams.Opex_var_reg(ismember(systemParamsReg.IndexID,'RNGA'),1,year) = systemParamsReg.Opex_var_reg(ismember(systemParamsReg.IndexID,'RNGA'),1,year).*setup.SC.Biomass_GasAddCost(year);
        systemParams.Opex_var_reg(ismember(systemParamsReg.IndexID,'RPET'),1,year) = systemParamsReg.Opex_var_reg(ismember(systemParamsReg.IndexID,'RPET'),1,year).*setup.SC.Biomass_OilAddCost(year);

        %%


        if setup.SC.feedin_cost>el_cost
            feedin_cost = 0.8*el_cost;
        else
            feedin_cost = setup.SC.feedin_cost;
        end

		results = RunSelfConsTransition(setup,setup.rootDir,setup.projName,systemParams,activeElements,modFiles,Reg,produsersTypes{i},costYear,setup.SC.share_of_SC(prosumersShareStep),el_cost,setup.SC.feedin_cost,setup.endHour);


        if sum(results.RPVO_EL)>0.3*sum(sum(systemParams.ValueLoad(1,:,1,1)))
            if prosumersShareStep<length(setup.SC.share_of_SC)
                prosumersShareStep=prosumersShareStep+1;
            end
        end
        if (costYear>setup.startYear) | setup.OvernightFlag
            systemParams.Instalations(1,systemParams.IndexYears==costYear,find(ismember(systemParams.IndexID,produsersNames{i})))=systemParamsReg.Instalations(1,systemParams.IndexYears==costYear,find(ismember(systemParams.IndexID,produsersNames{i})))+ceil(results.OPT_SIZE_RPVO-LowLimitPV)*(ceil(results.OPT_SIZE_RPVO-LowLimitPV)>0);

            systemParams.Instalations(1,systemParams.IndexYears==costYear,find(ismember(systemParams.IndexID,batteryNames{i})))=systemParamsReg.Instalations(1,systemParams.IndexYears==costYear,find(ismember(systemParams.IndexID,batteryNames{i})))+ceil(results.OPT_SIZE_SBAT-LowLimitSBAT)*(ceil(results.OPT_SIZE_SBAT-LowLimitSBAT)>0);

            systemParams.Instalations(1,systemParams.IndexYears==costYear,find(ismember(systemParams.IndexID,interfaceNames{i})))=systemParamsReg.Instalations(1,systemParams.IndexYears==costYear,find(ismember(systemParams.IndexID,interfaceNames{i})))+ceil(results.OPT_SIZE_IBAT-LowLimitIBAT)*(ceil(results.OPT_SIZE_IBAT-LowLimitIBAT)>0);


            for ii = 1:length(systemParams.IndexID)

                if systemParams.Active(ii)==1 & ~sum(ismember({'RPVO','SBAT','IBAT'},systemParams.IndexID{ii}))
                    try
                        eval(['instCap = results.OPT_SIZE_' systemParams.IndexID{ii} ';']);

                        systemParams.Instalations(:,find(systemParams.IndexYears==costYear),ii)=systemParams.Instalations(:,find(systemParams.IndexYears==costYear),ii)+max((instCap-shiftdim(systemParams.SizeLimits(ii,1,:),1))',0);

                    end
                end
            end
        end

        systemParams.Instalations(:,find(systemParams.IndexYears==costYear),ismember(systemParams.IndexID,'RHAR')) = sum(results.RHAR_FU);
        systemParams.Instalations(:,find(systemParams.IndexYears==costYear),ismember(systemParams.IndexID,'RPET')) = sum(results.RPET_FU);
        systemParams.Instalations(:,find(systemParams.IndexYears==costYear),ismember(systemParams.IndexID,'RNGA')) = sum(results.RNGA_FU);

        systemParams.Instalations(:,find(systemParams.IndexYears==costYear),ismember(systemParams.IndexID,'RBGA')) = sum(results.RBGA_FU);
        systemParams.Instalations(:,find(systemParams.IndexYears==costYear),ismember(systemParams.IndexID,'RBMW')) = sum(results.RBMW_FU);
        systemParams.Instalations(:,find(systemParams.IndexYears==costYear),ismember(systemParams.IndexID,'RWOO')) = sum(results.RWOO_FU);
        systemParams.Instalations(:,find(systemParams.IndexYears==costYear),ismember(systemParams.IndexID,'RWWO')) = sum(results.RWWO_FU);
        systemParams.RBGAFlows(find(systemParams.IndexYears==costYear),:) = results.RBGA_FU;


    end

    save([setup.rootDir filesep 'projects' filesep setup.projName filesep 'output' filesep 'Instalations_' produsersTypes{i} '_' num2str(Reg)],'systemParams');
end