function main_Self_Cons_Transition(setupData,systemData)
%
% FUNCTION main_Self_Cons_Transition(setupData, systemData)
%
% Runs the self-consumption analysis for a system with transition settings.
%
%
% INPUT:
%            setupData:    Structure with general settings and inputs for the model.
%            systemData:   Structure with system-specific input data.
%
% OUTPUT:
%            -
%
%Dmitrii Bogdanov
%last change 23.07.2025


endHour = 8760;

setupData.endHour = endHour;

produsersTypes=[{'RES'} {'COM'} {'IND'} ];

yearConsumptions = systemData.systemParams.SectorsCons;


numReg = length(systemData.systemParams.IndexNodes);

startReg = 1;

try setupData.SC_Regions;
    SC_Regions = setupData.SC_Regions;
catch
    setupData.SC_Regions = [startReg:numReg]';
    SC_Regions = setupData.SC_Regions;
end


if length(SC_Regions)==1
    for i = startReg:length(SC_Regions)

        Reg = SC_Regions(i);


        systemParamsReg{i} = systemData.systemParams;

        activeElements = systemData.activeElements;
        activeElements = ActiveComponentsProsumers(systemData.systemParams,activeElements);

        systemParamsReg{i}.IndexNodes = systemData.systemParams.IndexNodes(Reg:Reg);
        systemParamsReg{i}.IndexNumNodes = systemData.systemParams.IndexNumNodes(Reg:Reg);
        systemParamsReg{i}.SectorsCons = systemData.systemParams.SectorsCons(Reg:Reg,:);
        systemParamsReg{i}.Coords = systemData.systemParams.Coords(:,Reg:Reg);
        systemParamsReg{i}.Urbanisation = systemData.systemParams.Urbanisation(:,Reg:Reg);
        systemParamsReg{i}.SizeLimits = systemData.systemParams.SizeLimits(:,:,Reg:Reg);
        systemParamsReg{i}.ValueLoad = systemData.systemParams.ValueLoad(:,:,:,Reg:Reg);
        systemParamsReg{i}.Instalations = systemData.systemParams.Instalations(Reg:Reg,:,:);
        systemParamsReg{i}.ValueResource = systemData.systemParams.ValueResource(:,:,:,Reg:Reg);
        systemParamsReg{i}.ValueHydroDam = systemData.systemParams.ValueHydroDam(:,:,:,Reg:Reg);
        systemParamsReg{i}.Capex_reg = systemData.systemParams.Capex_reg(:,Reg:Reg,:);
        systemParamsReg{i}.Opex_fix_reg = systemData.systemParams.Opex_fix_reg(:,Reg:Reg,:);
        systemParamsReg{i}.Opex_var_reg = systemData.systemParams.Opex_var_reg(:,Reg:Reg,:);
        systemParamsReg{i}.ValueResourceTotal = systemData.systemParams.ValueResourceTotal(:,Reg:Reg);
        systemParamsReg{i}.DesalinationParams.regionPumpLong = systemData.systemParams.DesalinationParams.regionPumpLong(Reg:Reg,:);
        systemParamsReg{i}.DesalinationParams.regionPumpUp = systemData.systemParams.DesalinationParams.regionPumpUp(Reg:Reg,:);
        systemParamsReg{i}.AC_losses = systemData.systemParams.AC_losses(Reg:Reg,:);
        systemParamsReg{i}.TLlength = 0;

        systemParamsReg{i}.TemperatureAir = systemData.systemParams.TemperatureAir(:,Reg:Reg);
        systemParamsReg{i}.TemperatureWater = systemData.systemParams.TemperatureWater(:,Reg:Reg);
        systemParamsReg{i}.TemperatureGround = systemData.systemParams.TemperatureGround(:,Reg:Reg);


        setup{i} = setupData;
        setup{i}.SC.share_of_SC = setupData.SC.share_of_SC(Reg:Reg,:);
        setup{i}.SC.share_of_SC_RES = setupData.SC.share_of_SC_RES(Reg:Reg,:);
        setup{i}.SC.share_of_SC_COM = setupData.SC.share_of_SC_COM(Reg:Reg,:);
        setup{i}.SC.share_of_SC_IND = setupData.SC.share_of_SC_IND(Reg:Reg,:);
        setup{i}.SC.feedin_cost = setupData.SC.feedin_cost(Reg:Reg);

        setup{i}.SC.Biomass_WastesShare = setupData.SC.Biomass_WastesShare(Reg:Reg);
        setup{i}.SC.Biomass_BiomassShare = setupData.SC.Biomass_BiomassShare(Reg:Reg);
        setup{i}.SC.Biomass_BiogasShare = setupData.SC.Biomass_BiogasShare(Reg:Reg);
        setup{i}.SC.Biomass_WastesAddCost = setupData.SC.Biomass_WastesAddCost(:,Reg:Reg);
        setup{i}.SC.Biomass_BiomassAddCost = setupData.SC.Biomass_BiomassAddCost(:,Reg:Reg);
        setup{i}.SC.Biomass_BiogasAddCost = setupData.SC.Biomass_BiogasAddCost(:,Reg:Reg);
        setup{i}.SC.Biomass_GasAddCost = setupData.SC.Biomass_GasAddCost(:,Reg:Reg);
        setup{i}.SC.Biomass_OilAddCost = setupData.SC.Biomass_OilAddCost(:,Reg:Reg);

        setup{i}.Heat.shareOfIndivHeatSectors = setupData.Heat.shareOfIndivHeatSectors(Reg:Reg,:);
        setup{i}.Heat.shareOfDistrHeat = setupData.Heat.shareOfDistrHeat(Reg:Reg,:);
        setup{i}.Heat.shareOfLowHeatInd = setupData.Heat.shareOfLowHeatInd(Reg:Reg);
        setup{i}.Heat.shareOfHighHeatInd = setupData.Heat.shareOfHighHeatInd(Reg:Reg);

        setup{i}.groundHPshare_IH = setupData.groundHPshare_IH(Reg:Reg,:);
        setup{i}.airHPshare_IH = setupData.airHPshare_IH(Reg:Reg,:);
        setup{i}.waterHPshare_IH = setupData.waterHPshare_IH(Reg:Reg,:);
        setup{i}.groundHPshare_DH = setupData.groundHPshare_DH(Reg:Reg,:);
        setup{i}.airHPshare_DH = setupData.airHPshare_DH(Reg:Reg,:);
        setup{i}.waterHPshare_DH = setupData.waterHPshare_DH(Reg:Reg,:);

        PrepareSelfConsScenario(setup{i},Reg,systemParamsReg{i},activeElements,produsersTypes);




    end
else
    try setupData.ParallelOpt
    catch
        setupData.ParallelOpt=min(length(SC_Regions),feature('numcores')-1);
    end
    
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        parpool(setupData.ParallelOpt)
    end
    %parpool(min(length(SC_Regions)+1,feature('numcores')-1))
    parfor i = startReg:length(SC_Regions)

        pause(i/100);
        Reg = SC_Regions(i);


        systemParamsReg{i} = systemData.systemParams;

        activeElements = systemData.activeElements;
        activeElements = ActiveComponentsProsumers(systemData.systemParams,activeElements);


        systemParamsReg{i}.IndexNodes = systemData.systemParams.IndexNodes(Reg:Reg);
        systemParamsReg{i}.IndexNumNodes = systemData.systemParams.IndexNumNodes(Reg:Reg);
        systemParamsReg{i}.SectorsCons = systemData.systemParams.SectorsCons(Reg:Reg,:);
        systemParamsReg{i}.Coords = systemData.systemParams.Coords(:,Reg:Reg);
        systemParamsReg{i}.Urbanisation = systemData.systemParams.Urbanisation(:,Reg:Reg);
        systemParamsReg{i}.SizeLimits = systemData.systemParams.SizeLimits(:,:,Reg:Reg);
        systemParamsReg{i}.ValueLoad = systemData.systemParams.ValueLoad(:,:,:,Reg:Reg);
        systemParamsReg{i}.Instalations = systemData.systemParams.Instalations(Reg:Reg,:,:);
        systemParamsReg{i}.ValueResource = systemData.systemParams.ValueResource(:,:,:,Reg:Reg);
        systemParamsReg{i}.ValueHydroDam = systemData.systemParams.ValueHydroDam(:,:,:,Reg:Reg);
        systemParamsReg{i}.Capex_reg = systemData.systemParams.Capex_reg(:,Reg:Reg,:);
        systemParamsReg{i}.Opex_fix_reg = systemData.systemParams.Opex_fix_reg(:,Reg:Reg,:);
        systemParamsReg{i}.Opex_var_reg = systemData.systemParams.Opex_var_reg(:,Reg:Reg,:);
        systemParamsReg{i}.ValueResourceTotal = systemData.systemParams.ValueResourceTotal(:,Reg:Reg);
        systemParamsReg{i}.DesalinationParams.regionPumpLong = systemData.systemParams.DesalinationParams.regionPumpLong(Reg:Reg,:);
        systemParamsReg{i}.DesalinationParams.regionPumpUp = systemData.systemParams.DesalinationParams.regionPumpUp(Reg:Reg,:);
        systemParamsReg{i}.AC_losses = systemData.systemParams.AC_losses(Reg:Reg,:);
        systemParamsReg{i}.TLlength = 0;

        systemParamsReg{i}.TemperatureAir = systemData.systemParams.TemperatureAir(:,Reg:Reg);
        systemParamsReg{i}.TemperatureWater = systemData.systemParams.TemperatureWater(:,Reg:Reg);
        systemParamsReg{i}.TemperatureGround = systemData.systemParams.TemperatureGround(:,Reg:Reg);


        setup{i} = setupData;
        setup{i}.SC.share_of_SC = setupData.SC.share_of_SC(Reg:Reg,:);
        setup{i}.SC.share_of_SC_RES = setupData.SC.share_of_SC_RES(Reg:Reg,:);
        setup{i}.SC.share_of_SC_COM = setupData.SC.share_of_SC_COM(Reg:Reg,:);
        setup{i}.SC.share_of_SC_IND = setupData.SC.share_of_SC_IND(Reg:Reg,:);
        setup{i}.SC.feedin_cost = setupData.SC.feedin_cost(Reg:Reg);

        setup{i}.SC.Biomass_WastesShare = setupData.SC.Biomass_WastesShare(Reg:Reg);
        setup{i}.SC.Biomass_BiomassShare = setupData.SC.Biomass_BiomassShare(Reg:Reg);
        setup{i}.SC.Biomass_BiogasShare = setupData.SC.Biomass_BiogasShare(Reg:Reg);
        setup{i}.SC.Biomass_WastesAddCost = setupData.SC.Biomass_WastesAddCost(:,Reg:Reg);
        setup{i}.SC.Biomass_BiomassAddCost = setupData.SC.Biomass_BiomassAddCost(:,Reg:Reg);
        setup{i}.SC.Biomass_BiogasAddCost = setupData.SC.Biomass_BiogasAddCost(:,Reg:Reg);
        setup{i}.SC.Biomass_GasAddCost = setupData.SC.Biomass_GasAddCost(:,Reg:Reg);
        setup{i}.SC.Biomass_OilAddCost = setupData.SC.Biomass_OilAddCost(:,Reg:Reg);

        setup{i}.Heat.shareOfIndivHeatSectors = setupData.Heat.shareOfIndivHeatSectors(Reg:Reg,:);
        setup{i}.Heat.shareOfDistrHeat = setupData.Heat.shareOfDistrHeat(Reg:Reg,:);
        setup{i}.Heat.shareOfLowHeatInd = setupData.Heat.shareOfLowHeatInd(Reg:Reg);
        setup{i}.Heat.shareOfHighHeatInd = setupData.Heat.shareOfHighHeatInd(Reg:Reg);

        setup{i}.groundHPshare_IH = setupData.groundHPshare_IH(Reg:Reg,:);
        setup{i}.airHPshare_IH = setupData.airHPshare_IH(Reg:Reg,:);
        setup{i}.waterHPshare_IH = setupData.waterHPshare_IH(Reg:Reg,:);
        setup{i}.groundHPshare_DH = setupData.groundHPshare_DH(Reg:Reg,:);
        setup{i}.airHPshare_DH = setupData.airHPshare_DH(Reg:Reg,:);
        setup{i}.waterHPshare_DH = setupData.waterHPshare_DH(Reg:Reg,:);


        PrepareSelfConsScenario(setup{i},Reg,systemParamsReg{i},activeElements,produsersTypes);


    end
    delete(gcp('nocreate'));
end


for costYear = setupData.startYear:setupData.stepYear:setupData.endYear

    SCData.generation = zeros(8760,numReg);
    SCData.RPVO_EL = zeros(8760,numReg);
    SCData.OPT_SIZE_RPVO = zeros(1,numReg);


    SCData.ToGrid = zeros(8760,numReg);
    SCData.FromGrid = zeros(8760,numReg);

    SCData.RBGA_Cons = zeros(8760,numReg);
    SCData.RBMW_Cons = zeros(8760,numReg);
    SCData.RWWO_Cons = zeros(8760,numReg);
    SCData.RWOO_Cons = zeros(8760,numReg);

    SCData.RHAR_Cons = zeros(8760,numReg);
    SCData.RPET_Cons = zeros(8760,numReg);
    SCData.RNGA_Cons = zeros(8760,numReg);

    SCData.SC_demand = zeros(8760,numReg);

    SCData.SC_share = zeros(3,numReg);
    SCData.shareIndHeatElPros = zeros(3,numReg);

    for RegR = 1:numReg

        for i=1:length(produsersTypes)
            try
                load([setupData.rootDir filesep 'projects' filesep setupData.projName filesep 'output' filesep 'results_' num2str(RegR) '_' produsersTypes{i} '_' num2str(costYear) '.mat']);

                SCData.generation(:,RegR) = SCData.generation(:,RegR)+(results.RPVO_EL+results.SBAT_EL-results.EL_SBAT);
                SCData.RPVO_EL(:,RegR) = SCData.RPVO_EL(:,RegR)+(results.RPVO_EL);
                SCData.OPT_SIZE_RPVO(:,RegR) = SCData.OPT_SIZE_RPVO(:,RegR)+(results.OPT_SIZE_RPVO);


                SCData.ToGrid(:,RegR) = SCData.ToGrid(:,RegR)+results.EL_EXCESS;
                SCData.FromGrid(:,RegR) = SCData.FromGrid(:,RegR)+results.EL_GRID;

                SCData.RBGA_Cons(:,RegR) = SCData.RBGA_Cons(:,RegR)+results.RBGA_FU;
                SCData.RBMW_Cons(:,RegR) = SCData.RBMW_Cons(:,RegR)+results.RBMW_FU;
                SCData.RWWO_Cons(:,RegR) = SCData.RWWO_Cons(:,RegR)+results.RWWO_FU;
                SCData.RWOO_Cons(:,RegR) = SCData.RWOO_Cons(:,RegR)+results.RWOO_FU;

                SCData.RHAR_Cons(:,RegR) = SCData.RHAR_Cons(:,RegR)+results.RHAR_FU;
                SCData.RPET_Cons(:,RegR) = SCData.RPET_Cons(:,RegR)+results.RPET_FU;
                SCData.RNGA_Cons(:,RegR) = SCData.RNGA_Cons(:,RegR)+results.RNGA_FU;

                SCData.SC_demand(:,RegR) = SCData.SC_demand(:,RegR)+results.SC_demand;

                SCData.SC_share(i,RegR) = results.SC_share';
                SCData.shareIndHeatElPros(i,RegR) = results.shareIndHeatElPros';

            catch
                display(['No result for region ' num2str(RegR)])
            end

        end
    end

    %selfConsumption = generation;
    save([setupData.rootDir filesep 'projects' filesep setupData.projName filesep 'output' filesep 'SelfConsumptionData' '_' num2str(costYear) '.mat'],'SCData');

end
for costYear = setupData.startYear:setupData.stepYear:setupData.endYear
    for k=1:length(produsersTypes)
        try
            load([setupData.rootDir filesep 'projects' filesep setupData.projName filesep 'output' filesep 'results_' num2str(1) '_' produsersTypes{k} '_' num2str(costYear) '.mat']);
            namesFields = fieldnames(results);
        catch
            load([setupData.rootDir filesep 'projects' filesep setupData.projName filesep 'output'  filesep 'results_' num2str(47) '_' produsersTypes{k} '_' num2str(costYear) '.mat']);
            namesFields = fieldnames(results);
        end

        for i = 2:numReg
            try
                temp = load([setupData.rootDir filesep 'projects' filesep setupData.projName filesep 'output' filesep 'results_' num2str(i) '_' produsersTypes{k} '_' num2str(costYear) '.mat']);
                for j=1:length(namesFields)
                    eval(['results.' namesFields{j} '=[results.' namesFields{j} ' temp.results.' namesFields{j} '];'])
                end
            catch
                display(['No result for region ' num2str(i) ', put all zeros'])
                temp = load([setupData.rootDir filesep 'projects' filesep setupData.projName filesep 'output' filesep 'results_' num2str(1) '_' produsersTypes{k} '_' num2str(costYear) '.mat']);

                for j=1:length(namesFields)
                    eval(['results.' namesFields{j} '=[results.' namesFields{j} ' 0*temp.results.' namesFields{j} '];'])
                end
            end

        end
        save([setupData.rootDir filesep 'projects' filesep setupData.projName filesep 'output' filesep 'results_' setupData.projName '_' produsersTypes{k} '_' num2str(costYear) '.mat'],'results');
    end
end