function PrepareCountriesScenario(systemData,setup)
%
% retructures the input data according to given Country configuration,
% or for individual Regions simulations
%
% last change Dmitrii Bogdanov 07.11.2023


setup.CountriesType='byCountries';

setup.fixYear = 200;
startReg = 1;


try setup.projectSize;
catch
setup.projectSize = size(setup.Countries,1);
end



if setup.projectSize>1
    
    try setup.ParallelOpt
    catch
        setup.ParallelOpt=min(setup.projectSize,feature('numcores')-1);
    end

    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        parpool(setup.ParallelOpt)
    end

    parfor (numCon = startReg:size(setup.Countries,1),setup.ParallelOpt)


        if iscell(setup.Countries)
            Reg = setup.Countries{numCon}(1);
            Reg2 = setup.Countries{numCon}(end);

            Regions{numCon} = setup.Countries{numCon};

        else
            Reg = setup.Countries(numCon,1);
            Reg2 = setup.Countries(numCon,2);

            Regions{numCon} = [Reg:Reg2];

        end

        systemDataTemp{numCon} = systemData;
        systemDataTemp{numCon}.systemParams = systemData.systemParams;
        systemDataTemp{numCon}.systemParams.Regions = Regions{numCon};
        systemDataTemp{numCon}.systemParams.area = systemData.systemParams.area(Regions{numCon});
        systemDataTemp{numCon}.systemParams.population = systemData.systemParams.population(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.IndexNodes = systemData.systemParams.IndexNodes(Regions{numCon});
        systemDataTemp{numCon}.systemParams.IndexNumNodes = systemData.systemParams.IndexNumNodes(Regions{numCon});
        systemDataTemp{numCon}.systemParams.Coords = systemData.systemParams.Coords(:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.Urbanisation = systemData.systemParams.Urbanisation(:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.SizeLimits = systemData.systemParams.SizeLimits(:,:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.Instalations = systemData.systemParams.Instalations(Regions{numCon},:,:);
        systemDataTemp{numCon}.systemParams.ValueLoad = systemData.systemParams.ValueLoad(:,:,:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.ValueLoadTag = systemData.systemParams.ValueLoadTag(:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.ValueResource = systemData.systemParams.ValueResource(:,:,:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.ValueHydroDam = systemData.systemParams.ValueHydroDam(:,:,:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.ValueResourceTag = systemData.systemParams.ValueResourceTag(:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.ValueHydroDamTag = systemData.systemParams.ValueHydroDamTag(:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.Capex_reg = systemData.systemParams.Capex_reg(:,Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Opex_fix_reg = systemData.systemParams.Opex_fix_reg(:,Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Opex_var_reg = systemData.systemParams.Opex_var_reg(:,Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.ValueResourceTotal = systemData.systemParams.ValueResourceTotal(:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.DesalinationParams.regionPumpLong = systemData.systemParams.DesalinationParams.regionPumpLong(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.DesalinationParams.regionPumpUp = systemData.systemParams.DesalinationParams.regionPumpUp(Regions{numCon},:);


        systemDataTemp{numCon}.systemParams.Mobility.Demand = systemData.systemParams.Mobility.Demand(:,Regions{numCon},:);

        systemDataTemp{numCon}.systemParams.Mobility.LDV_Pass = systemData.systemParams.Mobility.LDV_Pass(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.BUS_Pass = systemData.systemParams.Mobility.BUS_Pass(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.W23_Pass = systemData.systemParams.Mobility.W23_Pass(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.LDV_km = systemData.systemParams.Mobility.LDV_km(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.BUS_km = systemData.systemParams.Mobility.BUS_km(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.W23_km = systemData.systemParams.Mobility.W23_km(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.LDV_PHEV_el = systemData.systemParams.Mobility.LDV_PHEV_el(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.BUS_PHEV_el = systemData.systemParams.Mobility.BUS_PHEV_el(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.W23_PHEV_el = systemData.systemParams.Mobility.W23_PHEV_el(Regions{numCon},:);

        systemDataTemp{numCon}.systemParams.Mobility.MDV_Tonne = systemData.systemParams.Mobility.MDV_Tonne(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.HDV_Tonne = systemData.systemParams.Mobility.HDV_Tonne(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.MDV_km = systemData.systemParams.Mobility.MDV_km(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.HDV_km = systemData.systemParams.Mobility.HDV_km(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.MDV_PHEV_el = systemData.systemParams.Mobility.MDV_PHEV_el(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.HDV_PHEV_el = systemData.systemParams.Mobility.HDV_PHEV_el(Regions{numCon},:);

        systemDataTemp{numCon}.systemParams.Mobility.Shares = systemData.systemParams.Mobility.Shares(Regions{numCon},:,:);

        systemDataTemp{numCon}.systemParams.EffCO2Scr_El = systemData.systemParams.EffCO2Scr_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffCO2Scr_He = systemData.systemParams.EffCO2Scr_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffMETH_CO = systemData.systemParams.EffMETH_CO(Regions{numCon},:);
        try
            systemDataTemp{numCon}.systemParams.EffFT_El = systemData.systemParams.EffFT_El(Regions{numCon},:);
        end
        systemDataTemp{numCon}.systemParams.EffFT_CO = systemData.systemParams.EffFT_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffFT_He = systemData.systemParams.EffFT_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffLNG_El = systemData.systemParams.EffLNG_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffLH2_El = systemData.systemParams.EffLH2_El(Regions{numCon},:);


        systemDataTemp{numCon}.systemParams.EffMeO_CO = systemData.systemParams.EffMeO_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffMeO_El = systemData.systemParams.EffMeO_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffMeO_He_out = systemData.systemParams.EffMeO_He_out(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffDME_CO = systemData.systemParams.EffDME_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffDME_El = systemData.systemParams.EffDME_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffDME_He_out = systemData.systemParams.EffDME_He_out(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffNH3_El = systemData.systemParams.EffNH3_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffNH3_He_out = systemData.systemParams.EffNH3_He_out(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffHyS_El = systemData.systemParams.EffHyS_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffMET_El = systemData.systemParams.EffMET_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffMET_He_out = systemData.systemParams.EffMET_He_out(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffWEL_He_out = systemData.systemParams.EffWEL_He_out(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffICB_Lime = systemData.systemParams.EffICB_Lime(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffICB_El = systemData.systemParams.EffICB_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffICB_He = systemData.systemParams.EffICB_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffICI_Lime = systemData.systemParams.EffICI_Lime(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffICI_El = systemData.systemParams.EffICI_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffICI_He = systemData.systemParams.EffICI_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISB_Ore = systemData.systemParams.EffISB_Ore(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISB_Scr = systemData.systemParams.EffISB_Scr(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISB_El = systemData.systemParams.EffISB_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISB_He = systemData.systemParams.EffISB_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISB_HC = systemData.systemParams.EffISB_HC(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISH_Ore = systemData.systemParams.EffISH_Ore(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISH_Scr = systemData.systemParams.EffISH_Scr(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISH_El = systemData.systemParams.EffISH_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISH_He = systemData.systemParams.EffISH_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISH_CC = systemData.systemParams.EffISH_CC(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISH_Hy = systemData.systemParams.EffISH_Hy(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISR_Scr = systemData.systemParams.EffISR_Scr(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISR_El = systemData.systemParams.EffISR_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISR_He = systemData.systemParams.EffISR_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISR_CC = systemData.systemParams.EffISR_CC(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISE_Ore = systemData.systemParams.EffISE_Ore(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISE_Scr = systemData.systemParams.EffISE_Scr(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISE_El = systemData.systemParams.EffISE_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISE_He = systemData.systemParams.EffISE_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISE_CC = systemData.systemParams.EffISE_CC(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAA_Ba = systemData.systemParams.EffIAA_Ba(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAA_Lime = systemData.systemParams.EffIAA_Lime(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAA_So = systemData.systemParams.EffIAA_So(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAA_He = systemData.systemParams.EffIAA_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAA_HeOut = systemData.systemParams.EffIAA_HeOut(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAM_Aa = systemData.systemParams.EffIAM_Aa(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAM_El = systemData.systemParams.EffIAM_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAM_HeOut = systemData.systemParams.EffIAM_HeOut(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAR_Al = systemData.systemParams.EffIAR_Al(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAR_El = systemData.systemParams.EffIAR_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAR_He = systemData.systemParams.EffIAR_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIPP_Wo = systemData.systemParams.EffIPP_Wo(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIPP_El = systemData.systemParams.EffIPP_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIPP_He = systemData.systemParams.EffIPP_He(Regions{numCon},:);

        systemDataTemp{numCon}.systemParams.EffSMC_CO = systemData.systemParams.EffSMC_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffBCS_CO = systemData.systemParams.EffBCS_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffCBC_CO = systemData.systemParams.EffCBC_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffCWC_CO = systemData.systemParams.EffCWC_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffFTB_El = systemData.systemParams.EffFTB_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffFTM_El = systemData.systemParams.EffFTM_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffPSC_CO = systemData.systemParams.EffPSC_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffPSC_El = systemData.systemParams.EffPSC_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffPSC_He = systemData.systemParams.EffPSC_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffPSP_CO = systemData.systemParams.EffPSP_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffPSP_El = systemData.systemParams.EffPSP_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffPSP_He = systemData.systemParams.EffPSP_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRLS_input = systemData.systemParams.EffRLS_input(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRLS_El = systemData.systemParams.EffRLS_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRSS_input = systemData.systemParams.EffRSS_input(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRSS_El = systemData.systemParams.EffRSS_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRMI_input = systemData.systemParams.EffRMI_input(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRMI_Wa = systemData.systemParams.EffRMI_Wa(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRMI_El = systemData.systemParams.EffRMI_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRMO_input = systemData.systemParams.EffRMO_input(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRMO_Wa = systemData.systemParams.EffRMO_Wa(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRMO_El = systemData.systemParams.EffRMO_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRME_input = systemData.systemParams.EffRME_input(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRME_IW = systemData.systemParams.EffRME_IW(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRME_MR = systemData.systemParams.EffRME_MR(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRME_CaCO3out = systemData.systemParams.EffRME_CaCO3out(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRME_MgCO3out = systemData.systemParams.EffRME_MgCO3out(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRME_SiO2out = systemData.systemParams.EffRME_SiO2out(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRME_El = systemData.systemParams.EffRME_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRME_He = systemData.systemParams.EffRME_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRSi_input = systemData.systemParams.EffRSi_input(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRSi_Hy = systemData.systemParams.EffRSi_Hy(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRSi_SiO2 = systemData.systemParams.EffRSi_SiO2(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRSi_SiCout = systemData.systemParams.EffRSi_SiCout(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRSi_El = systemData.systemParams.EffRSi_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRSi_He = systemData.systemParams.EffRSi_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRAD_Wa = systemData.systemParams.EffRAD_Wa(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRAD_El = systemData.systemParams.EffRAD_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffREW_MR = systemData.systemParams.EffREW_MR(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffREW_El = systemData.systemParams.EffREW_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRBC_Bio = systemData.systemParams.EffRBC_Bio(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRBC_Wa = systemData.systemParams.EffRBC_Wa(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRBC_El = systemData.systemParams.EffRBC_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRGE_El = systemData.systemParams.EffRGE_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRGE_He = systemData.systemParams.EffRGE_He(Regions{numCon},:);

        systemDataTemp{numCon}.systemParams.TemperatureAir = systemData.systemParams.TemperatureAir(:,Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.TemperatureWater = systemData.systemParams.TemperatureWater(:,Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.TemperatureGround = systemData.systemParams.TemperatureGround(:,Regions{numCon},:);


        systemDataTemp{numCon}.systemParams.EffHCS_CO = systemData.systemParams.EffHCS_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffGCS_CO = systemData.systemParams.EffGCS_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffSMC_CO = systemData.systemParams.EffSMC_CO(Regions{numCon},:);

        systemDataTemp{numCon}.systemParams.AC_losses = systemData.systemParams.AC_losses(Regions{numCon},:);

        systemDataTemp{numCon}.systemParams.SectorsCons = systemData.systemParams.SectorsCons(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.ElCostRES = systemData.systemParams.ElCostRES(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.ElCostCOM = systemData.systemParams.ElCostCOM(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.ElCostIND = systemData.systemParams.ElCostIND(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.name = num2str(Reg);

        try
            systemDataTemp{numCon}.systemParams.totProfilesAC = systemData.systemParams.totProfilesAC(:,Regions{numCon});
            systemDataTemp{numCon}.systemParams.totProfilesDC = systemData.systemParams.totProfilesDC(:,Regions{numCon});
        catch
        end


        setupTemp{numCon} = setup;

        setupTemp{numCon}.shareOfSNG = setup.shareOfSNG(Regions{numCon},:);

        try
            setupTemp{numCon}.miscParam = setup.miscParam(Regions{numCon});
        catch
            display(['Parameters ' 'miscParam' ' are not availiable'])
        end


        setupTemp{numCon}.SC.share_of_SC = setup.SC.share_of_SC(Regions{numCon},:);
        setupTemp{numCon}.SC.feedin_cost = setup.SC.feedin_cost(Regions{numCon});

        setupTemp{numCon}.SC.Biomass_WastesShare = setup.SC.Biomass_WastesShare(Regions{numCon});
        setupTemp{numCon}.SC.Biomass_BiomassShare = setup.SC.Biomass_BiomassShare(Regions{numCon});
        setupTemp{numCon}.SC.Biomass_BiogasShare = setup.SC.Biomass_BiogasShare(Regions{numCon});
        setupTemp{numCon}.SC.Biomass_WastesAddCost = setup.SC.Biomass_WastesAddCost(:,Regions{numCon});
        setupTemp{numCon}.SC.Biomass_BiomassAddCost = setup.SC.Biomass_BiomassAddCost(:,Regions{numCon});
        setupTemp{numCon}.SC.Biomass_BiogasAddCost = setup.SC.Biomass_BiogasAddCost(:,Regions{numCon});
        setupTemp{numCon}.SC.Biomass_GasAddCost = setup.SC.Biomass_GasAddCost(:,Regions{numCon});
        setupTemp{numCon}.SC.Biomass_OilAddCost = setup.SC.Biomass_OilAddCost(:,Regions{numCon});


        setupTemp{numCon}.Heat.shareOfIndivHeatSectors = setup.Heat.shareOfIndivHeatSectors(Regions{numCon},:);
        setupTemp{numCon}.Heat.shareOfDistrHeat = setup.Heat.shareOfDistrHeat(Regions{numCon},:);
        setupTemp{numCon}.Heat.shareOfLowHeatInd = setup.Heat.shareOfLowHeatInd(Regions{numCon});
        setupTemp{numCon}.Heat.shareOfHighHeatInd = setup.Heat.shareOfHighHeatInd(Regions{numCon});
        setupTemp{numCon}.Heat.DistrHeatEff = setup.Heat.DistrHeatEff(Regions{numCon},:);

        setupTemp{numCon}.DumpChargeShare_LDV = setup.DumpChargeShare_LDV(Regions{numCon},:);
        setupTemp{numCon}.V2GShare_LDV = setup.V2GShare_LDV(Regions{numCon},:);
        setupTemp{numCon}.DumpChargeShare_23W = setup.DumpChargeShare_23W(Regions{numCon},:);
        setupTemp{numCon}.V2GShare_23W = setup.V2GShare_23W(Regions{numCon},:);
        setupTemp{numCon}.DumpChargeShare_BUS = setup.DumpChargeShare_BUS(Regions{numCon},:);
        setupTemp{numCon}.V2GShare_BUS = setup.V2GShare_BUS(Regions{numCon},:);
        setupTemp{numCon}.DumpChargeShare_MDV = setup.DumpChargeShare_MDV(Regions{numCon},:);
        setupTemp{numCon}.V2GShare_MDV = setup.V2GShare_MDV(Regions{numCon},:);
        setupTemp{numCon}.DumpChargeShare_HDV = setup.DumpChargeShare_HDV(Regions{numCon},:);
        setupTemp{numCon}.V2GShare_HDV = setup.V2GShare_HDV(Regions{numCon},:);

        setupTemp{numCon}.importH2Cost = setup.importH2Cost(Regions{numCon},:);
        setupTemp{numCon}.importFTCost = setup.importFTCost(Regions{numCon},:);
        setupTemp{numCon}.importLNGCost = setup.importLNGCost(Regions{numCon},:);
        setupTemp{numCon}.importNH3Cost = setup.importNH3Cost(Regions{numCon},:);
        setupTemp{numCon}.importMeOHCost = setup.importMeOHCost(Regions{numCon},:);

        setupTemp{numCon}.TRAP_Eff = setup.TRAP_Eff(Regions{numCon},:);
        setupTemp{numCon}.TRAD_Eff = setup.TRAD_Eff(Regions{numCon},:);

        setupTemp{numCon}.TICH_CO = setup.TICH_CO(Regions{numCon},:);

        setupTemp{numCon}.groundHPshare_IH = setup.groundHPshare_IH(Regions{numCon},:);
        setupTemp{numCon}.airHPshare_IH = setup.airHPshare_IH(Regions{numCon},:);
        setupTemp{numCon}.waterHPshare_IH = setup.waterHPshare_IH(Regions{numCon},:);
        setupTemp{numCon}.groundHPshare_DH = setup.groundHPshare_DH(Regions{numCon},:);
        setupTemp{numCon}.airHPshare_DH = setup.airHPshare_DH(Regions{numCon},:);
        setupTemp{numCon}.waterHPshare_DH = setup.waterHPshare_DH(Regions{numCon},:);

        setupTemp{numCon}.AmmoniaDemand = setup.AmmoniaDemand(Regions{numCon},:);
        setupTemp{numCon}.GasToChem = setup.GasToChem(Regions{numCon});
        setupTemp{numCon}.OilToChem = setup.OilToChem(Regions{numCon});
        setupTemp{numCon}.CoalToChem = setup.CoalToChem(Regions{numCon});

        cross{numCon} = (ismember(systemData.systemParams.gridMat(:,1),Regions{numCon})&ismember(systemData.systemParams.gridMat(:,2),Regions{numCon}));

        if sum(cross{numCon})<1
            systemDataTemp{numCon}.systemParams.GridMax = 0;
        end
        systemDataTemp{numCon}.systemParams.gridMat = systemData.systemParams.gridMat(cross{numCon},:);
        systemDataTemp{numCon}.systemParams.gridMat = systemDataTemp{numCon}.systemParams.gridMat-Reg+1;
        systemDataTemp{numCon}.systemParams.TLlength = systemData.systemParams.TLlength(cross{numCon});
        systemDataTemp{numCon}.systemParams.TL_AC_LowLimits = systemData.systemParams.TL_AC_LowLimits(cross{numCon});
        systemDataTemp{numCon}.systemParams.TL_AC_UpLimits = systemData.systemParams.TL_AC_UpLimits(cross{numCon});
        systemDataTemp{numCon}.systemParams.TL_DC_LowLimits = systemData.systemParams.TL_DC_LowLimits(cross{numCon});
        systemDataTemp{numCon}.systemParams.TL_DC_UpLimits = systemData.systemParams.TL_DC_UpLimits(cross{numCon});
        systemDataTemp{numCon}.Reg1=Reg;
        systemDataTemp{numCon}.Reg2=Reg2;

        UReg{numCon} = unique([systemDataTemp{numCon}.systemParams.gridMat(:,2);systemDataTemp{numCon}.systemParams.gridMat(:,1)]);


        for i = 1:length(UReg{numCon})

            systemDataTemp{numCon}.systemParams.gridMat(systemDataTemp{numCon}.systemParams.gridMat==UReg{numCon}(i))=i;

        end
        try
            PrepareScenarioForTransition(systemDataTemp{numCon},setupTemp{numCon});

        catch ME
            warning('Problem with the %s region', num2str(Reg),'', ME.message);

        end
    end
else


    for numCon = startReg:size(setup.Countries,1)

        if iscell(setup.Countries)
            Reg = setup.Countries{numCon}(1);
            Reg2 = setup.Countries{numCon}(end);

            Regions{numCon} = setup.Countries{numCon};

        else
            Reg = setup.Countries(numCon,1);
            Reg2 = setup.Countries(numCon,2);

            Regions{numCon} = [Reg:Reg2];

        end

        systemDataTemp{numCon} = systemData;
        systemDataTemp{numCon}.systemParams = systemData.systemParams;
        systemDataTemp{numCon}.systemParams.Regions = Regions{numCon};
        systemDataTemp{numCon}.systemParams.area = systemData.systemParams.area(Regions{numCon});
        systemDataTemp{numCon}.systemParams.population = systemData.systemParams.population(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.IndexNodes = systemData.systemParams.IndexNodes(Regions{numCon});
        systemDataTemp{numCon}.systemParams.IndexNumNodes = systemData.systemParams.IndexNumNodes(Regions{numCon});
        systemDataTemp{numCon}.systemParams.Coords = systemData.systemParams.Coords(:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.Urbanisation = systemData.systemParams.Urbanisation(:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.SizeLimits = systemData.systemParams.SizeLimits(:,:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.Instalations = systemData.systemParams.Instalations(Regions{numCon},:,:);
        systemDataTemp{numCon}.systemParams.ValueLoad = systemData.systemParams.ValueLoad(:,:,:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.ValueLoadTag = systemData.systemParams.ValueLoadTag(:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.ValueResource = systemData.systemParams.ValueResource(:,:,:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.ValueHydroDam = systemData.systemParams.ValueHydroDam(:,:,:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.ValueResourceTag = systemData.systemParams.ValueResourceTag(:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.ValueHydroDamTag = systemData.systemParams.ValueHydroDamTag(:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.Capex_reg = systemData.systemParams.Capex_reg(:,Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Opex_fix_reg = systemData.systemParams.Opex_fix_reg(:,Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Opex_var_reg = systemData.systemParams.Opex_var_reg(:,Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.ValueResourceTotal = systemData.systemParams.ValueResourceTotal(:,Regions{numCon});
        systemDataTemp{numCon}.systemParams.DesalinationParams.regionPumpLong = systemData.systemParams.DesalinationParams.regionPumpLong(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.DesalinationParams.regionPumpUp = systemData.systemParams.DesalinationParams.regionPumpUp(Regions{numCon},:);

        %if setup.Mobility
        systemDataTemp{numCon}.systemParams.Mobility.Demand = systemData.systemParams.Mobility.Demand(:,Regions{numCon},:);

        systemDataTemp{numCon}.systemParams.Mobility.LDV_Pass = systemData.systemParams.Mobility.LDV_Pass(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.BUS_Pass = systemData.systemParams.Mobility.BUS_Pass(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.W23_Pass = systemData.systemParams.Mobility.W23_Pass(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.LDV_km = systemData.systemParams.Mobility.LDV_km(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.BUS_km = systemData.systemParams.Mobility.BUS_km(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.W23_km = systemData.systemParams.Mobility.W23_km(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.LDV_PHEV_el = systemData.systemParams.Mobility.LDV_PHEV_el(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.BUS_PHEV_el = systemData.systemParams.Mobility.BUS_PHEV_el(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.W23_PHEV_el = systemData.systemParams.Mobility.W23_PHEV_el(Regions{numCon},:);

        systemDataTemp{numCon}.systemParams.Mobility.MDV_Tonne = systemData.systemParams.Mobility.MDV_Tonne(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.HDV_Tonne = systemData.systemParams.Mobility.HDV_Tonne(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.MDV_km = systemData.systemParams.Mobility.MDV_km(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.HDV_km = systemData.systemParams.Mobility.HDV_km(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.MDV_PHEV_el = systemData.systemParams.Mobility.MDV_PHEV_el(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.Mobility.HDV_PHEV_el = systemData.systemParams.Mobility.HDV_PHEV_el(Regions{numCon},:);

        systemDataTemp{numCon}.systemParams.Mobility.Shares = systemData.systemParams.Mobility.Shares(Regions{numCon},:,:);
        %end
        systemDataTemp{numCon}.systemParams.EffCO2Scr_El = systemData.systemParams.EffCO2Scr_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffCO2Scr_He = systemData.systemParams.EffCO2Scr_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffMETH_CO = systemData.systemParams.EffMETH_CO(Regions{numCon},:);
        try
            systemDataTemp{numCon}.systemParams.EffFT_El = systemData.systemParams.EffFT_El(Regions{numCon},:);
        end
        systemDataTemp{numCon}.systemParams.EffFT_CO = systemData.systemParams.EffFT_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffFT_He = systemData.systemParams.EffFT_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffLNG_El = systemData.systemParams.EffLNG_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffLH2_El = systemData.systemParams.EffLH2_El(Regions{numCon},:);




        systemDataTemp{numCon}.systemParams.EffMeO_CO = systemData.systemParams.EffMeO_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffMeO_El = systemData.systemParams.EffMeO_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffMeO_He_out = systemData.systemParams.EffMeO_He_out(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffDME_CO = systemData.systemParams.EffDME_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffDME_El = systemData.systemParams.EffDME_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffDME_He_out = systemData.systemParams.EffDME_He_out(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffNH3_El = systemData.systemParams.EffNH3_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffNH3_He_out = systemData.systemParams.EffNH3_He_out(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffHyS_El = systemData.systemParams.EffHyS_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffMET_El = systemData.systemParams.EffMET_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffMET_He_out = systemData.systemParams.EffMET_He_out(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffWEL_He_out = systemData.systemParams.EffWEL_He_out(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffICB_Lime = systemData.systemParams.EffICB_Lime(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffICB_El = systemData.systemParams.EffICB_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffICB_He = systemData.systemParams.EffICB_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffICI_Lime = systemData.systemParams.EffICI_Lime(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffICI_El = systemData.systemParams.EffICI_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffICI_He = systemData.systemParams.EffICI_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISB_Ore = systemData.systemParams.EffISB_Ore(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISB_Scr = systemData.systemParams.EffISB_Scr(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISB_El = systemData.systemParams.EffISB_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISB_He = systemData.systemParams.EffISB_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISB_HC = systemData.systemParams.EffISB_HC(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISH_Ore = systemData.systemParams.EffISH_Ore(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISH_Scr = systemData.systemParams.EffISH_Scr(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISH_El = systemData.systemParams.EffISH_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISH_He = systemData.systemParams.EffISH_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISH_CC = systemData.systemParams.EffISH_CC(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISH_Hy = systemData.systemParams.EffISH_Hy(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISR_Scr = systemData.systemParams.EffISR_Scr(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISR_El = systemData.systemParams.EffISR_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISR_He = systemData.systemParams.EffISR_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISR_CC = systemData.systemParams.EffISR_CC(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISE_Ore = systemData.systemParams.EffISE_Ore(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISE_Scr = systemData.systemParams.EffISE_Scr(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISE_El = systemData.systemParams.EffISE_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISE_He = systemData.systemParams.EffISE_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffISE_CC = systemData.systemParams.EffISE_CC(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAA_Ba = systemData.systemParams.EffIAA_Ba(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAA_Lime = systemData.systemParams.EffIAA_Lime(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAA_So = systemData.systemParams.EffIAA_So(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAA_He = systemData.systemParams.EffIAA_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAA_HeOut = systemData.systemParams.EffIAA_HeOut(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAM_Aa = systemData.systemParams.EffIAM_Aa(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAM_El = systemData.systemParams.EffIAM_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAM_HeOut = systemData.systemParams.EffIAM_HeOut(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAR_Al = systemData.systemParams.EffIAR_Al(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAR_El = systemData.systemParams.EffIAR_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIAR_He = systemData.systemParams.EffIAR_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIPP_Wo = systemData.systemParams.EffIPP_Wo(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIPP_El = systemData.systemParams.EffIPP_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffIPP_He = systemData.systemParams.EffIPP_He(Regions{numCon},:);

        systemDataTemp{numCon}.systemParams.EffSMC_CO = systemData.systemParams.EffSMC_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffBCS_CO = systemData.systemParams.EffBCS_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffCBC_CO = systemData.systemParams.EffCBC_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffCWC_CO = systemData.systemParams.EffCWC_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffFTB_El = systemData.systemParams.EffFTB_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffFTM_El = systemData.systemParams.EffFTM_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffPSC_CO = systemData.systemParams.EffPSC_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffPSC_El = systemData.systemParams.EffPSC_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffPSC_He = systemData.systemParams.EffPSC_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffPSP_CO = systemData.systemParams.EffPSP_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffPSP_El = systemData.systemParams.EffPSP_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffPSP_He = systemData.systemParams.EffPSP_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRLS_input = systemData.systemParams.EffRLS_input(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRLS_El = systemData.systemParams.EffRLS_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRSS_input = systemData.systemParams.EffRSS_input(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRSS_El = systemData.systemParams.EffRSS_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRMI_input = systemData.systemParams.EffRMI_input(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRMI_Wa = systemData.systemParams.EffRMI_Wa(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRMI_El = systemData.systemParams.EffRMI_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRMO_input = systemData.systemParams.EffRMO_input(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRMO_Wa = systemData.systemParams.EffRMO_Wa(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRMO_El = systemData.systemParams.EffRMO_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRME_input = systemData.systemParams.EffRME_input(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRME_IW = systemData.systemParams.EffRME_IW(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRME_MR = systemData.systemParams.EffRME_MR(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRME_CaCO3out = systemData.systemParams.EffRME_CaCO3out(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRME_MgCO3out = systemData.systemParams.EffRME_MgCO3out(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRME_SiO2out = systemData.systemParams.EffRME_SiO2out(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRME_El = systemData.systemParams.EffRME_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRME_He = systemData.systemParams.EffRME_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRSi_input = systemData.systemParams.EffRSi_input(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRSi_Hy = systemData.systemParams.EffRSi_Hy(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRSi_SiO2 = systemData.systemParams.EffRSi_SiO2(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRSi_SiCout = systemData.systemParams.EffRSi_SiCout(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRSi_El = systemData.systemParams.EffRSi_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRSi_He = systemData.systemParams.EffRSi_He(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRAD_Wa = systemData.systemParams.EffRAD_Wa(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRAD_El = systemData.systemParams.EffRAD_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffREW_MR = systemData.systemParams.EffREW_MR(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffREW_El = systemData.systemParams.EffREW_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRBC_Bio = systemData.systemParams.EffRBC_Bio(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRBC_Wa = systemData.systemParams.EffRBC_Wa(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRBC_El = systemData.systemParams.EffRBC_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRGE_El = systemData.systemParams.EffRGE_El(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffRGE_He = systemData.systemParams.EffRGE_He(Regions{numCon},:);


        systemDataTemp{numCon}.systemParams.TemperatureAir = systemData.systemParams.TemperatureAir(:,Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.TemperatureWater = systemData.systemParams.TemperatureWater(:,Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.TemperatureGround = systemData.systemParams.TemperatureGround(:,Regions{numCon},:);


        systemDataTemp{numCon}.systemParams.EffHCS_CO = systemData.systemParams.EffHCS_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffGCS_CO = systemData.systemParams.EffGCS_CO(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.EffSMC_CO = systemData.systemParams.EffSMC_CO(Regions{numCon},:);

        systemDataTemp{numCon}.systemParams.AC_losses = systemData.systemParams.AC_losses(Regions{numCon},:);

        systemDataTemp{numCon}.systemParams.SectorsCons = systemData.systemParams.SectorsCons(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.ElCostRES = systemData.systemParams.ElCostRES(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.ElCostCOM = systemData.systemParams.ElCostCOM(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.ElCostIND = systemData.systemParams.ElCostIND(Regions{numCon},:);
        systemDataTemp{numCon}.systemParams.name = num2str(Reg);

        try
            systemDataTemp{numCon}.systemParams.totProfilesAC = systemData.systemParams.totProfilesAC(:,Regions{numCon});
            systemDataTemp{numCon}.systemParams.totProfilesDC = systemData.systemParams.totProfilesDC(:,Regions{numCon});
        catch
        end

        setupTemp{numCon} = setup;

        setupTemp{numCon}.shareOfSNG = setup.shareOfSNG(Regions{numCon},:);


        try
            setupTemp{numCon}.miscParam = setup.miscParam(Regions{numCon});
        catch
            display(['Parameters ' 'miscParam' ' are not availiable'])
        end


        setupTemp{numCon}.SC.share_of_SC = setup.SC.share_of_SC(Regions{numCon},:);
        setupTemp{numCon}.SC.feedin_cost = setup.SC.feedin_cost(Regions{numCon});

        setupTemp{numCon}.SC.Biomass_WastesShare = setup.SC.Biomass_WastesShare(Regions{numCon});
        setupTemp{numCon}.SC.Biomass_BiomassShare = setup.SC.Biomass_BiomassShare(Regions{numCon});
        setupTemp{numCon}.SC.Biomass_BiogasShare = setup.SC.Biomass_BiogasShare(Regions{numCon});
        setupTemp{numCon}.SC.Biomass_WastesAddCost = setup.SC.Biomass_WastesAddCost(:,Regions{numCon});
        setupTemp{numCon}.SC.Biomass_BiomassAddCost = setup.SC.Biomass_BiomassAddCost(:,Regions{numCon});
        setupTemp{numCon}.SC.Biomass_BiogasAddCost = setup.SC.Biomass_BiogasAddCost(:,Regions{numCon});
        setupTemp{numCon}.SC.Biomass_GasAddCost = setup.SC.Biomass_GasAddCost(:,Regions{numCon});
        setupTemp{numCon}.SC.Biomass_OilAddCost = setup.SC.Biomass_OilAddCost(:,Regions{numCon});


        setupTemp{numCon}.Heat.shareOfIndivHeatSectors = setup.Heat.shareOfIndivHeatSectors(Regions{numCon},:);
        setupTemp{numCon}.Heat.shareOfDistrHeat = setup.Heat.shareOfDistrHeat(Regions{numCon},:);
        setupTemp{numCon}.Heat.shareOfLowHeatInd = setup.Heat.shareOfLowHeatInd(Regions{numCon});
        setupTemp{numCon}.Heat.shareOfHighHeatInd = setup.Heat.shareOfHighHeatInd(Regions{numCon});
        setupTemp{numCon}.Heat.DistrHeatEff = setup.Heat.DistrHeatEff(Regions{numCon},:);

        setupTemp{numCon}.DumpChargeShare_LDV = setup.DumpChargeShare_LDV(Regions{numCon},:);
        setupTemp{numCon}.V2GShare_LDV = setup.V2GShare_LDV(Regions{numCon},:);
        setupTemp{numCon}.DumpChargeShare_23W = setup.DumpChargeShare_23W(Regions{numCon},:);
        setupTemp{numCon}.V2GShare_23W = setup.V2GShare_23W(Regions{numCon},:);
        setupTemp{numCon}.DumpChargeShare_BUS = setup.DumpChargeShare_BUS(Regions{numCon},:);
        setupTemp{numCon}.V2GShare_BUS = setup.V2GShare_BUS(Regions{numCon},:);
        setupTemp{numCon}.DumpChargeShare_MDV = setup.DumpChargeShare_MDV(Regions{numCon},:);
        setupTemp{numCon}.V2GShare_MDV = setup.V2GShare_MDV(Regions{numCon},:);
        setupTemp{numCon}.DumpChargeShare_HDV = setup.DumpChargeShare_HDV(Regions{numCon},:);
        setupTemp{numCon}.V2GShare_HDV = setup.V2GShare_HDV(Regions{numCon},:);


        setupTemp{numCon}.importH2Cost = setup.importH2Cost(Regions{numCon},:);
        setupTemp{numCon}.importLNGCost = setup.importLNGCost(Regions{numCon},:);
        setupTemp{numCon}.importFTCost = setup.importFTCost(Regions{numCon},:);
        setupTemp{numCon}.importNH3Cost = setup.importNH3Cost(Regions{numCon},:);
        setupTemp{numCon}.importMeOHCost = setup.importMeOHCost(Regions{numCon},:);

        setupTemp{numCon}.TRAP_Eff = setup.TRAP_Eff(Regions{numCon},:);
        setupTemp{numCon}.TRAD_Eff = setup.TRAD_Eff(Regions{numCon},:);

        setupTemp{numCon}.TICH_CO = setup.TICH_CO(Regions{numCon},:);

        setupTemp{numCon}.groundHPshare_IH = setup.groundHPshare_IH(Regions{numCon},:);
        setupTemp{numCon}.airHPshare_IH = setup.airHPshare_IH(Regions{numCon},:);
        setupTemp{numCon}.waterHPshare_IH = setup.waterHPshare_IH(Regions{numCon},:);
        setupTemp{numCon}.groundHPshare_DH = setup.groundHPshare_DH(Regions{numCon},:);
        setupTemp{numCon}.airHPshare_DH = setup.airHPshare_DH(Regions{numCon},:);
        setupTemp{numCon}.waterHPshare_DH = setup.waterHPshare_DH(Regions{numCon},:);

        setupTemp{numCon}.AmmoniaDemand = setup.AmmoniaDemand(Regions{numCon},:);
        setupTemp{numCon}.GasToChem = setup.GasToChem(Regions{numCon});
        setupTemp{numCon}.OilToChem = setup.OilToChem(Regions{numCon});
        setupTemp{numCon}.CoalToChem = setup.CoalToChem(Regions{numCon});

        cross{numCon} = (ismember(systemData.systemParams.gridMat(:,1),Regions{numCon})&ismember(systemData.systemParams.gridMat(:,2),Regions{numCon}));

        if sum(cross{numCon})<1
            systemDataTemp{numCon}.systemParams.GridMax = 0;
        end
        systemDataTemp{numCon}.systemParams.gridMat = systemData.systemParams.gridMat(cross{numCon},:);
        systemDataTemp{numCon}.systemParams.gridMat = systemDataTemp{numCon}.systemParams.gridMat-Reg+1;
        systemDataTemp{numCon}.systemParams.TLlength = systemData.systemParams.TLlength(cross{numCon});
        systemDataTemp{numCon}.systemParams.TL_AC_LowLimits = systemData.systemParams.TL_AC_LowLimits(cross{numCon});
        systemDataTemp{numCon}.systemParams.TL_AC_UpLimits = systemData.systemParams.TL_AC_UpLimits(cross{numCon});
        systemDataTemp{numCon}.systemParams.TL_DC_LowLimits = systemData.systemParams.TL_DC_LowLimits(cross{numCon});
        systemDataTemp{numCon}.systemParams.TL_DC_UpLimits = systemData.systemParams.TL_DC_UpLimits(cross{numCon});
        systemDataTemp{numCon}.Reg1=Reg;
        systemDataTemp{numCon}.Reg2=Reg2;

        UReg{numCon} = unique([systemDataTemp{numCon}.systemParams.gridMat(:,2);systemDataTemp{numCon}.systemParams.gridMat(:,1)]);


        for i = 1:length(UReg{numCon})

            systemDataTemp{numCon}.systemParams.gridMat(systemDataTemp{numCon}.systemParams.gridMat==UReg{numCon}(i))=i;

        end
        % try

            PrepareScenarioForTransition(systemDataTemp{numCon},setupTemp{numCon});
        % catch ME
        %     warning('Problem with the %s region', num2str(Reg),'', ME.message);
        % 
        % end
    end
end