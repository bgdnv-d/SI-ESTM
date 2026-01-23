function collectResults_SPE_R(setup,Countries,costYear,name)
%
% FUNCTION collectResults_SPE_R(setup, Countries, costYear, name)
%
%
% INPUT:
%            setup:      Structure with all necessary settings and input data.
%            Countries:  List of countries to include in the results.
%            costYear:   Year to which all costs should be adjusted.
%            name:       Label or identifier used for saving or organizing results.
%
% OUTPUT:
%            -
%
%Dmitrii Bogdanov
%last change 22.07.2025


pName = [setup.projType];

if setup.MacroMc
    pName = [pName '_Mc'];
end

if setup.MacroReg
    pName = [pName '_R'];
end

if setup.OvernightFlag
    pName = [pName '_Overnight'];
end

if setup.SC.Flag
    pName = [pName '_EL_SC'];
end

if setup.Heat.Flag
    pName = [pName '_HE'];
end

if setup.Mobility
    pName = [pName '_TR'];
end

if setup.IndustryFlag
    pName = [pName '_IND'];
end

if setup.GasFlag
    pName = [pName '_GAS'];
end

if setup.DesalinationFlag
    pName = [pName '_DES'];
end

if setup.Only.Flag
    pName = [pName '_Only'];
end
bbData = load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep 'simulation-input_' 'Base' '.mat']);

baseData = bbData;
activeElements = bbData.activeElements;
modFiles = bbData.modFiles;
datFiles = bbData.datFiles;
RegionsAll = [];


if strcmp(setup.projType, 'Regions')
    setup.Countries = [bbData.systemParams.IndexNumNodes, bbData.systemParams.IndexNumNodes];
    Countries = setup.Countries;
end

if iscell(setup.Countries)
    for i=1:size(Countries,1)
        RegionsAll = [RegionsAll [Countries{i}]];
    end
    RegionsAll = sort(RegionsAll);
else
    for i=1:size(Countries,1)
        RegionsAll = [RegionsAll [Countries(i,1):Countries(i,2)]];
    end
end


systemParams = baseData.systemParams;

systemParams.area = baseData.systemParams.area(RegionsAll);
systemParams.population = baseData.systemParams.population(RegionsAll,:);
systemParams.IndexNodes = baseData.systemParams.IndexNodes(RegionsAll);
systemParams.IndexNumNodes = baseData.systemParams.IndexNumNodes(RegionsAll);
systemParams.Coords = baseData.systemParams.Coords(:,RegionsAll);
systemParams.Urbanisation = baseData.systemParams.Urbanisation(:,RegionsAll);
systemParams.SizeLimits = baseData.systemParams.SizeLimits(:,:,RegionsAll);
systemParams.Instalations = baseData.systemParams.Instalations(RegionsAll,:,:);
systemParams.ValueLoad = baseData.systemParams.ValueLoad(:,:,:,RegionsAll);
systemParams.ValueLoadTag = baseData.systemParams.ValueLoadTag(:,RegionsAll);
systemParams.ValueResource = baseData.systemParams.ValueResource(:,:,:,RegionsAll);
systemParams.ValueHydroDam = baseData.systemParams.ValueHydroDam(:,:,:,RegionsAll);
systemParams.Capex_reg = baseData.systemParams.Capex_reg(:,RegionsAll,:);
systemParams.Opex_fix_reg = baseData.systemParams.Opex_fix_reg(:,RegionsAll,:);
systemParams.Opex_var_reg = baseData.systemParams.Opex_var_reg(:,RegionsAll,:);
systemParams.ValueResourceTotal = baseData.systemParams.ValueResourceTotal(:,RegionsAll);
systemParams.DesalinationParams.regionPumpLong = baseData.systemParams.DesalinationParams.regionPumpLong(RegionsAll,:);
systemParams.DesalinationParams.regionPumpUp = baseData.systemParams.DesalinationParams.regionPumpUp(RegionsAll,:);


%if setup.Mobility
systemParams.Mobility.Demand = baseData.systemParams.Mobility.Demand(:,RegionsAll,:);

systemParams.Mobility.LDV_Pass = baseData.systemParams.Mobility.LDV_Pass(RegionsAll,:);
systemParams.Mobility.BUS_Pass = baseData.systemParams.Mobility.BUS_Pass(RegionsAll,:);
systemParams.Mobility.W23_Pass = baseData.systemParams.Mobility.W23_Pass(RegionsAll,:);
systemParams.Mobility.LDV_km = baseData.systemParams.Mobility.LDV_km(RegionsAll,:);
systemParams.Mobility.BUS_km = baseData.systemParams.Mobility.BUS_km(RegionsAll,:);
systemParams.Mobility.W23_km = baseData.systemParams.Mobility.W23_km(RegionsAll,:);
systemParams.Mobility.LDV_PHEV_el = baseData.systemParams.Mobility.LDV_PHEV_el(RegionsAll,:);
systemParams.Mobility.BUS_PHEV_el = baseData.systemParams.Mobility.BUS_PHEV_el(RegionsAll,:);
systemParams.Mobility.W23_PHEV_el = baseData.systemParams.Mobility.W23_PHEV_el(RegionsAll,:);

systemParams.Mobility.MDV_Tonne = baseData.systemParams.Mobility.MDV_Tonne(RegionsAll,:);
systemParams.Mobility.HDV_Tonne = baseData.systemParams.Mobility.HDV_Tonne(RegionsAll,:);
systemParams.Mobility.MDV_km = baseData.systemParams.Mobility.MDV_km(RegionsAll,:);
systemParams.Mobility.HDV_km = baseData.systemParams.Mobility.HDV_km(RegionsAll,:);
systemParams.Mobility.MDV_PHEV_el = baseData.systemParams.Mobility.MDV_PHEV_el(RegionsAll,:);
systemParams.Mobility.HDV_PHEV_el = baseData.systemParams.Mobility.HDV_PHEV_el(RegionsAll,:);

systemParams.Mobility.Shares = baseData.systemParams.Mobility.Shares(RegionsAll,:,:);



systemParams.Mobility.Cons_Primary_S = zeros(length(RegionsAll),size(systemParams.IndexYears,2),size(systemParams.IndexIDM,1));
systemParams.Mobility.Cons_Secondary_S = zeros(length(RegionsAll),size(systemParams.IndexYears,2),5);




%end
systemParams.EffCO2Scr_El = baseData.systemParams.EffCO2Scr_El(RegionsAll,:);
systemParams.EffCO2Scr_He = baseData.systemParams.EffCO2Scr_He(RegionsAll,:);
systemParams.EffMETH_CO = baseData.systemParams.EffMETH_CO(RegionsAll,:);

systemParams.EffFT_CO = baseData.systemParams.EffFT_CO(RegionsAll,:);
systemParams.EffFT_He = baseData.systemParams.EffFT_He(RegionsAll,:);
systemParams.EffLNG_El = baseData.systemParams.EffLNG_El(RegionsAll,:);
systemParams.EffLH2_El = baseData.systemParams.EffLH2_El(RegionsAll,:);

%if setup.New
systemParams.EffMeO_CO = baseData.systemParams.EffMeO_CO(RegionsAll,:);
systemParams.EffMeO_El = baseData.systemParams.EffMeO_El(RegionsAll,:);
systemParams.EffMeO_He_out = baseData.systemParams.EffMeO_He_out(RegionsAll,:);
systemParams.EffDME_CO = baseData.systemParams.EffDME_CO(RegionsAll,:);
systemParams.EffDME_El = baseData.systemParams.EffDME_El(RegionsAll,:);
systemParams.EffDME_He_out = baseData.systemParams.EffDME_He_out(RegionsAll,:);
systemParams.EffNH3_El = baseData.systemParams.EffNH3_El(RegionsAll,:);
systemParams.EffNH3_He_out = baseData.systemParams.EffNH3_He_out(RegionsAll,:);
systemParams.EffHyS_El = baseData.systemParams.EffHyS_El(RegionsAll,:);
systemParams.EffMET_El = baseData.systemParams.EffMET_El(RegionsAll,:);
systemParams.EffMET_He_out = baseData.systemParams.EffMET_He_out(RegionsAll,:);
systemParams.EffWEL_He_out = baseData.systemParams.EffWEL_He_out(RegionsAll,:);
systemParams.EffICB_Lime = baseData.systemParams.EffICB_Lime(RegionsAll,:);
systemParams.EffICB_El = baseData.systemParams.EffICB_El(RegionsAll,:);
systemParams.EffICB_He = baseData.systemParams.EffICB_He(RegionsAll,:);
systemParams.EffICI_Lime = baseData.systemParams.EffICI_Lime(RegionsAll,:);
systemParams.EffICI_El = baseData.systemParams.EffICI_El(RegionsAll,:);
systemParams.EffICI_He = baseData.systemParams.EffICI_He(RegionsAll,:);
systemParams.EffISB_Ore = baseData.systemParams.EffISB_Ore(RegionsAll,:);
systemParams.EffISB_Scr = baseData.systemParams.EffISB_Scr(RegionsAll,:);
systemParams.EffISB_El = baseData.systemParams.EffISB_El(RegionsAll,:);
systemParams.EffISB_He = baseData.systemParams.EffISB_He(RegionsAll,:);
systemParams.EffISB_HC = baseData.systemParams.EffISB_HC(RegionsAll,:);
systemParams.EffISH_Ore = baseData.systemParams.EffISH_Ore(RegionsAll,:);
systemParams.EffISH_Scr = baseData.systemParams.EffISH_Scr(RegionsAll,:);
systemParams.EffISH_El = baseData.systemParams.EffISH_El(RegionsAll,:);
systemParams.EffISH_He = baseData.systemParams.EffISH_He(RegionsAll,:);
systemParams.EffISH_CC = baseData.systemParams.EffISH_CC(RegionsAll,:);
systemParams.EffISH_Hy = baseData.systemParams.EffISH_Hy(RegionsAll,:);
systemParams.EffISR_Scr = baseData.systemParams.EffISR_Scr(RegionsAll,:);
systemParams.EffISR_El = baseData.systemParams.EffISR_El(RegionsAll,:);
systemParams.EffISR_He = baseData.systemParams.EffISR_He(RegionsAll,:);
systemParams.EffISR_CC = baseData.systemParams.EffISR_CC(RegionsAll,:);
systemParams.EffISE_Ore = baseData.systemParams.EffISE_Ore(RegionsAll,:);
systemParams.EffISE_Scr = baseData.systemParams.EffISE_Scr(RegionsAll,:);
systemParams.EffISE_El = baseData.systemParams.EffISE_El(RegionsAll,:);
systemParams.EffISE_He = baseData.systemParams.EffISE_He(RegionsAll,:);
systemParams.EffISE_CC = baseData.systemParams.EffISE_CC(RegionsAll,:);
systemParams.EffIAA_Ba = baseData.systemParams.EffIAA_Ba(RegionsAll,:);
systemParams.EffIAA_Lime = baseData.systemParams.EffIAA_Lime(RegionsAll,:);
systemParams.EffIAA_So = baseData.systemParams.EffIAA_So(RegionsAll,:);
systemParams.EffIAA_He = baseData.systemParams.EffIAA_He(RegionsAll,:);
systemParams.EffIAA_HeOut = baseData.systemParams.EffIAA_HeOut(RegionsAll,:);
systemParams.EffIAM_Aa = baseData.systemParams.EffIAM_Aa(RegionsAll,:);
systemParams.EffIAM_El = baseData.systemParams.EffIAM_El(RegionsAll,:);
systemParams.EffIAM_HeOut = baseData.systemParams.EffIAM_HeOut(RegionsAll,:);
systemParams.EffIAR_Al = baseData.systemParams.EffIAR_Al(RegionsAll,:);
systemParams.EffIAR_El = baseData.systemParams.EffIAR_El(RegionsAll,:);
systemParams.EffIAR_He = baseData.systemParams.EffIAR_He(RegionsAll,:);
systemParams.EffIPP_Wo = baseData.systemParams.EffIPP_Wo(RegionsAll,:);
systemParams.EffIPP_El = baseData.systemParams.EffIPP_El(RegionsAll,:);
systemParams.EffIPP_He = baseData.systemParams.EffIPP_He(RegionsAll,:);
%end



systemParams.EffHCS_CO = baseData.systemParams.EffHCS_CO(RegionsAll,:);
systemParams.EffGCS_CO = baseData.systemParams.EffGCS_CO(RegionsAll,:);

systemParams.EffSMC_CO = baseData.systemParams.EffSMC_CO(RegionsAll,:);
systemParams.EffBCS_CO = baseData.systemParams.EffBCS_CO(RegionsAll,:);
systemParams.EffCBC_CO = baseData.systemParams.EffCBC_CO(RegionsAll,:);
systemParams.EffCWC_CO = baseData.systemParams.EffCWC_CO(RegionsAll,:);
systemParams.EffFTB_El = baseData.systemParams.EffFTB_El(RegionsAll,:);
systemParams.EffFTM_El = baseData.systemParams.EffFTM_El(RegionsAll,:);
systemParams.EffPSC_CO = baseData.systemParams.EffPSC_CO(RegionsAll,:);
systemParams.EffPSC_El = baseData.systemParams.EffPSC_El(RegionsAll,:);
systemParams.EffPSC_He = baseData.systemParams.EffPSC_He(RegionsAll,:);
systemParams.EffPSP_CO = baseData.systemParams.EffPSP_CO(RegionsAll,:);
systemParams.EffPSP_El = baseData.systemParams.EffPSP_El(RegionsAll,:);
systemParams.EffPSP_He = baseData.systemParams.EffPSP_He(RegionsAll,:);
systemParams.EffRLS_input = baseData.systemParams.EffRLS_input(RegionsAll,:);
systemParams.EffRLS_El = baseData.systemParams.EffRLS_El(RegionsAll,:);
systemParams.EffRSS_input = baseData.systemParams.EffRSS_input(RegionsAll,:);
systemParams.EffRSS_El = baseData.systemParams.EffRSS_El(RegionsAll,:);
systemParams.EffRMI_input = baseData.systemParams.EffRMI_input(RegionsAll,:);
systemParams.EffRMI_Wa = baseData.systemParams.EffRMI_Wa(RegionsAll,:);
systemParams.EffRMI_El = baseData.systemParams.EffRMI_El(RegionsAll,:);
systemParams.EffRMO_input = baseData.systemParams.EffRMO_input(RegionsAll,:);
systemParams.EffRMO_Wa = baseData.systemParams.EffRMO_Wa(RegionsAll,:);
systemParams.EffRMO_El = baseData.systemParams.EffRMO_El(RegionsAll,:);
systemParams.EffRME_input = baseData.systemParams.EffRME_input(RegionsAll,:);
systemParams.EffRME_IW = baseData.systemParams.EffRME_IW(RegionsAll,:);
systemParams.EffRME_MR = baseData.systemParams.EffRME_MR(RegionsAll,:);
systemParams.EffRME_CaCO3out = baseData.systemParams.EffRME_CaCO3out(RegionsAll,:);
systemParams.EffRME_MgCO3out = baseData.systemParams.EffRME_MgCO3out(RegionsAll,:);
systemParams.EffRME_SiO2out = baseData.systemParams.EffRME_SiO2out(RegionsAll,:);
systemParams.EffRME_El = baseData.systemParams.EffRME_El(RegionsAll,:);
systemParams.EffRME_He = baseData.systemParams.EffRME_He(RegionsAll,:);
systemParams.EffRSi_input = baseData.systemParams.EffRSi_input(RegionsAll,:);
systemParams.EffRSi_Hy = baseData.systemParams.EffRSi_Hy(RegionsAll,:);
systemParams.EffRSi_SiO2 = baseData.systemParams.EffRSi_SiO2(RegionsAll,:);
systemParams.EffRSi_SiCout = baseData.systemParams.EffRSi_SiCout(RegionsAll,:);
systemParams.EffRSi_El = baseData.systemParams.EffRSi_El(RegionsAll,:);
systemParams.EffRSi_He = baseData.systemParams.EffRSi_He(RegionsAll,:);
systemParams.EffRAD_Wa = baseData.systemParams.EffRAD_Wa(RegionsAll,:);
systemParams.EffRAD_El = baseData.systemParams.EffRAD_El(RegionsAll,:);
systemParams.EffREW_MR = baseData.systemParams.EffREW_MR(RegionsAll,:);
systemParams.EffREW_El = baseData.systemParams.EffREW_El(RegionsAll,:);
systemParams.EffRBC_Bio = baseData.systemParams.EffRBC_Bio(RegionsAll,:);
systemParams.EffRBC_Wa = baseData.systemParams.EffRBC_Wa(RegionsAll,:);
systemParams.EffRBC_El = baseData.systemParams.EffRBC_El(RegionsAll,:);
systemParams.EffRGE_El = baseData.systemParams.EffRGE_El(RegionsAll,:);
systemParams.EffRGE_He = baseData.systemParams.EffRGE_He(RegionsAll,:);

systemParams.AC_losses = baseData.systemParams.AC_losses(RegionsAll,:);

systemParams.SectorsCons = baseData.systemParams.SectorsCons(RegionsAll,:);
systemParams.ElCostRES = baseData.systemParams.ElCostRES(RegionsAll,:);
systemParams.ElCostCOM = baseData.systemParams.ElCostCOM(RegionsAll,:);
systemParams.ElCostIND = baseData.systemParams.ElCostIND(RegionsAll,:);
systemParams.name = name;

if ~strcmp(setup.projType,'Regions')
    try
        systemParams.totProfilesAC = baseData.systemParams.totProfilesAC(:,RegionsAll,:);
        systemParams.totProfilesDC = baseData.systemParams.totProfilesDC(:,RegionsAll,:);
    catch
    end

    cross = (ismember(baseData.systemParams.gridMat(:,1),RegionsAll)&ismember(baseData.systemParams.gridMat(:,2),RegionsAll));
    if sum(cross)<1
        systemParams.GridMax = 0;
    end
    systemParams.gridMat = baseData.systemParams.gridMat(cross,:);
    systemParams.gridMat = systemParams.gridMat-RegionsAll(1)+1;
    systemParams.TLlength = baseData.systemParams.TLlength(cross);
    systemParams.TL_AC_LowLimits = baseData.systemParams.TL_AC_LowLimits(cross);
    systemParams.TL_AC_UpLimits = baseData.systemParams.TL_AC_UpLimits(cross);
    systemParams.TL_DC_LowLimits = baseData.systemParams.TL_DC_LowLimits(cross);
    systemParams.TL_DC_UpLimits = baseData.systemParams.TL_DC_UpLimits(cross);
end
resBase = load([setup.rootDir filesep 'projects' filesep pName '_' num2str(costYear) filesep 'output' filesep 'matlab' filesep 'results_' num2str(1) '.mat']);

FN = fieldnames(resBase.results);

for i = 1:length(FN)

    resultsAll.(FN{i}) = zeros(size(resBase.results.(FN{i}),1),length(RegionsAll));

end

resultsAll.LINEneg_AC = zeros(8760,size(systemParams.gridMat,1));
resultsAll.LINEneg_DC = zeros(8760,size(systemParams.gridMat,1));
resultsAll.LINEpos_AC = zeros(8760,size(systemParams.gridMat,1));
resultsAll.LINEpos_DC = zeros(8760,size(systemParams.gridMat,1));
resultsAll.OPT_SIZE_TRTL = zeros(1,size(systemParams.gridMat,1));


for i = 1:length(Countries)
    Regions = Countries(i,1);
    if iscell(setup.Countries)
        Regions = setup.Countries{i};
    else
        Regions = Countries(i,1);
    end
    r = Regions(1);
    resBase = load([setup.rootDir filesep 'projects' filesep pName '_' num2str(costYear) filesep 'output' filesep 'matlab' filesep 'results_' num2str(r) '.mat']);
    SysBase = load([setup.rootDir filesep 'projects' filesep pName '_' num2str(costYear) filesep 'input-data' filesep 'simulation-input_' pName '_' num2str(costYear) '_'  num2str(r) '.mat']);
    %% SYSTEM PARAMETERS
    systemParams.area(Regions) = SysBase.systemParams.area;
    systemParams.population(Regions,:) = SysBase.systemParams.population;
    systemParams.IndexNodes(Regions) = SysBase.systemParams.IndexNodes;
    systemParams.IndexNumNodes(Regions) = SysBase.systemParams.IndexNumNodes;
    systemParams.Coords(:,Regions) = SysBase.systemParams.Coords;
    systemParams.Urbanisation(:,Regions) = SysBase.systemParams.Urbanisation;
    systemParams.SizeLimits(:,:,Regions) = SysBase.systemParams.SizeLimits;
    systemParams.Instalations(Regions,:,:) = SysBase.systemParams.Instalations;
    systemParams.ValueLoad(:,:,:,Regions) = SysBase.systemParams.ValueLoad;
    systemParams.ValueLoadTag(:,Regions) = SysBase.systemParams.ValueLoadTag;
    systemParams.ValueResource(:,:,:,Regions) = SysBase.systemParams.ValueResource;
    systemParams.ValueHydroDam (:,:,:,Regions)= SysBase.systemParams.ValueHydroDam;
    systemParams.Capex_reg(:,Regions,:) = SysBase.systemParams.Capex_reg;
    systemParams.Opex_fix_reg(:,Regions,:) = SysBase.systemParams.Opex_fix_reg;
    systemParams.Opex_var_reg(:,Regions,:) = SysBase.systemParams.Opex_var_reg;
    systemParams.ValueResourceTotal(:,Regions) = SysBase.systemParams.ValueResourceTotal;
    systemParams.DesalinationParams.regionPumpLong(Regions,:) = SysBase.systemParams.DesalinationParams.regionPumpLong;
    systemParams.DesalinationParams.regionPumpUp(Regions,:) = SysBase.systemParams.DesalinationParams.regionPumpUp;

    systemParams.Heat.shareOfDistrHeat(Regions,:)= SysBase.systemParams.Heat.shareOfDistrHeat;
    systemParams.Heat.shareOfHighHeatInd(Regions)= SysBase.systemParams.Heat.shareOfHighHeatInd;
    systemParams.Heat.shareOfLowHeatInd(Regions)= SysBase.systemParams.Heat.shareOfLowHeatInd;
    try
        systemParams.Heat.shareIndHeatElPros(Regions)= SysBase.systemParams.Heat.shareIndHeatElPros;

    catch
        systemParams.Heat.shareIndHeatElPros(Regions)= 0;
    end
    %if setup.Mobility
    systemParams.Mobility.Demand(:,Regions,:) = SysBase.systemParams.Mobility.Demand;

    systemParams.Mobility.LDV_Pass(Regions,:) = SysBase.systemParams.Mobility.LDV_Pass;
    systemParams.Mobility.BUS_Pass(Regions,:) = SysBase.systemParams.Mobility.BUS_Pass;
    systemParams.Mobility.W23_Pass(Regions,:) = SysBase.systemParams.Mobility.W23_Pass;
    systemParams.Mobility.LDV_km(Regions,:) = SysBase.systemParams.Mobility.LDV_km;
    systemParams.Mobility.BUS_km(Regions,:) = SysBase.systemParams.Mobility.BUS_km;
    systemParams.Mobility.W23_km(Regions,:) = SysBase.systemParams.Mobility.W23_km;
    systemParams.Mobility.LDV_PHEV_el(Regions,:) = SysBase.systemParams.Mobility.LDV_PHEV_el;
    systemParams.Mobility.BUS_PHEV_el(Regions,:) = SysBase.systemParams.Mobility.BUS_PHEV_el;
    systemParams.Mobility.W23_PHEV_el(Regions,:) = SysBase.systemParams.Mobility.W23_PHEV_el;

    systemParams.Mobility.MDV_Tonne(Regions,:) = SysBase.systemParams.Mobility.MDV_Tonne;
    systemParams.Mobility.HDV_Tonne(Regions,:) = SysBase.systemParams.Mobility.HDV_Tonne;
    systemParams.Mobility.MDV_km(Regions,:) = SysBase.systemParams.Mobility.MDV_km;
    systemParams.Mobility.HDV_km(Regions,:) = SysBase.systemParams.Mobility.HDV_km;
    systemParams.Mobility.MDV_PHEV_el(Regions,:) = SysBase.systemParams.Mobility.MDV_PHEV_el;
    systemParams.Mobility.HDV_PHEV_el(Regions,:) = SysBase.systemParams.Mobility.HDV_PHEV_el;

    systemParams.Mobility.Shares(Regions,:,:) = SysBase.systemParams.Mobility.Shares;

    systemParams.Mobility.Cons_Primary_S(Regions,:,:) = SysBase.systemParams.Mobility.Cons_Primary_S;
    systemParams.Mobility.Cons_Secondary_S(Regions,:,:) = SysBase.systemParams.Mobility.Cons_Secondary_S;
    systemParams.Mobility.Cons_Names_Secondary = SysBase.systemParams.Mobility.Cons_Names_Secondary;
    %end
    systemParams.EffCO2Scr_El(Regions,:) = SysBase.systemParams.EffCO2Scr_El;
    systemParams.EffCO2Scr_He(Regions,:) = SysBase.systemParams.EffCO2Scr_He;
    systemParams.EffMETH_CO(Regions,:) = SysBase.systemParams.EffMETH_CO;

    systemParams.EffFT_CO(Regions,:) = SysBase.systemParams.EffFT_CO;
    systemParams.EffFT_He(Regions,:) = SysBase.systemParams.EffFT_He;
    systemParams.EffLNG_El(Regions,:) = SysBase.systemParams.EffLNG_El;
    systemParams.EffLH2_El(Regions,:) = SysBase.systemParams.EffLH2_El;

    %if setup.New
    systemParams.EffMeO_CO(Regions,:) = SysBase.systemParams.EffMeO_CO;
    systemParams.EffMeO_El(Regions,:) = SysBase.systemParams.EffMeO_El;
    systemParams.EffMeO_He_out(Regions,:) = SysBase.systemParams.EffMeO_He_out;
    systemParams.EffDME_CO(Regions,:) = SysBase.systemParams.EffDME_CO;
    systemParams.EffDME_El(Regions,:) = SysBase.systemParams.EffDME_El;
    systemParams.EffDME_He_out(Regions,:) = SysBase.systemParams.EffDME_He_out;
    systemParams.EffNH3_El(Regions,:) = SysBase.systemParams.EffNH3_El;
    systemParams.EffNH3_He_out(Regions,:) = SysBase.systemParams.EffNH3_He_out;
    systemParams.EffHyS_El(Regions,:) = SysBase.systemParams.EffHyS_El;
    systemParams.EffMET_El(Regions,:) = SysBase.systemParams.EffMET_El;
    systemParams.EffMET_He_out(Regions,:) = SysBase.systemParams.EffMET_He_out;
    systemParams.EffWEL_He_out(Regions,:) = SysBase.systemParams.EffWEL_He_out;
    systemParams.EffICB_Lime(Regions,:) = SysBase.systemParams.EffICB_Lime;
    systemParams.EffICB_El(Regions,:) = SysBase.systemParams.EffICB_El;
    systemParams.EffICB_He(Regions,:) = SysBase.systemParams.EffICB_He;
    systemParams.EffICI_Lime(Regions,:) = SysBase.systemParams.EffICI_Lime;
    systemParams.EffICI_El(Regions,:) = SysBase.systemParams.EffICI_El;
    systemParams.EffICI_He(Regions,:) = SysBase.systemParams.EffICI_He;
    systemParams.EffISB_Ore(Regions,:) = SysBase.systemParams.EffISB_Ore;
    systemParams.EffISB_Scr(Regions,:) = SysBase.systemParams.EffISB_Scr;
    systemParams.EffISB_El(Regions,:) = SysBase.systemParams.EffISB_El;
    systemParams.EffISB_He(Regions,:) = SysBase.systemParams.EffISB_He;
    systemParams.EffISB_HC(Regions,:) = SysBase.systemParams.EffISB_HC;
    systemParams.EffISH_Ore(Regions,:) = SysBase.systemParams.EffISH_Ore;
    systemParams.EffISH_Scr(Regions,:) = SysBase.systemParams.EffISH_Scr;
    systemParams.EffISH_El(Regions,:) = SysBase.systemParams.EffISH_El;
    systemParams.EffISH_He(Regions,:) = SysBase.systemParams.EffISH_He;
    systemParams.EffISH_CC(Regions,:) = SysBase.systemParams.EffISH_CC;
    systemParams.EffISH_Hy(Regions,:) = SysBase.systemParams.EffISH_Hy;
    systemParams.EffISR_Scr(Regions,:) = SysBase.systemParams.EffISR_Scr;
    systemParams.EffISR_El(Regions,:) = SysBase.systemParams.EffISR_El;
    systemParams.EffISR_He(Regions,:) = SysBase.systemParams.EffISR_He;
    systemParams.EffISR_CC(Regions,:) = SysBase.systemParams.EffISR_CC;
    systemParams.EffISE_Ore(Regions,:) = SysBase.systemParams.EffISE_Ore;
    systemParams.EffISE_Scr(Regions,:) = SysBase.systemParams.EffISE_Scr;
    systemParams.EffISE_El(Regions,:) = SysBase.systemParams.EffISE_El;
    systemParams.EffISE_He(Regions,:) = SysBase.systemParams.EffISE_He;
    systemParams.EffISE_CC(Regions,:) = SysBase.systemParams.EffISE_CC;
    systemParams.EffIAA_Ba(Regions,:) = SysBase.systemParams.EffIAA_Ba;
    systemParams.EffIAA_Lime(Regions,:) = SysBase.systemParams.EffIAA_Lime;
    systemParams.EffIAA_So(Regions,:) = SysBase.systemParams.EffIAA_So;
    systemParams.EffIAA_He(Regions,:) = SysBase.systemParams.EffIAA_He;
    systemParams.EffIAA_HeOut(Regions,:) = SysBase.systemParams.EffIAA_HeOut;
    systemParams.EffIAM_Aa(Regions,:) = SysBase.systemParams.EffIAM_Aa;
    systemParams.EffIAM_El(Regions,:) = SysBase.systemParams.EffIAM_El;
    systemParams.EffIAM_HeOut(Regions,:) = SysBase.systemParams.EffIAM_HeOut;
    systemParams.EffIAR_Al(Regions,:) = SysBase.systemParams.EffIAR_Al;
    systemParams.EffIAR_El(Regions,:) = SysBase.systemParams.EffIAR_El;
    systemParams.EffIAR_He(Regions,:) = SysBase.systemParams.EffIAR_He;
    systemParams.EffIPP_Wo(Regions,:) = SysBase.systemParams.EffIPP_Wo;
    systemParams.EffIPP_El(Regions,:) = SysBase.systemParams.EffIPP_El;
    systemParams.EffIPP_He(Regions,:) = SysBase.systemParams.EffIPP_He;
    %end

    systemParams.EffHCS_CO(Regions,:) = SysBase.systemParams.EffHCS_CO;
    systemParams.EffGCS_CO(Regions,:) = SysBase.systemParams.EffGCS_CO;

    systemParams.AC_losses(Regions,:) = SysBase.systemParams.AC_losses;

    systemParams.SectorsCons(Regions,:) = SysBase.systemParams.SectorsCons;
    systemParams.ElCostRES(Regions,:) = SysBase.systemParams.ElCostRES;
    systemParams.ElCostCOM(Regions,:) = SysBase.systemParams.ElCostCOM;
    systemParams.ElCostIND(Regions,:) = SysBase.systemParams.ElCostIND;
    systemParams.name = name;

    Regions
    %% RESULTS MAIN
    for ii = 1:length(FN)
        if ~ismember(FN{ii},{'OPT_SIZE_TRTL','OPT_SIZE_THAO','LINEneg_AC','LINEpos_AC','LINEneg_DC','LINEpos_DC','SolvParams'})
            FN{ii};
            resultsAll.(FN{ii})(:,Regions) = resBase.results.(FN{ii});
        end
    end

    %% RESULTS GRID intra

    cross = (ismember(baseData.systemParams.gridMat(:,1),Regions)&ismember(baseData.systemParams.gridMat(:,2),Regions));
    if sum(cross)
        Regions
        resultsAll.OPT_SIZE_TRTL(cross) = resBase.results.OPT_SIZE_TRTL;
        resultsAll.OPT_SIZE_THAO(cross) = resBase.results.OPT_SIZE_THAO;

        resultsAll.LINEneg_AC(:,cross) = resBase.results.LINEneg_AC;
        resultsAll.LINEneg_DC(:,cross) = resBase.results.LINEneg_DC;
        resultsAll.LINEpos_AC(:,cross) = resBase.results.LINEpos_AC;
        resultsAll.LINEpos_DC(:,cross) = resBase.results.LINEpos_DC;
    end




end


%% Grid el prep

effCS = systemParams.TransmissionEfficiencies(ismember(systemParams.IndexIDTr,'TRCS'));
%HVDC
for nodeNo=1:max(max(systemParams.gridMat))
    neg  = find(systemParams.gridMat(:,1) == nodeNo);
    pos  = find(systemParams.gridMat(:,2) == nodeNo);

    transPower_DC=zeros(8760,1);

    for k=1:length(neg)
        effTL = 1-(1-systemParams.TransmissionEfficiencies(ismember(systemParams.IndexIDTr,'TRTL')))*systemParams.TLlength(neg(k))/1000;
        transPower_DC = transPower_DC + resultsAll.LINEneg_DC(:,neg(k))*effCS*effTL - resultsAll.LINEpos_DC(:,neg(k));
    end


    for k=1:length(pos)
        effTL = 1-(1-systemParams.TransmissionEfficiencies(ismember(systemParams.IndexIDTr,'TRTL')))*systemParams.TLlength(pos(k))/1000;
        transPower_DC = transPower_DC + resultsAll.LINEpos_DC(:,pos(k))*effCS*effTL - resultsAll.LINEneg_DC(:,pos(k));
    end

    resultsAll.GRID_DC(:,nodeNo) = transPower_DC;

end

%HVAC
for nodeNo=1:max(max(systemParams.gridMat))
    neg  = find(systemParams.gridMat(:,1) == nodeNo);
    pos  = find(systemParams.gridMat(:,2) == nodeNo);

    transPower_AC=zeros(8760,1);

    for k=1:length(neg)
        effTL = 1-(1-systemParams.TransmissionEfficiencies(ismember(systemParams.IndexIDTr,'THAO')))*systemParams.TLlength(neg(k))/1000;
        transPower_AC = transPower_AC + resultsAll.LINEneg_AC(:,neg(k))*effTL - resultsAll.LINEpos_AC(:,neg(k));
    end


    for k=1:length(pos)
        effTL = 1-(1-systemParams.TransmissionEfficiencies(ismember(systemParams.IndexIDTr,'THAO')))*systemParams.TLlength(pos(k))/1000;
        transPower_AC = transPower_AC + resultsAll.LINEpos_AC(:,pos(k))*effTL - resultsAll.LINEneg_AC(:,pos(k));
    end

    resultsAll.GRID_AC(:,nodeNo) = transPower_AC;

end

qq = 0

results = resultsAll;

save([setup.rootDir filesep 'projects' filesep pName '_' num2str(costYear) filesep 'output' filesep 'matlab' filesep 'results_' name '.mat'],'results','-v7.3');
save([setup.rootDir filesep 'projects' filesep pName '_' num2str(costYear) filesep 'output' filesep 'matlab' filesep 'sysParams_' pName '_' num2str(costYear) '_' name '.mat'],'systemParams');
save([setup.rootDir filesep 'projects' filesep pName '_' num2str(costYear) filesep 'input-data' filesep 'simulation-input_' pName '_' num2str(costYear) '_' name '.mat'],'systemParams','activeElements', 'modFiles','datFiles');















