function writeMisc(fileName,setup,struct,activeElements,transCapa,fossilLimit,fossilGasLimit,storePeriod,costYear)
%
% FUNCTION writeMisc(fileName, setup, struct, activeElements, transCapa, fossilLimit, fossilGasLimit, storePeriod, costYear)
%
% Writes miscellaneous data for the scenario to a file.
%
% INPUT:
%            fileName: Desired filename without extension
%            setup: Structure that contains all necessary settings and data for processing
%            struct: Parameter structure obtained from Excel sheets
%            activeElements: Structure with active elements in the system
%            transCapa: Transmission capacity data
%            fossilLimit: Limit on total fossil fuel use
%            fossilGasLimit: Limit on fossil gas use
%            storePeriod: Storage period
%            costYear: Year to which costs are adjusted
%
% OUTPUT:
%            -
%
%Dmitrii Bogdanov
%last change 24.07.2025


% input data in kW, kWh, kg for CO2 and tonnes of product for industry, m3
% per day for desalination
try setup.units;
catch
    setup.units = 'kW';
end
switch setup.units
    case 'kW'
        m_factor = 1;
    case 'MW'
        m_factor = 0.001;
    case 'GW'
        m_factor = 0.000001;
end

WaterConsumption = 0.000466; %water consumption m3 per 1kwh of methanation imput

% feed-in time series for feed-in resources


WriteParamFileBegin(fileName);

Vpump_el_demand = struct.Opex_var(find(ismember(struct.IndexID,'WVPU')),(struct.IndexYears==costYear));
Hpump_el_demand = struct.Opex_var(find(ismember(struct.IndexID,'WHPU')),(struct.IndexYears==costYear));


WriteParams(fileName,'Vdist',round(struct.DesalinationParams.regionPumpUp(:,struct.IndexYears==costYear)/1000,3),[1:length(struct.IndexNodes)]')
WriteParams(fileName,'Hdist',round(struct.DesalinationParams.regionPumpLong(:,struct.IndexYears==costYear)/1000,3),[1:length(struct.IndexNodes)]')

WriteParams(fileName,'RO_El_cons',round(struct.Efficiency(find(ismember(struct.IndexID,'WROD'))).*struct.Opex_var_reg(find(ismember(struct.IndexID,'WROD')),:,struct.IndexYears==costYear)',4),[1:length(struct.IndexNodes)]')
WriteParams(fileName,'MSS_El_cons',round(struct.Efficiency(find(ismember(struct.IndexID,'WMSS'))).*struct.Opex_var_reg(find(ismember(struct.IndexID,'WMSS')),:,struct.IndexYears==costYear)',4),[1:length(struct.IndexNodes)]')
WriteParams(fileName,'MSC_El_cons',round(struct.Efficiency(find(ismember(struct.IndexID,'WMSC'))).*struct.Opex_var_reg(find(ismember(struct.IndexID,'WMSC')),:,struct.IndexYears==costYear)',4),[1:length(struct.IndexNodes)]')
WriteParams(fileName,'MDS_El_cons',round(struct.Efficiency(find(ismember(struct.IndexID,'WMDS'))).*struct.Opex_var_reg(find(ismember(struct.IndexID,'WMDS')),:,struct.IndexYears==costYear)',4),[1:length(struct.IndexNodes)]')
WriteParams(fileName,'MDC_El_cons',round(struct.Efficiency(find(ismember(struct.IndexID,'WMDC'))).*struct.Opex_var_reg(find(ismember(struct.IndexID,'WMDC')),:,struct.IndexYears==costYear)',4),[1:length(struct.IndexNodes)]')

WriteParams(fileName,'MSS_Heat_cons',round(struct.Opex_var(find(ismember(struct.IndexID,'WMSS')),struct.IndexYears==costYear).*struct.Opex_var_reg(find(ismember(struct.IndexID,'WMSS')),:,struct.IndexYears==costYear)',4),[1:length(struct.IndexNodes)]')
WriteParams(fileName,'MSC_Gas_cons',round(struct.Opex_var(find(ismember(struct.IndexID,'WMSC')),struct.IndexYears==costYear).*struct.Opex_var_reg(find(ismember(struct.IndexID,'WMSC')),:,struct.IndexYears==costYear)',4),[1:length(struct.IndexNodes)]')
WriteParams(fileName,'MDS_Heat_cons',round(struct.Opex_var(find(ismember(struct.IndexID,'WMDS')),struct.IndexYears==costYear).*struct.Opex_var_reg(find(ismember(struct.IndexID,'WMDS')),:,struct.IndexYears==costYear)',4),[1:length(struct.IndexNodes)]')
WriteParams(fileName,'MDC_Gas_cons',round(struct.Opex_var(find(ismember(struct.IndexID,'WMDC')),struct.IndexYears==costYear).*struct.Opex_var_reg(find(ismember(struct.IndexID,'WMDC')),:,struct.IndexYears==costYear)',4),[1:length(struct.IndexNodes)]')

WriteParams(fileName,'Hdemand',round(Hpump_el_demand*1000,3))
WriteParams(fileName,'Vdemand',round(Vpump_el_demand*1000,3))

WriteParams(fileName,'water_coef',WaterConsumption)

WriteParams(fileName,'transCapa',transCapa)
WriteParams(fileName,'fossilLimit',fossilLimit)
WriteParams(fileName,'fossilGasLimit',fossilGasLimit',[1:length(struct.IndexNodes)]')
WriteParams(fileName,'storePeriod_HDAM',storePeriod)

% WriteParams(fileName,'AvF_Oil',round(struct.AvF_Oil,3)',[1:length(struct.IndexNodes)]')
% WriteParams(fileName,'AvF_Coal',round(struct.AvF_Coal,3)',[1:length(struct.IndexNodes)]')
% WriteParams(fileName,'AvF_Gas',round(struct.AvF_Gas,3)',[1:length(struct.IndexNodes)]')
% WriteParams(fileName,'AvF_Nuc',round(struct.AvF_Nuc,3)',[1:length(struct.IndexNodes)]')

WriteParams(fileName,'feedInEfficiency',round(struct.FeedInEfficiencies(ismember(struct.IndexIDE,activeElements.labels.elFeedIn))',3),activeElements.labels.elFeedIn)

WriteParams(fileName,'urbanPop',round(struct.Urbanisation',3),[1:length(struct.IndexNodes)]')

WriteParams(fileName,'shareOfREN',round(struct.REN_Share_Lim(:,struct.IndexYears==costYear),3),[1:length(struct.IndexNodes)]')
WriteParams(fileName,'shareOfREN_Area',round(struct.REN_Share_Lim_Area(struct.IndexYears==costYear),3))

struct.TL_DC_LowLimits(struct.TL_DC_LowLimits>struct.TL_DC_UpLimits) = struct.TL_DC_UpLimits(struct.TL_DC_LowLimits>struct.TL_DC_UpLimits);
WriteParams(fileName,'TransmissionLines_DC_LowerLim',round(struct.TL_DC_LowLimits*m_factor,3),[1:length(struct.TL_DC_LowLimits)])

WriteParams(fileName,'TransmissionLines_DC_UpperLim',round(struct.TL_DC_UpLimits*m_factor,3),[1:length(struct.TL_DC_UpLimits)])

struct.TL_AC_LowLimits(struct.TL_AC_LowLimits>struct.TL_AC_UpLimits) = struct.TL_AC_UpLimits(struct.TL_AC_LowLimits>struct.TL_AC_UpLimits);
WriteParams(fileName,'TransmissionLines_AC_LowerLim',round(struct.TL_AC_LowLimits*m_factor,3),[1:length(struct.TL_AC_LowLimits)])

WriteParams(fileName,'TransmissionLines_AC_UpperLim',round(struct.TL_AC_UpLimits*m_factor,3),[1:length(struct.TL_AC_UpLimits)])

WriteParams(fileName,'GridMax',struct.GridMax);
% CO2 pricing
WriteParams(fileName,'fossilCO2Cost',round(struct.fossilCO2Cost(struct.IndexYears==costYear),3));
%DAC
WriteParams(fileName,'DAC_El_cons',round(struct.EffCO2Scr_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'DAC_He_cons',round(struct.EffCO2Scr_He(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
%SNG
WriteParams(fileName,'SNG_CO2_cons',round(struct.EffMETH_CO(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');


% FT
try
    WriteParams(fileName,'FT_El_cons',round(struct.EffFT_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
end
WriteParams(fileName,'FT_CO2_cons',round(struct.EffFT_CO(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'FT_He_prod',round(struct.EffFT_He(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
% Bio FT
WriteParams(fileName,'BioFT_El_cons',round(struct.EffFTB_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
% Methanol FT
WriteParams(fileName,'MethanolFT_El_cons',round(struct.EffFTM_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

% LNG-LH2
WriteParams(fileName,'LNG_El_cons',round(struct.EffLNG_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'LH2_El_cons',round(struct.EffLH2_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');



WriteParams(fileName,'MeO_CO2_cons',round(struct.EffMeO_CO(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'MeO_El_cons',round(struct.EffMeO_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'MeO_He_prod',round(struct.EffMeO_He_out(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'DME_CO2_cons',round(struct.EffDME_CO(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'DME_El_cons',round(struct.EffDME_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'DME_He_prod',round(struct.EffDME_He_out(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'NH3_El_cons',round(struct.EffNH3_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'NH3_He_prod',round(struct.EffNH3_He_out(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'HyS_El_cons',round(struct.EffHyS_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'SNG_El_cons',round(struct.EffMET_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'SNG_He_prod',round(struct.EffMET_He_out(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'WEL_He_prod',round(struct.EffWEL_He_out(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'ICB_Lime_cons',round(struct.EffICB_Lime(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ICB_El_cons',round(struct.EffICB_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ICB_He_cons',round(struct.EffICB_He(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'ICI_Lime_cons',round(struct.EffICI_Lime(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ICI_El_cons',round(struct.EffICI_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ICI_He_cons',round(struct.EffICI_He(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'ISB_Ore_cons',round(struct.EffISB_Ore(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ISB_El_cons',round(struct.EffISB_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ISB_He_cons',round(struct.EffISB_He(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ISB_HC_cons',round(struct.EffISB_HC(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'ISH_Ore_cons',round(struct.EffISH_Ore(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ISH_El_cons',round(struct.EffISH_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ISH_He_cons',round(struct.EffISH_He(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ISH_Hy_cons',round(struct.EffISH_Hy(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ISH_CC_cons',round(struct.EffISH_CC(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'ISR_Scrap_cons',round(struct.EffISR_Scr(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ISR_El_cons',round(struct.EffISR_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ISR_He_cons',round(struct.EffISR_He(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ISR_CC_cons',round(struct.EffISR_CC(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'ISE_Ore_cons',round(struct.EffISE_Ore(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ISE_El_cons',round(struct.EffISE_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ISE_He_cons',round(struct.EffISE_He(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ISE_CC_cons',round(struct.EffISE_CC(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'IAA_Ba_cons',round(struct.EffIAA_Ba(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'IAA_So_cons',round(struct.EffIAA_So(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'IAA_He_cons',round(struct.EffIAA_He(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'IAA_He_prod',round(struct.EffIAA_HeOut(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'IAM_Aa_cons',round(struct.EffIAM_Aa(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'IAM_El_cons',round(struct.EffIAM_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'IAM_He_prod',round(struct.EffIAM_HeOut(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'IAR_Al_cons',round(struct.EffIAR_Al(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'IAR_El_cons',round(struct.EffIAR_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'IAR_He_cons',round(struct.EffIAR_He(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'IPP_Wo_cons',round(struct.EffIPP_Wo(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'IPP_El_cons',round(struct.EffIPP_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'IPP_He_cons',round(struct.EffIPP_He(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

%CCS

WriteParams(fileName,'shareOfSNG',round(struct.shareOfSNG(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'shareOfDistrHeat',setup.Heat.shareOfDistrHeat(:,struct.IndexYears==costYear),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'shareOfLowHeatInd',setup.Heat.shareOfLowHeatInd',[1:length(struct.IndexNodes)]');
WriteParams(fileName,'shareOfHighHeatInd',setup.Heat.shareOfHighHeatInd',[1:length(struct.IndexNodes)]');
WriteParams(fileName,'DistrHeatEff',round(setup.Heat.DistrHeatEff(:,struct.IndexYears==costYear),3),[1:length(struct.IndexNodes)]');


WriteParams(fileName,'shareOfProsWOO',setup.SC.Biomass_BiomassShare',[1:length(struct.IndexNodes)]');
WriteParams(fileName,'shareOfProsBGA',setup.SC.Biomass_BiogasShare',[1:length(struct.IndexNodes)]');

WriteParams(fileName,'DumpChargeShare_LDV',setup.DumpChargeShare_LDV(:,struct.IndexYears==costYear),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'V2GShare_LDV',setup.V2GShare_LDV(:,struct.IndexYears==costYear),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'DumpChargeShare_23W',setup.DumpChargeShare_23W(:,struct.IndexYears==costYear),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'V2GShare_23W',setup.V2GShare_23W(:,struct.IndexYears==costYear),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'DumpChargeShare_BUS',setup.DumpChargeShare_BUS(:,struct.IndexYears==costYear),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'V2GShare_BUS',setup.V2GShare_BUS(:,struct.IndexYears==costYear),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'DumpChargeShare_MDV',setup.DumpChargeShare_MDV(:,struct.IndexYears==costYear),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'V2GShare_MDV',setup.V2GShare_MDV(:,struct.IndexYears==costYear),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'DumpChargeShare_HDV',setup.DumpChargeShare_HDV(:,struct.IndexYears==costYear),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'V2GShare_HDV',setup.V2GShare_HDV(:,struct.IndexYears==costYear),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'ProsumersPV_EL',round(sum(struct.SCData.RPVO_EL,1)*m_factor,3)',[1:length(struct.IndexNodes)]');
WriteParams(fileName,'ProsumersPV_cap',round(sum(struct.SCData.OPT_SIZE_RPVO,1)*m_factor,3)',[1:length(struct.IndexNodes)]');

try setup.fixTech.Flag;
catch
    setup.fixTech.Flag = 0;
end

if (setup.fixTech.Flag)*(costYear>setup.startYear)

    for tt = 1:length(setup.fixTech.Tech)
        temp = (struct.prevStepCap.(setup.fixTech.Tech{tt}))';
        temp = floor(temp);
        WriteParams(fileName,['minimum_' setup.fixTech.Tech{tt}],floor(temp*m_factor),[1:length(struct.IndexNodes)]');
        if strcmp(setup.fixTech.Tech{tt},'RWIN')
            temp = struct.prevStepCap.RWIO';
            temp = floor(temp);
            WriteParams(fileName,'minimum_RWIO',floor(temp*m_factor),[1:length(struct.IndexNodes)]');
        end
    end
end

try
    WriteParams(fileName,'importH2Cost',round(setup.importH2Cost(:,struct.IndexYears==costYear),3),[1:length(struct.IndexNodes)]');
    WriteParams(fileName,'importLNGCost',round(setup.importLNGCost(:,struct.IndexYears==costYear),3),[1:length(struct.IndexNodes)]');
    WriteParams(fileName,'importFTCost',round(setup.importFTCost(:,struct.IndexYears==costYear),3),[1:length(struct.IndexNodes)]');
end
try
    WriteParams(fileName,'importMeOHCost',round(setup.importMeOHCost(:,struct.IndexYears==costYear),3),[1:length(struct.IndexNodes)]');
    WriteParams(fileName,'importNH3Cost',round(setup.importNH3Cost(:,struct.IndexYears==costYear),3),[1:length(struct.IndexNodes)]');
end




WriteParams(fileName,'THCS_CCeff',round(struct.EffHCS_CO(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TGCS_CCeff',round(struct.EffGCS_CO(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TSMC_CCeff',round(struct.EffSMC_CO(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'TBCS_CCeff',round(struct.EffBCS_CO(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TCBC_CCeff',round(struct.EffCBC_CO(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TCWC_CCeff',round(struct.EffCWC_CO(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TPSC_CCeff',round(struct.EffPSC_CO(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TPSP_CCeff',round(struct.EffPSP_CO(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');


WriteParams(fileName,'TPSC_El_cons',round(struct.EffPSC_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TPSC_He_cons',round(struct.EffPSC_He(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TPSP_El_cons',round(struct.EffPSP_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TPSP_He_cons',round(struct.EffPSP_He(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'TRSL_input',round(struct.EffRLS_input(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRSL_El_cons',round(struct.EffRLS_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRSS_input',round(struct.EffRSS_input(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRSS_El_cons',round(struct.EffRSS_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'TRMI_input',round(struct.EffRMI_input(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRMI_Wa_cons',round(struct.EffRMI_Wa(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRMI_El_cons',round(struct.EffRMI_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRMO_input',round(struct.EffRMO_input(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRMO_Wa_cons',round(struct.EffRMO_Wa(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRMO_El_cons',round(struct.EffRMO_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'TRME_input',round(struct.EffRME_input(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRME_El_cons',round(struct.EffRME_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRME_He_cons',round(struct.EffRME_He(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'TRSi_input',round(struct.EffRSi_input(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRSi_Hy_cons',round(struct.EffRSi_Hy(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRSi_El_cons',round(struct.EffRSi_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRSi_He_cons',round(struct.EffRSi_He(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');


WriteParams(fileName,'TRAD_Wa_cons',round(struct.EffRAD_Wa(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRAD_El_cons',round(struct.EffRAD_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'TREW_El_cons',round(struct.EffREW_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'TRBC_Bio_cons',round(struct.EffRBC_Bio(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRBC_Wa_cons',round(struct.EffRBC_Wa(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRBC_El_cons',round(struct.EffRBC_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'TRGE_El_cons',round(struct.EffRGE_El(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRGE_He_cons',round(struct.EffRGE_He(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'TRAP_eff',round(setup.TRAP_Eff(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'TRAD_eff',round(setup.TRAD_Eff(:,struct.IndexYears==costYear),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'TandD_Loss',round(1./(1-struct.AC_losses(:,struct.IndexYears==costYear)/100),4),[1:length(struct.IndexNodes)]');

WriteParams(fileName,'TICH_CO' ,repmat(round(setup.TICH_CO(:,struct.IndexYears==costYear)/setup.endHour,4),1,setup.endHour)',[1:setup.endHour],[1:length(struct.IndexNodes)]);

WriteParamFileEnd(fileName);