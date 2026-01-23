function [results] = prepareResults(setup,namesAndInd,values,totalCost)
%
% FUNCTION prepareResults(setup, namesAndInd, values, totalCost)
%
% Prepares and structures result data from the model run.
%
%
% INPUT:
%            setup:          Structure that contains all necessary settings and data for processing.
%            namesAndInd:    Cell array with names and indices for each result entry.
%            values:         Corresponding values for the result entries.
%            totalCost:      Total cost value to be included in the result structure.
%
% OUTPUT:
%            results:        Structure containing organized results and cost information.
%
%Dmitrii Bogdanov
%last change 23.07.2025


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

for varNo=1:length(namesAndInd)
    varName = regexp(namesAndInd{varNo},'\<[A-z0-9]+(?=\()','match');
    varIndices = regexp(namesAndInd{varNo},'[\w+,]+\w+(?=\))','match');
    if isempty(varIndices); varIndices = regexp(namesAndInd{varNo},'\w+(?=\))','match');end
    varIndicesSplit = regexp(varIndices,',','split');
    try
        varIndicesSplit = varIndicesSplit{1};
    catch
        display(num2str(varNo));
    end

    for k=1:length(varIndicesSplit)
        check(k) = isempty(regexp(varIndicesSplit{k},'[a-zA-Z]+','match'));
    end

    switch varName{1}

        %% Base system
        % power


        case 'gen_electTransf'
            eval(['results.' varIndicesSplit{3} '_EL(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'gen_electFeedIn'
            eval(['results.' varIndicesSplit{3} '_EL(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'gen_Hydro'
            eval(['results.' varIndicesSplit{3} '_EL(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);

        case 'instCapacity_electTransf'
            eval(['results.OPT_SIZE_' varIndicesSplit{2} '(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);
        case 'instCapacity_electFeedIn'
            eval(['results.OPT_SIZE_' varIndicesSplit{2} '(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);
        case 'instCapacity_Hydro'
            eval(['results.OPT_SIZE_' varIndicesSplit{2} '(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);
        case 'instCapacityDam'
            eval(['results.OPT_SIZE_' varIndicesSplit{2} '(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);

            % Storage
        case 'chargePower'
            eval(['results.EL_' varIndicesSplit{3} '(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'dischargePower'
            eval(['results.' varIndicesSplit{3} '_EL(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'storageState'
            eval(['results.SoC_' varIndicesSplit{3} '(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);

        case 'storageCapacity'
            eval(['results.OPT_SIZE_' varIndicesSplit{2} '(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);
        case 'storageInterface'
            eval(['results.OPT_SIZE_' varIndicesSplit{2} '(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);

        case 'storageDam'
            eval(['results.OPT_SIZE_ST_' varIndicesSplit{2} '(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);
        case 'storageStateDam'
            eval(['results.SoC_' varIndicesSplit{3} '(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'chargePowerDam'
            eval(['results.Ch_' varIndicesSplit{3} '(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);

        case 'dischargePowerDam'
            eval(['results.DCh_' varIndicesSplit{3} '(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);

        case 'geoheat_out'
            results.RGEO_HE(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'excessGeoHeat'
            results.GHE_EXCESS(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;


            % Grids
        case 'fromgrid'
            results.EL_GRID(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'instCapacitiyTransmLines_DC'
            results.OPT_SIZE_TRTL(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;
        case 'instCapacitiyTransmLines_AC'
            results.OPT_SIZE_THAO(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;
        case 'transPower_DC'
            results.GRID_DC(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'transPower_AC'
            results.GRID_AC(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'line_DC'
            switch varIndicesSplit{3}
                case 'pos'
                    results.LINEpos_DC(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
                case 'neg'
                    results.LINEneg_DC(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
            end
        case 'line_AC'
            switch varIndicesSplit{3}
                case 'pos'
                    results.LINEpos_AC(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
                case 'neg'
                    results.LINEneg_AC(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
            end

            % CO2 exhaust
        case 'fossilCO2_RNGA'
            results.RNGA_CO(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'fossilCO2_RPET'
            results.RPET_CO(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'fossilCO2_RHAR'
            results.RHAR_CO(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'fossilCO2'
            results.Fossil_CO(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

            %heat
        case 'gen_heatTransf'
            eval(['results.' varIndicesSplit{3} '_HE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'gen_heatTransf_EP'
            eval(['results.' varIndicesSplit{3} '_HE_EP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'gen_heatTransf_HP'
            eval(['results.' varIndicesSplit{3} '_HE_HP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);

        case 'gen_chpElTransf'
            eval(['results.' varIndicesSplit{3} '_EL(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'gen_chpHeTransf'
            eval(['results.' varIndicesSplit{3} '_HE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'gen_heatFeedIn'
            eval(['results.' varIndicesSplit{3} '_HE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'gen_heatFeedIn_EP'
            eval(['results.' varIndicesSplit{3} '_HE_EP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'gen_heatFeedIn_HP'
            eval(['results.' varIndicesSplit{3} '_HE_HP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
            %Cap heat
        case 'instCapacity_heatTransf'
            eval(['results.OPT_SIZE_' varIndicesSplit{2} '(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);
        case 'instCapacity_heatTransf_EP'
            eval(['results.OPT_SIZE_' varIndicesSplit{2} '_EP(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);
        case 'instCapacity_heatTransf_HP'
            eval(['results.OPT_SIZE_' varIndicesSplit{2} '_HP(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);
        case 'instCapacity_chpTransf'
            eval(['results.OPT_SIZE_' varIndicesSplit{2} '(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);
        case 'instCapacity_heatFeedIn'
            eval(['results.OPT_SIZE_' varIndicesSplit{2} '(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);
        case 'instCapacity_heatFeedIn_EP'
            eval(['results.OPT_SIZE_' varIndicesSplit{2} '_EP(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);
        case 'instCapacity_heatFeedIn_HP'
            eval(['results.OPT_SIZE_' varIndicesSplit{2} '_HP(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);

            % Fuels prod
        case 'gen_ptgTransf'
            eval(['results.' varIndicesSplit{3} '_GA(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'instCapacity_ptgTransf'
            eval(['results.OPT_SIZE_' varIndicesSplit{2} '(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);
            % Fuels inflow
        case 'inputFuel'
            eval(['results.' varIndicesSplit{3} '_FU(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'instCapacity_Resource'
            eval(['results.OPT_SIZE_' varIndicesSplit{2} '(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);
            % NG flows
        case 'importGasToOCGT'
            results.TOCG_GAS_FOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importGasToCCGT'
            results.TCCG_GAS_FOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importGasToTGCS'
            results.TGCS_GAS_FOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importGasToTICM'
            results.TICM_GAS_FOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importGasToCHP'
            results.TCNG_GAS_FOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importGasToHeat'
            results.TDNG_GAS_FOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importGasToIndHeat'
            results.LHIN_GAS_FOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importGasToIndivHeat'
            results.Indiv_GAS_FOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importGasToInd'
            results.LIGA_GAS_FOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importGasToDes'
            results.WDES_GAS_FOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'demandProsGasFos'
            results.Pros_GAS_FOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importGasToLNG'
            results.TLNG_GAS_FOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importGasToStRef'
            results.TSMR_GAS_FOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importGasToTSMR'
            results.TSMR_GAS_FOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importGasToTSMC'
            results.TSMC_GAS_FOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

            %SNG flows
        case 'SNGGasToOCGT'
            results.TOCG_GAS_REN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'SNGGasToCCGT'
            results.TCCG_GAS_REN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'SNGGasToTGCS'
            results.TGCS_GAS_REN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'SNGGasToTICM'
            results.TICM_GAS_REN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'SNGGasToCHP'
            results.TCNG_GAS_REN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'SNGGasToHeat'
            results.TDNG_GAS_REN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'SNGGasToIndHeat'
            results.LHIN_GAS_REN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'SNGGasToInd'
            results.LIGA_GAS_REN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'SNGGasToDes'
            results.WDES_GAS_REN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'demandProsGasRen'
            results.Pros_GAS_REN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'SNGGasToLNG'
            results.TLNG_GAS_REN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

            %Coal flows
        case 'CoalToPower'
            results.RHAR_THPP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'CoalToTHCS'
            results.RHAR_THCS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'CoalToCHP'
            results.RHAR_TCCO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'CoalToHeat'
            results.RHAR_TDCO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'CoalToIndHeat'
            results.RHAR_LHIN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

            %Oil flows
        case 'OilToPower'
            results.RPET_TICG(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'OilToTICM'
            results.RPET_TICM(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'OilToCCGT'
            results.RPET_TCCG(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'OilToOCGT'
            results.RPET_TOCG(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'OilToTGCS'
            results.RPET_TGCS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'OilToTCNG'
            results.RPET_TCNG(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'OilToCHP'
            results.RPET_TCOI(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'OilToHeat'
            results.RPET_TDOI(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'OilToIndHeat'
            results.RPET_LHIN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'FOSdiesel'
            results.RPET_DI(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'FOSkerosene'
            results.RPET_KE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
            %RWOO flows
        case 'BioToPower'
            results.RWOO_TBPP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BioToPowerCCS'
            results.RWOO_TBCS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BioToCHP'
            results.RWOO_TCBP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BioToCHPCCS'
            results.RWOO_TCBC(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BioToHeat'
            results.RWOO_TDBP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BioToIndHeat'
            results.RWOO_LHIN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BioToCook'
            results.RWOO_COOK(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BioToBDS'
            results.RWOO_RBDS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BioToBiocharCDR'
            results.RWOO_TRBC(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
            %RWWO flows
        case 'MBBioToPower'
            results.RWWO_TBPP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MBBioToPowerCCS'
            results.RWWO_TBCS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MBBioToCHP'
            results.RWWO_TCBP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MBBioToCHPCCS'
            results.RWWO_TCBC(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MBBioToHeat'
            results.RWWO_TDBP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MBBioToIndHeat'
            results.RWWO_LHIN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MBBioToCook'
            results.RWWO_COOK(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MBBioToBDS'
            results.RWWO_RBDS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
            %RBMW
        case 'WastesToCHPCCS'
            results.RBMW_TCWC(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'WastesToCHP'
            results.RBMW_TMSW(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

            %RBDS
        case 'BIOdiesel'
            results.RBDS_DI(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BIOkerosene'
            results.RBDS_KE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
            %RBGA
        case 'BGAToCHP'
            results.RBGA_TCHP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BGAToIndHeat'
            results.RBGA_LHIN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'BGAToIndCook'
            results.RBGA_COOK(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BGAToBME'
            results.RBGA_TBGU(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'demandProsBGA'
            results.RBGA_THBG(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;


        case 'unusedBiomass'
            results.UnusedBiomass(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;

            % e-fuels import

        case 'importH2'
            results.importH2(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importH2_regas'
            results.importH2_regas(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importH2_liq'
            results.importH2_liq(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importLNG'
            results.importLNG(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importLNG_regas'
            results.importLNG_regas(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importLNG_liq'
            results.importLNG_liq(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importFT'
            results.importFT(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importFT_kerosene'
            results.importFT_kerosene(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importFT_diesel'
            results.importFT_diesel(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importEL'
            results.importEL(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importMeOH'
            results.importMeOH(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'importNH3'
            results.importNH3(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

            %Power to heat
        case 'powerToTDHR'
            results.EL_TDHR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'powerToTDHP'
            results.EL_TDHP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'powerToTHHR'
            results.EL_THHR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'powerToTHHP'
            results.EL_THHP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'powerToTDGE'
            results.EL_TDGE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

            % Power to Gas

        case 'excess_El'
            results.EL_EXCESS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'excess_He'
            results.HE_EXCESS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'excess_HeLo'
            results.HE_EXCESS_Local(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'excess_HeLo_EP'
            results.HE_EXCESS_Local_EP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'excess_HeLo_HP'
            results.HE_EXCESS_Local_HP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'fromgrid'
            results.EL_GRID(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'fromgridToEl'
            results.EL_GRIDforEl(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'fromgridToHe'
            results.EL_GRIDforHe(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'electrolysis'
            results.EL_TWEL(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'HefromFT'
            results.TFTU_HE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'HefromMET'
            results.TMET_HE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'HefromWEL'
            results.TWEL_HE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'HefromNH3'
            results.TNH3_HE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'HefromMEO'
            results.TMEO_HE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'co2Scrubbing_El'
            results.EL_TCOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'co2Scrubbing_He'
            results.HE_TCOS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'NH_LINH'
            results.NH_LINH(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MeOtoMobility'
            results.MeOtoMobility(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'MeOtoMethanolFT'
            results.ME_TFTM(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'ME_LIME'
            results.ME_LIME(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'NH3toMobility'
            results.NH3toMobility(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'h2toMET'
            results.HY_TMET(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'h2toFT'
            results.HY_TFTU(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'h2toLH2'
            results.HY_TLH2(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'h2toMobility'
            results.HY_TRSP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'h2toMEO'
            results.HY_TMeO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'h2toNH3'
            results.HY_TNH3(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'h2toTISH'
            results.HY_TISH(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'h2toCCGT'
            results.HY_TCCG(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'h2toOCGT'
            results.HY_TOCG(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'h2toTGCS'
            results.HY_TGCS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'h2toTCNG'
            results.HY_TCNG(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'h2toIndHeat'
            results.HY_LHIN(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'h2toTICM'
            results.HY_TICM(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'h2toTRSi'
            results.HY_TRSi(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'CO2toMET'
            results.CO_TMET(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'CO2toFT'
            results.CO_TFTU(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'CO2toMEO'
            results.CO_TMeO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'CO2toCCS'
            results.CO_LCCS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'CO2fromSMC'
            results.TSMC_CO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'CO2fromGCS'
            results.TGCS_CO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'CO2fromHCS'
            results.THCS_CO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'CO2fromCBC'
            results.TCBC_CO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'CO2fromBCS'
            results.TBCS_CO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'CO2fromCWC'
            results.TCWC_CO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;


        case 'FTkeroseneProd'
            results.TFTU_KE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'FTnaphtaProd'
            results.TFTU_NA(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'FTdieselProd'
            results.TFTU_DI(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BioFTkeroseneProd'
            results.TFTB_KE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BioFTnaphtaProd'
            results.TFTB_NA(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BioFTdieselProd'
            results.TFTB_DI(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'methanolFTkeroseneProd'
            results.TFTM_KE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'methanolFTnaphtaProd'
            results.TFTM_NA(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'methanolFTdieselProd'
            results.TFTM_DI(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'BIOdiesel'
            results.BIOdiesel(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BIOkerosene'
            results.BIOkerosene(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'FTdiesel'
            results.FTdiesel(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'FTkerosene'
            results.FTkerosene(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BioFTdiesel'
            results.BioFTdiesel(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BioFTkerosene'
            results.BioFTkerosene(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'methanolFTdiesel'
            results.methanolFTdiesel(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'methanolFTkerosene'
            results.methanolFTkerosene(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'GasToLNG'
            results.GA_TLNG(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'EltoFT'
            results.EL_TFTU(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'EltoLNG'
            results.EL_TLNG(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'EltoLH2'
            results.EL_TLH2(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'EltoMET'
            results.EL_TMET(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'EltoMEO'
            results.EL_TMeO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'EltoNH3'
            results.EL_TNH3(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'EltoHyS'
            results.EL_THyS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'EltoBioFT'
            results.EL_TFTB(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'EltoMethanolFT'
            results.EL_TFTM(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'heatToPower'
            results.HE_TSTU(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'powerToHeat'
            results.EL_HE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'biomethane'
            results.RBME(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'geoheat'
            results.SoC_RGEO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'LossToHeat'
            results.Loss_HE(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

            %% Transport

        case 'iceFu_Dem'
            eval(['results.FU_' varIndicesSplit{3} '(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'bevEl_Dem'
            eval(['results.EL_' varIndicesSplit{3} '(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'fceFu_Dem'
            eval(['results.HY_' varIndicesSplit{3} '(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'phvFu_Dem'
            eval(['results.FU_' varIndicesSplit{3} '(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);
        case 'phvEl_Dem'
            eval(['results.EL_' varIndicesSplit{3} '(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;']);

        case 'elToLDV'
            results.elToLDV(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'LDV_chRate'
            results.LDV_chRate(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;
        case 'elToW23'
            results.elToW23(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'W23_chRate'
            results.W23_chRate(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;
        case 'elToBus'
            results.elToBus(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'BUS_chRate'
            results.BUS_chRate(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;
        case 'elToMDV'
            results.elToMDV(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MDV_chRate'
            results.MDV_chRate(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;
        case 'elToHDV'
            results.elToHDV(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'HDV_chRate'
            results.HDV_chRate(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;

        case 'elToLDV_V2G'
            results.elToLDV_V2G(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToW23_V2G'
            results.elToW23_V2G(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToBus_V2G'
            results.elToBus_V2G(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToMDV_V2G'
            results.elToMDV_V2G(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToHDV_V2G'
            results.elToHDV_V2G(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'elToLDV_G2V'
            results.elToLDV_G2V(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToW23_G2V'
            results.elToW23_G2V(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToBus_G2V'
            results.elToBus_G2V(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToMDV_G2V'
            results.elToMDV_G2V(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToHDV_G2V'
            results.elToHDV_G2V(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'elToLDV_DumpCh'
            results.elToLDV_DumpCh(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToW23_DumpCh'
            results.elToW23_DumpCh(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToBus_DumpCh'
            results.elToBus_DumpCh(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToMDV_DumpCh'
            results.elToMDV_DumpCh(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToHDV_DumpCh'
            results.elToHDV_DumpCh(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'elToLDV_SmartCh'
            results.elToLDV_SmartCh(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToW23_SmartCh'
            results.elToW23_SmartCh(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToBus_SmartCh'
            results.elToBus_SmartCh(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToMDV_SmartCh'
            results.elToMDV_SmartCh(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToHDV_SmartCh'
            results.elToHDV_SmartCh(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'elToLDV_DumpDisch'
            results.elToLDV_DumpDisch(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToW23_DumpDisch'
            results.elToW23_DumpDisch(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToBus_DumpDisch'
            results.elToBus_DumpDisch(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToMDV_DumpDisch'
            results.elToMDV_DumpDisch(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToHDV_DumpDisch'
            results.elToHDV_DumpDisch(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'elToLDV_SmartDisch'
            results.elToLDV_SmartDisch(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToW23_SmartDisch'
            results.elToW23_SmartDisch(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToBus_SmartDisch'
            results.elToBus_SmartDisch(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToMDV_SmartDisch'
            results.elToMDV_SmartDisch(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToHDV_SmartDisch'
            results.elToHDV_SmartDisch(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'MRPF_Dem'
            results.MRPF_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MRPE_Dem'
            results.MRPE_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MRFF_Dem'
            results.MRFF_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MRFE_Dem'
            results.MRFE_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'MMPF_Dem'
            results.MMPF_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MMPE_Dem'
            results.MMPE_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MMPH_Dem'
            results.MMPH_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MMPG_Dem'
            results.MMPG_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MMPA_Dem'
            results.MMPA_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MMPM_Dem'
            results.MMPM_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MMFF_Dem'
            results.MMFF_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MMFE_Dem'
            results.MMFE_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MMFH_Dem'
            results.MMFH_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MMFG_Dem'
            results.MMFG_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MMFA_Dem'
            results.MMFA_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MMFM_Dem'
            results.MMFM_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToMP'
            results.elToMP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToMF'
            results.elToMF(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MP_chRate'
            results.MP_chRate(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;
        case 'MF_chRate'
            results.MP_chRate(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;


        case 'MAPF_Dem'
            results.MAPF_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MAPE_Dem'
            results.MAPE_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MAPH_Dem'
            results.MAPH_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MAFF_Dem'
            results.MAFF_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MAFE_Dem'
            results.MAFE_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MAFH_Dem'
            results.MAFH_Dem(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToAP'
            results.elToAP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'elToAF'
            results.elToAF(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'AP_chRate'
            results.AP_chRate(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;
        case 'AF_chRate'
            results.AF_chRate(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;

        case 'EltoMobility'
            results.EltoMobility(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'diesToMobility'
            results.diesToMobility(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'kersToMobility'
            results.kersToMobility(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'LNGtoMobility'
            results.LNGtoMobility(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'LH2toMobility'
            results.LH2toMobility(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

            %%% Desalination

            % Water production
        case 'waterDesalination'
            results.Desalination(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'waterRODesalination'
            results.WROD_Desalination(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'waterMSSDesalination'
            results.WMSS_Desalination(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'waterMSCDesalination'
            results.WMSC_Desalination(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'waterMDSDesalination'
            results.WMDS_Desalination(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'waterMDCDesalination'
            results.WMDC_Desalination(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'desalination_Gas_demand'
            results.GA_WDES(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'desalination_Gas_demand_MSC'
            results.GA_WMSC(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'desalination_Gas_demand_MDC'
            results.GA_WMDC(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'desalination_Heat_demand'
            results.HE_WDES(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'desalination_Heat_demand_MSS'
            results.HE_WMSS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'desalination_Heat_demand_MDS'
            results.HE_WMDS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'desalination_El_demand'
            results.EL_WDES(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'desalination_El_demand_RO'
            results.EL_WROD(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'desalination_El_demand_MSS'
            results.EL_WMSS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'desalination_El_demand_MSC'
            results.EL_WMSC(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'desalination_El_demand_MDS'
            results.EL_WMDS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'desalination_El_demand_MDC'
            results.EL_WMDC(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'Vpump_demand'
            results.EL_WVPU(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'Hpump_demand'
            results.EL_WHPU(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'TotalDesalinationElDemand'
            results.TotDesDem(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'desalination_El_production'
            results.WDES_EL(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'desalination_El_production_MSC'
            results.WMSC_EL(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'desalination_El_production_MDC'
            results.WMDC_EL(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'waterStorageState'
            results.SoC_SWAT(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'instCapacitiy_Desalination'
            eval(['results.OPT_SIZE_' varIndicesSplit{2} '(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;']);
            %% industry
        case 'TICB_PR'
            results.TICB_PR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TICI_PR'
            results.TICI_PR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TISB_PR'
            results.TISB_PR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TISH_PR'
            results.TISH_PR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TISR_PR'
            results.TISR_PR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TISE_PR'
            results.TISE_PR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TIAA_PR'
            results.TIAA_PR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TIAM_PR'
            results.TIAM_PR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TIAR_PR'
            results.TIAR_PR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TIPP_PR'
            results.TIPP_PR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'EL_TICB'
            results.EL_TICB(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'EL_TICI'
            results.EL_TICI(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'EL_TISB'
            results.EL_TISB(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'EL_TISH'
            results.EL_TISH(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'EL_TISR'
            results.EL_TISR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'EL_TISE'
            results.EL_TISE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'EL_TIAM'
            results.EL_TIAM(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'EL_TIAR'
            results.EL_TIAR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'EL_TIPP'
            results.EL_TIPP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'HE_TICB'
            results.HE_TICB(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'HE_TICI'
            results.HE_TICI(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'HE_TISB'
            results.HE_TISB(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'HE_TISH'
            results.HE_TISH(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'HE_TISR'
            results.HE_TISR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'HE_TISE'
            results.HE_TISE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'HE_TIAA'
            results.HE_TIAA(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'HE_TIAR'
            results.HE_TIAR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'HE_TIPP'
            results.HE_TIPP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'HC_TISB'
            results.HC_TISB(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'CC_TISH'
            results.CC_TISH(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'CC_TISR'
            results.CC_TISR(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'CC_TISE'
            results.CC_TISE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'TIAM_HE'
            results.TIAM_HE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TIAR_HE'
            results.TIAR_HE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'El_to_TPSC'
            results.EL_TPSC(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'El_to_TPSP'
            results.EL_TPSP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'He_to_TPSC'
            results.HE_TPSC(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'He_to_TPSP'
            results.HE_TPSP(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'Wa_to_TRAD'
            results.Wa_to_TRAD(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'Wa_to_TRMI'
            results.Wa_to_TRMI(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'Wa_to_TRMO'
            results.Wa_to_TRMO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'Wa_to_TRBC'
            results.Wa_to_TRBC(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'El_to_TRSS'
            results.EL_TRSS(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'El_to_TRSL'
            results.EL_TRSL(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'El_to_TRMI'
            results.EL_TRMI(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'El_to_TRMO'
            results.EL_TRMO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'El_to_TRME'
            results.EL_TRME(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'El_to_TRSi'
            results.EL_TRSi(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'El_to_TRAD'
            results.EL_TRAD(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'El_to_TREW'
            results.EL_TREW(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'El_to_TRBC'
            results.EL_TRBC(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'El_to_TRGE'
            results.EL_TRGE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'He_to_TRME'
            results.HE_TRME(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'He_to_TRSi'
            results.HE_TRSi(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'He_to_TRGE'
            results.HE_TRGE(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TICB_CO'
            results.TICB_CO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TICI_CO'
            results.TICI_CO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TISB_CO'
            results.TISB_CO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TISR_CO'
            results.TISR_CO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TISH_CO'
            results.TISH_CO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TISE_CO'
            results.TISE_CO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TIAA_CO'
            results.TIAA_CO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TIAR_CO'
            results.TIAR_CO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TIAM_CO'
            results.TIAM_CO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'TIPP_CO'
            results.TIPP_CO(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'CO2EmissionsfromInd'
            results.CO2EmissionsfromInd(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'EltoInd'
            results.EltoInd(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'HighHetoInd'
            results.HighHetoInd(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MidHetoInd'
            results.MidHetoInd(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'MidHefromInd'
            results.MidHefromInd(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'LowHefromInd'
            results.LowHefromInd(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'TIPP_WWtot'
            results.TIPP_WWtot(str2num(varIndicesSplit{1})) = values(varNo)/m_factor;

        case 'addElPlus'
            results.addElPlus(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'addElMinus'
            results.addElMinus(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
        case 'addElSoC'
            results.addElSoC(str2num(varIndicesSplit{1}),str2num(varIndicesSplit{2})) = values(varNo)/m_factor;

        case 'extraHydroRoR'
            results.extraHydroRoR(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo)/m_factor;
            %%
            %%Capacities


    end
    clear varName varIndices  varIndicesSplit check;
end

results.totalCost = totalCost;

