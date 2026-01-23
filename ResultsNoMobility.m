function results = results_WO_Mobility(results, numNodes, endHour)
%
% FUNCTION results_WO_Mobility(results, numNodes, endHour)
%
% Removes mobility-related results from the dataset.
%
%
% INPUT:
%            results:     Structure containing the original model results.
%            numNodes:    Number of nodes or regions in the model.
%            endHour:     Last hour of the simulation time period.
%
% OUTPUT:
%            results:     Updated results structure without mobility related data.
%
%Dmitrii Bogdanov
%last change 24.07.2025


%% Transport


results.FU_MRBI = zeros(endHour,numNodes);
results.FU_MRBP = zeros(endHour,numNodes);
results.FU_MRHI = zeros(endHour,numNodes);
results.FU_MRHP = zeros(endHour,numNodes);
results.FU_MRLI = zeros(endHour,numNodes);
results.FU_MRLP = zeros(endHour,numNodes);
results.FU_MRMI = zeros(endHour,numNodes);
results.FU_MRMP = zeros(endHour,numNodes);
results.FU_MRWI = zeros(endHour,numNodes);
results.FU_MRWP = zeros(endHour,numNodes);

results.EL_MRBB = zeros(endHour,numNodes);
results.EL_MRBP = zeros(endHour,numNodes);
results.EL_MRHB = zeros(endHour,numNodes);
results.EL_MRHP = zeros(endHour,numNodes);
results.EL_MRLB = zeros(endHour,numNodes);
results.EL_MRLP = zeros(endHour,numNodes);
results.EL_MRMB = zeros(endHour,numNodes);
results.EL_MRMP = zeros(endHour,numNodes);
results.EL_MRWB = zeros(endHour,numNodes);
results.EL_MRWP = zeros(endHour,numNodes);

results.HY_MRBF = zeros(endHour,numNodes);
results.HY_MRHF = zeros(endHour,numNodes);
results.HY_MRLF = zeros(endHour,numNodes);
results.HY_MRMF = zeros(endHour,numNodes);
results.HY_MRWF = zeros(endHour,numNodes);

results.elToLDV = zeros(endHour,numNodes);
results.LDV_chRate = zeros(1,numNodes);
results.elToW23 = zeros(endHour,numNodes);
results.W23_chRate = zeros(1,numNodes);
results.elToBus = zeros(endHour,numNodes);
results.BUS_chRate = zeros(1,numNodes);
results.elToMDV = zeros(endHour,numNodes);
results.MDV_chRate = zeros(1,numNodes);
results.elToHDV = zeros(endHour,numNodes);
results.HDV_chRate = zeros(1,numNodes);

results.MRPF_Dem = zeros(endHour,numNodes);
results.MRPE_Dem = zeros(endHour,numNodes);
results.MRFF_Dem = zeros(endHour,numNodes);
results.MRFE_Dem = zeros(endHour,numNodes);

results.MMPF_Dem = zeros(endHour,numNodes);
results.MMPE_Dem = zeros(endHour,numNodes);
results.MMPH_Dem = zeros(endHour,numNodes);
results.MMPG_Dem = zeros(endHour,numNodes);
results.MMPA_Dem = zeros(endHour,numNodes);
results.MMPM_Dem = zeros(endHour,numNodes);
results.MMFF_Dem = zeros(endHour,numNodes);
results.MMFE_Dem = zeros(endHour,numNodes);
results.MMFH_Dem = zeros(endHour,numNodes);
results.MMFG_Dem = zeros(endHour,numNodes);
results.MMFA_Dem = zeros(endHour,numNodes);
results.MMFM_Dem = zeros(endHour,numNodes);
results.elToMP = zeros(endHour,numNodes);
results.elToMF = zeros(endHour,numNodes);
results.MP_chRate = zeros(1,numNodes);
results.MP_chRate = zeros(1,numNodes);


results.MAPF_Dem = zeros(endHour,numNodes);
results.MAPE_Dem = zeros(endHour,numNodes);
results.MAPH_Dem = zeros(endHour,numNodes);
results.MAFF_Dem = zeros(endHour,numNodes);
results.MAFE_Dem = zeros(endHour,numNodes);
results.MAFH_Dem = zeros(endHour,numNodes);
results.elToAP = zeros(endHour,numNodes);
results.elToAF = zeros(endHour,numNodes);
results.AP_chRate = zeros(1,numNodes);
results.AF_chRate = zeros(1,numNodes);

results.EltoMobility = zeros(endHour,numNodes);
results.diesToMobility = zeros(endHour,numNodes);
results.kersToMobility = zeros(endHour,numNodes);
results.LNGtoMobility = zeros(endHour,numNodes);
results.LH2toMobility = zeros(endHour,numNodes);
