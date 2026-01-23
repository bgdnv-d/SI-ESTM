function writePhysicalData(fileName,setup,systemParams,activeElements,costYear)

% write physical parameters to gmpl conform data file
%
%Dmitrii Bogdanov
%last change 24.07.2025


allTransf = logical(activeElements.activeElTransformer + activeElements.activeCHPTransformer + activeElements.activeHeatTransformer  + activeElements.activeGasTransformer);
% reshaping of data structure
for h=1:length(systemParams.IndexFuels)
    TMPefficiencyTransf(1,[1:length(systemParams.IndexID(allTransf))],h) = systemParams.EtaTrans(allTransf,h)';
end
efficiencyTransf = repmat(TMPefficiencyTransf,[length(systemParams.IndexNodes) 1 1]);
efficiencyTransf(isnan(efficiencyTransf)) = 0;


WriteParamFileBegin(fileName);

% storage params
WriteParams(fileName,'efficiencyIn',round(systemParams.EtaStorage(ismember(systemParams.IndexIDS,activeElements.labels.storage),1),3),activeElements.labels.storage)
WriteParams(fileName,'efficiencyOut',round(1./systemParams.EtaStorage(ismember(systemParams.IndexIDS,activeElements.labels.storage),2),3),activeElements.labels.storage)
WriteParams(fileName,'energyPowerRatioIn',round(systemParams.EnPoRatioStorage(ismember(systemParams.IndexIDS,activeElements.labels.storage),1),3),activeElements.labels.storage)
WriteParams(fileName,'energyPowerRatioOut',round(systemParams.EnPoRatioStorage(ismember(systemParams.IndexIDS,activeElements.labels.storage),2),3),activeElements.labels.storage)

WriteParams(fileName,'selfDischarge',round(systemParams.SelfDischargeStorage(ismember(systemParams.IndexIDS,activeElements.labels.storage),1),4),activeElements.labels.storage)


% power plant (transformer) params
WriteParams(fileName,'efficiencyTransf',round(efficiencyTransf,3),[1:length(systemParams.IndexNodes)],[systemParams.IndexID(allTransf)],systemParams.IndexFuels');
WriteParams(fileName,'rampability',round(systemParams.RampTrans(allTransf),3),[systemParams.IndexID(allTransf)]);

try setup.ModelType;
catch
    setup.ModelType = 'Main';
end
%if strcmp(setup.ModelType,'Main')
switch setup.ModelType
    case 'Main'
        if setup.Mobility
            WriteParams(fileName,'transpConsPrim',shiftdim(round(systemParams.Mobility.Cons_Primary_S(:,systemParams.IndexYears==costYear,:),4),2)',[1:length(systemParams.IndexNodes)],systemParams.Mobility.Cons_Names_Primary)
            WriteParams(fileName,'transpConsSecond',shiftdim(round(systemParams.Mobility.Cons_Secondary_S(:,systemParams.IndexYears==costYear,:),4),2)',[1:length(systemParams.IndexNodes)],systemParams.Mobility.Cons_Names_Secondary)

            WriteParams(fileName,'transpTypeShares',shiftdim(round(systemParams.Mobility.SharesTot(:,systemParams.IndexYears==costYear,ismember(systemParams.IndexID,systemParams.IndexIDM)),4),2)',[1:length(systemParams.IndexNodes)],systemParams.IndexIDM);

            WriteParams(fileName,'shareFuMRLP',round(1-systemParams.Mobility.LDV_PHEV_el(:,systemParams.IndexYears==costYear),4),[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'shareFuMRWP',round(1-systemParams.Mobility.W23_PHEV_el(:,systemParams.IndexYears==costYear),4),[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'shareFuMRBP',round(1-systemParams.Mobility.BUS_PHEV_el(:,systemParams.IndexYears==costYear),4),[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'shareFuMRMP',round(1-systemParams.Mobility.MDV_PHEV_el(:,systemParams.IndexYears==costYear),4),[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'shareFuMRHP',round(1-systemParams.Mobility.HDV_PHEV_el(:,systemParams.IndexYears==costYear),4),[1:length(systemParams.IndexNodes)]);

            WriteParams(fileName,'shareElMRLP',round(systemParams.Mobility.LDV_PHEV_el(:,systemParams.IndexYears==costYear),4),[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'shareElMRWP',round(systemParams.Mobility.W23_PHEV_el(:,systemParams.IndexYears==costYear),4),[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'shareElMRBP',round(systemParams.Mobility.BUS_PHEV_el(:,systemParams.IndexYears==costYear),4),[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'shareElMRMP',round(systemParams.Mobility.MDV_PHEV_el(:,systemParams.IndexYears==costYear),4),[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'shareElMRHP',round(systemParams.Mobility.HDV_PHEV_el(:,systemParams.IndexYears==costYear),4),[1:length(systemParams.IndexNodes)]);
        end

        WriteParams(fileName,'TDHP_COP',round(systemParams.HP.COP_DH(:,:,systemParams.IndexYears==costYear),4),[1:size(systemParams.HP.COP_DH(:,:,systemParams.IndexYears==costYear),1)],[1:size(systemParams.HP.COP_DH(:,:,systemParams.IndexYears==costYear),2)]);
    otherwise
        WriteParams(fileName,'THHP_COP',round(systemParams.HP.COP_IH(:,:,systemParams.IndexYears==costYear),4),[1:size(systemParams.HP.COP_IH(:,:,systemParams.IndexYears==costYear),1)],[1:size(systemParams.HP.COP_DH(:,:,systemParams.IndexYears==costYear),2)]);
end


WriteParams(fileName,'GasEmissions',systemParams.GasEmissionsMain);
WriteParams(fileName,'GasEmissionsAdditional',systemParams.GasEmissionsAdditional);
WriteParams(fileName,'CoalEmissions',systemParams.CoalEmissionsMain);
WriteParams(fileName,'CoalEmissionsAdditional',systemParams.CoalEmissionsAdditional);
WriteParams(fileName,'OilEmissions',systemParams.OilEmissionsMain);
WriteParams(fileName,'OilEmissionsAdditional',systemParams.OilEmissionsAdditional);
WriteParams(fileName,'BiomassEmissions',systemParams.BiomassEmissionsMain);
WriteParams(fileName,'BiomassEmissionsAdditional',systemParams.BiomassEmissionsAdditional);
WriteParams(fileName,'WastesEmissions',systemParams.WastesEmissionsMain);
WriteParams(fileName,'WastesEmissionsAdditional',systemParams.WastesEmissionsAdditional);

WriteParams(fileName,'ICB_Emissions',systemParams.TICB_EmissionsMain);
WriteParams(fileName,'ICB_EmissionsAdditional',systemParams.TICB_EmissionsAdditional);
WriteParams(fileName,'ICI_Emissions',systemParams.TICI_EmissionsMain);
WriteParams(fileName,'ICI_EmissionsAdditional',systemParams.TICI_EmissionsAdditional);
WriteParams(fileName,'ISB_Emissions',systemParams.TISB_EmissionsMain);
WriteParams(fileName,'ISB_EmissionsAdditional',systemParams.TISB_EmissionsAdditional);
WriteParams(fileName,'ISR_Emissions',systemParams.TISR_EmissionsMain);
WriteParams(fileName,'ISR_EmissionsAdditional',systemParams.TISR_EmissionsAdditional);
WriteParams(fileName,'ISE_Emissions',systemParams.TISE_EmissionsMain);
WriteParams(fileName,'ISE_EmissionsAdditional',systemParams.TISE_EmissionsAdditional);
WriteParams(fileName,'ISH_Emissions',systemParams.TISH_EmissionsMain);
WriteParams(fileName,'ISH_EmissionsAdditional',systemParams.TISH_EmissionsAdditional);
WriteParams(fileName,'IAA_Emissions',systemParams.TIAA_EmissionsMain);
WriteParams(fileName,'IAA_EmissionsAdditional',systemParams.TIAA_EmissionsAdditional);
WriteParams(fileName,'IAM_Emissions',systemParams.TIAM_EmissionsMain);
WriteParams(fileName,'IAM_EmissionsAdditional',systemParams.TIAM_EmissionsAdditional);
WriteParams(fileName,'IAR_Emissions',systemParams.TIAR_EmissionsMain);
WriteParams(fileName,'IAR_EmissionsAdditional',systemParams.TIAR_EmissionsAdditional);
WriteParams(fileName,'IPP_Emissions',systemParams.TIPP_EmissionsMain);
WriteParams(fileName,'IPP_EmissionsAdditional',systemParams.TIPP_EmissionsAdditional);
WriteParamFileEnd(fileName);