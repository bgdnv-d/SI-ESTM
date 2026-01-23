function activeElements = writeCostFileCSP(fileName,setup,struct,activeElements,costYear)
% writes gnuMathProg conform data files
%	INPUT:
%		fileName - desired filename without extension
%		struct - parameter structure obtained from excel sheets
%		activeElements - structure with active elements in the system
%
%Dmitrii Bogdanov
%last change 24.07.2025


% electrical transformer
capexElTransformerActive = repmat((struct.Capex(activeElements.activeElTransformer,(struct.IndexYears==costYear)).*CapitalRecoveryFactor(struct.WACC,struct.Lifetime(activeElements.activeElTransformer,(struct.IndexYears==costYear))))',length(struct.IndexNodes),1).*struct.Capex_reg(activeElements.activeElTransformer,:,(struct.IndexYears==costYear))';
opex_fixElTransformerActive = repmat(struct.Opex_fix(activeElements.activeElTransformer,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_fix_reg(activeElements.activeElTransformer,:,(struct.IndexYears==costYear))';
opex_varElTransformerActive = repmat(struct.Opex_var(activeElements.activeElTransformer,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_var_reg(activeElements.activeElTransformer,:,(struct.IndexYears==costYear))';

%% change OGCT and ICG cost for infeasibility suppression
try setup.flexStart;
catch
    try
        setup.flexStart = setup.flex2015;
    catch
        setup.flexStart = 0;
    end
end



% heat transformer
capexHeatTransformerActive = repmat((struct.Capex(activeElements.activeHeatTransformer,(struct.IndexYears==costYear)).*CapitalRecoveryFactor(struct.WACC,struct.Lifetime(activeElements.activeHeatTransformer,(struct.IndexYears==costYear))))',length(struct.IndexNodes),1).*struct.Capex_reg(activeElements.activeHeatTransformer,:,(struct.IndexYears==costYear))';
opex_fixHeatTransformerActive = repmat(struct.Opex_fix(activeElements.activeHeatTransformer,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_fix_reg(activeElements.activeHeatTransformer,:,(struct.IndexYears==costYear))';
opex_varHeatTransformerActive = repmat(struct.Opex_var(activeElements.activeHeatTransformer,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_var_reg(activeElements.activeHeatTransformer,:,(struct.IndexYears==costYear))';


% CHP transformer

capexCHPTransformerActive = repmat((struct.Capex(activeElements.activeCHPTransformer,(struct.IndexYears==costYear)).*CapitalRecoveryFactor(struct.WACC,struct.Lifetime(activeElements.activeCHPTransformer,(struct.IndexYears==costYear))))',length(struct.IndexNodes),1).*struct.Capex_reg(activeElements.activeCHPTransformer,:,(struct.IndexYears==costYear))';
opex_fixCHPTransformerActive = repmat(struct.Opex_fix(activeElements.activeCHPTransformer,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_fix_reg(activeElements.activeCHPTransformer,:,(struct.IndexYears==costYear))';
opex_varCHPTransformerActive = repmat(struct.Opex_var(activeElements.activeCHPTransformer,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_var_reg(activeElements.activeCHPTransformer,:,(struct.IndexYears==costYear))';


% gas transformer
capexGasTransformerActive = repmat((struct.Capex(activeElements.activeGasTransformer,(struct.IndexYears==costYear)).*CapitalRecoveryFactor(struct.WACC,struct.Lifetime(activeElements.activeGasTransformer,(struct.IndexYears==costYear))))',length(struct.IndexNodes),1).*struct.Capex_reg(activeElements.activeGasTransformer,:,(struct.IndexYears==costYear))';
opex_fixGasTransformerActive = repmat(struct.Opex_fix(activeElements.activeGasTransformer,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_fix_reg(activeElements.activeGasTransformer,:,(struct.IndexYears==costYear))';
opex_varGasTransformerActive = repmat(struct.Opex_var(activeElements.activeGasTransformer,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_var_reg(activeElements.activeGasTransformer,:,(struct.IndexYears==costYear))';

% storage
capexStorageActive = repmat((struct.Capex(activeElements.activeStorage,(struct.IndexYears==costYear)).*CapitalRecoveryFactor(struct.WACC,struct.Lifetime(activeElements.activeStorage,(struct.IndexYears==costYear))))',length(struct.IndexNodes),1).*struct.Capex_reg(activeElements.activeStorage,:,(struct.IndexYears==costYear))';
opex_fixStorageActive = repmat(struct.Opex_fix(activeElements.activeStorage,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_fix_reg(activeElements.activeStorage,:,(struct.IndexYears==costYear))';
opex_varStorageActive = repmat(struct.Opex_var(activeElements.activeStorage,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_var_reg(activeElements.activeStorage,:,(struct.IndexYears==costYear))';
try
    % storage inteface
    capexStorageIntActive = repmat((struct.Capex(activeElements.activeStorageInterface,(struct.IndexYears==costYear)).*CapitalRecoveryFactor(struct.WACC,struct.Lifetime(activeElements.activeStorageInterface,(struct.IndexYears==costYear))))',length(struct.IndexNodes),1).*struct.Capex_reg(activeElements.activeStorageInterface,:,(struct.IndexYears==costYear))';
    opex_fixStorageIntActive = repmat(struct.Opex_fix(activeElements.activeStorageInterface,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_fix_reg(activeElements.activeStorageInterface,:,(struct.IndexYears==costYear))';
    opex_varStorageIntActive = repmat(struct.Opex_var(activeElements.activeStorageInterface,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_var_reg(activeElements.activeStorageInterface,:,(struct.IndexYears==costYear))';
catch
end
% resource
capexResourceActive = repmat((struct.Capex(activeElements.activeResource,(struct.IndexYears==costYear)).*CapitalRecoveryFactor(struct.WACC,struct.Lifetime(activeElements.activeResource,(struct.IndexYears==costYear))))',length(struct.IndexNodes),1).*struct.Capex_reg(activeElements.activeResource,:,(struct.IndexYears==costYear))';
opex_fixResourceActive = repmat(struct.Opex_fix(activeElements.activeResource,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_fix_reg(activeElements.activeResource,:,(struct.IndexYears==costYear))';
opex_varResourceActive = repmat(struct.Opex_var(activeElements.activeResource,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_var_reg(activeElements.activeResource,:,(struct.IndexYears==costYear))';

% electrical feed-in resource
capexElectricalFeedInActive = repmat((struct.Capex(activeElements.activeElFeedIn,(struct.IndexYears==costYear)).*CapitalRecoveryFactor(struct.WACC,struct.Lifetime(activeElements.activeElFeedIn,(struct.IndexYears==costYear))))',length(struct.IndexNodes),1).*struct.Capex_reg(activeElements.activeElFeedIn,:,(struct.IndexYears==costYear))';
opex_fixElectricalFeedInActive = repmat(struct.Opex_fix(activeElements.activeElFeedIn,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_fix_reg(activeElements.activeElFeedIn,:,(struct.IndexYears==costYear))';
opex_varElectricalFeedInActive = repmat(struct.Opex_var(activeElements.activeElFeedIn,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_var_reg(activeElements.activeElFeedIn,:,(struct.IndexYears==costYear))';

% heat feed-in resource
capexHeatFeedInActive = repmat((struct.Capex(activeElements.activeHeatFeedIn,(struct.IndexYears==costYear)).*CapitalRecoveryFactor(struct.WACC,struct.Lifetime(activeElements.activeHeatFeedIn,(struct.IndexYears==costYear))))',length(struct.IndexNodes),1).*struct.Capex_reg(activeElements.activeHeatFeedIn,:,(struct.IndexYears==costYear))';
opex_fixHeatFeedInActive = repmat(struct.Opex_fix(activeElements.activeHeatFeedIn,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_fix_reg(activeElements.activeHeatFeedIn,:,(struct.IndexYears==costYear))';
opex_varHeatFeedInActive = repmat(struct.Opex_var(activeElements.activeHeatFeedIn,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_var_reg(activeElements.activeHeatFeedIn,:,(struct.IndexYears==costYear))';

% hydro resource (with dams)
capexHydroActive = repmat((struct.Capex(activeElements.activeHydro,(struct.IndexYears==costYear)).*CapitalRecoveryFactor(struct.WACC,struct.Lifetime(activeElements.activeHydro,(struct.IndexYears==costYear))))',length(struct.IndexNodes),1).*struct.Capex_reg(activeElements.activeHydro,:,(struct.IndexYears==costYear))';
opex_fixHydroActive = repmat(struct.Opex_fix(activeElements.activeHydro,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_fix_reg(activeElements.activeHydro,:,(struct.IndexYears==costYear))';
opex_varHydroActive = repmat(struct.Opex_var(activeElements.activeHydro,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_var_reg(activeElements.activeHydro,:,(struct.IndexYears==costYear))';


% transmission
TLlengthIndex = max(1,length(struct.TLlength));
capexTransmission = repmat((struct.Capex(activeElements.activeTransmission,(struct.IndexYears==costYear)).*CapitalRecoveryFactor(struct.WACC,struct.Lifetime(activeElements.activeTransmission,(struct.IndexYears==costYear))))',TLlengthIndex,1);
opex_fixTransmission = repmat(struct.Opex_fix(activeElements.activeTransmission,(struct.IndexYears==costYear))',TLlengthIndex,1);
opex_varTransmission = repmat(struct.Opex_var(activeElements.activeTransmission,(struct.IndexYears==costYear))',TLlengthIndex,1);


% Desalination
capexDesalinationActive = repmat((struct.Capex(activeElements.activeDesalination,(struct.IndexYears==costYear)).*CapitalRecoveryFactor(struct.WACC,struct.Lifetime(activeElements.activeDesalination,(struct.IndexYears==costYear))))',length(struct.IndexNodes),1).*struct.Capex_reg(activeElements.activeDesalination,:,(struct.IndexYears==costYear))';
opex_fixDesalinationActive = repmat(struct.Opex_fix(activeElements.activeDesalination,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_fix_reg(activeElements.activeDesalination,:,(struct.IndexYears==costYear))';
opex_varDesalinationActive = repmat(struct.Opex_var(activeElements.activeDesalination,(struct.IndexYears==costYear))',length(struct.IndexNodes),1).*struct.Opex_var_reg(activeElements.activeDesalination,:,(struct.IndexYears==costYear))';

%% speciall WACC for coal and Nuclear

try
    capexElTransformerActive(:,ismember(activeElements.labels.elTransformer,'THPP')) = repmat((struct.Capex(ismember(struct.IndexID,'THPP'),(struct.IndexYears==costYear)).*CapitalRecoveryFactor(struct.WACC_Spec,struct.Lifetime(ismember(struct.IndexID,'THPP'),(struct.IndexYears==costYear))))',length(struct.IndexNodes),1).*struct.Capex_reg(ismember(struct.IndexID,'THPP'),:,(struct.IndexYears==costYear))';
    capexElTransformerActive(:,ismember(activeElements.labels.elTransformer,'TNUC')) = repmat((struct.Capex(ismember(struct.IndexID,'TNUC'),(struct.IndexYears==costYear)).*CapitalRecoveryFactor(struct.WACC_Spec,struct.Lifetime(ismember(struct.IndexID,'TNUC'),(struct.IndexYears==costYear))))',length(struct.IndexNodes),1).*struct.Capex_reg(ismember(struct.IndexID,'TNUC'),:,(struct.IndexYears==costYear))';
    capexCHPTransformerActive(:,ismember(activeElements.labels.chpTransformer,'TCCO')) = repmat((struct.Capex(ismember(struct.IndexID,'TCCO'),(struct.IndexYears==costYear)).*CapitalRecoveryFactor(struct.WACC_Spec,struct.Lifetime(ismember(struct.IndexID,'TCCO'),(struct.IndexYears==costYear))))',length(struct.IndexNodes),1).*struct.Capex_reg(ismember(struct.IndexID,'TCCO'),:,(struct.IndexYears==costYear))';

    capexHeatTransformerActive(:,ismember(activeElements.labels.heatTransformer,'TDCO')) = repmat((struct.Capex(ismember(struct.IndexID,'TDCO'),(struct.IndexYears==costYear)).*CapitalRecoveryFactor(struct.WACC_Spec,struct.Lifetime(ismember(struct.IndexID,'TDCO'),(struct.IndexYears==costYear))))',length(struct.IndexNodes),1).*struct.Capex_reg(ismember(struct.IndexID,'TDCO'),:,(struct.IndexYears==costYear))';

catch


end

if costYear == setup.startYear
    if setup.flexStart

        capexElTransformerActive(:,find(ismember(activeElements.labels.elTransformer,'TOCG'))) = 1999;
        capexElTransformerActive(:,find(ismember(activeElements.labels.elTransformer,'TICG'))) = 1999;
    end
    capexHeatTransformerActive(:,find(ismember(activeElements.labels.heatTransformer,'THOI'))) = 999;
end



% write financial parameter file
WriteParamFileBegin(fileName);

WriteParams(fileName,'Capex_electTransf',round(capexElTransformerActive,3),[1:length(struct.IndexNodes)],activeElements.labels.elTransformer)
WriteParams(fileName,'Opex_fix_electTransf',round(opex_fixElTransformerActive,3),[1:length(struct.IndexNodes)],activeElements.labels.elTransformer)
WriteParams(fileName,'Opex_var_electTransf',round(opex_varElTransformerActive*1000,3),[1:length(struct.IndexNodes)],activeElements.labels.elTransformer)

WriteParams(fileName,'Capex_heatTransf',round(capexHeatTransformerActive,3),[1:length(struct.IndexNodes)],activeElements.labels.heatTransformer)
WriteParams(fileName,'Opex_fix_heatTransf',round(opex_fixHeatTransformerActive,3),[1:length(struct.IndexNodes)],activeElements.labels.heatTransformer)
WriteParams(fileName,'Opex_var_heatTransf',round(opex_varHeatTransformerActive*1000,3),[1:length(struct.IndexNodes)],activeElements.labels.heatTransformer)

WriteParams(fileName,'Capex_chpTransf',round(capexCHPTransformerActive,3),[1:length(struct.IndexNodes)],activeElements.labels.chpTransformer)
WriteParams(fileName,'Opex_fix_chpTransf',round(opex_fixCHPTransformerActive,3),[1:length(struct.IndexNodes)],activeElements.labels.chpTransformer)
WriteParams(fileName,'Opex_var_chpTransf',round(opex_varCHPTransformerActive*1000,3),[1:length(struct.IndexNodes)],activeElements.labels.chpTransformer)

WriteParams(fileName,'Capex_ptgTransf',round(capexGasTransformerActive,3),[1:length(struct.IndexNodes)],activeElements.labels.gasTransformer)
WriteParams(fileName,'Opex_fix_ptgTransf',round(opex_fixGasTransformerActive,3),[1:length(struct.IndexNodes)],activeElements.labels.gasTransformer)
WriteParams(fileName,'Opex_var_ptgTransf',round(opex_varGasTransformerActive*1000,3),[1:length(struct.IndexNodes)],activeElements.labels.gasTransformer)

WriteParams(fileName,'Capex_electFeedIn',round(capexElectricalFeedInActive,3),[1:length(struct.IndexNodes)],activeElements.labels.elFeedIn)
WriteParams(fileName,'Opex_fix_electFeedIn',round(opex_fixElectricalFeedInActive,3),[1:length(struct.IndexNodes)],activeElements.labels.elFeedIn)
WriteParams(fileName,'Opex_var_electFeedIn',round(opex_varElectricalFeedInActive*1000,3),[1:length(struct.IndexNodes)],activeElements.labels.elFeedIn)

WriteParams(fileName,'Capex_heatFeedIn',round(capexHeatFeedInActive,3),[1:length(struct.IndexNodes)],activeElements.labels.heatFeedIn)
WriteParams(fileName,'Opex_fix_heatFeedIn',round(opex_fixHeatFeedInActive,3),[1:length(struct.IndexNodes)],activeElements.labels.heatFeedIn)
WriteParams(fileName,'Opex_var_heatFeedIn',round(opex_varHeatFeedInActive*1000,3),[1:length(struct.IndexNodes)],activeElements.labels.heatFeedIn)

WriteParams(fileName,'Capex_hydro',round(capexHydroActive,3),[1:length(struct.IndexNodes)],activeElements.labels.hydro)
WriteParams(fileName,'Opex_fix_hydro',round(opex_fixHydroActive,3),[1:length(struct.IndexNodes)],activeElements.labels.hydro)
WriteParams(fileName,'Opex_var_hydro',round(opex_varHydroActive*1000,3),[1:length(struct.IndexNodes)],activeElements.labels.hydro)

WriteParams(fileName,'Capex_storage',round(capexStorageActive,3),[1:length(struct.IndexNodes)],activeElements.labels.storage)
WriteParams(fileName,'Opex_fix_storage',round(opex_fixStorageActive,3),[1:length(struct.IndexNodes)],activeElements.labels.storage)
WriteParams(fileName,'Opex_var_storage',round(opex_varStorageActive*1000,3),[1:length(struct.IndexNodes)],activeElements.labels.storage)
try
    WriteParams(fileName,'Capex_storint',round(capexStorageIntActive,3),[1:length(struct.IndexNodes)],activeElements.labels.storageInterface)
    WriteParams(fileName,'Opex_fix_storint',round(opex_fixStorageIntActive,3),[1:length(struct.IndexNodes)],activeElements.labels.storageInterface)
    WriteParams(fileName,'Opex_var_storint',round(opex_varStorageIntActive*1000,3),[1:length(struct.IndexNodes)],activeElements.labels.storageInterface)
catch
end
WriteParams(fileName,'Capex_resource',round(capexResourceActive,3),[1:length(struct.IndexNodes)],activeElements.labels.resource)
WriteParams(fileName,'Opex_fix_resource',round(opex_fixResourceActive,3),[1:length(struct.IndexNodes)],activeElements.labels.resource)
WriteParams(fileName,'Opex_var_resource',round(opex_varResourceActive*1000,3),[1:length(struct.IndexNodes)],activeElements.labels.resource)

WriteParams(fileName,'Capex_desal',round(capexDesalinationActive,3),[1:length(struct.IndexNodes)],activeElements.labels.desalination)
WriteParams(fileName,'Opex_fix_desal',round(opex_fixDesalinationActive,3),[1:length(struct.IndexNodes)],activeElements.labels.desalination)
WriteParams(fileName,'Opex_var_desal',opex_varDesalinationActive,[1:length(struct.IndexNodes)],activeElements.labels.desalination)


WriteParams(fileName,'Capex_TransmLine',round(capexTransmission,3),[1:TLlengthIndex],activeElements.labels.transmission)
WriteParams(fileName,'Opex_fix_TransmLine',round(opex_fixTransmission,3),[1:TLlengthIndex],activeElements.labels.transmission)
WriteParams(fileName,'Opex_var_TransmLine',round(opex_varTransmission,3),[1:TLlengthIndex],activeElements.labels.transmission)

WriteParamFileEnd(fileName);
