function output = GetSystemFromExcelTransition(setup,files)
% FUNCTION GetSystemFromExcelTransition(setup,files)
% Interface for loading system data from excel files to system parameters structure
% INPUTS:
% 			setup: 	model setup structure containing path to folder containing excel files
% 			files: 	struct, contains names of excel files
%  				files.financial: string, e.g., 'Modelparams financial.xlsx'
%  				files.limits: string, e.g.,setup 'Modelparams limits.xlsx'
%  				files.physical: string, e.g., 'Modelparams physical.xlsx'
%  			files.geographical: string, e.g., 'Modelparams geographical.xlsx'
% OUTPUTS:
% 			output: system parameters data structure
% Dmitrii Bogdanov
% Refactored for modern MATLAB performance (xlsread replaced)

path = [setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep];
nr_of_hours=8760; %changing this number to values other than 8760 might produce unexpected results

% Read Excel Files -------------------------------------------------------
disp('Reading Excel Files');

%% Physical
[numPT,txtPT,~] = local_xlsread([path,files.physical],'Transformer');
[~,txtPR,rawPR] = local_xlsread([path,files.physical],'Resource');
[~,txtPE,rawPE] = local_xlsread([path,files.physical],'Set-up Efficiency');

try
    [~,~,rawPEfY] = local_xlsread([path,files.physical],'TransformerYear');
catch
    disp('No Transformer Year Efficiency data');
end

try
    [~,txtPD,~] = local_xlsread([path,files.physical],'Desalination');
    [~,~,rawPD_V] = local_xlsread([path,files.physical],'Desalination_D_V');
    [~,~,rawPD_H] = local_xlsread([path,files.physical],'Desalination_D_H');
catch
    disp('No desalination data');
end

[~,txtPH,rawPH] = local_xlsread([path,files.physical],'Hydro');
[numPS,txtPS,~] = local_xlsread([path,files.physical],'Storage');
[~,txtPL,rawPL] = local_xlsread([path,files.physical],'Load');
[~,~,rawPSC] = local_xlsread([path,files.physical],'Demand_Sectors');
[numPTr,txtPTr,~] = local_xlsread([path,files.physical],'Transmission');

try
    [~,~,rawPTr_L] = local_xlsread([path,files.physical],'AC_Losses');
catch
    disp('No AC losses data');
end

try
    [numPC,~,~] = local_xlsread([path,files.physical],'CO2');
catch
    disp('No CO2 emissions data');
end

%if setup.Mobility
[numPM,txtPM,~] = local_xlsread([path,files.physical],'Mobility');
%end
[numPY,~,~] = local_xlsread([path,files.physical],'Year');

%% Limits
[~,txtLT,rawLT] = local_xlsread([path,files.limits],'Transformer');
[~,txtLR,rawLR] = local_xlsread([path,files.limits],'Resource');
try
    [~,txtLD,rawLD] = local_xlsread([path,files.limits],'Desalination');
catch
    disp('No desalination data');
end
[~,txtLM,rawLM] = local_xlsread([path,files.limits],'Mobility');
[~,txtLH,rawLH] = local_xlsread([path,files.limits],'Hydro');
[~,txtLS,rawLS] = local_xlsread([path,files.limits],'Storage');
[~,txtLL,rawLL] = local_xlsread([path,files.limits],'Load');
[~,txtLTr,rawLTr] = local_xlsread([path,files.limits],'Transmission');

%% Geographical
[numG,txtG,~] = local_xlsread([path,files.geographical],'Node');
[~,~,rawGP] = local_xlsread([path,files.geographical],'Population');
[~,txtGT,rawGT] = local_xlsread([path,files.geographical],'Temperature');

%% Instalation
[~,~,rawTI_T] = local_xlsread([path,files.instalation],'Transformer');
[~,~,rawTI_R] = local_xlsread([path,files.instalation],'Resource');
[~,~,rawTI_M] = local_xlsread([path,files.instalation],'Mobility');
[~,~,rawTI_H] = local_xlsread([path,files.instalation],'Hydro');
[~,~,rawTI_D] = local_xlsread([path,files.instalation],'Desalination');
[~,~,rawTI_S] = local_xlsread([path,files.instalation],'Storage');
[~,~,rawTI_L] = local_xlsread([path,files.instalation],'Load');
[~,~,rawTI_Tr] = local_xlsread([path,files.instalation],'Transmission');

%% Financial
[numFT,~,~] = local_xlsread([path,files.financial],'Transformer');
[numFR,~,~] = local_xlsread([path,files.financial],'Resource');
[numFH,~,~] = local_xlsread([path,files.financial],'Hydro');
[numFS,~,~] = local_xlsread([path,files.financial],'Storage');
[numFL,~,~] = local_xlsread([path,files.financial],'Load');
[numFTr,~,~] = local_xlsread([path,files.financial],'Transmission');
[~,~,rawFSC_res] = local_xlsread([path,files.financial],'ElCostRES');
[~,~,rawFSC_com] = local_xlsread([path,files.financial],'ElCostCOM');
[~,~,rawFSC_ind] = local_xlsread([path,files.financial],'ElCostIND');
[numFW,~,~] = local_xlsread([path,files.financial],'WACC');
[numFC,~,~] = local_xlsread([path,files.financial],'CO2');
[numFM,~,~] = local_xlsread([path,files.financial],'Mobility');
[numFD,~,~] = local_xlsread([path,files.financial],'Desalination');

%if setup.Mobility
%% Financial - Mobility
[numMRoP,~,~] = local_xlsread([path,files.mobility],'Road p-km');
[numMRoT,~,~] = local_xlsread([path,files.mobility],'Road t-km');
[numMRaP,~,~] = local_xlsread([path,files.mobility],'Rail p-km');
[numMRaT,~,~] = local_xlsread([path,files.mobility],'Rail t-km');
[numMMP,~,~] = local_xlsread([path,files.mobility],'Marine p-km');
[numMMT,~,~] = local_xlsread([path,files.mobility],'Marine t-km');
[numMAP,~,~] = local_xlsread([path,files.mobility],'Aviation p-km');
[numMAT,~,~] = local_xlsread([path,files.mobility],'Aviation t-km');

[numMRL,~,~] = local_xlsread([path,files.mobility],'Road LDV Share');
[numMRW,~,~] = local_xlsread([path,files.mobility],'Road 23W Share');
[numMRB,~,~] = local_xlsread([path,files.mobility],'Road BUS Share');
[numMRM,~,~] = local_xlsread([path,files.mobility],'Road MDV Share ');
[numMRH,~,~] = local_xlsread([path,files.mobility],'Road HDV Share');

[numMRLI,~,~] = local_xlsread([path,files.mobility],'LDV Share ICE');
[numMRLB,~,~] = local_xlsread([path,files.mobility],'LDV Share BEV');
[numMRLP,~,~] = local_xlsread([path,files.mobility],'LDV Share PHEV');
[numMRLF,~,~] = local_xlsread([path,files.mobility],'LDV Share FCEV');

[numMRWI,~,~] = local_xlsread([path,files.mobility],'23W Share ICE');
[numMRWB,~,~] = local_xlsread([path,files.mobility],'23W Share BEV');
[numMRWP,~,~] = local_xlsread([path,files.mobility],'23W Share PHEV');
[numMRWF,~,~] = local_xlsread([path,files.mobility],'23W Share FCEV');

[numMRBI,~,~] = local_xlsread([path,files.mobility],'BUS Share ICE');
[numMRBB,~,~] = local_xlsread([path,files.mobility],'BUS Share BEV');
[numMRBP,~,~] = local_xlsread([path,files.mobility],'BUS Share PHEV');
[numMRBF,~,~] = local_xlsread([path,files.mobility],'BUS Share FCEV');

[numMRMI,~,~] = local_xlsread([path,files.mobility],'MDV Share ICE');
[numMRMB,~,~] = local_xlsread([path,files.mobility],'MDV Share BEV');
[numMRMP,~,~] = local_xlsread([path,files.mobility],'MDV Share PHEV');
[numMRMF,~,~] = local_xlsread([path,files.mobility],'MDV Share FCEV');

[numMRHI,~,~] = local_xlsread([path,files.mobility],'HDV Share ICE');
[numMRHB,~,~] = local_xlsread([path,files.mobility],'HDV Share BEV');
[numMRHP,~,~] = local_xlsread([path,files.mobility],'HDV Share PHEV');
[numMRHF,~,~] = local_xlsread([path,files.mobility],'HDV Share FCEV');

[numMRLp,~,~] = local_xlsread([path,files.mobility],'LDV Pass');
[numMRWp,~,~] = local_xlsread([path,files.mobility],'23W Pass');
[numMRBp,~,~] = local_xlsread([path,files.mobility],'BUS Pass');
[numMRMt,~,~] = local_xlsread([path,files.mobility],'MDV Tonne');
[numMRHt,~,~] = local_xlsread([path,files.mobility],'HDV Tonne');

[numMRLkm,~,~] = local_xlsread([path,files.mobility],'Annual km Per LDV');
[numMRWkm,~,~] = local_xlsread([path,files.mobility],'Annual km Per 23W');
[numMRBkm,~,~] = local_xlsread([path,files.mobility],'Annual km Per BUS');
[numMRMkm,~,~] = local_xlsread([path,files.mobility],'Annual km Per MDV');
[numMRHkm,~,~] = local_xlsread([path,files.mobility],'Annual km Per HDV');

[numMRLel,~,~] = local_xlsread([path,files.mobility],'LVD Share Of Elec for PHEV');
[numMRWel,~,~] = local_xlsread([path,files.mobility],'23W Share Of Elec for PHEV');
[numMRBel,~,~] = local_xlsread([path,files.mobility],'BUS Share Of Elec for PHEV');
[numMRMel,~,~] = local_xlsread([path,files.mobility],'MDV Share Of Elec for PHEV');
[numMRHel,~,~] = local_xlsread([path,files.mobility],'HDV Share Of Elec for PHEV');

[numMRPE,~,~] = local_xlsread([path,files.mobility],'Rail Pass Elec Share');
[numMRFE,~,~] = local_xlsread([path,files.mobility],'Rail Tonne Elec Share');
[numMRPF,~,~] = local_xlsread([path,files.mobility],'Rail Pass Liqu Share');
[numMRFF,~,~] = local_xlsread([path,files.mobility],'Rail Tonne Liqu Share');

[numMMPF,~,~] = local_xlsread([path,files.mobility],'Marine Pass Liqu Share');
[numMMPE,~,~] = local_xlsread([path,files.mobility],'Marine Pass Elec Share');
[numMMPH,~,~] = local_xlsread([path,files.mobility],'Marine Pass H2 Share');
[numMMPG,~,~] = local_xlsread([path,files.mobility],'Marine Pass LNG Share');
[numMMPA,~,~] = local_xlsread([path,files.mobility],'Marine Pass NH3 Share');
[numMMPM,~,~] = local_xlsread([path,files.mobility],'Marine Pass MeOH Share');

[numMMFF,~,~] = local_xlsread([path,files.mobility],'Marine Tonne Liqu Share');
[numMMFE,~,~] = local_xlsread([path,files.mobility],'Marine Tonne Elec Share');
[numMMFH,~,~] = local_xlsread([path,files.mobility],'Marine Tonne H2 Share');
[numMMFG,~,~] = local_xlsread([path,files.mobility],'Marine Tonne LNG Share');
[numMMFA,~,~] = local_xlsread([path,files.mobility],'Marine Tonne NH3 Share');
[numMMFM,~,~] = local_xlsread([path,files.mobility],'Marine Tonne MeOH Share');

[numMAPF,~,~] = local_xlsread([path,files.mobility],'Aviation Pass Liqu Share');
[numMAPE,~,~] = local_xlsread([path,files.mobility],'Aviation Pass Elec Share');
[numMAPH,~,~] = local_xlsread([path,files.mobility],'Aviation Pass H2 Share');

[numMAFF,~,~] = local_xlsread([path,files.mobility],'Aviation Tonne Liqu Share');
[numMAFE,~,~] = local_xlsread([path,files.mobility],'Aviation Tonne Elec Share');
[numMAFH,~,~] = local_xlsread([path,files.mobility],'Aviation Tonne H2 Share');

%% Get Year indeses (temporarry)
IndexYears=numFT(1,1:find(numFT(1,1:end) == max(numFT(1,1:end))));

zz=size(IndexYears,2); endzz=zz;
for ii=1:zz
    if isnan(IndexYears(1,ii))
        endzz=ii-1;
    end
end
IndexYears=IndexYears(1:endzz);
output.IndexYears=IndexYears;

% Loop for regions using local_xlsread
for ii = length(IndexYears):-1:1
    try
        [numFT_reg(:,:,ii)] = local_xlsread([path,files.financial],['TransformerRegions_' num2str(IndexYears(ii))]);
    catch
        [numFT_reg(:,:,ii)] = local_xlsread([path,files.financial],['TransformerRegions_' num2str(2015)]);
    end
end

for ii = length(IndexYears):-1:1
    try
        [numFR_reg(:,:,ii)] = local_xlsread([path,files.financial],['ResourceRegions_' num2str(IndexYears(ii))]);
    catch
        [numFR_reg(:,:,ii)] = local_xlsread([path,files.financial],['ResourceRegions_' num2str(2015)]);
    end
end

for ii = length(IndexYears):-1:1
    try
        [numFH_reg(:,:,ii)] = local_xlsread([path,files.financial],['HydroRegions_' num2str(IndexYears(ii))]);
    catch
        [numFH_reg(:,:,ii)] = local_xlsread([path,files.financial],['HydroRegions_' num2str(2015)]);
    end
end

for ii = length(IndexYears):-1:1
    try
        [numFS_reg(:,:,ii)] = local_xlsread([path,files.financial],['StorageRegions_' num2str(IndexYears(ii))]);
    catch
        [numFS_reg(:,:,ii)] = local_xlsread([path,files.financial],['StorageRegions_' num2str(2015)]);
    end
end

for ii = length(IndexYears):-1:1
    try
        [numFM_reg(:,:,ii)] = local_xlsread([path,files.financial],['MobilityRegions_' num2str(IndexYears(ii))]);
    catch
        [numFM_reg(:,:,ii)] = local_xlsread([path,files.financial],['MobilityRegions_' num2str(2015)]);
    end
end

for ii = length(IndexYears):-1:1
    try
        [numFD_reg(:,:,ii)] = local_xlsread([path,files.financial],['DesalinationRegions_' num2str(IndexYears(ii))]);
    catch
        [numFD_reg(:,:,ii)] = local_xlsread([path,files.financial],['DesalinationRegions_' num2str(2015)]);
    end
end
disp('Done Reading')

%%%%%%%%%%%%%% INDEXES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IndexYears -------------------------------------------------------------
IndexYears=numFT(1,1:find(numFT(1,1:end) == max(numFT(1,1:end))));
IndexYears=IndexYears(IndexYears<=setup.endYear);
zz=size(IndexYears,2); endzz=zz;
for ii=1:zz
    if isnan(IndexYears(1,ii))
        endzz=ii-1;
    end
end
IndexYears=IndexYears(1:endzz);
output.IndexYears=IndexYears;

% IndexNodes ------------------------------------------------------------
nr_of_nodes=1+max(numG(1:end,1))-min(numG(1:end,1));
IndexNodes=txtG(2:nr_of_nodes+1,2);
IndexNumNodes=numG(1:nr_of_nodes,1);
output.IndexNumNodes=IndexNumNodes;
output.IndexNodes=IndexNodes;

% IndexID's -------------------------------------------------------------
IndexIDT=txtPT(4:end,1);
zz=size(IndexIDT,1); endzz=zz;
for ii=1:zz
    if strcmp(IndexIDT,'')
        endzz=ii-1;
    end
end
IndexIDT=IndexIDT(1:endzz);
output.IndexIDT=IndexIDT;

IndexIDR=txtPR(2:end,1);
zz=size(IndexIDR,1); endzz=zz;
for ii=1:zz
    if strcmp(IndexIDR,'')
        endzz=ii-1;
    end
end
IndexIDR=IndexIDR(1:endzz);
output.IndexIDR=IndexIDR;

IndexIDE=txtPE(2:end,1);
zz=size(IndexIDE,1); endzz=zz;
for ii=1:zz
    if strcmp(IndexIDE,'')
        endzz=ii-1;
    end
end
IndexIDE=IndexIDE(1:endzz);
output.IndexIDE=IndexIDE;

IndexIDH=txtPH(2:end,1);
zz=size(IndexIDH,1); endzz=zz;
for ii=1:zz
    if strcmp(IndexIDH,'')
        endzz=ii-1;
    end
end
IndexIDH=IndexIDH(1:endzz);
output.IndexIDH=IndexIDH;

IndexIDS=txtPS(2:end,1);
zz=size(IndexIDS,1); endzz=zz;
for ii=1:zz
    if strcmp(IndexIDS,'')
        endzz=ii-1;
    end
end
IndexIDS=IndexIDS(1:endzz);
output.IndexIDS=IndexIDS;

IndexIDM=txtLM(3:end,1);
zz=size(IndexIDM,1); endzz=zz;
for ii=1:zz
    if strcmp(IndexIDM,'')
        endzz=ii-1;
    end
end
IndexIDM=IndexIDM(1:endzz);
output.IndexIDM=IndexIDM;

IndexIDD=txtPD(2:end,1);
zz=size(IndexIDD,1); endzz=zz;
for ii=1:zz
    if strcmp(IndexIDD,'')
        endzz=ii-1;
    end
end
IndexIDD=IndexIDD(1:endzz);
output.IndexIDD=IndexIDD;

IndexIDL=txtPL(2:end,1);
zz=size(IndexIDL,1); endzz=zz;
for ii=1:zz
    if strcmp(IndexIDL,'')
        endzz=ii-1;
    end
end
IndexIDL=IndexIDL(1:endzz);
output.IndexIDL=IndexIDL;

IndexIDTr=txtPTr(2:end,1);
zz=size(IndexIDTr,1); endzz=zz;
for ii=1:zz
    if strcmp(IndexIDTr,'')
        endzz=ii-1;
    end
end
IndexIDTr=IndexIDTr(1:endzz);
output.IndexIDTr=IndexIDTr;

IndexID=[IndexIDT; IndexIDS; IndexIDR; IndexIDH; IndexIDM; IndexIDD; IndexIDL; IndexIDTr];
output.IndexID=IndexID;


% IndexFuels ------------------------------------------------------------
IndexFuels=txtPT(1,3:13);
zz=size(IndexFuels,2); endzz=zz;
for ii=1:zz
    if strcmp(IndexFuels(ii),'')
        endzz=ii-1;
    end
end
IndexFuels=IndexFuels(1:endzz);
output.IndexFuels=IndexFuels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Year ------------------------------------------------------------------
output.year=numPY(1,1);

% Coords ----------------------------------------------------------------
output.Coords=numG(1:nr_of_nodes,4:6)';  %[longitude latitude timezone]

% Urbanisation level ----------------------------------------------------
output.Urbanisation = numG(1:nr_of_nodes,7)';
% Time zome level ----------------------------------------------------
output.timeZone = numG(1:nr_of_nodes,6)';
% Area ----------------------------------------------------
output.area = numG(1:nr_of_nodes,9)';
% Population ----------------------------------------------------
output.population = cell2mat(rawGP(2:1+nr_of_nodes,3:length(IndexYears)+2));

% WACC ------------------------------------------------------------------
output.WACC=numFW(1,1);
output.WACC_SC=numFW(1,2:4);
output.WACC_Spec=numFW(1,5);
% CO2 ------------------------------------------------------------------
output.fossilCO2Cost=numFC(2,1:end);
% Active? ---------------------------------------------------------------
output.Active=[~strcmp(txtLT(3:size(IndexIDT,1)+2,3),'-');...
    ~strcmp(txtLS(3:size(IndexIDS,1)+2,3),'-');...
    ~strcmp(txtLR(3:size(IndexIDR,1)+2,3),'-');...
    ~strcmp(txtLH(3:size(IndexIDH,1)+2,3),'-');...
    ~strcmp(txtLM(3:size(IndexIDM,1)+2,3),'-');...
    ~strcmp(txtLD(3:size(IndexIDD,1)+2,3),'-');...
    ~strcmp(txtLL(3:size(IndexIDL,1)+2,3),'-');...
    ~strcmp(txtLTr(3:size(IndexIDTr,1)+2,3),'-')];

output.Type = [(txtLT(3:size(IndexIDT,1)+2,3));...
    (txtLS(3:size(IndexIDS,1)+2,3));...
    (txtLR(3:size(IndexIDR,1)+2,3));...
    (txtLH(3:size(IndexIDH,1)+2,3));...
    (txtLM(3:size(IndexIDM,1)+2,3));...
    (txtLD(3:size(IndexIDD,1)+2,3));...
    (txtLL(3:size(IndexIDL,1)+2,3));...
    (txtLTr(3:size(IndexIDTr,1)+2,3))];

% SizeLimits ------------------------------------------------------------
SizeLimitsraw=[rawLT(3:size(IndexIDT,1)+2,4:3+nr_of_nodes);...
    rawLS(3:size(IndexIDS,1)+2,4:3+nr_of_nodes);...
    rawLR(3:size(IndexIDR,1)+2,4:3+nr_of_nodes);...
    rawLH(3:size(IndexIDH,1)+2,4:3+nr_of_nodes);...
    rawLM(3:size(IndexIDM,1)+2,4:3+nr_of_nodes);...
    rawLD(3:size(IndexIDD,1)+2,4:3+nr_of_nodes);...
    rawLL(3:size(IndexIDL,1)+2,4:3+nr_of_nodes);...
    rawLTr(3:size(IndexIDTr,1)+2,4:3+nr_of_nodes)];
if size(SizeLimitsraw,1)~=size(IndexID,1)
    warning('Number of Technologies and given limits on Technologies is different')
end
for node=1:nr_of_nodes
    for ii=1:size(SizeLimitsraw,1)
        if ischar(SizeLimitsraw{ii,node})
            [m0 , ~]=regexp(SizeLimitsraw{ii,node},'-','split');

            if strcmp(m0(2),'')
                SizeLimits(ii,:,node)=[str2num(m0{1}) Inf];
            else
                SizeLimits(ii,:,node)=[str2num(m0{1}) str2num(m0{2})];
            end
        else
            if ismissing(SizeLimitsraw{ii,node})
                SizeLimits(ii,:,node)=[0 Inf];
            else
                if SizeLimitsraw{ii,node}<0
                    SizeLimits(ii,:,node)=[0 -SizeLimitsraw{ii,node}];
                else
                    SizeLimits(ii,:,node)=[SizeLimitsraw{ii,node} SizeLimitsraw{ii,node}];
                end
            end
        end
    end
end
output.SizeLimits=SizeLimits;

% Capex, Opex -----------------------------------------------------------
Finraw=[numFT(2:size(IndexIDT,1)*5+1,1:size(IndexYears,2)+0);...
    numFS(2:size(IndexIDS,1)*5+1,1:size(IndexYears,2)+0);...
    numFR(2:size(IndexIDR,1)*5+1,1:size(IndexYears,2)+0);...
    numFH(2:size(IndexIDH,1)*5+1,1:size(IndexYears,2)+0);...
    numFM(2:size(IndexIDM,1)*5+1,1:size(IndexYears,2)+0);...
    numFD(2:size(IndexIDD,1)*5+1,1:size(IndexYears,2)+0);...
    numFL(2:size(IndexIDL,1)*5+1,1:size(IndexYears,2)+0);...
    numFTr(2:size(IndexIDTr,1)*5+1,1:size(IndexYears,2)+0)];

output.Capex=Finraw(1:5:size(Finraw,1)-4,:);
output.Opex_fix=Finraw(2:5:size(Finraw,1)-3,:);
output.Opex_var=Finraw(3:5:size(Finraw,1)-2,:);
output.Lifetime=Finraw(4:5:size(Finraw,1)-1,:);
output.Efficiency=Finraw(5:5:size(Finraw,1)-0,:);

for ii = length(IndexYears):-1:1
    Finraw2=[numFT_reg(1:size(IndexIDT,1)*3+0,1:size(IndexNodes,1)+0,ii);...
        numFS_reg(1:size(IndexIDS,1)*3+0,1:size(IndexNodes,1)+0,ii);...
        numFR_reg(1:size(IndexIDR,1)*3+0,1:size(IndexNodes,1)+0,ii);...
        numFH_reg(1:size(IndexIDH,1)*3+0,1:size(IndexNodes,1)+0,ii);...
        numFM_reg(1:size(IndexIDM,1)*3+0,1:size(IndexNodes,1)+0,ii);...
        numFD_reg(1:size(IndexIDD,1)*3+0,1:size(IndexNodes,1)+0,ii);...
        ones(size(IndexIDL,1)*3,size(IndexNodes,1));...
        ones(size(IndexIDTr,1)*3,size(IndexNodes,1))];

    output.Capex_reg(:,:,ii)=Finraw2(1:3:size(Finraw2,1)-2,:);
    output.Opex_fix_reg(:,:,ii)=Finraw2(2:3:size(Finraw2,1)-1,:);
    output.Opex_var_reg(:,:,ii)=Finraw2(3:3:size(Finraw2,1),:);
end


% Rampability -----------------------------------------------------------
output.RampTrans=numPT(:,12);

% Efficiency Transformers -----------------------------------------------
EtaTrans=numPT(:,1:11);
EtaTrans(~strcmp(txtPT(4:end,3:12),''))=-9999;
output.EtaTrans=EtaTrans;
output.EtaFormulaTrans=strrep(txtPT(4:end,3:12),',','.');

% Efficiency Storages ---------------------------------------------------
EtaStorage=numPS(:,3:4);
EtaStorage(~strcmp(txtPS(2:end,5:6),''))=-9999;
output.EtaStorage=EtaStorage;
output.EtaFormulaStorage=strrep(txtPS(2:end,5:6),',','.');

% Energy-Power Ratio Storages -------------------------------------------
output.EnPoRatioStorage=numPS(:,1:2);

% Self Discharge Rate Storages ------------------------------------------
SelfDischargeStorage=numPS(:,5);
SelfDischargeStorage(~strcmp(txtPS(2:end,7),''))=-9999;
output.SelfDischargeStorage=SelfDischargeStorage;
output.SelfDischargeFormulaStorage=strrep(txtPS(2:end,7),',','.');

% Loads -----------------------------------------------------------------
ltype=txtPL(2:end,3);
for ii=1:size(rawPL,1)-1
    for jj=1:size(rawPR,2)-3
        temp=rawPL{ii+1,jj+3};
        if ischar(temp)
            temp=str2double(temp);
        end
        numPL2(ii,jj)=temp;
    end
end

for ii=1:size(ltype,1)
    switch ltype{ii}
        case 'constant'
            ValueLoad(ii,1:nr_of_hours,1:2,:)=reshape(ones(nr_of_hours*2,1)*numPL2(ii,1:nr_of_nodes),nr_of_hours,2,[]);
        case 'daily'
            % FEHLT NOCH
        case 'monthly'
            % FEHLT NOCH
        case 'hourly'
            for nn=1:nr_of_nodes
                tag=rawPL{ii+1,nn+3};
                if exist([path filesep tag '.mat'])==2
                    % tag is a file name
                    temp=load([path filesep tag '.mat']);
                    if isstruct(temp)
                        ss=fieldnames(temp);
                        for iii=1:size(ss,1)
                            if numel(temp.(ss{iii}))==nr_of_hours
                                temp=temp.(ss{iii});
                            end
                        end
                    end
                else
                    if exist(['GetFromDatabase_',tag])
                        %tag must be a database name
                        temp=eval(['GetFromDatabase_',tag,'(output.Coords(:,',int2str(nn),'),output.year);']);
                    else
                        disp([tag,' is neither a valid file name nor does the database function exist'])
                        temp=zeros(nr_of_hours,1);
                    end
                end

                ValueLoad(ii,:,1,nn)=temp;
                ValueLoad(ii,:,2,nn)=temp;
                ValueLoadTag{ii,nn}=tag;
                %end
            end
        case 'hourlyA'
            for nn=1:nr_of_nodes
                tag=rawPL{ii+1,nn+3};
                if exist([path filesep tag '_' num2str(setup.startYear) '.mat'])==2
                    % tag is a file name
                    temp=load([path filesep tag '_' num2str(setup.startYear) '.mat']);
                    if isstruct(temp)
                        ss=fieldnames(temp);
                        for iii=1:size(ss)
                            disp([path filesep tag '_' num2str(setup.startYear) '.mat'])
                            if numel(temp.(ss{iii}))==nr_of_hours
                                temp=temp.(ss{iii});
                            end
                        end
                    end
                else
                    if exist(['GetFromDatabase_',tag])
                        %tag must be a database name
                        temp=eval(['GetFromDatabase_',tag,'(output.Coords(:,',int2str(nn),'),output.year);']);
                    else
                        disp([tag,' is neither a valid file name nor does the database function exist'])
                        temp=zeros(nr_of_hours,1);
                    end
                end

                ValueLoad(ii,:,1,nn)=temp;
                ValueLoad(ii,:,2,nn)=temp;
                ValueLoadTag{ii,nn}=tag;
            end
    end
end
output.ValueLoadTag=ValueLoadTag;
output.ValueLoad=ValueLoad;

% Resources -------------------------------------------------------------
rowTypes=txtPR(2:end,3);
ValueResourceTag = cell(size(rowTypes,1),nr_of_nodes);
for ii=1:size(rawPR,1)-1
    for jj=1:size(rawPR,2)-3
        temp=rawPR{ii+1,jj+3};
        if ischar(temp)
            temp=str2double(temp);
        end
        numPR2(ii,jj)=temp;
    end
end
for ii=1:size(rowTypes,1)
    switch rowTypes{ii}
        case 'total'
            ValueResource(ii,1:nr_of_hours,1,:)=zeros(nr_of_hours,1,nr_of_nodes);
            ValueResource(ii,1:nr_of_hours,2,:)=Inf(nr_of_hours,1,nr_of_nodes);
            ValueResourceTotal(ii,:)=(numPR2(ii,1:nr_of_nodes));
        case 'constant'
            ValueResource(ii,1:nr_of_hours,1:2,:)=reshape(ones(nr_of_hours*2,1)*numPR2(ii,1:nr_of_nodes),nr_of_hours,2,[]);
            ValueResourceTotal(ii,:)=(8760*numPR2(ii,1:nr_of_nodes));
        case 'infinite'
            ValueResource(ii,1:nr_of_hours,1,:)=zeros(nr_of_hours,1,nr_of_nodes);
            ValueResource(ii,1:nr_of_hours,2,:)=Inf(nr_of_hours,1,nr_of_nodes);
            ValueResourceTotal(ii,:)=1e14*ones(1,nr_of_nodes);
        case 'limited'
            ValueResource(ii,1:nr_of_hours,1:2,:)=reshape(ones(nr_of_hours*2,1)*numPR2(ii,1:nr_of_nodes),nr_of_hours,2,[]);
            ValueResource(ii,1:nr_of_hours,1,:)=zeros(nr_of_hours,1,nr_of_nodes);
            ValueResourceTotal(ii,:)=(8760*numPR2(ii,1:nr_of_nodes));
        case 'hourly'
            for nn=1:nr_of_nodes
                tag=rawPR{ii+1,nn+3};
                if exist([path filesep tag '.mat'],'file')==2
                    % tag is a file name
                    temp=load([path filesep tag '.mat']);
                    if isstruct(temp)
                        ss=fieldnames(temp);
                        for iii=1:size(ss,1)
                            if numel(temp.(ss{iii}))==nr_of_hours
                                temp=temp.(ss{iii});
                            end
                        end
                    end
                else
                    if exist(['GetFromDatabase_',tag])
                        %tag must be a database name
                        temp=eval(['GetFromDatabase_',tag,'(output.Coords(:,',int2str(nn),'),output.year);']);
                    else
                        disp([tag,' is neither a valid file name nor does the database function exist'])
                        temp=zeros(1,nr_of_hours);
                    end
                end
                if strcmp(IndexIDR(ii),'RWIN')||strcmp(IndexIDR(ii),'ROWI')||strcmp(IndexIDR(ii),'RWIO')
                    ValueResource(ii,:,1,nn)=temp./95.*98;
                    ValueResource(ii,:,2,nn)=temp./95.*98;
                    ValueResourceTotal(ii,nn)=sum(temp./95.*98);
                    output.WindCF = 0.98;
                else
                    ValueResource(ii,:,1,nn)=temp;
                    ValueResource(ii,:,2,nn)=temp;
                    ValueResourceTotal(ii,nn)=sum(temp);
                end
            end
        case 'hourly limited'
            for nn=1:nr_of_nodes
                tag=rawPR{ii+1,nn+3};
                if exist(tag)==2
                    % tag is a file name
                    temp=load([tag '.mat']);
                    if isstruct(temp)
                        ss=fieldnames(temp);
                        for iii=1:size(ss)
                            if numel(temp.(ss{iii}))==nr_of_hours
                                temp=temp.(ss{iii});
                            end
                        end
                    end
                else
                    if exist(['GetFromDatabase_',tag])
                        %tag must be a database name
                        temp=eval(['GetFromDatabase_',tag,'(output.Coords(:,',int2str(nn),'),output.year);']);
                    else
                        disp([tag,' is neither a valid file name nor does the database function exist'])
                        temp=zeros(1,nr_of_hours);
                    end
                end
                ValueResource(ii,:,1,nn)=zeros(nr_of_hours,1,nr_of_nodes);
                ValueResource(ii,:,2,nn)=temp;
                ValueResourceTotal(ii,nn)=sum(temp);
            end
        case 'hourlyA'
            for nn=1:nr_of_nodes
                tag=rawPR{ii+1,nn+3};
                if exist([path filesep tag '_' num2str(setup.startYear) '.mat'])==2
                    % tag is a file name
                    temp=load([path filesep tag '_' num2str(setup.startYear) '.mat']);
                    if isstruct(temp)
                        ss=fieldnames(temp);
                        for iii=1:size(ss)
                            if numel(temp.(ss{iii}))==nr_of_hours
                                temp=temp.(ss{iii});
                            end
                        end
                    end
                else
                    if exist(['GetFromDatabase_',tag])
                        %tag must be a database name
                        temp=eval(['GetFromDatabase_',tag,'(output.Coords(:,',int2str(nn),'),output.year);']);
                    else
                        disp([tag,' is neither a valid file name nor does the database function exist'])
                        temp=zeros(nr_of_hours,1);
                    end
                end
                ValueResource(ii,:,1,nn)=temp;
                ValueResource(ii,:,2,nn)=temp;
                ValueResourceTotal(ii,nn)=sum(temp);
                ValueResourceTag{ii,nn}=tag;
            end
    end
end
output.ValueResource=ValueResource;
output.ValueResourceTotal=ValueResourceTotal;
output.ValueResourceTag=ValueResourceTag;
% Set-up Efficiency -------------------------------------------------------------
for ii=1:size(rawPE,1)-1
    temp=rawPE{ii+1,3};
    if ischar(temp)
        temp=str2double(temp);
    end
    numPE2(ii)=temp;
end
output.FeedInEfficiencies=numPE2;
% Electricity transmission Efficiency -------------------------------------------------------------
output.TransmissionLinesEfficiencies=numPTr(1);
output.TransmissionConvertionEfficiencies=numPTr(2);

output.TransmissionEfficiencies=numPTr(1:length(output.IndexIDTr));
% Hydro -------------------------------------------------------------
rowTypes=txtPH(2:end,3);
ValueHydroDamTag = cell(size(rowTypes,1),nr_of_nodes);
for ii=1:size(rawPH,1)-1
    for jj=1:size(rawPH,2)-3
        temp=rawPH{ii+1,jj+3};
        if ischar(temp)
            temp=str2double(temp);
        end
        numPH2(ii,jj)=temp;
    end
end
for ii=1:size(rowTypes,1)
    switch rowTypes{ii}
        case 'constant'
            ValueHydroDam(ii,1:nr_of_hours,1:2,:)=reshape(ones(nr_of_hours*2,1)*numPH2(ii,1:nr_of_nodes),nr_of_hours,2,[]);
        case 'infinite'
            ValueHydroDam(ii,1:nr_of_hours,1,:)=zeros(nr_of_hours,1,nr_of_nodes);
            ValueHydroDam(ii,1:nr_of_hours,2,:)=Inf(nr_of_hours,1,nr_of_nodes);
        case 'limited'
            ValueHydroDam(ii,1:nr_of_hours,1:2,:)=reshape(ones(nr_of_hours*2,1)*numPH2(ii,1:nr_of_nodes),nr_of_hours,2,[]);
            ValueHydroDam(ii,1:nr_of_hours,1,:)=zeros(nr_of_hours,1,nr_of_nodes);
        case 'hourly'
            for nn=1:nr_of_nodes
                tag=rawPH{ii+1,nn+3};
                if exist([path filesep tag '.mat'])==2
                    % tag is a file name
                    temp=load([path filesep tag '.mat']);
                    if isstruct(temp)
                        ss=fieldnames(temp);
                        for iii=1:size(ss,1)
                            if numel(temp.(ss{iii}))==nr_of_hours
                                temp=temp.(ss{iii});
                            end
                        end
                    end
                else
                    if exist(['GetFromDatabase_',tag])
                        %tag must be a database name
                        temp=eval(['GetFromDatabase_',tag,'(output.Coords(:,',int2str(nn),'),output.year);']);
                    else
                        disp([tag,' is neither a valid file name nor does the database function exist'])
                        temp=zeros(1,nr_of_hours);
                    end
                end
                ValueHydroDam(ii,:,1,nn)=temp;
                ValueHydroDam(ii,:,2,nn)=temp;
            end
        case 'hourly limited'
            for nn=1:nr_of_nodes
                tag=rawPH{ii+1,nn+3};
                if exist(tag)==2
                    % tag is a file name
                    temp=load(tag);
                    if isstruct(temp)
                        ss=fieldnames(temp);
                        for iii=1:size(ss)
                            if numel(temp.(ss{iii}))==nr_of_hours
                                temp=temp.(ss{iii});
                            end
                        end
                    end
                else
                    if exist(['GetFromDatabase_',tag])
                        %tag must be a database name
                        temp=eval(['GetFromDatabase_',tag,'(output.Coords(:,',int2str(nn),'),output.year);']);
                    else
                        disp([tag,' is neither a valid file name nor does the database function exist'])
                        temp=zeros(1,nr_of_hours);
                    end
                end
                ValueHydroDam(ii,:,1,nn)=zeros(nr_of_hours,1,nr_of_nodes);
                ValueHydroDam(ii,:,2,nn)=temp;
            end
        case 'hourlyA'
            for nn=1:nr_of_nodes
                tag=rawPH{ii+1,nn+3};
                if exist([path filesep tag '_' num2str(setup.startYear) '.mat'])==2
                    % tag is a file name
                    temp=load([path filesep tag '_' num2str(setup.startYear) '.mat']);
                    if isstruct(temp)
                        ss=fieldnames(temp);
                        for iii=1:size(ss)
                            if numel(temp.(ss{iii}))==nr_of_hours
                                temp=temp.(ss{iii});
                            end
                        end
                    end
                else
                    if exist(['GetFromDatabase_',tag])
                        %tag must be a database name
                        temp=eval(['GetFromDatabase_',tag,'(output.Coords(:,',int2str(nn),'),output.year);']);
                    else
                        disp([tag,' is neither a valid file name nor does the database function exist'])
                        temp=zeros(nr_of_hours,1);
                    end
                end
                ValueHydroDam(ii,:,1,nn)=temp;
                ValueHydroDam(ii,:,2,nn)=temp;
                ValueHydroDamTag{ii,nn}=tag;
            end
    end
end
output.ValueHydroDam = ValueHydroDam;
output.ValueHydroDamTag = ValueHydroDamTag;

%% Temperature Data Processing ---------------------------------------------
% 1. Extract Metadata and Data
rowNames = txtGT(2:4, 1);
rowTypes = txtGT(2:4, 2);

% 2. Parse Numeric Data (Vectorized)
rawDataBlock = rawGT(2:end, 3:end);

% 3. Pre-allocation
numRows = size(rowTypes, 1);
Temperature = zeros(numRows, nr_of_hours, nr_of_nodes);

% 4. Main Processing Loop
for i = 1:numRows
    currentType = rowTypes{i};
    currentName = rowNames{i};
    
    switch lower(currentType) % Handle case-insensitivity
        
        case {'constant', 'limited'}
            % Spread the single value across all hours
            tempData = repmat(rawDataBlock(i, 1:nr_of_nodes), nr_of_hours, 1);
            Temperature(i, :, :) = cell2mat(tempData);
        case 'infinite'
            Temperature(i, :, :) = Inf(nr_of_hours, nr_of_nodes);
        case 'hourly'
            for n = 1:nr_of_nodes
                tag = rawDataBlock{i, n}; % Get filename or DB tag
                fullFilePath = fullfile(path, [tag '.mat']);
                
                if exist(fullFilePath, 'file') == 2
                    % Load from MAT file
                    tempVec = loadVariableFromFile(fullFilePath, nr_of_hours);
                elseif exist(['GetFromDatabase_' tag], 'file')
                    % Load from Database Function (Replaces eval)
                    dbFunc = str2func(['GetFromDatabase_' tag]);
                    % Assuming output.Coords and output.year exist in workspace
                    tempVec = dbFunc(output.Coords(:, n), output.year);
                else
                    warning('%s is neither a valid file nor a database function.', tag);
                    tempVec = zeros(1, nr_of_hours);
                end
                Temperature(i, :, n) = tempVec;
            end
            
        case 'hourly limited'
            for n = 1:nr_of_nodes
                tag = rawDataBlock{i, n};
                if exist(tag, 'file') == 2
                    tempVec = loadVariableFromFile(tag, nr_of_hours);
                else
                    tempVec = zeros(1, nr_of_hours); % Default fallback
                end
                Temperature(i, :, n) = tempVec;
            end
    end
    output.(currentName) = shiftdim(Temperature(i, :, :), 1);
end

function data = loadVariableFromFile(filename, targetSize)
    % Helper to load a .mat file and find the variable matching the target size
    loadedStruct = load(filename);
    fields = fieldnames(loadedStruct);
    data = [];
    % Search for the field that matches the number of hours
    for k = 1:numel(fields)
        if numel(loadedStruct.(fields{k})) == targetSize
            data = loadedStruct.(fields{k});
            return;
        end
    end
    % If no matching field is found, return empty or handle error
    if isempty(data)
       error('No variable of length %d found in %s', targetSize, filename);
    end
end


%% Instalation -----------------------------------------------------------------
Instraw=[rawTI_T(2:size(IndexIDT,1)*1+1,1:3);...
    rawTI_S(2:size(IndexIDS,1)*1+1,1:3);...
    rawTI_R(2:size(IndexIDR,1)*1+1,1:3);...
    rawTI_H(2:size(IndexIDH,1)*1+1,1:3);...
    rawTI_M(2:size(IndexIDM,1)*1+1,1:3);...
    rawTI_D(2:size(IndexIDD,1)*1+1,1:3);...
    rawTI_L(2:size(IndexIDL,1)*1+1,1:3);...
    rawTI_Tr(2:size(IndexIDTr,1)*1+1,1:3)];
installations = zeros(nr_of_nodes,length(IndexYears), length(Instraw));
for ii=1:size(Instraw,1)
    if ~strcmp(Instraw(ii,3),'-')
        [~,~,t3] = local_xlsread([path filesep Instraw{ii,3}]);
        try
            installations(:,:,ii) = cell2mat(t3(2:2+nr_of_nodes-1,3:3+length(IndexYears)-1));
        catch
            disp([Instraw{ii,3} ': data for last years is missing'])
            installations(:,:,ii) = [cell2mat(t3(2:2+nr_of_nodes-1,3:end)) repmat(cell2mat(t3(2:2+nr_of_nodes-1,end)),1,length(IndexYears)-size(t3(2:2+nr_of_nodes-1,3:end),2))];
        end
    else
        installations(:,:,ii) = zeros(nr_of_nodes,length(IndexYears));
    end
end
installations(isnan(installations))=0;
output.Instalations=installations;

%% Specific efficiencies
i = 2;
output.EffCO2Scr_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffCO2Scr_He = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffMETH_CO = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;

output.EffFT_CO = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffFT_He = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffFT_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffLNG_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffLH2_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffGCS_CO = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffHCS_CO = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;


output.EffMeO_CO = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffMeO_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffMeO_He_out = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffDME_CO = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffDME_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffDME_He_out = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffNH3_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffNH3_He_out = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffHyS_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffMET_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffMET_He_out = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffWEL_He_out = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffICB_Lime = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffICB_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffICB_He = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffICI_Lime = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffICI_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffICI_He = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISB_Ore = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISB_Scr = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISB_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISB_He = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISB_HC = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISH_Ore = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISH_Scr = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISH_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISH_He = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISH_CC = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISH_Hy = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISR_Scr = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISR_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISR_He = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISR_CC = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISE_Ore = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISE_Scr = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISE_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISE_He = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffISE_CC = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffIAA_Ba = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffIAA_Lime = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffIAA_So = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffIAA_He = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffIAA_HeOut = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffIAM_Aa = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffIAM_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffIAM_HeOut = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffIAR_Al = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffIAR_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffIAR_He = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);i= i+1;
output.EffIPP_Wo = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffIPP_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffIPP_He = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffSMC_CO = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffBCS_CO = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffCBC_CO = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffCWC_CO = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffFTB_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i = i+1;
output.EffFTM_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffPSC_CO = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffPSC_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffPSC_He = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffPSP_CO = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffPSP_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffPSP_He = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRLS_input = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRLS_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRSS_input = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRSS_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;

output.EffRMI_input = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRMI_Wa = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRMI_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRMO_input = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRMO_Wa = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRMO_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;

output.EffRME_input = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRME_IW = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRME_MR = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRME_CaCO3out = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRME_MgCO3out = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRME_SiO2out = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRME_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRME_He = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;

output.EffRSi_input = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRSi_Hy = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRSi_SiO2 = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRSi_SiCout = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRSi_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRSi_He = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;

output.EffRAD_Wa = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRAD_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;

output.EffREW_MR = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffREW_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;

output.EffRBC_Bio = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRBC_Wa = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRBC_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;

output.EffRGE_El = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1); i= i+1;
output.EffRGE_He = repmat(cell2mat(rawPEfY(i,4:length(IndexYears)+3)),nr_of_nodes,1);

%% CO2 emissions numbers
i=1;
output.OilEmissionsMain = numPC(i,1);
output.OilEmissionsAdditional = numPC(i,2); i=i+1;
output.GasEmissionsMain = numPC(i,1);
output.GasEmissionsAdditional = numPC(i,2); i=i+1;
output.CoalEmissionsMain = numPC(i,1);
output.CoalEmissionsAdditional = numPC(i,2); i=i+1;
output.BiomassEmissionsMain = numPC(i,1);
output.BiomassEmissionsAdditional = numPC(i,2); i=i+1;
output.WastesEmissionsMain = numPC(i,1);
output.WastesEmissionsAdditional = numPC(i,2); i=i+1;


output.TICB_EmissionsMain = numPC(i,1);
output.TICB_EmissionsAdditional = numPC(i,2); i=i+1;
output.TICI_EmissionsMain = numPC(i,1);
output.TICI_EmissionsAdditional = numPC(i,2); i=i+1;
output.TISB_EmissionsMain = numPC(i,1);
output.TISB_EmissionsAdditional = numPC(i,2); i=i+1;
output.TISH_EmissionsMain = numPC(i,1);
output.TISH_EmissionsAdditional = numPC(i,2); i=i+1;
output.TISR_EmissionsMain = numPC(i,1);
output.TISR_EmissionsAdditional = numPC(i,2); i=i+1;
output.TISE_EmissionsMain = numPC(i,1);
output.TISE_EmissionsAdditional = numPC(i,2); i=i+1;
output.TIAA_EmissionsMain = numPC(i,1);
output.TIAA_EmissionsAdditional = numPC(i,2); i=i+1;
output.TIAM_EmissionsMain = numPC(i,1);
output.TIAM_EmissionsAdditional = numPC(i,2); i=i+1;
output.TIAR_EmissionsMain = numPC(i,1);
output.TIAR_EmissionsAdditional = numPC(i,2); i=i+1;
output.TIPP_EmissionsMain = numPC(i,1);
output.TIPP_EmissionsAdditional = numPC(i,2);

%% AC losses -------------------------------------
AC_losses = cell2mat(rawPTr_L(2:2+nr_of_nodes-1,3:3+length(IndexYears)-1));

output.AC_losses = AC_losses;

output.GridMax = 10^20;
%% Desalination vertical pumping
output.DesalinationParams.regionPumpLong = cell2mat(rawPD_V(2:1+nr_of_nodes,3:length(IndexYears)+2));

%% Desalination horizontal pumping
output.DesalinationParams.regionPumpUp = cell2mat(rawPD_H(2:1+nr_of_nodes,3:length(IndexYears)+2));

%% Sectors consumption
output.SectorsCons = cell2mat(rawPSC(2:1+nr_of_nodes,3:3+2));

%% Electricity costs by regions,years
output.ElCostRES = cell2mat(rawFSC_res(2:1+nr_of_nodes,2:length(IndexYears)+1));
output.ElCostCOM = cell2mat(rawFSC_com(2:1+nr_of_nodes,2:length(IndexYears)+1));
output.ElCostIND = cell2mat(rawFSC_ind(2:1+nr_of_nodes,2:length(IndexYears)+1));

%if setup.Mobility
% Demands
output.Mobility.Demand(1,1:nr_of_nodes,1:length(IndexYears)) = numMRoP(2:1+nr_of_nodes,3:length(IndexYears)+2).*numMRL(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in p*km
output.Mobility.Demand(2,1:nr_of_nodes,1:length(IndexYears)) = numMRoP(2:1+nr_of_nodes,3:length(IndexYears)+2).*numMRW(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in p*km
output.Mobility.Demand(3,1:nr_of_nodes,1:length(IndexYears)) = numMRoP(2:1+nr_of_nodes,3:length(IndexYears)+2).*numMRB(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in p*km
output.Mobility.Demand(4,1:nr_of_nodes,1:length(IndexYears)) = numMRoT(2:1+nr_of_nodes,3:length(IndexYears)+2).*numMRM(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in t*km
output.Mobility.Demand(5,1:nr_of_nodes,1:length(IndexYears)) = numMRoT(2:1+nr_of_nodes,3:length(IndexYears)+2).*numMRH(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in t*km
output.Mobility.Demand(6,1:nr_of_nodes,1:length(IndexYears)) = numMRaP(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in p*km
output.Mobility.Demand(7,1:nr_of_nodes,1:length(IndexYears)) = numMRaT(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in t*km
output.Mobility.Demand(8,1:nr_of_nodes,1:length(IndexYears)) = numMMP(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in p*km
output.Mobility.Demand(9,1:nr_of_nodes,1:length(IndexYears)) = numMMT(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in t*km
output.Mobility.Demand(10,1:nr_of_nodes,1:length(IndexYears)) = numMAP(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in p*km
output.Mobility.Demand(11,1:nr_of_nodes,1:length(IndexYears)) = numMAT(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in t*km

output.Instalations(:,:,ismember(IndexID,'LMLL')) = numMRoP(2:1+nr_of_nodes,3:length(IndexYears)+2).*numMRL(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in p*km
output.Instalations(:,:,ismember(IndexID,'LMLW')) = numMRoP(2:1+nr_of_nodes,3:length(IndexYears)+2).*numMRW(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in p*km
output.Instalations(:,:,ismember(IndexID,'LMLB')) = numMRoP(2:1+nr_of_nodes,3:length(IndexYears)+2).*numMRB(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in p*km
output.Instalations(:,:,ismember(IndexID,'LMLM')) = numMRoT(2:1+nr_of_nodes,3:length(IndexYears)+2).*numMRM(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in t*km
output.Instalations(:,:,ismember(IndexID,'LMLH')) = numMRoT(2:1+nr_of_nodes,3:length(IndexYears)+2).*numMRH(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in t*km
output.Instalations(:,:,ismember(IndexID,'LMRP')) = numMRaP(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in p*km
output.Instalations(:,:,ismember(IndexID,'LMRF')) = numMRaT(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in t*km
output.Instalations(:,:,ismember(IndexID,'LMMP')) = numMMP(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in p*km
output.Instalations(:,:,ismember(IndexID,'LMMF')) = numMMT(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in t*km
output.Instalations(:,:,ismember(IndexID,'LMAP')) = numMAP(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in p*km
output.Instalations(:,:,ismember(IndexID,'LMAF')) = numMAT(2:1+nr_of_nodes,3:length(IndexYears)+2)*10^6; %in t*km

% tonne/passengers per car
output.Mobility.LDV_Pass(1:nr_of_nodes,1:length(IndexYears)) = numMRLp(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.W23_Pass(1:nr_of_nodes,1:length(IndexYears)) = numMRWp(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.BUS_Pass(1:nr_of_nodes,1:length(IndexYears)) = numMRBp(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.MDV_Tonne(1:nr_of_nodes,1:length(IndexYears)) = numMRMt(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.HDV_Tonne(1:nr_of_nodes,1:length(IndexYears)) = numMRHt(2:1+nr_of_nodes,3:length(IndexYears)+2);
% km per car
output.Mobility.LDV_km(1:nr_of_nodes,1:length(IndexYears)) = numMRLkm(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.W23_km(1:nr_of_nodes,1:length(IndexYears)) = numMRWkm(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.BUS_km(1:nr_of_nodes,1:length(IndexYears)) = numMRBkm(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.MDV_km(1:nr_of_nodes,1:length(IndexYears)) = numMRMkm(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.HDV_km(1:nr_of_nodes,1:length(IndexYears)) = numMRHkm(2:1+nr_of_nodes,3:length(IndexYears)+2);
% PHEV el share
output.Mobility.LDV_PHEV_el(1:nr_of_nodes,1:length(IndexYears)) = numMRLel(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.W23_PHEV_el(1:nr_of_nodes,1:length(IndexYears)) = numMRWel(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.BUS_PHEV_el(1:nr_of_nodes,1:length(IndexYears)) = numMRBel(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.MDV_PHEV_el(1:nr_of_nodes,1:length(IndexYears)) = numMRMel(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.HDV_PHEV_el(1:nr_of_nodes,1:length(IndexYears)) = numMRHel(2:1+nr_of_nodes,3:length(IndexYears)+2);

% Road transp shares
output.Mobility.LDV_Shares(1,1:nr_of_nodes,1:length(IndexYears)) = numMRLI(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.LDV_Shares(2,1:nr_of_nodes,1:length(IndexYears)) = numMRLB(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.LDV_Shares(3,1:nr_of_nodes,1:length(IndexYears)) = numMRLF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.LDV_Shares(4,1:nr_of_nodes,1:length(IndexYears)) = numMRLP(2:1+nr_of_nodes,3:length(IndexYears)+2);

output.Mobility.W23_Shares(1,1:nr_of_nodes,1:length(IndexYears)) = numMRWI(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.W23_Shares(2,1:nr_of_nodes,1:length(IndexYears)) = numMRWB(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.W23_Shares(3,1:nr_of_nodes,1:length(IndexYears)) = numMRWF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.W23_Shares(4,1:nr_of_nodes,1:length(IndexYears)) = numMRWP(2:1+nr_of_nodes,3:length(IndexYears)+2);

output.Mobility.BUS_Shares(1,1:nr_of_nodes,1:length(IndexYears)) = numMRBI(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.BUS_Shares(2,1:nr_of_nodes,1:length(IndexYears)) = numMRBB(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.BUS_Shares(3,1:nr_of_nodes,1:length(IndexYears)) = numMRBF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.BUS_Shares(4,1:nr_of_nodes,1:length(IndexYears)) = numMRBP(2:1+nr_of_nodes,3:length(IndexYears)+2);

output.Mobility.MDV_Shares(1,1:nr_of_nodes,1:length(IndexYears)) = numMRMI(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.MDV_Shares(2,1:nr_of_nodes,1:length(IndexYears)) = numMRMB(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.MDV_Shares(3,1:nr_of_nodes,1:length(IndexYears)) = numMRMF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.MDV_Shares(4,1:nr_of_nodes,1:length(IndexYears)) = numMRMP(2:1+nr_of_nodes,3:length(IndexYears)+2);

output.Mobility.HDV_Shares(1,1:nr_of_nodes,1:length(IndexYears)) = numMRHI(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.HDV_Shares(2,1:nr_of_nodes,1:length(IndexYears)) = numMRHB(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.HDV_Shares(3,1:nr_of_nodes,1:length(IndexYears)) = numMRHF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.HDV_Shares(4,1:nr_of_nodes,1:length(IndexYears)) = numMRHP(2:1+nr_of_nodes,3:length(IndexYears)+2);
% Rail shares
output.Mobility.Rail_Shares(1,1:nr_of_nodes,1:length(IndexYears)) = numMRPF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Rail_Shares(2,1:nr_of_nodes,1:length(IndexYears)) = numMRPE(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Rail_Shares(3,1:nr_of_nodes,1:length(IndexYears)) = numMRFF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Rail_Shares(4,1:nr_of_nodes,1:length(IndexYears)) = numMRFE(2:1+nr_of_nodes,3:length(IndexYears)+2);
% Marine shares
output.Mobility.Marine_Shares(1,1:nr_of_nodes,1:length(IndexYears)) = numMMPF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Marine_Shares(2,1:nr_of_nodes,1:length(IndexYears)) = numMMPE(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Marine_Shares(3,1:nr_of_nodes,1:length(IndexYears)) = numMMPH(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Marine_Shares(4,1:nr_of_nodes,1:length(IndexYears)) = numMMPG(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Marine_Shares(5,1:nr_of_nodes,1:length(IndexYears)) = numMMPA(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Marine_Shares(6,1:nr_of_nodes,1:length(IndexYears)) = numMMPM(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Marine_Shares(7,1:nr_of_nodes,1:length(IndexYears)) = numMMFF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Marine_Shares(8,1:nr_of_nodes,1:length(IndexYears)) = numMMFE(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Marine_Shares(9,1:nr_of_nodes,1:length(IndexYears)) = numMMFH(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Marine_Shares(10,1:nr_of_nodes,1:length(IndexYears)) = numMMFG(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Marine_Shares(11,1:nr_of_nodes,1:length(IndexYears)) = numMMFA(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Marine_Shares(12,1:nr_of_nodes,1:length(IndexYears)) = numMMFM(2:1+nr_of_nodes,3:length(IndexYears)+2);
% Avia shares
output.Mobility.Aviation_Shares(1,1:nr_of_nodes,1:length(IndexYears)) = numMAPF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Aviation_Shares(2,1:nr_of_nodes,1:length(IndexYears)) = numMAPE(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Aviation_Shares(3,1:nr_of_nodes,1:length(IndexYears)) = numMAPH(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Aviation_Shares(4,1:nr_of_nodes,1:length(IndexYears)) = numMAFF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Aviation_Shares(5,1:nr_of_nodes,1:length(IndexYears)) = numMAFE(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Aviation_Shares(6,1:nr_of_nodes,1:length(IndexYears)) = numMAFH(2:1+nr_of_nodes,3:length(IndexYears)+2);

% All Shares
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRLI')) = numMRLI(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRLB')) = numMRLB(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRLF')) = numMRLF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRLP')) = numMRLP(2:1+nr_of_nodes,3:length(IndexYears)+2);

output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRWI')) = numMRWI(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRWB')) = numMRWB(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRWF')) = numMRWF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRWP')) = numMRWP(2:1+nr_of_nodes,3:length(IndexYears)+2);

output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRBI')) = numMRBI(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRBB')) = numMRBB(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRBF')) = numMRBF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRBP')) = numMRBP(2:1+nr_of_nodes,3:length(IndexYears)+2);

output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRMI')) = numMRMI(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRMB')) = numMRMB(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRMF')) = numMRMF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRMP')) = numMRMP(2:1+nr_of_nodes,3:length(IndexYears)+2);

output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRHI')) = numMRHI(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRHB')) = numMRHB(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRHF')) = numMRHF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRHP')) = numMRHP(2:1+nr_of_nodes,3:length(IndexYears)+2);
% Rail shares
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRPF')) = numMRPF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRPE')) = numMRPE(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRFF')) = numMRFF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MRFE')) = numMRFE(2:1+nr_of_nodes,3:length(IndexYears)+2);
%Marine shares
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MMPF')) = numMMPF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MMPE')) = numMMPE(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MMPH')) = numMMPH(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MMPG')) = numMMPG(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MMPA')) = numMMPA(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MMPM')) = numMMPM(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MMFF')) = numMMFF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MMFE')) = numMMFE(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MMFH')) = numMMFH(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MMFG')) = numMMFG(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MMFA')) = numMMFA(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MMFM')) = numMMFM(2:1+nr_of_nodes,3:length(IndexYears)+2);
% Avia shares
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MAPF')) = numMAPF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MAPE')) = numMAPE(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MAPH')) = numMAPH(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MAFF')) = numMAFF(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MAFE')) = numMAFE(2:1+nr_of_nodes,3:length(IndexYears)+2);
output.Mobility.Shares(1:nr_of_nodes,1:length(IndexYears),ismember(IndexID,'MAFH')) = numMAFH(2:1+nr_of_nodes,3:length(IndexYears)+2);

% Consumptions
output.Mobility.Cons_Names = txtPM(2:length(IndexIDM)+5+1,1);
output.Mobility.Cons_Units = txtPM(2:length(IndexIDM)+5+1,3);
output.Mobility.Cons_Type = txtPM(2:length(IndexIDM)+5+1,5);
output.Mobility.Lifetime = numPM(2:length(IndexIDM)+5+1,1);
output.Mobility.Cons_Values = numPM(2:length(IndexIDM)+5+1,3:length(IndexYears)+2);

[numRC,~,~] = local_xlsread([path,setup.files.ramping],'Additional financial');
output.rampingCost = numRC;
end

% -------------------------------------------------------------------------
% LOCAL HELPER FUNCTION TO REPLACE XLSREAD
% -------------------------------------------------------------------------
function [num, txt, raw] = local_xlsread(filename, sheetname)
    % Wraps readcell to mimic the behavior of xlsread:
    % 1. Returns raw cell array (raw)
    % 2. Returns numeric matrix with NaNs for non-numeric (num)
    % 3. Returns text cell array with empty strings for numbers (txt)
    % 4. Attempts to trim the 'num' output to match typical xlsread behavior 
    %    (stripping outer empty/text headers to find the data block)

    if nargin < 2 || isempty(sheetname)
        try
            sheets = sheetnames(filename);
        catch ME
            error('Error obtaining sheet names from %s: %s', filename, ME.message);
        end
        if isempty(sheets)
            error('No sheets found in %s.', filename);
        end
        sheetname = sheets(1); % first sheet (string)
    end
        
    try
        % Read all data using readcell
        % 'TreatAsMissing' can be used if needed, but default is usually fine
        raw = readcell(filename, 'Sheet', sheetname);
    catch ME
        % Handle errors (e.g., sheet not found)
        error('Error reading %s in %s: %s', sheetname, filename, ME.message);
    end

    [rows, cols] = size(raw);
    
    % Initialize outputs
    num_full = nan(rows, cols);
    txt = cell(rows, cols);
    txt(:) = {''}; % Initialize with empty char vectors
    
    % Iterate once to populate num and txt
    for r = 1:rows
        for c = 1:cols
            val = raw{r,c};
            
            % Check for missing/empty
            if ismissing(val)
                continue; 
            end
            
            % Populate NUM
            if isnumeric(val)
                num_full(r,c) = double(val);
            end
            
            % Populate TXT
            % xlsread returns text cell array where numbers are empty strings
            if ischar(val) || isstring(val)
                txt{r,c} = char(val);
            end
        end
    end
    
    % Mimic xlsread trimming for 'num'
    % xlsread typically strips leading/trailing rows/cols that are non-numeric
    row_has_data = any(~isnan(num_full), 2);
    col_has_data = any(~isnan(num_full), 1);
    
    if any(row_has_data) && any(col_has_data)
        first_row = find(row_has_data, 1, 'first');
        last_row = find(row_has_data, 1, 'last');
        first_col = find(col_has_data, 1, 'first');
        last_col = find(col_has_data, 1, 'last');
        
        num = num_full(first_row:last_row, first_col:last_col);
    else
        num = [];
    end
    
    % Note: xlsread 'txt' output behavior is often the same size as 'raw' 
    % or slightly trimmed, but keeping it same size as 'raw' (with empty strings)
    % is safer for preserving indices in mixed-data workflows.
end
