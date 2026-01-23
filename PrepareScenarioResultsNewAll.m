function calc = PrepareScenarioResultsNewAll(setup,costYear,baseData,Reg1,Reg2,type,varargin)
%
% FUNCTION prepareScenarioResultsNewAll(setup, costYear, baseData, Reg1, Reg2, type, varargin)
%
% Prepares scenario result data for all regions and technologies.
%
% This function processes data for a given cost year and scenario type, combining information
% from two regions and the base dataset. Optional input can include region groupings or filters.
%
% INPUT:
%            setup:      Structure that contains all necessary settings and data for processing.
%            costYear:   Year to which all costs are adjusted.
%            baseData:   Structure with base scenario data.
%            Reg1:       Name or index of the first region.
%            Reg2:       Name or index of the second region.
%            type:       ADD!!!!!!!!!!!!!!!!!!!!
%            varargin:   ADD!!!!!!!!!!!!!!!!!!!!
%
% OUTPUT:
%            calc:       ADD!!!!!!!!!!!!!!!!!!!!
%
%Dmitrii Bogdanov
%last change 23.07.2025


try
    setup.Big;
catch
    setup.Big = 0;
end

Regs = setup.Regions;

% last change 31.10.2022
name = num2str(Reg1);
nameAdd = [];
scriptVersion = 5.0; % historical !!!
if nargin == 7
    nameAdd = ['_' varargin{1}];
    name = varargin{1};
end

rootDir = setup.rootDir;
projType = setup.projType;
OvernightFlag = setup.OvernightFlag;
SCFlag = setup.SC.Flag;
GasFlag = setup.GasFlag;
DesalinationFlag = setup.DesalinationFlag;
Capex_Rehab = setup. Capex_Rehab;
Lifetime_Rehab = setup.Lifetime_Rehab;
Only = setup.Only;
%% Get project name


try
    setup.MacroMc;
catch
    setup.MacroMc = 0;
end

try
    setup.MacroReg;
catch
    setup.MacroReg = 0;
end


pName = [projType];

if setup.MacroMc
    pName = [pName '_Mc'];
end

if setup.MacroReg
    pName = [pName '_R'];
end

if OvernightFlag
    pName = [pName '_Overnight'];
end

if SCFlag
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


if GasFlag
    pName = [pName '_GAS'];
end

if DesalinationFlag
    pName = [pName '_DES'];
end

if Only.Flag
    pName = [pName '_Only'];
end

pName = [pName '_' num2str(costYear)];


output.name = pName;

%% Load data
% results file
try
    load([rootDir filesep 'projects' filesep pName filesep 'output' filesep 'results_'  name]);
catch
    load([rootDir filesep 'projects' filesep pName filesep 'output' filesep 'results']);
end

% simulation input for the scenario
try
    load([rootDir filesep 'projects' filesep pName filesep 'input-data' filesep 'simulation-input_' pName '_' name]);

catch
    load([rootDir filesep 'projects' filesep pName filesep 'input-data' filesep 'simulation-input_' pName]);

end
% original input data
if setup.MacroMc
    baseData = load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep 'Macro' filesep 'simulation-input_' 'Base' '.mat']);
else
    if setup.Big ==1
        baseData = load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep 'simulation-input_' 'Base' nameAdd '.mat']);
    else
        baseData = load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep 'simulation-input_' 'Base' '.mat']);
    end
end


%% Temp changes for versions compatibility

try results.HY_TRSi;
catch
    results.HY_TRSi=zeros(size(results.RWIN_EL));
    warning(['HY_TRSi is fixed!'])
end

try results.TPSP_GA;
catch
    results.TPSP_GA=zeros(size(results.RWIN_EL));
    warning(['TPSP_GA is fixed!'])
end

try results.TBCS_CO;
catch
    results.TBCS_CO=zeros(size(results.RWIN_EL));
    warning(['TBCS_CO is fixed!'])
end

try results.TCBC_CO;
catch
    results.TCBC_CO=zeros(size(results.RWIN_EL));
    warning(['TCBC_CO is fixed!'])
end

try results.TPSC_GA;
catch
    results.TPSC_GA=zeros(size(results.RWIN_EL));
    warning(['TPSC_GA is fixed!'])
end

try results.TCWC_CO;
catch
    results.TCWC_CO=zeros(size(results.RWIN_EL));
    warning(['TCWC_CO is fixed!'])
end

try results.THCS_CO;
catch
    results.THCS_CO=zeros(size(results.RWIN_EL));
    warning(['THCS_CO is fixed!'])
end

try results.TGCS_CO;
catch
    results.TGCS_CO=zeros(size(results.RWIN_EL));
    warning(['TGCS_CO is fixed!'])
end

try results.TSMC_CO;
catch
    results.TSMC_CO=zeros(size(results.RWIN_EL));
    warning(['TSMC_CO is fixed!'])
end

try results.EL_TPSP;
catch
    results.EL_TPSP=zeros(size(results.RWIN_EL));
    warning(['EL_TPSP is fixed!'])
end


try results.EL_TPSC;
catch
    results.EL_TPSC=zeros(size(results.RWIN_EL));
    warning(['EL_TPSC is fixed!'])
end

try results.EL_TRSS;
catch
    results.EL_TRSS=zeros(size(results.RWIN_EL));
    results.EL_TRSL=zeros(size(results.RWIN_EL));
    results.EL_TRMI=zeros(size(results.RWIN_EL));
    results.EL_TRMO=zeros(size(results.RWIN_EL));
    results.EL_TRME=zeros(size(results.RWIN_EL));
    results.EL_TRSi=zeros(size(results.RWIN_EL));
    results.EL_TRAD=zeros(size(results.RWIN_EL));
    results.EL_TREW=zeros(size(results.RWIN_EL));
    results.EL_TRBC=zeros(size(results.RWIN_EL));
    results.EL_TRGE=zeros(size(results.RWIN_EL));

    results.TRSS_GA=zeros(size(results.RWIN_EL));
    results.TRSL_GA=zeros(size(results.RWIN_EL));
    results.TRMI_GA=zeros(size(results.RWIN_EL));
    results.TRMO_GA=zeros(size(results.RWIN_EL));
    results.TRME_GA=zeros(size(results.RWIN_EL));
    results.TRSi_GA=zeros(size(results.RWIN_EL));
    results.TRAD_GA=zeros(size(results.RWIN_EL));
    results.TREW_GA=zeros(size(results.RWIN_EL));
    results.TRBC_GA=zeros(size(results.RWIN_EL));
    results.TRGE_GA=zeros(size(results.RWIN_EL));
end

try results.RPBO_EL;
catch
    results.RPBO_EL=zeros(size(results.RWIN_EL));
    results.RPBA_EL=zeros(size(results.RWIN_EL));
    results.RPBV_EL=zeros(size(results.RWIN_EL));
    results.RPVF_EL=zeros(size(results.RWIN_EL));
end

try results.TBCS_EL;
catch
    results.TBCS_EL=zeros(size(results.RWIN_EL));
    results.TCBC_EL=zeros(size(results.RWIN_EL));
    results.TCWC_EL=zeros(size(results.RWIN_EL));
end

try results.RWOO_TBCS;
catch
    results.RWOO_TBCS=zeros(size(results.RWIN_EL));
    results.RWOO_TCBC=zeros(size(results.RWIN_EL));    
    results.RWWO_TBCS=zeros(size(results.RWIN_EL));
    results.RWWO_TCBC=zeros(size(results.RWIN_EL));
end

try results.HY_TRSi;
catch
    results.HY_TRSi=zeros(size(results.RWIN_EL));
    warning(['HY_TRSi is fixed!'])
end
try results.MMFA_Dem;
catch
    results.MMFA_Dem=zeros(size(results.MMFG_Dem));
    results.MMPA_Dem=zeros(size(results.MMFG_Dem));
    warning(['Ammonia for marine is fixed!'])
end

try results.MMFM_Dem;
catch
    results.MMFM_Dem=zeros(size(results.MMFG_Dem));
    results.MMPM_Dem=zeros(size(results.MMFG_Dem));
    warning(['Methanol for marine is fixed!'])
end
try results.RDSH_HE;
catch
    results.RDSH_HE=zeros(size(results.RWIN_EL));
    results.OPT_SIZE_RDSH=zeros(size(results.OPT_SIZE_RWIN));
    warning(['RDSH_HE is fixed!'])
end

try results.importEL;
catch
    results.importEL=zeros(size(results.RWIN_EL));
    setup.importELCost= 1000*ones(size(setup.importH2Cost));
    warning(['importEL is fixed!'])
end

try results.HY_LHIN;
catch
    results.HY_LHIN=zeros(size(results.RWIN_EL));
    warning(['HY_LHIN is fixed!'])
end

try results.RWOO_LHIN;
catch
    results.RWOO_LHIN=zeros(size(results.RWIN_EL));
    warning(['RWOO_LHIN is fixed!'])
end

try results.RWWO_LHIN;
catch
    results.RWWO_LHIN=zeros(size(results.RWIN_EL));
    warning(['RWWO_LHIN is fixed!'])
end
try results.RBGA_LHIN;
catch
    results.RBGA_LHIN=zeros(size(results.RWIN_EL));
    warning(['RBGA_LHIN is fixed!'])
end

try results.RHAR_LHIN;
catch
    results.RHAR_LHIN=zeros(size(results.RWIN_EL));
    warning(['RHAR_LHIN is fixed!'])
end

try results.RPET_LHIN;
catch
    results.RPET_LHIN=zeros(size(results.RWIN_EL));
    warning(['RPET_LHIN is fixed!'])
end

try results.RPET_LHIN;
catch
    results.RPET_LHIN=zeros(size(results.RWIN_EL));
    warning(['RPET_LHIN is fixed!'])
end

try results.EL_TFTU;

catch
    results.EL_TFTU=zeros(size(results.RWIN_EL));
    warning(['EL_TFTU is fixed!'])
end

try results.LHIN_GAS_FOS;

catch
    results.LHIN_GAS_FOS = zeros(size(results.RWIN_EL));
    warning(['LHIN_GAS_FOS is fixed!'])
end

try results.LHIN_GAS_REN;

catch
    results.LHIN_GAS_REN = zeros(size(results.RWIN_EL));
    warning(['LHIN_GAS_REN is fixed!'])
end
try results.OPT_SIZE_RPBO;
catch
    results.OPT_SIZE_RPBO = zeros(size(results.OPT_SIZE_RWIN));
    results.OPT_SIZE_RPBA = zeros(size(results.OPT_SIZE_RWIN));
    results.OPT_SIZE_RPBV = zeros(size(results.OPT_SIZE_RWIN));
    results.OPT_SIZE_RPVF = zeros(size(results.OPT_SIZE_RWIN));
    results.OPT_SIZE_TBCS = zeros(size(results.OPT_SIZE_RWIN));
    results.OPT_SIZE_TCBC = zeros(size(results.OPT_SIZE_RWIN));
    results.OPT_SIZE_TCWC = zeros(size(results.OPT_SIZE_RWIN));
    
end

try results.OPT_SIZE_TICM;

catch
    results.OPT_SIZE_TICM = zeros(size(results.OPT_SIZE_RWIN));
    results.TICM_EL = zeros(size(results.RWIN_EL));
    results.TICM_GAS_FOS = zeros(size(results.RWIN_EL));
    results.TICM_GAS_REN = zeros(size(results.RWIN_EL));
    results.RPET_TICM = zeros(size(results.RWIN_EL));
    results.HY_TICM = zeros(size(results.RWIN_EL));
    warning(['TICM is fixed!'])
end

try results.OPT_SIZE_TGCS;

catch
    results.OPT_SIZE_TGCS = zeros(size(results.OPT_SIZE_RWIN));
    results.TGCS_EL = zeros(size(results.RWIN_EL));
    results.TGCS_GAS_FOS = zeros(size(results.RWIN_EL));
    results.TGCS_GAS_REN = zeros(size(results.RWIN_EL));
    results.RPET_TGCS = zeros(size(results.RWIN_EL));
    results.HY_TGCS = zeros(size(results.RWIN_EL));
    warning(['TGCS is fixed!'])
end

try results.OPT_SIZE_THCS;

catch
    results.OPT_SIZE_THCS = zeros(size(results.OPT_SIZE_RWIN));
    results.THCS_EL = zeros(size(results.RWIN_EL));
    results.RHAR_THCS = zeros(size(results.RWIN_EL));
    warning(['THCS is fixed!'])
end

try results.RPET_TCCG;

catch
    results.RPET_TCCG = zeros(size(results.RWIN_EL));
    results.RPET_TOCG = zeros(size(results.RWIN_EL));
    results.RPET_TCNG = zeros(size(results.RWIN_EL));
    results.HY_TCCG = zeros(size(results.RWIN_EL));
    results.HY_TOCG = zeros(size(results.RWIN_EL));
    results.HY_TCNG = zeros(size(results.RWIN_EL));
    warning(['multifuel gas turbines are fixed!'])
end

try
    results.importFT;
catch
    results.importFT = zeros(size(results.RWIN_EL));
    results.importFT_diesel = zeros(size(results.RWIN_EL));
    results.importFT_kerosene = zeros(size(results.RWIN_EL));
    results.importH2 = zeros(size(results.RWIN_EL));
    results.importH2_liq = zeros(size(results.RWIN_EL));
    results.importH2_regas = zeros(size(results.RWIN_EL));
    results.importLNG = zeros(size(results.RWIN_EL));
    results.importLNG_liq = zeros(size(results.RWIN_EL));
    results.importLNG_regas = zeros(size(results.RWIN_EL));

end


try
    results.MMPA_Dem;
catch
    results.MMPA_Dem = ones(size(results.MMFF_Dem)).*repmat((systemParams.Instalations(:,systemParams.IndexYears==2020,ismember(systemParams.IndexID,'LMMP')).*systemParams.Mobility.SharesTot(:,systemParams.IndexYears==2020,ismember(systemParams.IndexID,'MMPA')).*systemParams.Mobility.Cons_Primary_S(:,systemParams.IndexYears==2020,ismember(systemParams.Mobility.Cons_Names_Primary,'MMPA')))'/8760,8760,1);
end

try
    results.MMPM_Dem;
catch
    results.MMPM_Dem = ones(size(results.MMFF_Dem)).*repmat((systemParams.Instalations(:,systemParams.IndexYears==2020,ismember(systemParams.IndexID,'LMMP')).*systemParams.Mobility.SharesTot(:,systemParams.IndexYears==2020,ismember(systemParams.IndexID,'MMPM')).*systemParams.Mobility.Cons_Primary_S(:,systemParams.IndexYears==2020,ismember(systemParams.Mobility.Cons_Names_Primary,'MMPM')))'/8760,8760,1);
end

try results.exportH2
catch
    results.exportH2 = zeros(size(results.RWIN_EL));
end
try results.exportLNG
catch
    results.exportLNG = zeros(size(results.RWIN_EL));
end
try results.exportFT
catch
    results.exportFT = zeros(size(results.RWIN_EL));
end
try results.exportFT_kerosene
catch
    results.exportFT_kerosene = zeros(size(results.RWIN_EL));
end
try results.exportFT_diesel
catch
    results.exportFT_diesel = zeros(size(results.RWIN_EL));
end
try results.exportNH3
catch
    results.exportNH3 = zeros(size(results.RWIN_EL));
end
try results.exportMeOH
catch
    results.exportMeOH = zeros(size(results.RWIN_EL));
end
try results.exportEL
catch
    results.exportEL = zeros(size(results.RWIN_EL));
end


for ii = 1:length(results.OPT_SIZE_SBAT)

    results.OPT_SIZE_SBAT(ii) = max(results.OPT_SIZE_SBAT(ii),systemParams.SizeLimits(ismember(systemParams.IndexID,'SBAT'),1,ii));
    results.OPT_SIZE_IBAT(ii) = max(results.OPT_SIZE_IBAT(ii),systemParams.SizeLimits(ismember(systemParams.IndexID,'IBAT'),1,ii));

end

%% fix coal
if sum(sum(results.RHAR_FU)) < sum(sum(results.RHAR_LHIN+results.RHAR_TCCO+results.RHAR_TDCO+results.RHAR_THCS+results.RHAR_THPP+results.HC_TISB))
    results.RHAR_FU=results.RHAR_LHIN+results.RHAR_TCCO+results.RHAR_TDCO+results.RHAR_THCS+results.RHAR_THPP+results.HC_TISB;
end


%% add other PED
try
    nn = xlsread([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep 'PED_remaining_GAS.xlsx']);
    results.RNGA_OTHER = repmat((nn(2:size(systemParams.IndexNodes)+1,2+find(ismember(systemParams.IndexYears,costYear)))/8760),1,8760)';
    results.RNGA_FU = results.RNGA_FU + results.RNGA_OTHER;
catch
    results.RNGA_OTHER = zeros(size(results.RNGA_FU));
end
results.RHAR_CO = results.RHAR_CO + results.RNGA_OTHER * systemParams.GasEmissionsMain;

try
    nn = xlsread([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep 'PED_remaining_OIL.xlsx']);
    results.RPET_OTHER = repmat((nn(2:size(systemParams.IndexNodes)+1,2+find(ismember(systemParams.IndexYears,costYear)))/8760),1,8760)';
    results.RPET_FU = results.RPET_FU + results.RPET_OTHER;
catch
    results.RPET_OTHER = zeros(size(results.RPET_FU));
end
results.RPET_CO = results.RPET_CO + results.RPET_OTHER * systemParams.OilEmissionsMain;

try
    nn = xlsread([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep 'PED_remaining_COAL.xlsx']);
    results.RHAR_OTHER = repmat((nn(2:size(systemParams.IndexNodes)+1,2+find(ismember(systemParams.IndexYears,costYear)))/8760),1,8760)';
    results.RHAR_FU = results.RHAR_FU + results.RHAR_OTHER;
catch
    results.RHAR_OTHER = zeros(size(results.RHAR_FU));
end
results.RHAR_CO = results.RHAR_CO + results.RHAR_OTHER * systemParams.CoalEmissionsMain;

%% load prosumers sector results
if SCFlag
    if setup.MacroMc
        resultsSC = ReadProsResTransMacro(setup,systemParams,Reg1,Reg2,costYear);
    else
        if setup.Big
            resultsSC = ReadProsResTrans(setup,systemParams,Reg1,Reg2,costYear,name);
        else
            resultsSC = ReadProsResTrans(setup,systemParams,Reg1,Reg2,costYear);
        end
    end



else

    resultsSC.OPT_SIZE_IBAC = zeros(1,length(systemParams.IndexNodes));
    resultsSC.OPT_SIZE_IBAI = zeros(1,length(systemParams.IndexNodes));
    resultsSC.OPT_SIZE_IBAR = zeros(1,length(systemParams.IndexNodes));
    resultsSC.OPT_SIZE_IHHS = zeros(1,length(systemParams.IndexNodes));
    resultsSC.OPT_SIZE_SHHS = zeros(1,length(systemParams.IndexNodes));
    resultsSC.OPT_SIZE_SBAC = zeros(1,length(systemParams.IndexNodes));
    resultsSC.OPT_SIZE_SBAI = zeros(1,length(systemParams.IndexNodes));
    resultsSC.OPT_SIZE_SBAR = zeros(1,length(systemParams.IndexNodes));

    resultsSC.OPT_SIZE_RPVC = zeros(1,length(systemParams.IndexNodes));
    resultsSC.OPT_SIZE_RPVI = zeros(1,length(systemParams.IndexNodes));
    resultsSC.OPT_SIZE_RPVR = zeros(1,length(systemParams.IndexNodes));
    resultsSC.OPT_SIZE_RRSH = zeros(1,length(systemParams.IndexNodes));

    resultsSC.OPT_SIZE_THBG = zeros(1,length(systemParams.IndexNodes));
    resultsSC.OPT_SIZE_THBP = zeros(1,length(systemParams.IndexNodes));
    resultsSC.OPT_SIZE_THHP = zeros(1,length(systemParams.IndexNodes));
    resultsSC.OPT_SIZE_THHR = zeros(1,length(systemParams.IndexNodes));
    resultsSC.OPT_SIZE_THNG = zeros(1,length(systemParams.IndexNodes));
    resultsSC.OPT_SIZE_THOI = zeros(1,length(systemParams.IndexNodes));

    resultsSC.Cap_IBAC = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));
    resultsSC.Cap_IBAI = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));
    resultsSC.Cap_IBAR = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));
    resultsSC.Cap_IHHS = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));
    resultsSC.Cap_SHHS = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));
    resultsSC.Cap_SBAC = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));
    resultsSC.Cap_SBAI = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));
    resultsSC.Cap_SBAR = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));

    resultsSC.Cap_RPVC = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));
    resultsSC.Cap_RPVI = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));
    resultsSC.Cap_RPVR = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));
    resultsSC.Cap_RRSH = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));

    resultsSC.Cap_THBG = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));
    resultsSC.Cap_THBP = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));
    resultsSC.Cap_THHP = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));
    resultsSC.Cap_THHR = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));
    resultsSC.Cap_THNG = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));
    resultsSC.Cap_THOI = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears));

    resultsSC.EL_EXCESS = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.EL_GRID = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.EL_SHHS = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.EL_SBAC = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.EL_SBAI = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.EL_SBAR = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.EL_THHP = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.EL_THHR = zeros(8760,length(systemParams.IndexNodes));

    resultsSC.SHHS_EL = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.SBAC_EL = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.SBAI_EL = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.SBAR_EL = zeros(8760,length(systemParams.IndexNodes));

    resultsSC.THBG_HE = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.THBP_HE = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.THHP_HE = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.THHR_HE = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.THNG_HE = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.THOI_HE = zeros(8760,length(systemParams.IndexNodes));

    resultsSC.RBGA_FU = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.RHAR_FU = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.RNGA_FU = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.RPET_FU = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.RWOO_FU = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.RWWO_FU = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.RRSH_HE = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.RPVC_EL = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.RPVI_EL = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.RPVR_EL = zeros(8760,length(systemParams.IndexNodes));

    resultsSC.HE_EXCESS_Local = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.SoC_SHHS = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.SoC_SBAC = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.SoC_SBAI = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.SoC_SBAR = zeros(8760,length(systemParams.IndexNodes));

    resultsSC.RES.DemP = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.RES.Dem = zeros(1,length(systemParams.IndexNodes));
    resultsSC.RES.EL_EXCESS = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.RES.EL_GRID = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.RES.EL_SBAT = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.RES.OPT_SIZE_IBAT = zeros(length(systemParams.IndexNodes));
    resultsSC.RES.OPT_SIZE_SBAT = zeros(length(systemParams.IndexNodes));
    resultsSC.RES.OPT_SIZE_RPVO = zeros(length(systemParams.IndexNodes));
    resultsSC.RES.RPVO_EL = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.RES.SBAT_EL = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.RES.SC_demand = zeros(8760,length(systemParams.IndexNodes));

    resultsSC.COM.DemP = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.COM.Dem = zeros(1,length(systemParams.IndexNodes));
    resultsSC.COM.EL_EXCESS = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.COM.EL_GRID = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.COM.EL_SBAT = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.COM.OPT_SIZE_IBAT = zeros(length(systemParams.IndexNodes));
    resultsSC.COM.OPT_SIZE_SBAT = zeros(length(systemParams.IndexNodes));
    resultsSC.COM.OPT_SIZE_RPVO = zeros(length(systemParams.IndexNodes));
    resultsSC.COM.RPVO_EL = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.COM.SBAT_EL = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.COM.SC_demand = zeros(8760,length(systemParams.IndexNodes));

    resultsSC.IND.DemP = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.IND.Dem = zeros(1,length(systemParams.IndexNodes));
    resultsSC.IND.EL_EXCESS = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.IND.EL_GRID = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.IND.EL_SBAT = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.IND.OPT_SIZE_IBAT = zeros(length(systemParams.IndexNodes));
    resultsSC.IND.OPT_SIZE_SBAT = zeros(length(systemParams.IndexNodes));
    resultsSC.IND.OPT_SIZE_RPVO = zeros(length(systemParams.IndexNodes));
    resultsSC.IND.RPVO_EL = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.IND.SBAT_EL = zeros(8760,length(systemParams.IndexNodes));
    resultsSC.IND.SC_demand = zeros(8760,length(systemParams.IndexNodes));

end
%% fix oil
try
    if sum(results.RPET_FU,'all') < sum(results.RPET_DI+results.RPET_KE+results.RPET_LHIN+results.RPET_TCCG+results.RPET_TCNG+results.RPET_TCOI+results.RPET_TDOI+results.RPET_TGCS+results.RPET_TICG+results.RPET_TICM+results.RPET_TOCG+resultsSC.RPET_FU,'all')
        results.RPET_FU = (results.RPET_DI+results.RPET_KE+results.RPET_LHIN+results.RPET_TCCG+results.RPET_TCNG+results.RPET_TCOI+results.RPET_TDOI+results.RPET_TGCS+results.RPET_TICG+results.RPET_TICM+results.RPET_TOCG+resultsSC.RPET_FU);
    end
end

%% Basic information
% Number of regions
regNumb = size(systemParams.IndexNodes,1);
output.regNumb = regNumb;

yearNumb = find(systemParams.IndexYears==costYear);

%% Data preparation (for versions compatibility)




% if no grid
try results.GRID_AC;

catch
    results.GRID_AC = zeros(size(results.RWIN_EL));
    results.GRID_DC = zeros(size(results.RWIN_EL));
end


%demand computations
demand = shiftdim(sum(systemParams.ValueLoad(1,:,1,:)),3);



allActiveExLoad = (activeElements.activeAllTransformer|activeElements.activeHydro|activeElements.activeElFeedIn|activeElements.activeHeatFeedIn|activeElements.activeResource|activeElements.activeStorage|activeElements.activeStorageInterface|activeElements.activeDesalination|activeElements.activeTransmission);
allActiveExLoad = (activeElements.activeLHeatTransformer|activeElements.activeLHeatFeedIn|activeElements.activeAllTransformer|activeElements.activeHydro|activeElements.activeElFeedIn|activeElements.activeHeatFeedIn|activeElements.activeResource|activeElements.activeStorage|activeElements.activeStorageInterface|activeElements.activeDesalination|activeElements.activeTransmission|ismember(systemParams.IndexID,{'SHHS','IHHS'}));
allActiveExLoadLabels = systemParams.IndexID(allActiveExLoad);

if costYear==2015
    instalMask = (repmat(systemParams.IndexYears,size(systemParams.IndexID,1),1)<=costYear);
else
    instalMask = ((repmat(systemParams.IndexYears,size(systemParams.IndexID,1),1)<=costYear).*(repmat(systemParams.IndexYears,size(systemParams.IndexID,1),1)>=(costYear-systemParams.Lifetime)));

end

%% fix Hydro Instalations

%Reconstructed hydro capacities
RoR_LL = sum(shiftdim(systemParams.Instalations(:,1:find(systemParams.IndexYears==costYear-Lifetime_Rehab)-1,find(ismember(systemParams.IndexID,'RRRI'))),1));
Dam_LL = sum(shiftdim(systemParams.Instalations(:,1:find(systemParams.IndexYears==costYear-Lifetime_Rehab)-1,find(ismember(systemParams.IndexID,'HDAM'))),1));
PHS_LL = sum(shiftdim(systemParams.Instalations(:,1:find(systemParams.IndexYears==costYear-Lifetime_Rehab)-1,find(ismember(systemParams.IndexID,'SPHS'))),1));

% get capacities
actCount = 0;
for i=1:length(allActiveExLoad)
    if allActiveExLoad(i) & ~activeElements.activeTransmission(i)
        actCount = actCount + 1;
        fString = ['OPT_SIZE_' systemParams.IndexID{i}];
        if isfield(results,fString)

            capacityCumul(:,actCount) = getfield(results,fString)';
            if setup.OvernightFlag
                for reg = 1:size(systemParams.IndexNodes,1)
                    if strcmp(systemParams.IndexID{i},'RRRI')
                        capacityCumul(reg,actCount) = capacityCumul(reg,actCount)-RoR_LL(reg);
                    end
                    if strcmp(systemParams.IndexID{i},'HDAM')
                        capacityCumul(reg,actCount) = capacityCumul(reg,actCount)-Dam_LL(reg);
                    end
                    if strcmp(systemParams.IndexID{i},'SPHS')
                        capacityCumul(reg,actCount) = capacityCumul(reg,actCount)-PHS_LL(reg);
                    end
                    capacity(reg,:,actCount) = zeros(1,length(systemParams.IndexYears)); % (regions,years,tech)
                    capacity(reg,systemParams.IndexYears==costYear,actCount) = capacityCumul(reg,actCount);
                end
            else
                for reg = 1:size(systemParams.IndexNodes,1)
                    capacity(reg,:,actCount) = systemParams.Instalations(reg,:,i).*instalMask(i,:); % (regions,years,tech)
                    if strcmp(systemParams.IndexID{i},'RRRI')
                        capacityCumul(reg,actCount) = capacityCumul(reg,actCount)-RoR_LL(reg);
                        capacity(reg,:,actCount) = capacity(reg,:,actCount).*(systemParams.IndexYears>costYear-Lifetime_Rehab-1); % (regions,years,tech)
                    end
                    if strcmp(systemParams.IndexID{i},'HDAM')
                        capacityCumul(reg,actCount) = capacityCumul(reg,actCount)-Dam_LL(reg);
                        capacity(reg,:,actCount) = capacity(reg,:,actCount).*(systemParams.IndexYears>costYear-Lifetime_Rehab-1); % (regions,years,tech)
                    end
                    if strcmp(systemParams.IndexID{i},'SPHS')
                        capacityCumul(reg,actCount) = capacityCumul(reg,actCount)-PHS_LL(reg);
                        capacity(reg,:,actCount) = capacity(reg,:,actCount).*(systemParams.IndexYears>costYear-Lifetime_Rehab-1); % (regions,years,tech)
                    end

                    capacity(reg,systemParams.IndexYears==costYear,actCount) = capacityCumul(reg,actCount) - sum(capacity(reg,systemParams.IndexYears<costYear,actCount));
                end
            end
        else
            capacity(:,:,actCount) = zeros(size(systemParams.IndexNodes,1),size(systemParams.IndexYears,2),1);
        end
    elseif activeElements.activeTransmission(i)
        actCount = actCount + 1;
        fString = ['OPT_SIZE_' systemParams.IndexID{i}];
        if isfield(results,fString)
            if prod(size(getfield(results,fString)))==0
                systemParams.IndexID{i}
                capacityCumul(1,actCount) = 0;
            else
                capacityCumul(1,actCount) = sum(getfield(results,fString)');
            end
            capacity(1,systemParams.IndexYears==costYear,actCount) = capacityCumul(1,actCount);
        else
            capacity(:,:,actCount) = zeros(size(systemParams.IndexNodes,1),size(systemParams.IndexYears,2),1);
            if strcmp(systemParams.IndexID{i},'TRCS')
                capacityCumul(:,actCount) = max(abs(results.GRID_DC))';
                capacity(:,systemParams.IndexYears==costYear,actCount) = capacityCumul(:,actCount);
            end

        end

    end
end


if strcmp(type, 'Temp')

    for reg = 1:size(systemParams.IndexNodes,1)
        for i = 1:size(capacity,3)
            temp = sum(capacity(reg,:,i),2);
            capacity(reg,:,i) = capacity(reg,:,i)*0;
            capacity(reg,systemParams.IndexYears==costYear,i) = temp;
        end
    end

end

if OvernightFlag

end

capacity(capacity<0)=0;

if Only.Flag
    if SCFlag
        shareOfSector = Only.Power./(Only.Power+Only.Gas+Only.Desalination);
    elseif GasFlag
        shareOfSector = Only.Gas./(Only.Power+Only.Gas+Only.Desalination);
    elseif DesalinationFlag
        shareOfSector = Only.Desalination./(Only.Power+Only.Gas+Only.Desalination);
    end
    RoR_LL = RoR_LL*shareOfSector;
    Dam_LL = Dam_LL*shareOfSector;
    PHS_LL = PHS_LL*shareOfSector;
end

% get capex
CRFvalue = [];
for i=systemParams.IndexYears
    CRFvalue = [CRFvalue, CapitalRecoveryFactor(systemParams.WACC,systemParams.Lifetime(allActiveExLoad,systemParams.IndexYears==i))];
end
CRFvalue = permute(repmat(CRFvalue,1,1,size(systemParams.IndexNodes,1)),[3 2 1]);
if ndims(systemParams.Capex_reg)==2
    systemParams.Capex_reg = repmat(systemParams.Capex_reg,1,1,length(systemParams.IndexYears));
    systemParams.Opex_fix_reg = repmat(systemParams.Opex_fix_reg,1,1,length(systemParams.IndexYears));
    systemParams.Opex_var_reg = repmat(systemParams.Opex_var_reg,1,1,length(systemParams.IndexYears));
end
%% Add new Uranium costs
warning('Adds new Uranium costs, check if values are up to date')
systemParams.Opex_var(find(ismember(systemParams.IndexID,'RURA')),:) = repmat([0.0026],1,length(systemParams.IndexYears));

%% Add Crude oil refinery cost
warning('Adds new oil refinery cost, check if values are up to date')
LCOR = 7.18/1000;%[EUR/kWh]
systemParams.Opex_var(find(ismember(systemParams.IndexID,'RPET')),:) = systemParams.Opex_var(find(ismember(systemParams.IndexID,'RPET')),:) + LCOR;

%% Recalculated Excess heat


capexCRF = (permute(repmat(systemParams.Capex(allActiveExLoad,:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(allActiveExLoad,:,:),[2 3 1]).* CRFvalue);
capex = (permute(repmat(systemParams.Capex(allActiveExLoad,:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(allActiveExLoad,:,:),[2 3 1]));
opex_fix = (permute(repmat(systemParams.Opex_fix(allActiveExLoad,:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_fix_reg(allActiveExLoad,:,:),[2 3 1]));
opex_var = (permute(repmat(systemParams.Opex_var(allActiveExLoad,:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_var_reg(allActiveExLoad,:,:),[2 3 1]));


totalCapex = sum(sum(capacity.*capex,2),3)...
    +Dam_LL'.*Capex_Rehab...
    +RoR_LL'.*Capex_Rehab;

totalCapexNew = sum(sum(capacity(:,find(systemParams.IndexYears==costYear),:).*capex(:,find(systemParams.IndexYears==costYear),:),2),3)...
    ;
indivHeatCapexNew = capex(:,find(systemParams.IndexYears==costYear),ismember(allActiveExLoadLabels,'RRSH')).*(resultsSC.Cap_RRSH(:,find(systemParams.IndexYears==costYear))) +...
    capex(:,find(systemParams.IndexYears==costYear),ismember(allActiveExLoadLabels,'THHR')).*(resultsSC.Cap_THHR(:,find(systemParams.IndexYears==costYear))) +...
    capex(:,find(systemParams.IndexYears==costYear),ismember(allActiveExLoadLabels,'THHP')).*(resultsSC.Cap_THHP(:,find(systemParams.IndexYears==costYear))) +...
    capex(:,find(systemParams.IndexYears==costYear),ismember(allActiveExLoadLabels,'THNG')).*(resultsSC.Cap_THNG(:,find(systemParams.IndexYears==costYear))) +...
    capex(:,find(systemParams.IndexYears==costYear),ismember(allActiveExLoadLabels,'THOI')).*(resultsSC.Cap_THOI(:,find(systemParams.IndexYears==costYear))) +...%    opex_capex(:,:,ismember(allActiveExLoadLabels,'THCO')).*resultsSC.Cap_THCO +opex_var(:,:,ismember(allActiveExLoadLabels,'THCO')).*sum(resultsSC.THCO_HE).*(resultsSC.Cap_THCO/resultsSC.OPT_SIZE_THCO) +...
    capex(:,find(systemParams.IndexYears==costYear),ismember(allActiveExLoadLabels,'THBP')).*(resultsSC.Cap_THBP(:,find(systemParams.IndexYears==costYear))) +...
    capex(:,find(systemParams.IndexYears==costYear),ismember(allActiveExLoadLabels,'THBG')).*(resultsSC.Cap_THBG(:,find(systemParams.IndexYears==costYear)));

sfields = fieldnames(results);
for m=1:length(sfields)
    temp = regexp(sfields{m},'\w+_EL$','match');
    if ~isempty(temp) && isempty(regexp(temp{:},'^S\w+','match'))
        temp_sum(:,m) = sum(getfield(results,temp{:}));
    end
end
prod_electricity_sys = sum(temp_sum,2);


for k=1:length(allActiveExLoadLabels)
    compOutput(:,k) = zeros(size(systemParams.IndexNodes));
    for l=1:length(sfields)
        temp = regexp(sfields{l},['^' allActiveExLoadLabels{k} '_\w+'],'match');
        if ~isempty(temp)
            compOutput(:,k) = compOutput(:,k) + sum(getfield(results,temp{:}))';
        end
    end
end

% compOutput by years
capSharesByYears = capacity./repmat(sum(capacity,2),1,length(systemParams.IndexYears),1);
capSharesByYears(isnan(capSharesByYears)) = 0;

totalOpex = sum(sum(capacity.*opex_fix,2),3) + sum(sum(capSharesByYears.*permute(repmat(compOutput,1,1,length(systemParams.IndexYears)),[1 3 2]).*opex_var,2),3);

opex_capex = capexCRF+opex_fix;
opex_capexRehab_RoR = repmat((Capex_Rehab.* CapitalRecoveryFactor(systemParams.WACC,Lifetime_Rehab)),regNumb,length(systemParams.IndexYears),1).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'RRRI'),:,:),[2 3 1])+opex_fix(:,:,find(ismember(allActiveExLoadLabels,'RRRI')));
opex_capexRehab_Dam = repmat((Capex_Rehab.* CapitalRecoveryFactor(systemParams.WACC,Lifetime_Rehab)),regNumb,length(systemParams.IndexYears),1).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'HDAM'),:,:),[2 3 1])+opex_fix(:,:,find(ismember(allActiveExLoadLabels,'HDAM')));
opex_capexRehab_PHS = repmat((Capex_Rehab/8.* CapitalRecoveryFactor(systemParams.WACC,Lifetime_Rehab)),regNumb,length(systemParams.IndexYears),1).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'SPHS'),:,:),[2 3 1])+opex_fix(:,:,find(ismember(allActiveExLoadLabels,'SPHS')));


%% Grid calculations
% This formulation includes cost of converter station
if isfield(activeElements,'activeTransmission')
    capex_TL_DC = systemParams.Capex((ismember(systemParams.IndexID,'TRTL')),1);
    opex_fix_TL_DC = systemParams.Opex_fix((ismember(systemParams.IndexID,'TRTL')),1);
    opex_var_TL_DC = systemParams.Opex_var((ismember(systemParams.IndexID,'TRTL')),1);
    capex_TL_AC = systemParams.Capex((ismember(systemParams.IndexID,'THAO')),1);
    opex_fix_TL_AC = systemParams.Opex_fix((ismember(systemParams.IndexID,'THAO')),1);
    opex_var_TL_AC = systemParams.Opex_var((ismember(systemParams.IndexID,'THAO')),1);
    capex_CS = systemParams.Capex((ismember(systemParams.IndexID,'TRCS')),1);
    opex_fix_CS = systemParams.Opex_fix((ismember(systemParams.IndexID,'TRCS')),1);
else
    results.GRID_AC = zeros(size(results.RWIN_EL));
    results.GRID_DC = zeros(size(results.RWIN_EL));
end



% calculating share of each node at total transmission power
importPower = zeros(size(results.GRID_AC));
exportPower = zeros(size(results.GRID_AC));

importPower((results.GRID_AC + results.GRID_DC) >0) = (results.GRID_AC((results.GRID_AC + results.GRID_DC)>0) + results.GRID_DC((results.GRID_AC + results.GRID_DC)>0));
exportPower((results.GRID_AC + results.GRID_DC) <0) = (results.GRID_AC((results.GRID_AC + results.GRID_DC)<0) + results.GRID_DC((results.GRID_AC + results.GRID_DC)<0));

shareImport = sum(importPower,1) / sum(sum(importPower));
shareExport = sum(exportPower,1) / sum(sum(exportPower));
z=0.5;
Share = (z * shareExport + (1-z) * shareImport)';

for rr = Regs
    try
        qq = load([setup.rootDir filesep 'projects\Base\input-data' filesep baseData.systemParams.ValueLoadTag{1,rr} '_' num2str(costYear) '.mat']);
    catch
        qq = load([setup.rootDir filesep 'projects\Base\input-data' filesep baseData.systemParams.ValueLoadTag{1,rr} '.mat']);
    end
    profileName = (fieldnames(qq));
    baseData.systemParams.ValueLoad(1,:,1,rr) = qq.(profileName{1});

end
demand_orig = shiftdim(shiftdim(baseData.systemParams.ValueLoad(1,:,1,Regs),3),2)./repmat((sum(shiftdim(shiftdim(baseData.systemParams.ValueLoad(1,:,1,Regs),3),2),1)./baseData.systemParams.Instalations(Regs,systemParams.IndexYears==costYear,find(ismember(systemParams.IndexID,'LELE')))'),8760,1);
demand_orig(isnan(demand_orig))=0;

costYear

%demand calculation
genSC = resultsSC.RPVR_EL+resultsSC.RPVC_EL+resultsSC.RPVI_EL;
batSC = resultsSC.SBAR_EL+resultsSC.SBAC_EL+resultsSC.SBAI_EL-(resultsSC.EL_SBAR+resultsSC.EL_SBAC+resultsSC.EL_SBAI);
elDemSC = (resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand);
heDemSC = (resultsSC.EL_THHP+resultsSC.EL_THHR);
fromGridSC = (resultsSC.EL_GRID);
toGridSC = (resultsSC.EL_EXCESS);

el_cons_localpowerFromGrid_P=-((genSC+batSC-toGridSC-elDemSC).*((genSC+batSC-toGridSC)<elDemSC))';
el_cons_localpowerFromGrid=-sum((genSC+batSC-toGridSC-elDemSC).*((genSC+batSC-toGridSC)<elDemSC))';
el_cons_localheatFromGrid=-sum(((genSC+batSC-toGridSC-elDemSC).*((genSC+batSC-toGridSC-elDemSC)>0)-heDemSC).*((genSC+batSC-toGridSC-elDemSC).*((genSC+batSC-toGridSC-elDemSC)>0)<heDemSC))'+...
    sum(heDemSC.*((genSC+batSC-toGridSC)<elDemSC))';

el_cons_localheatFromGrid_P=-((genSC+batSC-toGridSC-elDemSC-heDemSC).*((genSC+batSC-toGridSC-elDemSC)>0).*((genSC+batSC-toGridSC-elDemSC-heDemSC)<0))'+...
    (heDemSC.*((genSC+batSC-toGridSC)<elDemSC))';
el_cons_localheatFromGrid=-sum((genSC+batSC-toGridSC-elDemSC-heDemSC).*((genSC+batSC-toGridSC-elDemSC)>0).*((genSC+batSC-toGridSC-elDemSC-heDemSC)<0))'+...
    sum(heDemSC.*((genSC+batSC-toGridSC)<elDemSC))';


resultsSC.GRID_EL_EX = -((demand_orig - resultsSC.RES.SC_demand-resultsSC.COM.SC_demand-resultsSC.IND.SC_demand)-resultsSC.EL_EXCESS + el_cons_localpowerFromGrid_P'+el_cons_localheatFromGrid_P');
resultsSC.GRID_EL_EX(resultsSC.GRID_EL_EX<0)=0;
((demand_orig - resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand)-resultsSC.EL_EXCESS + el_cons_localpowerFromGrid_P'+el_cons_localheatFromGrid_P')*(1-systemParams.AC_losses(1,find(systemParams.IndexYears==costYear))/100);

demand_base = shiftdim(sum(systemParams.ValueLoad(1,:,1,:)),3);
demand_imports = sum(results.GRID_AC,1)' + sum(results.GRID_DC,1)';
demand_sys = shiftdim(sum(systemParams.ValueLoad(1,:,1,:)),3) - sum(results.GRID_AC,1)' - sum(results.GRID_DC,1)' + sum(resultsSC.EL_EXCESS,1)' - sum(resultsSC.GRID_EL_EX)';
demand_des = sum(results.EL_WDES+results.EL_WHPU+results.EL_WVPU,1)'.*(1./(1-systemParams.AC_losses(:,systemParams.IndexYears==costYear)/100));
demand_tot = demand_sys + resultsSC.RES.Dem' + resultsSC.COM.Dem' + resultsSC.IND.Dem';
demand_tot = demand_sys + sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand+resultsSC.EL_THHP+resultsSC.EL_THHR,1)'-el_cons_localpowerFromGrid-el_cons_localheatFromGrid;
demand_origTot = demand+sum((resultsSC.RES.Dem+resultsSC.COM.Dem+resultsSC.IND.Dem+sum(resultsSC.EL_EXCESS-resultsSC.GRID_EL_EX,1)),1)';
demand_origTot = demand+sum((resultsSC.RES.Dem+resultsSC.COM.Dem+resultsSC.IND.Dem+sum(resultsSC.EL_EXCESS,1)),1)';
demand_origTot2 = demand_base+sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand,1)'+sum(resultsSC.EL_EXCESS-resultsSC.GRID_EL_EX,1)'-(el_cons_localpowerFromGrid+el_cons_localheatFromGrid)./((100-systemParams.AC_losses(:,yearNumb))/100);
demand_origTot = demand_sys+sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand,1)';
demand_origFinal = baseData.systemParams.Instalations(Regs,yearNumb,ismember(systemParams.IndexID,'LELE'));

demand_El_Power_tot = demand_base+sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand,1)'+sum(resultsSC.EL_EXCESS-resultsSC.GRID_EL_EX,1)'-(el_cons_localpowerFromGrid+el_cons_localheatFromGrid)./((100-systemParams.AC_losses(:,yearNumb))/100);
demand_El_Power_sys = demand_base+sum(resultsSC.EL_EXCESS-resultsSC.GRID_EL_EX,1)'-(el_cons_localheatFromGrid)./((100-systemParams.AC_losses(:,yearNumb))/100);
demand_El_Power_pros = sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand,1)'-(el_cons_localpowerFromGrid)./((100-systemParams.AC_losses(:,yearNumb))/100);

LocalGridsLoss = ((sum(demand_orig)')-(sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand))'-(sum(resultsSC.EL_EXCESS,1))'+(el_cons_localpowerFromGrid+el_cons_localheatFromGrid))*(1/(1-systemParams.AC_losses(1,find(systemParams.IndexYears==costYear))/100)-1);

TandDLossMain = (demand_origTot+(-el_cons_localpowerFromGrid-el_cons_localheatFromGrid)./((100-systemParams.AC_losses(:,yearNumb))/100)-demand_origFinal);
TandDLossOther = sum(results.EL_WDES+results.EL_WHPU+results.EL_WVPU+results.EltoMobility+results.EL_TDHR+results.EL_TDHP+results.EL_TDGE+results.EL_TICB+results.EL_TICI+results.EL_TISB+results.EL_TISH+results.EL_TISR+results.EL_TISE+results.EL_TIAM+results.EL_TIAR+results.EL_TIPP,1)'.*(-1+1./(1-systemParams.AC_losses(:,systemParams.IndexYears==costYear)/100));

shareOfSynGas=(sum(results.TMET_GA,1)./sum(results.TMET_GA+results.TBGU_GA+results.RNGA_FU,1))';
shareOfSynGas(isnan(shareOfSynGas))=0;

shareOfBMEonly=(sum(results.TBGU_GA,1)./sum(results.TMET_GA+results.TBGU_GA+results.RNGA_FU,1))';
shareOfBMEonly(isnan(shareOfBMEonly))=0;



%Share of SNG and biomethane used for Industry
shareOfGasInd = round(sum(results.LIGA_GAS_REN,1))'./round(sum(results.TMET_GA+results.TBGU_GA+results.importLNG_regas,1))';%*GasFlag;
shareOfGasInd(isnan(shareOfGasInd))=0;shareOfGasInd(isinf(shareOfGasInd))=1;
shareOfGasInd(shareOfGasInd<0)=0;
shareOfGasInd(shareOfGasInd>1)=1;
shareOfGasDes = round(sum(results.WDES_GAS_REN,1))'./round(sum(results.TMET_GA+results.TBGU_GA+results.importLNG_regas,1))';
shareOfGasDes(isnan(shareOfGasDes))=0;shareOfGasDes(isinf(shareOfGasDes))=1;
shareOfGasDH = round(sum(results.TDNG_GAS_REN,1))'./round(sum(results.TMET_GA+results.TBGU_GA+results.importLNG_regas,1))';
shareOfGasDH(isnan(shareOfGasDH))=0;shareOfGasDH(isinf(shareOfGasDH))=1;
shareOfGasIH = round(sum(results.Pros_GAS_REN,1))'./round(sum(results.TMET_GA+results.TBGU_GA+results.importLNG_regas,1))';
shareOfGasIH(isnan(shareOfGasIH))=0;shareOfGasIH(isinf(shareOfGasIH))=1;
shareOfGasINH = round(sum(results.LHIN_GAS_REN,1))'./round(sum(results.TMET_GA+results.TBGU_GA+results.importLNG_regas,1))';
shareOfGasINH(isnan(shareOfGasINH))=0;shareOfGasINH(isinf(shareOfGasINH))=1;
shareOfGasTrans = round(sum(results.TLNG_GAS_REN,1))'./round(sum(results.TMET_GA+results.TBGU_GA+results.importLNG_regas,1))';
shareOfGasTrans(isnan(shareOfGasTrans))=0;shareOfGasTrans(isinf(shareOfGasTrans))=1;
shareOfGasEl = round(sum(results.TCCG_GAS_REN+results.TOCG_GAS_REN+results.TCNG_GAS_REN+results.TGCS_GAS_REN+results.TICM_GAS_REN,1))'./round(sum(results.TMET_GA+results.TBGU_GA+results.importLNG_regas,1)+0.001)';
shareOfGasEl(isnan(shareOfGasEl))=0;shareOfGasEl(isinf(shareOfGasEl))=1;
shareOfGasExport = round(sum(results.exportLNG,1))'./round(sum(results.TMET_GA+results.TBGU_GA+results.importLNG_regas,1))';
shareOfGasExport(isnan(shareOfGasExport))=0;shareOfGasExport(isinf(shareOfGasExport))=1;

shareOfGasStor = shareOfGasEl; % old code


shareOfGasIndBME = sum(results.TBGU_GA,1)'./(sum(results.TBGU_GA+results.TMET_GA+results.importLNG_regas+0.001,1))';
shareOfGasIndBME(isnan(shareOfGasIndBME))=0;
shareOfGasIndSNG = sum(results.TMET_GA,1)'./(sum(results.TBGU_GA+results.TMET_GA+results.importLNG_regas+0.001,1))';
shareOfGasIndSNG(isnan(shareOfGasIndSNG))=0;
shareOfGasIndSNGimp = sum(results.importLNG_regas,1)'./(sum(results.TBGU_GA+results.TMET_GA+results.importLNG_regas+0.001,1))';
shareOfGasIndSNGimp(isnan(shareOfGasIndSNGimp))=0;

shareOfHyGreen = sum(results.TWEL_GA,1)'./(sum(results.TWEL_GA+results.TSMR_GA+results.importH2_regas+0.001,1))';
shareOfHyGreen(isnan(shareOfHyGreen))=0;
shareOfHyBlack = sum(results.TSMR_GA,1)'./(sum(results.TWEL_GA+results.TSMR_GA+results.importH2_regas+0.001,1))';
shareOfHyBlack(isnan(shareOfHyBlack))=0;
shareOfHyImpor = sum(results.importH2_regas,1)'./(sum(results.TWEL_GA+results.TSMR_GA+results.importH2_regas+0.001,1))';
shareOfHyImpor(isnan(shareOfHyImpor))=0;

shareOfBGA_Upgr = sum(results.RBGA_TBGU,1)'./(sum(results.RBGA_FU,1))';
shareOfBGA_Upgr(isnan(shareOfBGA_Upgr))=0;
shareOfBGA_DirEl = sum(results.RBGA_TCHP,1)'./(sum(results.RBGA_FU,1))';
shareOfBGA_DirEl(isnan(shareOfBGA_DirEl))=0;
shareOfBGA_DirIH = sum(resultsSC.RBGA_FU,1)'./(sum(results.RBGA_FU,1))';
shareOfBGA_DirIH(isnan(shareOfBGA_DirIH))=0;

shareOfBGA_Trans = shareOfBGA_Upgr.*shareOfGasTrans;
shareOfBGA_Ind = shareOfBGA_Upgr.*shareOfGasInd;
shareOfBGA_Des = shareOfBGA_Upgr.*shareOfGasDes;
shareOfBGA_DH = shareOfBGA_Upgr.*shareOfGasDH;
shareOfBGA_INH = shareOfBGA_Upgr.*shareOfGasINH;
shareOfBGA_IH = shareOfBGA_DirIH + shareOfBGA_Upgr.*shareOfGasIH;
shareOfBGA_El = shareOfBGA_DirEl + shareOfBGA_Upgr.*shareOfGasEl;

if ~(setup.Heat.Flag|setup.Mobility)

    shareOfBGA_El = ones(size(shareOfBGA_El));
end


%Share of NGA in GT  fuel
shareOfGT_NGA = sum(results.TGCS_GAS_FOS+results.TOCG_GAS_FOS+results.TCCG_GAS_FOS+results.TCNG_GAS_FOS,1)'./sum(results.TGCS_GAS_FOS+results.TOCG_GAS_FOS+results.TCCG_GAS_FOS+results.TCNG_GAS_FOS+results.TCNG_GAS_REN+results.TOCG_GAS_REN+results.TCCG_GAS_REN+results.TGCS_GAS_REN+results.RPET_TGCS+results.RPET_TOCG+results.RPET_TCCG+results.RPET_TCNG+results.HY_TGCS+results.HY_TOCG+results.HY_TCCG+results.HY_TCNG,1)';
shareOfGT_NGA(isnan(shareOfGT_NGA))=0;

shareOfGT_SNG = sum(results.TOCG_GAS_REN+results.TCCG_GAS_REN+results.TCNG_GAS_REN+results.TGCS_GAS_REN,1)'.*(1-shareOfGasIndBME)./sum(results.TGCS_GAS_FOS+results.TOCG_GAS_FOS+results.TCCG_GAS_FOS+results.TCNG_GAS_FOS+results.TCNG_GAS_REN+results.TOCG_GAS_REN+results.TCCG_GAS_REN+results.TGCS_GAS_REN+results.RPET_TGCS+results.RPET_TOCG+results.RPET_TCCG+results.RPET_TCNG+results.HY_TGCS+results.HY_TOCG+results.HY_TCCG+results.HY_TCNG,1)';
shareOfGT_SNG(isnan(shareOfGT_SNG))=0;

shareOfGT_BME = sum(results.TOCG_GAS_REN+results.TCCG_GAS_REN+results.TCNG_GAS_REN+results.TGCS_GAS_REN,1)'.*(shareOfGasIndBME)./sum(results.TGCS_GAS_FOS+results.TOCG_GAS_FOS+results.TCCG_GAS_FOS+results.TCNG_GAS_FOS+results.TCNG_GAS_REN+results.TOCG_GAS_REN+results.TCCG_GAS_REN+results.TGCS_GAS_REN+results.RPET_TGCS+results.RPET_TOCG+results.RPET_TCCG+results.RPET_TCNG+results.HY_TGCS+results.HY_TOCG+results.HY_TCCG+results.HY_TCNG,1)';
shareOfGT_BME(isnan(shareOfGT_BME))=0;

shareOfGT_SNGimp = sum(results.TOCG_GAS_REN+results.TCCG_GAS_REN+results.TCNG_GAS_REN+results.TGCS_GAS_REN,1)'.*(shareOfGasIndSNGimp)./sum(results.TGCS_GAS_FOS+results.TOCG_GAS_FOS+results.TCCG_GAS_FOS+results.TCNG_GAS_FOS+results.TCNG_GAS_REN+results.TOCG_GAS_REN+results.TCCG_GAS_REN+results.TGCS_GAS_REN+results.RPET_TGCS+results.RPET_TOCG+results.RPET_TCCG+results.RPET_TCNG+results.HY_TGCS+results.HY_TOCG+results.HY_TCCG+results.HY_TCNG,1)';
shareOfGT_SNGimp(isnan(shareOfGT_SNGimp))=0;

shareOfGT_Oil = sum(results.RPET_TGCS+results.RPET_TOCG+results.RPET_TCCG+results.RPET_TCNG,1)'./sum(results.TGCS_GAS_FOS+results.TOCG_GAS_FOS+results.TCCG_GAS_FOS+results.TCNG_GAS_FOS+results.TCNG_GAS_REN+results.TOCG_GAS_REN+results.TCCG_GAS_REN+results.TGCS_GAS_REN+results.RPET_TGCS+results.RPET_TOCG+results.RPET_TCCG+results.RPET_TCNG+results.HY_TGCS+results.HY_TOCG+results.HY_TCCG+results.HY_TCNG,1)';
shareOfGT_Oil(isnan(shareOfGT_Oil))=0;

shareOfGT_HY = sum(results.HY_TGCS+results.HY_TOCG+results.HY_TCCG+results.HY_TCNG,1)'./sum(results.TGCS_GAS_FOS+results.TOCG_GAS_FOS+results.TCCG_GAS_FOS+results.TCNG_GAS_FOS+results.TCNG_GAS_REN+results.TOCG_GAS_REN+results.TCCG_GAS_REN+results.TGCS_GAS_REN+results.RPET_TGCS+results.RPET_TOCG+results.RPET_TCCG+results.RPET_TCNG+results.HY_TGCS+results.HY_TOCG+results.HY_TCCG+results.HY_TCNG,1)';
shareOfGT_HY(isnan(shareOfGT_HY))=0;

shareOfGT_HYgreen = sum(results.HY_TGCS+results.HY_TOCG+results.HY_TCCG+results.HY_TCNG,1)'.*(shareOfHyGreen)./sum(results.TGCS_GAS_FOS+results.TOCG_GAS_FOS+results.TCCG_GAS_FOS+results.TCNG_GAS_FOS+results.TCNG_GAS_REN+results.TOCG_GAS_REN+results.TCCG_GAS_REN+results.TGCS_GAS_REN+results.RPET_TGCS+results.RPET_TOCG+results.RPET_TCCG+results.RPET_TCNG+results.HY_TGCS+results.HY_TOCG+results.HY_TCCG+results.HY_TCNG,1)';
shareOfGT_HYgreen(isnan(shareOfGT_HYgreen))=0;

shareOfGT_HYblack = sum(results.HY_TGCS+results.HY_TOCG+results.HY_TCCG+results.HY_TCNG,1)'.*(shareOfHyBlack)./sum(results.TGCS_GAS_FOS+results.TOCG_GAS_FOS+results.TCCG_GAS_FOS+results.TCNG_GAS_FOS+results.TCNG_GAS_REN+results.TOCG_GAS_REN+results.TCCG_GAS_REN+results.TGCS_GAS_REN+results.RPET_TGCS+results.RPET_TOCG+results.RPET_TCCG+results.RPET_TCNG+results.HY_TGCS+results.HY_TOCG+results.HY_TCCG+results.HY_TCNG,1)';
shareOfGT_HYblack(isnan(shareOfGT_HYblack))=0;

shareOfGT_HYimpor = sum(results.HY_TGCS+results.HY_TOCG+results.HY_TCCG+results.HY_TCNG,1)'.*(shareOfHyImpor)./sum(results.TGCS_GAS_FOS+results.TOCG_GAS_FOS+results.TCCG_GAS_FOS+results.TCNG_GAS_FOS+results.TCNG_GAS_REN+results.TOCG_GAS_REN+results.TCCG_GAS_REN+results.TGCS_GAS_REN+results.RPET_TGCS+results.RPET_TOCG+results.RPET_TCCG+results.RPET_TCNG+results.HY_TGCS+results.HY_TOCG+results.HY_TCCG+results.HY_TCNG,1)';
shareOfGT_HYimpor(isnan(shareOfGT_HYimpor))=0;

shareOfICM_NGA = sum(results.TICM_GAS_FOS,1)'./sum(results.TICM_GAS_FOS+results.TICM_GAS_REN+results.RPET_TICM+results.HY_TICM,1)';
shareOfICM_NGA(isnan(shareOfICM_NGA))=0;

shareOfICM_SNG = sum(results.TICM_GAS_REN,1)'.*(1-shareOfGasIndBME)./sum(results.TICM_GAS_FOS+results.TICM_GAS_REN+results.RPET_TICM+results.HY_TICM,1)';
shareOfICM_SNG(isnan(shareOfICM_SNG))=0;

shareOfICM_BME = sum(results.TICM_GAS_REN,1)'.*(shareOfGasIndBME)./sum(results.TICM_GAS_FOS+results.TICM_GAS_REN+results.RPET_TICM+results.HY_TICM,1)';
shareOfICM_BME(isnan(shareOfICM_BME))=0;

shareOfICM_SNGimp = sum(results.TICM_GAS_REN,1)'.*(shareOfGasIndSNGimp)./sum(results.TICM_GAS_FOS+results.TICM_GAS_REN+results.RPET_TICM+results.HY_TICM,1)';
shareOfICM_SNGimp(isnan(shareOfICM_SNGimp))=0;

shareOfICM_Oil = sum(results.RPET_TICM,1)'./sum(results.TICM_GAS_FOS+results.TICM_GAS_REN+results.RPET_TICM+results.HY_TICM,1)';
shareOfICM_Oil(isnan(shareOfICM_Oil))=0;

shareOfICM_HY = sum(results.HY_TICM,1)'./sum(results.TICM_GAS_FOS+results.TICM_GAS_REN+results.RPET_TICM+results.HY_TICM,1)';
shareOfICM_HY(isnan(shareOfICM_HY))=0;

shareOfICM_HYgreen = sum(results.HY_TICM,1)'.*(shareOfHyGreen)./sum(results.TICM_GAS_FOS+results.TICM_GAS_REN+results.RPET_TICM+results.HY_TICM,1)';
shareOfICM_HYgreen(isnan(shareOfICM_HYgreen))=0;

shareOfICM_HYblack = sum(results.HY_TICM,1)'.*(shareOfHyBlack)./sum(results.TICM_GAS_FOS+results.TICM_GAS_REN+results.RPET_TICM+results.HY_TICM,1)';
shareOfICM_HYblack(isnan(shareOfICM_HYblack))=0;

shareOfICM_HYimpor = sum(results.HY_TICM,1)'.*(shareOfHyImpor)./sum(results.TICM_GAS_FOS+results.TICM_GAS_REN+results.RPET_TICM+results.HY_TICM,1)';
shareOfICM_HYimpor(isnan(shareOfICM_HYimpor))=0;
%Share of NGA in HHB  fuel




% share of Biogas converted to biomthane
shareOfgasBGA = (sum(results.RBGA_TBGU,1))'./(sum(results.RBGA_TCHP+results.RBGA_COOK+resultsSC.RBGA_FU,1)+sum(results.RBGA_TBGU,1))';
shareOfgasBGA(isnan(shareOfgasBGA))=0;

shareOfNGAInd = sum(results.LIGA_GAS_FOS,1)'./sum(results.RNGA_FU,1)';
shareOfNGAInd(isnan(shareOfNGAInd))=0;
shareOfNGAHy = sum(results.TSMR_GAS_FOS+results.TSMC_GAS_FOS,1)'./sum(results.RNGA_FU,1)';
shareOfNGAHy(isnan(shareOfNGAHy))=0;
shareOfNGAHySMR = sum(results.TSMR_GAS_FOS,1)'./sum(results.RNGA_FU,1)';
shareOfNGAHySMR(isnan(shareOfNGAHySMR))=0;
shareOfNGAHySMC = sum(results.TSMC_GAS_FOS,1)'./sum(results.RNGA_FU,1)';
shareOfNGAHySMC(isnan(shareOfNGAHySMC))=0;
shareOfNGATrans = sum(results.TLNG_GAS_FOS,1)'./sum(results.RNGA_FU,1)';
shareOfNGATrans(isnan(shareOfNGATrans))=0;
shareOfNGAGen = sum(results.RNGA_FU-results.LIGA_GAS_FOS-results.TSMR_GAS_FOS-results.TLNG_GAS_FOS,1)'./sum(results.RNGA_FU,1)';
shareOfNGAGen(isnan(shareOfNGAGen))=0;

shareOfPETTrans = sum(results.RPET_DI+results.RPET_KE,1)'./sum(results.RPET_FU,1)';
shareOfPETTrans(isnan(shareOfPETTrans))=0;


shareOfBDSTrans = sum(results.RBDS_DI+results.RBDS_KE,1)'./sum(results.RBDS_FU,1)';
shareOfBDSTrans(isnan(shareOfBDSTrans))=0;


shareOfElCoalGHG = sum(results.RHAR_THPP+results.RHAR_THCS+results.RHAR_TCCO,1)'./sum(results.RHAR_FU,1)';
shareOfElCoalGHG(isnan(shareOfElCoalGHG))=0;shareOfElCoalGHG(shareOfElCoalGHG>1)=1;shareOfElCoalGHG(shareOfElCoalGHG<0)=0;
shareOfElOilGHG = sum(results.RPET_TICG+results.RPET_TICM+results.RPET_TCOI+results.RPET_TCCG+results.RPET_TOCG+results.RPET_TGCS+results.RPET_TCNG,1)'./sum(results.RPET_FU,1)';
shareOfElOilGHG(isnan(shareOfElOilGHG))=0;shareOfElOilGHG(shareOfElOilGHG>1)=1;shareOfElOilGHG(shareOfElOilGHG<0)=0;
shareOfElGasGHG = sum(results.TCCG_GAS_FOS+results.TOCG_GAS_FOS+results.TCNG_GAS_FOS+results.TGCS_GAS_FOS+results.TICM_GAS_FOS,1)'./sum(results.RNGA_FU,1)';
shareOfElGasGHG(isnan(shareOfElGasGHG))=0;shareOfElGasGHG(shareOfElGasGHG>1)=1;shareOfElGasGHG(shareOfElGasGHG<0)=0;

share.H2_REN = sum(results.TWEL_GA+results.importH2_regas,1)'./sum(results.TWEL_GA+results.TSMR_GA+results.TSMC_GA+results.importH2_regas,1)';
share.H2_REN(isinf(share.H2_REN))=1;share.H2_REN(isnan(share.H2_REN))=0;share.H2_REN((share.H2_REN>1))=1;share.H2_REN((share.H2_REN<0))=0;

share.H2toCH4 = sum(results.HY_TMET,1)'./sum(results.TWEL_GA+results.TSMR_GA+results.TSMC_GA+results.importH2_regas,1)';
share.H2toCH4(isinf(share.H2toCH4))=1;share.H2toCH4(isnan(share.H2toCH4))=0;share.H2toCH4((share.H2toCH4>1))=1;share.H2toCH4((share.H2toCH4<0))=0;
share.H2toElDir = round(sum(results.HY_TCCG+results.HY_TOCG+results.HY_TCNG+results.HY_TICM,1))'./round(sum(results.TWEL_GA+results.TSMR_GA+results.TSMC_GA+results.importH2_regas,1))';
share.H2toElDir(isinf(share.H2toElDir)) = 1; share.H2toElDir(isnan(share.H2toElDir)) = 0;share.H2toElDir((share.H2toElDir>1))=1;share.H2toElDir((share.H2toElDir<0))=0;
share.H2toEL = share.H2toElDir + share.H2toCH4.*shareOfGasEl;
share.H2toLNG = share.H2toCH4.*shareOfGasTrans;
share.H2toLNG(isinf(share.H2toLNG))=1;share.H2toLNG(isnan(share.H2toLNG))=0;share.H2toLNG((share.H2toLNG>1))=1;share.H2toLNG((share.H2toLNG<0))=0;
share.H2toLH2 = sum(results.HY_TLH2,1)'./sum(results.TWEL_GA+results.TSMR_GA+results.TSMC_GA+results.importH2_regas,1)';
share.H2toLH2(isinf(share.H2toLH2))=1;share.H2toLH2(isnan(share.H2toLH2))=0;share.H2toLH2((share.H2toLH2>1))=1;share.H2toLH2((share.H2toLH2<0))=0;
share.H2toFTU = sum(results.HY_TFTU,1)'./sum(results.TWEL_GA+results.TSMR_GA+results.TSMC_GA+results.importH2_regas,1)';
share.H2toFTU(isinf(share.H2toFTU))=1;share.H2toFTU(isnan(share.H2toFTU))=0;share.H2toFTU((share.H2toFTU>1))=1;share.H2toFTU((share.H2toFTU<0))=0;
share.H2toM = sum(results.HY_TRSP,1)'./sum(results.TWEL_GA+results.TSMR_GA+results.TSMC_GA+results.importH2_regas,1)';
share.H2toM(isinf(share.H2toM))=1;share.H2toM(isnan(share.H2toM))=0;share.H2toM((share.H2toM>1))=1;share.H2toM((share.H2toM<0))=0;
share.H2toNH3 = sum(results.HY_TNH3,1)'./sum(results.TWEL_GA+results.TSMR_GA+results.TSMC_GA+results.importH2_regas,1)';
share.H2toNH3(isinf(share.H2toNH3))=1;share.H2toNH3(isnan(share.H2toNH3))=0;share.H2toNH3((share.H2toNH3>1))=1;share.H2toNH3((share.H2toNH3<0))=0;
share.H2toMeO = sum(results.HY_TMeO,1)'./sum(results.TWEL_GA+results.TSMR_GA+results.TSMC_GA+results.importH2_regas,1)';
share.H2toMeO(isinf(share.H2toMeO))=1;share.H2toMeO(isnan(share.H2toMeO))=0;share.H2toMeO((share.H2toMeO>1))=1;share.H2toMeO((share.H2toMeO<0))=0;
share.H2toISH = sum(results.HY_TISH,1)'./sum(results.TWEL_GA+results.TSMR_GA+results.TSMC_GA+results.importH2_regas,1)';
share.H2toISH(isinf(share.H2toISH))=1;share.H2toISH(isnan(share.H2toISH))=0;share.H2toISH((share.H2toISH>1))=1;share.H2toISH((share.H2toISH<0))=0;
share.H2toINH = sum(results.HY_LHIN,1)'./sum(results.TWEL_GA+results.TSMR_GA+results.TSMC_GA+results.importH2_regas,1)';
share.H2toINH(isinf(share.H2toINH))=1;share.H2toINH(isnan(share.H2toINH))=0;share.H2toINH((share.H2toINH>1))=1;share.H2toINH((share.H2toINH<0))=0;
share.H2toExport = sum(results.importH2,1)'./sum(results.TWEL_GA+results.TSMR_GA+results.TSMC_GA+results.importH2_regas,1)';
share.H2toExport(isinf(share.H2toExport))=1;share.H2toExport(isnan(share.H2toExport))=0;share.H2toExport((share.H2toExport>1))=1;share.H2toExport((share.H2toExport<0))=0;
share.H2toCDR = sum(results.HY_TRSi,1)'./sum(results.TWEL_GA+results.TSMR_GA+results.TSMC_GA+results.importH2_regas,1)';
share.H2toCDR(isinf(share.H2toCDR))=1;share.H2toCDR(isnan(share.H2toCDR))=0;share.H2toCDR((share.H2toCDR>1))=1;share.H2toCDR((share.H2toCDR<0))=0;

share.CO2toCH4 = sum(results.CO_TMET,1)'./sum(results.TCOS_GA+results.TPSP_GA+results.TBCS_CO+results.TCBC_CO+results.TPSC_GA+results.TCWC_CO+results.THCS_CO+results.TGCS_CO+results.TSMC_CO,1)';
share.CO2toCH4(isinf(share.CO2toCH4))=1;share.CO2toCH4(isnan(share.CO2toCH4))=0;share.CO2toCH4((share.CO2toCH4>1))=1;share.CO2toCH4((share.CO2toCH4<0))=0;
share.CO2toLNG = share.CO2toCH4.*shareOfGasTrans;
share.CO2toLNG(isinf(share.CO2toLNG))=1;share.CO2toLNG(isnan(share.CO2toLNG))=0;share.CO2toLNG((share.CO2toLNG>1))=1;share.CO2toLNG((share.CO2toLNG<0))=0;
share.CO2toFTU = sum(results.CO_TFTU,1)'./sum(results.TCOS_GA+results.TPSP_GA+results.TBCS_CO+results.TCBC_CO+results.TPSC_GA+results.TCWC_CO+results.THCS_CO+results.TGCS_CO+results.TSMC_CO,1)';
share.CO2toFTU(isinf(share.CO2toFTU))=1;share.CO2toFTU(isnan(share.CO2toFTU))=0;share.CO2toFTU((share.CO2toFTU>1))=1;share.CO2toFTU((share.CO2toFTU<0))=0;
share.CO2toMeO = sum(results.CO_TMeO,1)'./sum(results.TCOS_GA+results.TPSP_GA+results.TBCS_CO+results.TCBC_CO+results.TPSC_GA+results.TCWC_CO+results.THCS_CO+results.TGCS_CO+results.TSMC_CO,1)';
share.CO2toMeO(isinf(share.CO2toMeO))=1;share.CO2toMeO(isnan(share.CO2toMeO))=0;share.CO2toMeO((share.CO2toMeO>1))=1;share.CO2toMeO((share.CO2toMeO<0))=0;
share.CO2toCDR = sum(results.CO_LCCS,1)'./sum(results.TCOS_GA+results.TPSP_GA+results.TBCS_CO+results.TCBC_CO+results.TPSC_GA+results.TCWC_CO+results.THCS_CO+results.TGCS_CO+results.TSMC_CO,1)';
share.CO2toCDR(isinf(share.CO2toCDR))=1;share.CO2toCDR(isnan(share.CO2toCDR))=0;share.CO2toCDR((share.CO2toCDR>1))=1;share.CO2toCDR((share.CO2toCDR<0))=0;

share.MeOtoTransp = repmat(sum(results.MeOtoMobility,'all')'./sum(results.MeOtoMobility+results.ME_LIME+results.exportMeOH,'all'),size(systemParams.IndexNodes,1),1);
share.MeOtoTransp(isinf(share.MeOtoTransp))=1;share.MeOtoTransp(isnan(share.MeOtoTransp))=0;share.MeOtoTransp((share.MeOtoTransp>1))=1;share.MeOtoTransp((share.MeOtoTransp<0))=0;
share.MeOtoIndustry = sum(results.ME_LIME,1)'./sum(results.MeOtoMobility+results.ME_LIME+results.exportMeOH,1)';
share.MeOtoIndustry = repmat(sum(results.ME_LIME,'all')'./sum(results.MeOtoMobility+results.ME_LIME+results.exportMeOH,'all'),size(systemParams.IndexNodes,1),1);
share.MeOtoIndustry(isinf(share.MeOtoIndustry))=1;share.MeOtoIndustry(isnan(share.MeOtoIndustry))=0;share.MeOtoIndustry((share.MeOtoIndustry>1))=1;share.MeOtoIndustry((share.MeOtoIndustry<0))=0;
share.MeOtoExport = sum(results.exportMeOH,1)'./sum(results.MeOtoMobility+results.ME_LIME+results.exportMeOH,1)';
share.MeOtoExport = repmat(sum(results.exportMeOH,'all')'./sum(results.MeOtoMobility+results.ME_LIME+results.exportMeOH,'all'),size(systemParams.IndexNodes,1),1);
share.MeOtoExport(isinf(share.MeOtoExport))=1;share.MeOtoExport(isnan(share.MeOtoExport))=0;share.MeOtoExport((share.MeOtoExport>1))=1;share.MeOtoExport((share.MeOtoExport<0))=0;

share.NH3toTransp = sum(results.NH3toMobility,1)'./sum(results.NH3toMobility+results.NH_LINH+results.exportNH3,1)';
share.NH3toTransp = repmat(sum(results.NH3toMobility,'all')'./sum(results.NH3toMobility+results.NH_LINH+results.exportNH3,'all'),size(systemParams.IndexNodes,1),1);
share.NH3toTransp(isinf(share.NH3toTransp))=1;share.NH3toTransp(isnan(share.NH3toTransp))=0;share.NH3toTransp((share.NH3toTransp>1))=1;share.NH3toTransp((share.NH3toTransp<0))=0;
share.NH3toIndustry = sum(results.NH_LINH,1)'./sum(results.NH3toMobility+results.NH_LINH+results.exportNH3,1)';
share.NH3toIndustry = repmat(sum(results.NH_LINH,'all')'./sum(results.NH3toMobility+results.NH_LINH+results.exportNH3,'all'),size(systemParams.IndexNodes,1),1);
share.NH3toIndustry(isinf(share.NH3toIndustry))=1;share.NH3toIndustry(isnan(share.NH3toIndustry))=0;share.NH3toIndustry((share.NH3toIndustry>1))=1;share.NH3toIndustry((share.NH3toIndustry<0))=0;
share.NH3toExport = sum(results.exportNH3,1)'./sum(results.NH3toMobility+results.NH_LINH+results.exportNH3,1)';
share.NH3toExport = repmat(sum(results.exportNH3,'all')'./sum(results.NH3toMobility+results.NH_LINH+results.exportNH3,'all'),size(systemParams.IndexNodes,1),1);
share.NH3toExport(isinf(share.NH3toExport))=1;share.NH3toExport(isnan(share.NH3toExport))=0;share.NH3toExport((share.NH3toExport>1))=1;share.NH3toExport((share.NH3toExport<0))=0;

try results.FTdiesel;
catch
    results.FTdiesel = results.TFTU_DI - results.exportFT_diesel;
    results.FTkerosene = results.TFTU_KE - results.exportFT_kerosene;
end


share.FTtoTransp = sum(results.FTdiesel+results.FTkerosene,'all')'./sum(results.TFTU_GA,'all')';
share.FTtoTransp(isinf(share.FTtoTransp))=1;share.FTtoTransp(isnan(share.FTtoTransp))=0;share.FTtoTransp((share.FTtoTransp>1))=1;share.FTtoTransp((share.FTtoTransp<0))=0;
share.FTtoIndustry = sum(results.TFTU_GA-results.FTdiesel-results.FTkerosene-results.exportFT,'all')'./sum(results.TFTU_GA,'all')';
share.FTtoIndustry(isinf(share.FTtoIndustry))=1;share.FTtoIndustry(isnan(share.FTtoIndustry))=0;share.FTtoIndustry((share.FTtoIndustry>1))=1;share.FTtoIndustry((share.FTtoIndustry<0))=0;
share.FTtoExport = sum(results.exportFT,'all')'./sum(results.TFTU_GA,'all')';
share.FTtoExport(isinf(share.FTtoExport))=1;share.FTtoExport(isnan(share.FTtoExport))=0;share.FTtoExport((share.FTtoExport>1))=1;share.FTtoExport((share.FTtoExport<0))=0;


share.NGAtoInd = sum(results.LIGA_GAS_FOS,1)'./sum(results.RNGA_FU,1)';
share.NGAtoInd(isnan(share.NGAtoInd))=0;
share.NGAtoHy = sum(results.TSMR_GAS_FOS+results.TSMC_GAS_FOS,1)'./sum(results.RNGA_FU,1)';
share.NGAtoHy(isnan(share.NGAtoHy))=0;
share.NGAtoLNG = sum(results.TLNG_GAS_FOS,1)'./sum(results.RNGA_FU,1)';
share.NGAtoLNG(isnan(share.NGAtoLNG))=0;
share.NGAtoEl = sum(results.TCCG_GAS_FOS+results.TOCG_GAS_FOS+results.TCNG_GAS_FOS+results.TGCS_GAS_FOS+results.TICM_GAS_FOS,1)'./sum(results.RNGA_FU,1)';
share.NGAtoEl(isnan(share.NGAtoEl))=0;
share.NGAtoHe = sum(results.TDNG_GAS_FOS+results.Pros_GAS_FOS,1)'./sum(results.RNGA_FU,1)';
share.NGAtoHe(isnan(share.NGAtoHe))=0;

share.COALtoEl = sum(results.RHAR_THPP+results.RHAR_THCS+results.RHAR_TCCO,1)'./sum(results.RHAR_FU,1)';
share.COALtoEl(isnan(share.COALtoEl))=0;
share.COALtoHe = sum(results.RHAR_TDCO,1)'./sum(results.RHAR_FU,1)';
share.COALtoHe(isnan(share.COALtoHe))=0;

share.PETtoTr = sum(results.RPET_DI+results.RPET_KE,1)'./sum(results.RPET_FU,1)';
share.PETtoTr(isnan(share.PETtoTr))=0;
share.PETtoEl = sum(results.RPET_TICG+results.RPET_TICM+results.RPET_TCOI+results.RPET_TCCG+results.RPET_TOCG+results.RPET_TGCS+results.RPET_TCNG,1)'./sum(results.RPET_FU,1)';
share.PETtoEl(isnan(share.PETtoEl))=0;
share.PETtoHe = sum(results.RPET_TDOI+resultsSC.RPET_FU,1)'./sum(results.RPET_FU,1)';
share.PETtoHe(isnan(share.PETtoHe))=0;

share.TICMfromHY = sum(results.HY_TICM,1)'./sum(results.RPET_TICM + results.TICM_GAS_FOS + results.TICM_GAS_REN + results.HY_TICM,1)';
share.TICMfromPET = sum(results.RPET_TICM,1)'./sum(results.RPET_TICM + results.TICM_GAS_FOS + results.TICM_GAS_REN + results.HY_TICM,1)';
share.TICMfromGA = sum(results.TICM_GAS_FOS + results.TICM_GAS_REN,1)'./sum(results.RPET_TICM + results.TICM_GAS_FOS + results.TICM_GAS_REN + results.HY_TICM,1)';
share.TICMfromHY(isnan(share.TICMfromHY))=0;
share.TICMfromPET(isnan(share.TICMfromPET))=0;
share.TICMfromGA(isnan(share.TICMfromGA))=0;

share.TCCGfromHY = sum(results.HY_TCCG,1)'./sum(results.RPET_TCCG + results.TCCG_GAS_FOS + results.TCCG_GAS_REN + results.HY_TCCG,1)';
share.TCCGfromPET = sum(results.RPET_TCCG,1)'./sum(results.RPET_TCCG + results.TCCG_GAS_FOS + results.TCCG_GAS_REN + results.HY_TCCG,1)';
share.TCCGfromGA = sum(results.TCCG_GAS_FOS + results.TCCG_GAS_REN,1)'./sum(results.RPET_TCCG + results.TCCG_GAS_FOS + results.TCCG_GAS_REN + results.HY_TCCG,1)';
share.TCCGfromHY(isnan(share.TCCGfromHY))=0;
share.TCCGfromPET(isnan(share.TCCGfromPET))=0;
share.TCCGfromGA(isnan(share.TCCGfromGA))=0;

share.TOCGfromHY = sum(results.HY_TOCG,1)'./sum(results.RPET_TOCG + results.TOCG_GAS_FOS + results.TOCG_GAS_REN + results.HY_TOCG,1)';
share.TOCGfromPET = sum(results.RPET_TOCG,1)'./sum(results.RPET_TOCG + results.TOCG_GAS_FOS + results.TOCG_GAS_REN + results.HY_TOCG,1)';
share.TOCGfromGA = sum(results.TOCG_GAS_FOS + results.TOCG_GAS_REN,1)'./sum(results.RPET_TOCG + results.TOCG_GAS_FOS + results.TOCG_GAS_REN + results.HY_TOCG,1)';
share.TOCGfromHY(isnan(share.TOCGfromHY))=0;
share.TOCGfromPET(isnan(share.TOCGfromPET))=0;
share.TOCGfromGA(isnan(share.TOCGfromGA))=0;

share.TGCSfromHY = sum(results.HY_TGCS,1)'./sum(results.RPET_TGCS + results.TGCS_GAS_FOS + results.TGCS_GAS_REN + results.HY_TGCS,1)';
share.TGCSfromPET = sum(results.RPET_TGCS,1)'./sum(results.RPET_TGCS + results.TGCS_GAS_FOS + results.TGCS_GAS_REN + results.HY_TGCS,1)';
share.TGCSfromGA = sum(results.TGCS_GAS_FOS + results.TGCS_GAS_REN,1)'./sum(results.RPET_TGCS + results.TGCS_GAS_FOS + results.TGCS_GAS_REN + results.HY_TGCS,1)';
share.TGCSfromHY(isnan(share.TGCSfromHY))=0;
share.TGCSfromPET(isnan(share.TGCSfromPET))=0;
share.TGCSfromGA(isnan(share.TGCSfromGA))=0;

share.TCNGfromHY = sum(results.HY_TCNG,1)'./sum(results.RPET_TCNG + results.TCNG_GAS_FOS + results.TCNG_GAS_REN + results.HY_TCNG,1)';
share.TCNGfromPET = sum(results.RPET_TCNG,1)'./sum(results.RPET_TCNG + results.TCNG_GAS_FOS + results.TCNG_GAS_REN + results.HY_TCNG,1)';
share.TCNGfromGA = sum(results.TCNG_GAS_FOS + results.TCNG_GAS_REN,1)'./sum(results.RPET_TCNG + results.TCNG_GAS_FOS + results.TCNG_GAS_REN + results.HY_TCNG,1)';
share.TCNGfromHY(isnan(share.TCNGfromHY))=0;
share.TCNGfromPET(isnan(share.TCNGfromPET))=0;
share.TCNGfromGA(isnan(share.TCNGfromGA))=0;

% split CSP and TSTU
results.OPT_SIZE_RCSP(results.OPT_SIZE_RCSP<0) = 0;
results.OPT_SIZE_TSTU(results.OPT_SIZE_TSTU<0) = 0;



share.CSPtoEl = (sum(results.HE_TSTU)./sum(results.RCSP_HE))';
share.CSPtoEl(share.CSPtoEl > 1) = 1;
share.CSPtoHe = 1 - share.CSPtoEl;

share.STUtoHtP = (results.OPT_SIZE_TSTU-results.OPT_SIZE_RCSP*systemParams.EtaTrans(ismember(systemParams.IndexIDT,'TSTU'),ismember(systemParams.IndexFuels,'FHOT')))'./(results.OPT_SIZE_TSTU+0.0001)';
share.STUtoHtP(share.STUtoHtP<0)=0;
share.STUtoCSP = 1 - share.STUtoHtP;



%% test

shareOfSynLigToTransp = (sum(results.TFTU_DI+results.TFTU_KE,1)'./sum(results.TFTU_DI+results.TFTU_KE,1)');
shareOfSynLigToTransp(isnan(shareOfSynLigToTransp)) = 0;

demand_transp = sum(results.EltoMobility,1)'.*(1./(1-systemParams.AC_losses(:,systemParams.IndexYears==costYear)/100))+sum(results.EL_TLH2,1)'+sum(results.EL_TFTU,1)'+sum(results.EL_TLNG,1)'+sum(results.EL_TMET,1)'.*shareOfGasTrans +sum(results.EL_TMeO,1)'.*share.MeOtoTransp+sum(results.EL_TNH3,1)'.*share.NH3toTransp +sum(results.EL_TWEL,1)'.*(share.H2toM+share.H2toFTU*share.FTtoTransp+share.H2toLNG+share.H2toLH2+share.H2toNH3.*share.NH3toTransp+share.H2toMeO.*share.MeOtoTransp)+sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*(share.CO2toFTU+share.CO2toLNG+share.CO2toMeO.*share.MeOtoTransp);

demand_export = sum(results.exportEL,1)'+sum(results.EL_TMET,1)'.*shareOfGasExport+sum(results.EL_TMeO,1)'.*share.MeOtoExport+sum(results.EL_TNH3,1)'.*share.NH3toExport +sum(results.EL_TWEL,1)'.*(share.H2toExport + share.H2toCH4.*shareOfGasExport +share.H2toFTU*share.FTtoExport+share.H2toNH3.*share.NH3toExport+share.H2toMeO.*share.MeOtoExport)+sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*(share.CO2toFTU*share.FTtoExport+share.CO2toCH4.*shareOfGasExport+share.CO2toMeO.*share.MeOtoExport);
heatProdElCons = (sum(results.EL_TDHR+results.EL_TDHP+results.EL_TDGE,1)'.*(1./(1-systemParams.AC_losses(:,systemParams.IndexYears==costYear)/100))+sum(results.EL_TMET,1)'.*(shareOfGasDH+shareOfGasIH+shareOfGasINH)+sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*share.CO2toCH4.*(shareOfGasDH+shareOfGasIH+shareOfGasINH)+sum(results.EL_TWEL,1)'.*(share.H2toCH4.*(shareOfGasDH+shareOfGasIH+shareOfGasINH)+share.H2toINH))*setup.Heat.Flag;
gasProdElCons = (sum(results.EL_TMET,1)'.*(shareOfGasInd)+sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*share.CO2toCH4.*(shareOfGasInd)+sum(results.EL_TWEL,1)'.*share.H2toCH4.*(shareOfGasInd)).*shareOfGasInd;%.*GasFlag;

demand_CDR = sum(results.EL_TRSS+results.EL_TRSL+results.EL_TRMI+results.EL_TRMO+results.EL_TRME+results.EL_TRSi+results.EL_TRAD+results.EL_TREW+results.EL_TRBC+results.EL_TRGE,1)'+...
    sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*(share.CO2toCDR)+...
    sum(results.EL_TWEL,1)'.*(share.H2toCDR);

[sum(results.EL_TRSS+results.EL_TRSL+results.EL_TRMI+results.EL_TRMO+results.EL_TRME+results.EL_TRSi+results.EL_TRAD+results.EL_TREW+results.EL_TRBC+results.EL_TRGE,1)',...
    sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*(share.CO2toCDR),...
    sum(results.EL_TWEL,1)'.*(share.H2toCDR)]./demand_CDR;

demand_El_Heat = heatProdElCons+sum(resultsSC.EL_THHP+resultsSC.EL_THHR)'-el_cons_localheatFromGrid+el_cons_localheatFromGrid./((100-systemParams.AC_losses(:,yearNumb))/100)*setup.Heat.Flag;

demand_industry = (sum(results.EL_TICB,1)'+sum(results.EL_TICI,1)'+sum(results.EL_TISB,1)'+sum(results.EL_TISH,1)'+sum(results.EL_TISR,1)'+sum(results.EL_TISE,1)'+sum(results.EL_TIAM,1)'+sum(results.EL_TIAR,1)'+sum(results.EL_TIPP,1)').*(1./(1-systemParams.AC_losses(:,systemParams.IndexYears==costYear)/100))+sum(results.EL_TMeO,1)'.*share.MeOtoIndustry+sum(results.EL_TNH3,1)'.*share.NH3toIndustry+sum(results.EL_TWEL,1)'.*(share.H2toNH3.*share.NH3toIndustry+share.H2toMeO.*share.MeOtoIndustry)+sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*(share.CO2toMeO.*share.MeOtoIndustry);
heat_demand_industry = sum(results.HE_TICB,1)'+sum(results.HE_TICI,1)'+sum(results.HE_TISB,1)'+sum(results.HE_TISH,1)'+sum(results.HE_TISR,1)'+sum(results.HE_TISE,1)'+sum(results.HE_TIAA,1)'+sum(results.EL_TIAR,1)'+sum(results.HE_TIPP,1)';

stor_electricityOUT_tot =  (sum(results.SBAT_EL,1)'+sum(results.SPHS_EL,1)'+sum(results.SACA_EL,1)'+sum(results.TCCG_EL,1)'+sum(results.TOCG_EL,1)'+sum(results.TSTU_EL,1)');

results.GRID = results.GRID_AC + results.GRID_DC;



storLosses_sys = (sum(results.EL_THyS,1)'+sum(results.EL_SACA,1)'-sum(results.SACA_EL,1)'+sum(results.EL_SBAT,1)'-sum(results.SBAT_EL,1)'+sum(results.EL_SPHS,1)'-sum(results.SPHS_EL,1)'+((sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*share.CO2toCH4+sum(results.EL_TWEL,1)'.*share.H2toCH4).*shareOfGasEl + sum(results.EL_TMET,1)'.*shareOfGasEl - sum(results.TCCG_EL+results.TOCG_EL+results.TCNG_EL+results.TICM_EL,1)'.*(shareOfGT_SNG)) +(sum(results.EL_TWEL,1)'.*share.H2toElDir-sum(results.TCCG_EL+results.TOCG_EL+results.TCNG_EL+results.TICM_EL,1)'.*(shareOfGT_HYgreen)) + (sum(results.EL_TDHR-results.TSTU_EL,1)'.*(sum(results.EL_TDHR-results.TSTU_EL,1)'>0))*(~setup.Heat.Flag))+...
    sum(results.EL_SLD2,1)'-sum(results.SLD2_EL,1)'+sum(results.EL_SW22,1)'-sum(results.SW22_EL,1)'+sum(results.EL_SBU2,1)'-sum(results.SBU2_EL,1)'+sum(results.EL_SMD2,1)'-sum(results.SMD2_EL,1)'+sum(results.EL_SHD2,1)'-sum(results.SHD2_EL,1)'+...
    (sum(results.EL_SLD2,1)'+sum(results.EL_SW22,1)'+sum(results.EL_SBU2,1)'+sum(results.EL_SMD2,1)'+sum(results.EL_SHD2,1)')./((100-systemParams.AC_losses(:,yearNumb))/100);


%%



prod_electricity_SC = sum(resultsSC.RES.RPVO_EL,1)'+ sum(resultsSC.COM.RPVO_EL,1)'+ sum(resultsSC.IND.RPVO_EL,1)';
storLosses_SC = sum(resultsSC.RES.EL_SBAT,1)'+ sum(resultsSC.COM.EL_SBAT,1)'+ sum(resultsSC.IND.EL_SBAT,1)'-sum(resultsSC.RES.SBAT_EL,1)'-sum(resultsSC.COM.SBAT_EL,1)'-sum(resultsSC.IND.SBAT_EL,1)';

resultsSC.Dem = resultsSC.RES.Dem + resultsSC.COM.Dem + resultsSC.IND.Dem;

RES.CapexCRF_RPVR = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'RPVR'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'RPVR'),:,:),[2 3 1]).*repmat(CapitalRecoveryFactor(systemParams.WACC_SC(1),systemParams.Lifetime(ismember(systemParams.IndexID,'RPVR'),:))',length(systemParams.IndexNodes),1,1));
COM.CapexCRF_RPVC = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'RPVC'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'RPVC'),:,:),[2 3 1]).*repmat(CapitalRecoveryFactor(systemParams.WACC_SC(2),systemParams.Lifetime(ismember(systemParams.IndexID,'RPVC'),:))',length(systemParams.IndexNodes),1,1));
IND.CapexCRF_RPVI = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'RPVI'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'RPVI'),:,:),[2 3 1]).*repmat(CapitalRecoveryFactor(systemParams.WACC_SC(3),systemParams.Lifetime(ismember(systemParams.IndexID,'RPVI'),:))',length(systemParams.IndexNodes),1,1));
RES.CapexCRF_SBAT = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'SBAR'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'SBAR'),:,:),[2 3 1]).*repmat(CapitalRecoveryFactor(systemParams.WACC_SC(1),systemParams.Lifetime(ismember(systemParams.IndexID,'SBAR'),:))',length(systemParams.IndexNodes),1,1));
COM.CapexCRF_SBAT = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'SBAC'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'SBAC'),:,:),[2 3 1]).*repmat(CapitalRecoveryFactor(systemParams.WACC_SC(2),systemParams.Lifetime(ismember(systemParams.IndexID,'SBAC'),:))',length(systemParams.IndexNodes),1,1));
IND.CapexCRF_SBAT = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'SBAI'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'SBAI'),:,:),[2 3 1]).*repmat(CapitalRecoveryFactor(systemParams.WACC_SC(3),systemParams.Lifetime(ismember(systemParams.IndexID,'SBAI'),:))',length(systemParams.IndexNodes),1,1));
RES.CapexCRF_IBAT = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'IBAR'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'IBAR'),:,:),[2 3 1]).*repmat(CapitalRecoveryFactor(systemParams.WACC_SC(1),systemParams.Lifetime(ismember(systemParams.IndexID,'IBAR'),:))',length(systemParams.IndexNodes),1,1));
COM.CapexCRF_IBAT = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'IBAC'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'IBAC'),:,:),[2 3 1]).*repmat(CapitalRecoveryFactor(systemParams.WACC_SC(2),systemParams.Lifetime(ismember(systemParams.IndexID,'IBAC'),:))',length(systemParams.IndexNodes),1,1));
IND.CapexCRF_IBAT = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'IBAI'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'IBAI'),:,:),[2 3 1]).*repmat(CapitalRecoveryFactor(systemParams.WACC_SC(3),systemParams.Lifetime(ismember(systemParams.IndexID,'IBAI'),:))',length(systemParams.IndexNodes),1,1));

RES.Capex_RPVR = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'RPVR'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'RPVR'),:,:),[2 3 1]));
COM.Capex_RPVC = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'RPVC'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'RPVC'),:,:),[2 3 1]));
IND.Capex_RPVI = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'RPVI'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'RPVI'),:,:),[2 3 1]));
RES.Capex_SBAT = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'SBAR'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'SBAR'),:,:),[2 3 1]));
COM.Capex_SBAT = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'SBAC'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'SBAC'),:,:),[2 3 1]));
IND.Capex_SBAT = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'SBAI'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'SBAI'),:,:),[2 3 1]));
RES.Capex_IBAT = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'IBAR'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'IBAR'),:,:),[2 3 1]));
COM.Capex_IBAT = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'IBAC'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'IBAC'),:,:),[2 3 1]));
IND.Capex_IBAT = (permute(repmat(systemParams.Capex(ismember(systemParams.IndexID,'IBAI'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Capex_reg(ismember(systemParams.IndexID,'IBAI'),:,:),[2 3 1]));

RES.Opex_fix_RPVR = (permute(repmat(systemParams.Opex_fix(ismember(systemParams.IndexID,'RPVR'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_fix_reg(ismember(systemParams.IndexID,'RPVR'),:,:),[2 3 1]));
COM.Opex_fix_RPVC = (permute(repmat(systemParams.Opex_fix(ismember(systemParams.IndexID,'RPVC'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_fix_reg(ismember(systemParams.IndexID,'RPVC'),:,:),[2 3 1]));
IND.Opex_fix_RPVI = (permute(repmat(systemParams.Opex_fix(ismember(systemParams.IndexID,'RPVI'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_fix_reg(ismember(systemParams.IndexID,'RPVI'),:,:),[2 3 1]));
RES.Opex_fix_SBAT = (permute(repmat(systemParams.Opex_fix(ismember(systemParams.IndexID,'SBAR'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_fix_reg(ismember(systemParams.IndexID,'SBAR'),:,:),[2 3 1]));
COM.Opex_fix_SBAT = (permute(repmat(systemParams.Opex_fix(ismember(systemParams.IndexID,'SBAC'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_fix_reg(ismember(systemParams.IndexID,'SBAC'),:,:),[2 3 1]));
IND.Opex_fix_SBAT = (permute(repmat(systemParams.Opex_fix(ismember(systemParams.IndexID,'SBAI'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_fix_reg(ismember(systemParams.IndexID,'SBAI'),:,:),[2 3 1]));
RES.Opex_fix_IBAT = (permute(repmat(systemParams.Opex_fix(ismember(systemParams.IndexID,'IBAR'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_fix_reg(ismember(systemParams.IndexID,'IBAR'),:,:),[2 3 1]));
COM.Opex_fix_IBAT = (permute(repmat(systemParams.Opex_fix(ismember(systemParams.IndexID,'IBAC'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_fix_reg(ismember(systemParams.IndexID,'IBAC'),:,:),[2 3 1]));
IND.Opex_fix_IBAT = (permute(repmat(systemParams.Opex_fix(ismember(systemParams.IndexID,'IBAI'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_fix_reg(ismember(systemParams.IndexID,'IBAI'),:,:),[2 3 1]));

RES.Opex_var_RPVR = (permute(repmat(systemParams.Opex_var(ismember(systemParams.IndexID,'RPVO'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_var_reg(ismember(systemParams.IndexID,'RPVR'),:,:),[2 3 1]));
COM.Opex_var_RPVC = (permute(repmat(systemParams.Opex_var(ismember(systemParams.IndexID,'RPVC'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_var_reg(ismember(systemParams.IndexID,'RPVC'),:,:),[2 3 1]));
IND.Opex_var_RPVI = (permute(repmat(systemParams.Opex_var(ismember(systemParams.IndexID,'RPVO'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_var_reg(ismember(systemParams.IndexID,'RPVO'),:,:),[2 3 1]));
RES.Opex_var_SBAT = (permute(repmat(systemParams.Opex_var(ismember(systemParams.IndexID,'SBAR'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_var_reg(ismember(systemParams.IndexID,'SBAR'),:,:),[2 3 1]));
COM.Opex_var_SBAT = (permute(repmat(systemParams.Opex_var(ismember(systemParams.IndexID,'SBAC'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_var_reg(ismember(systemParams.IndexID,'SBAC'),:,:),[2 3 1]));
IND.Opex_var_SBAT = (permute(repmat(systemParams.Opex_var(ismember(systemParams.IndexID,'SBAI'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_var_reg(ismember(systemParams.IndexID,'SBAI'),:,:),[2 3 1]));
RES.Opex_var_IBAT = (permute(repmat(systemParams.Opex_var(ismember(systemParams.IndexID,'IBAR'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_var_reg(ismember(systemParams.IndexID,'IBAR'),:,:),[2 3 1]));
COM.Opex_var_IBAT = (permute(repmat(systemParams.Opex_var(ismember(systemParams.IndexID,'IBAC'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_var_reg(ismember(systemParams.IndexID,'IBAC'),:,:),[2 3 1]));
IND.Opex_var_IBAT = (permute(repmat(systemParams.Opex_var(ismember(systemParams.IndexID,'IBAI'),:),1,1,size(systemParams.IndexNodes,1)),[3 2 1]).*permute(systemParams.Opex_var_reg(ismember(systemParams.IndexID,'IBAI'),:,:),[2 3 1]));

genShareRPVR = resultsSC.Cap_RPVR./repmat(resultsSC.OPT_SIZE_RPVR',1,size(systemParams.IndexYears,2));
genShareRPVR(isnan(genShareRPVR))=0;
genShareRPVC = resultsSC.Cap_RPVC./repmat(resultsSC.OPT_SIZE_RPVC',1,size(systemParams.IndexYears,2));
genShareRPVC(isnan(genShareRPVC))=0;
genShareRPVI = resultsSC.Cap_RPVI./repmat(resultsSC.OPT_SIZE_RPVI',1,size(systemParams.IndexYears,2));
genShareRPVI(isnan(genShareRPVI))=0;

genShareSBAR = resultsSC.Cap_SBAR./repmat(resultsSC.OPT_SIZE_SBAR',1,size(systemParams.IndexYears,2));
genShareSBAR(isnan(genShareSBAR))=0;
genShareSBAC = resultsSC.Cap_SBAC./repmat(resultsSC.OPT_SIZE_SBAC',1,size(systemParams.IndexYears,2));
genShareSBAC(isnan(genShareSBAC))=0;
genShareSBAI = resultsSC.Cap_SBAI./repmat(resultsSC.OPT_SIZE_SBAI',1,size(systemParams.IndexYears,2));
genShareSBAI(isnan(genShareSBAI))=0;


results.RES.PV_Annual = sum(resultsSC.Cap_RPVR.*(RES.CapexCRF_RPVR+RES.Opex_fix_RPVR)+(repmat(sum(resultsSC.RPVR_EL)',1,size(systemParams.IndexYears,2)).*genShareRPVR).*RES.Opex_var_RPVR,2);
results.COM.PV_Annual = sum(resultsSC.Cap_RPVC.*(COM.CapexCRF_RPVC+COM.Opex_fix_RPVC)+(repmat(sum(resultsSC.RPVC_EL)',1,size(systemParams.IndexYears,2)).*genShareRPVC).*COM.Opex_var_RPVC,2);
results.IND.PV_Annual = sum(resultsSC.Cap_RPVI.*(IND.CapexCRF_RPVI+IND.Opex_fix_RPVI)+(repmat(sum(resultsSC.RPVI_EL)',1,size(systemParams.IndexYears,2)).*genShareRPVI).*IND.Opex_var_RPVI,2);
results.RES.Bat_Annual = sum(resultsSC.Cap_SBAR.*(RES.CapexCRF_SBAT+RES.Opex_fix_SBAT)+resultsSC.Cap_IBAR.*(RES.CapexCRF_IBAT+RES.Opex_fix_IBAT)+(repmat(sum(resultsSC.SBAR_EL)',1,size(systemParams.IndexYears,2)).*genShareSBAR).*RES.Opex_var_SBAT,2);
results.COM.Bat_Annual = sum(resultsSC.Cap_SBAC.*(COM.CapexCRF_SBAT+COM.Opex_fix_SBAT)+resultsSC.Cap_IBAC.*(COM.CapexCRF_IBAT+COM.Opex_fix_IBAT)+(repmat(sum(resultsSC.SBAC_EL)',1,size(systemParams.IndexYears,2)).*genShareSBAC).*COM.Opex_var_SBAT,2);
results.IND.Bat_Annual = sum(resultsSC.Cap_SBAI.*(IND.CapexCRF_SBAT+IND.Opex_fix_SBAT)+resultsSC.Cap_IBAI.*(IND.CapexCRF_IBAT+IND.Opex_fix_IBAT)+(repmat(sum(resultsSC.SBAI_EL)',1,size(systemParams.IndexYears,2)).*genShareSBAI).*IND.Opex_var_SBAT,2);

results.RES.Annual=results.RES.PV_Annual+results.RES.Bat_Annual;
results.COM.Annual=results.COM.PV_Annual+results.COM.Bat_Annual;
results.IND.Annual=results.IND.PV_Annual+results.IND.Bat_Annual;

results.RES.PV_Capex = sum(resultsSC.Cap_RPVR.*(RES.Capex_RPVR),2);
results.COM.PV_Capex = sum(resultsSC.Cap_RPVC.*(COM.Capex_RPVC),2);
results.IND.PV_Capex = sum(resultsSC.Cap_RPVI.*(IND.Capex_RPVI),2);
results.RES.Bat_Capex = sum(resultsSC.Cap_SBAR.*(RES.Capex_SBAT),2)+sum(resultsSC.Cap_IBAR.*(RES.Capex_IBAT),2);
results.COM.Bat_Capex = sum(resultsSC.Cap_SBAC.*(COM.Capex_SBAT),2)+sum(resultsSC.Cap_IBAC.*(COM.Capex_IBAT),2);
results.IND.Bat_Capex = sum(resultsSC.Cap_SBAI.*(IND.Capex_SBAT),2)+sum(resultsSC.Cap_IBAI.*(IND.Capex_IBAT),2);

results.RES.Capex=results.RES.PV_Capex+results.RES.Bat_Capex;
results.COM.Capex=results.COM.PV_Capex+results.COM.Bat_Capex;
results.IND.Capex=results.IND.PV_Capex+results.IND.Bat_Capex;

results.RES.CapexNew=sum(resultsSC.Cap_RPVR(:,find(systemParams.IndexYears==costYear)).*(RES.Capex_RPVR(:,find(systemParams.IndexYears==costYear))),2)+sum(resultsSC.Cap_SBAR(:,find(systemParams.IndexYears==costYear)).*(RES.Capex_SBAT(:,find(systemParams.IndexYears==costYear))),2)+sum(resultsSC.Cap_IBAR(:,find(systemParams.IndexYears==costYear)).*(RES.Capex_IBAT(:,find(systemParams.IndexYears==costYear))),2);
results.COM.CapexNew=sum(resultsSC.Cap_RPVC(:,find(systemParams.IndexYears==costYear)).*(COM.Capex_RPVC(:,find(systemParams.IndexYears==costYear))),2)+sum(resultsSC.Cap_SBAC(:,find(systemParams.IndexYears==costYear)).*(COM.Capex_SBAT(:,find(systemParams.IndexYears==costYear))),2)+sum(resultsSC.Cap_IBAC(:,find(systemParams.IndexYears==costYear)).*(COM.Capex_IBAT(:,find(systemParams.IndexYears==costYear))),2);
results.IND.CapexNew=sum(resultsSC.Cap_RPVI(:,find(systemParams.IndexYears==costYear)).*(IND.Capex_RPVI(:,find(systemParams.IndexYears==costYear))),2)+sum(resultsSC.Cap_SBAI(:,find(systemParams.IndexYears==costYear)).*(IND.Capex_SBAT(:,find(systemParams.IndexYears==costYear))),2)+sum(resultsSC.Cap_IBAI(:,find(systemParams.IndexYears==costYear)).*(IND.Capex_IBAT(:,find(systemParams.IndexYears==costYear))),2);

results.RES.PV_Opex = sum(resultsSC.Cap_RPVR.*(RES.Opex_fix_RPVR)+(repmat(sum(resultsSC.RPVR_EL)',1,size(systemParams.IndexYears,2)).*genShareRPVR).*RES.Opex_var_RPVR,2);
results.COM.PV_Opex = sum(resultsSC.Cap_RPVC.*(COM.Opex_fix_RPVC)+(repmat(sum(resultsSC.RPVC_EL)',1,size(systemParams.IndexYears,2)).*genShareRPVC).*COM.Opex_var_RPVC,2);
results.IND.PV_Opex = sum(resultsSC.Cap_RPVI.*(IND.Opex_fix_RPVI)+(repmat(sum(resultsSC.RPVI_EL)',1,size(systemParams.IndexYears,2)).*genShareRPVI).*IND.Opex_var_RPVI,2);

results.RES.Bat_Opex = sum(resultsSC.Cap_SBAR.*(RES.Opex_fix_SBAT)+resultsSC.Cap_IBAR.*(RES.Opex_fix_IBAT)+(repmat(sum(resultsSC.SBAR_EL)',1,size(systemParams.IndexYears,2)).*genShareSBAR).*RES.Opex_var_SBAT,2);
results.COM.Bat_Opex = sum(resultsSC.Cap_SBAC.*(COM.Opex_fix_SBAT)+resultsSC.Cap_IBAC.*(COM.Opex_fix_IBAT)+(repmat(sum(resultsSC.SBAC_EL)',1,size(systemParams.IndexYears,2)).*genShareSBAC).*COM.Opex_var_SBAT,2);
results.IND.Bat_Opex = sum(resultsSC.Cap_SBAI.*(IND.Opex_fix_SBAT)+resultsSC.Cap_IBAI.*(IND.Opex_fix_IBAT)+(repmat(sum(resultsSC.SBAI_EL)',1,size(systemParams.IndexYears,2)).*genShareSBAR).*IND.Opex_var_SBAT,2);

results.RES.Opex=results.RES.PV_Opex+results.RES.Bat_Opex;
results.COM.Opex=results.COM.PV_Opex+results.COM.Bat_Opex;
results.IND.Opex=results.IND.PV_Opex+results.IND.Bat_Opex;

results.RES.LCOE = (results.RES.PV_Annual)./sum(resultsSC.RPVR_EL,1)';
results.COM.LCOE = (results.COM.PV_Annual)./sum(resultsSC.RPVC_EL,1)';
results.IND.LCOE = (results.IND.PV_Annual)./sum(resultsSC.RPVI_EL,1)';

primaryCost_SC = results.RES.PV_Annual+results.COM.PV_Annual+results.IND.PV_Annual;
storageCost_SC = results.RES.Bat_Annual+results.COM.Bat_Annual+results.IND.Bat_Annual;

primaryLCOEreal_SC = (primaryCost_SC)./sum((resultsSC.RPVI_EL+resultsSC.RPVC_EL+resultsSC.RPVR_EL),1)';
primaryLCOEreal_SC = primaryLCOEreal_SC.*round((resultsSC.Dem')/10^6)./round((resultsSC.Dem')/10^6);
primaryLCOEreal_SC(isnan(primaryLCOEreal_SC)) =0;
primaryLCOEreal_SC_average = sum(primaryCost_SC)./sum(sum((resultsSC.RPVR_EL+resultsSC.RPVC_EL+resultsSC.RPVI_EL),1)');
primaryLCOEreal_SC_average(isnan(primaryLCOEreal_SC_average)) =0;

primaryLCOE_SC = (primaryLCOEreal_SC.*(prod_electricity_SC-storLosses_SC))./(sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand+resultsSC.EL_THHP+resultsSC.EL_THHR,1)'-el_cons_localpowerFromGrid-el_cons_localheatFromGrid);
primaryLCOE_SC(isnan(primaryLCOE_SC)) =0;
primaryLCOE_SC_average = sum(primaryLCOEreal_SC.*(prod_electricity_SC-storLosses_SC))./sum((sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand+resultsSC.EL_THHP+resultsSC.EL_THHR,1)'-el_cons_localpowerFromGrid-el_cons_localheatFromGrid));
primaryLCOE_SC_average(isnan(primaryLCOE_SC_average)) =0;

LCOS_SC = (storageCost_SC + primaryLCOEreal_SC.*storLosses_SC)./(sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand+resultsSC.EL_THHP+resultsSC.EL_THHR,1)'-el_cons_localpowerFromGrid-el_cons_localheatFromGrid);
LCOS_SC = LCOS_SC.*round((resultsSC.Dem)'/10^6)./round((resultsSC.Dem)'/10^6);
LCOS_SC(isnan(LCOS_SC)) =0;
LCOS_SC_average = sum(storageCost_SC + primaryLCOEreal_SC.*storLosses_SC)./sum((sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand+resultsSC.EL_THHP+resultsSC.EL_THHR,1)'-el_cons_localpowerFromGrid-el_cons_localheatFromGrid));
LCOS_SC_average(isnan(LCOS_SC_average)) =0;


Transfer_SC =  - 0.02*sum(resultsSC.EL_EXCESS,1)'./(sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand+resultsSC.EL_THHP+resultsSC.EL_THHR,1)'-el_cons_localpowerFromGrid-el_cons_localheatFromGrid);
Transfer_SC(isnan(Transfer_SC)) =0;
Transfer_SC_average =  sum(- 0.02*sum(resultsSC.EL_EXCESS,1)')./sum((sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand+resultsSC.EL_THHP+resultsSC.EL_THHR,1)'-el_cons_localpowerFromGrid-el_cons_localheatFromGrid));
Transfer_SC_average(isnan(Transfer_SC_average))=0;

SC_LCOE = ((results.RES.PV_Annual+results.COM.PV_Annual+results.IND.PV_Annual + results.RES.Bat_Annual+results.COM.Bat_Annual+results.IND.Bat_Annual - 0.02*sum(resultsSC.EL_EXCESS,1)').*((results.RES.PV_Annual+results.COM.PV_Annual+results.IND.PV_Annual + results.RES.Bat_Annual+results.COM.Bat_Annual+results.IND.Bat_Annual - 0.02*sum(resultsSC.EL_EXCESS,1)')>0))./(sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand+resultsSC.EL_THHP+resultsSC.EL_THHR,1)'-el_cons_localpowerFromGrid-el_cons_localheatFromGrid);
SC_LCOE = SC_LCOE.*((sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand+resultsSC.EL_THHP+resultsSC.EL_THHR,1)'-el_cons_localpowerFromGrid-el_cons_localheatFromGrid)>10^3);
SC_LCOE(isnan(SC_LCOE)) =0;
SC_LCOE_average = sum(((results.RES.PV_Annual+results.COM.PV_Annual+results.IND.PV_Annual + results.RES.Bat_Annual+results.COM.Bat_Annual+results.IND.Bat_Annual - 0.02*sum(resultsSC.EL_EXCESS,1)').*((results.RES.PV_Annual+results.COM.PV_Annual+results.IND.PV_Annual + results.RES.Bat_Annual+results.COM.Bat_Annual+results.IND.Bat_Annual - 0.02*sum(resultsSC.EL_EXCESS,1)')>0)))./sum((sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand+resultsSC.EL_THHP+resultsSC.EL_THHR,1)'-el_cons_localpowerFromGrid-el_cons_localheatFromGrid));
SC_LCOE_average(isnan(SC_LCOE_average)) =0;
SC_LCOE_average = SC_LCOE_average.*(sum(sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand+resultsSC.EL_THHP+resultsSC.EL_THHR,1)'-el_cons_localpowerFromGrid-el_cons_localheatFromGrid)>10^3);



%% Cost calc

%% Grid
if costYear == 2025
    qq=1;
end

exportsFlows_DC = results.GRID_DC;
exportsFlows_DC(exportsFlows_DC>0) = 0;
importsFlows_DC = results.GRID_DC;
importsFlows_DC(importsFlows_DC<0) = 0;

exportsFlows_AC = results.GRID_AC;
exportsFlows_AC(exportsFlows_AC>0) = 0;
importsFlows_AC = results.GRID_AC;
importsFlows_AC(importsFlows_AC<0) = 0;

for i = 1:length(systemParams.IndexNodes)
    temp = repmat(importsFlows_DC(:,i),1,length(systemParams.IndexNodes)).*exportsFlows_DC./repmat(sum(exportsFlows_DC,2),1,length(systemParams.IndexNodes));
    temp(isnan(temp))=0;
    imports_DC(i,:) = sum(temp,1); %imported to reg i
    temp = repmat(importsFlows_AC(:,i),1,length(systemParams.IndexNodes)).*exportsFlows_AC./repmat(sum(exportsFlows_AC,2),1,length(systemParams.IndexNodes));
    temp(isnan(temp))=0;
    imports_AC(i,:) = sum(temp,1);

    temp = repmat(exportsFlows_DC(:,i),1,length(systemParams.IndexNodes)).*importsFlows_DC./repmat(sum(importsFlows_DC,2),1,length(systemParams.IndexNodes));
    temp(isnan(temp))=0;
    exports_DC(i,:) = sum(temp,1); %exported from reg i
    temp = repmat(exportsFlows_AC(:,i),1,length(systemParams.IndexNodes)).*importsFlows_AC./repmat(sum(importsFlows_AC,2),1,length(systemParams.IndexNodes));
    temp(isnan(temp))=0;
    exports_AC(i,:) = sum(temp,1);



end

try results.OPT_SIZE_TRTL;

    % losses split between exporters
    qq = results.GRID_DC; qq(results.GRID_DC>0) = 0;
    transmLosses_sys_DC = sum(sum(-results.GRID_DC)) * sum(qq,1)'/sum(sum(qq));
    qq = results.GRID_AC; qq(results.GRID_AC>0) = 0;
    transmLosses_sys_AC =sum(sum(-results.GRID_AC)) * sum(qq,1)'/sum(sum(qq));
    % losses split between net exporters

    transmLosses_sys_DC(isnan(transmLosses_sys_DC)) = 0;
    transmLosses_sys_AC(isnan(transmLosses_sys_AC)) = 0;
    transmLosses_sys_DC((transmLosses_sys_DC<0)) = 0;
    transmLosses_sys_AC((transmLosses_sys_AC<0)) = 0;
    transmLosses_sys = transmLosses_sys_DC + transmLosses_sys_AC;

    demand_sys = demand_sys - transmLosses_sys;
    demand_tot = demand_tot - transmLosses_sys;
    demand_origTot = demand_origTot - transmLosses_sys;

    allDem_sys = demand_sys+demand_des+gasProdElCons+heatProdElCons+demand_transp+demand_industry+demand_export+demand_CDR;
    allDem_sys(allDem_sys<=0) = 0;
    allDem_tot = demand_tot+demand_des+gasProdElCons+heatProdElCons+demand_transp+demand_industry+demand_export+demand_CDR;
    allDem_tot(allDem_tot<=0) = 0;

    transmCapex = systemParams.Capex(ismember(systemParams.IndexID,'TRTL'),yearNumb)*sum(results.OPT_SIZE_TRTL'.*systemParams.TLlength) + ...
        systemParams.Capex(ismember(systemParams.IndexID,'THAO'),yearNumb)*sum(results.OPT_SIZE_THAO'.*systemParams.TLlength) + ...
        systemParams.Capex(ismember(systemParams.IndexID,'TRCS'),yearNumb)*sum(max(abs(results.GRID_DC)));
    transmCapexNew = systemParams.Capex(ismember(systemParams.IndexID,'TRTL'),yearNumb)*sum((results.OPT_SIZE_TRTL'-systemParams.TL_DC_LowLimits).*systemParams.TLlength) + ...
        systemParams.Capex(ismember(systemParams.IndexID,'THAO'),yearNumb)*sum((results.OPT_SIZE_THAO'-systemParams.TL_AC_LowLimits).*systemParams.TLlength) + ...
        systemParams.Capex(ismember(systemParams.IndexID,'TRCS'),yearNumb)*sum(max(abs(results.GRID_DC)));
    transmCapexLinesNew = systemParams.Capex(ismember(systemParams.IndexID,'TRTL'),yearNumb)*sum((results.OPT_SIZE_TRTL'-systemParams.TL_DC_LowLimits).*systemParams.TLlength) + ...
        systemParams.Capex(ismember(systemParams.IndexID,'THAO'),yearNumb)*sum((results.OPT_SIZE_THAO'-systemParams.TL_AC_LowLimits).*systemParams.TLlength) + ...
        0;
    transmOpex = systemParams.Opex_fix(ismember(systemParams.IndexID,'TRTL'),yearNumb)*sum(results.OPT_SIZE_TRTL'.*systemParams.TLlength) + ...
        systemParams.Opex_fix(ismember(systemParams.IndexID,'THAO'),yearNumb)*sum(results.OPT_SIZE_THAO'.*systemParams.TLlength) + ...
        systemParams.Opex_fix(ismember(systemParams.IndexID,'TRCS'),yearNumb)*sum(max(abs(results.GRID_DC)));

    capexTL_DC_CRF = systemParams.Capex(ismember(systemParams.IndexID,'TRTL'),yearNumb)  .* CapitalRecoveryFactor(systemParams.WACC,systemParams.Lifetime(ismember(systemParams.IndexID,'TRTL')));
    capexTL_AC_CRF = systemParams.Capex(ismember(systemParams.IndexID,'THAO'),yearNumb)  .* CapitalRecoveryFactor(systemParams.WACC,systemParams.Lifetime(ismember(systemParams.IndexID,'THAO')));
    capexCS_CRF = systemParams.Capex(ismember(systemParams.IndexID,'TRCS'),yearNumb)  .* CapitalRecoveryFactor(systemParams.WACC,systemParams.Lifetime(ismember(systemParams.IndexID,'TRCS')));
    capexTL_DC = systemParams.Capex(ismember(systemParams.IndexID,'TRTL'),yearNumb);
    capexTL_AC = systemParams.Capex(ismember(systemParams.IndexID,'THAO'),yearNumb);
    capexCS = systemParams.Capex(ismember(systemParams.IndexID,'TRCS'),yearNumb);
    opexTL_DC_fix = systemParams.Opex_fix(ismember(systemParams.IndexID,'TRTL'),yearNumb);
    opexTL_AC_fix = systemParams.Opex_fix(ismember(systemParams.IndexID,'THAO'),yearNumb);
    opexCS_fix = systemParams.Opex_fix(ismember(systemParams.IndexID,'TRCS'),yearNumb);
    %opexCS_fix=1.8;
    transmAnCost = (capexTL_AC_CRF+opexTL_AC_fix)*sum(results.OPT_SIZE_THAO'.*systemParams.TLlength) +...
        (capexTL_DC_CRF+opexTL_DC_fix)*sum(results.OPT_SIZE_TRTL'.*systemParams.TLlength) +...
        (capexCS_CRF+opexCS_fix)*sum(results.OPT_SIZE_TRTL);
    linesCapex = (capexTL_AC)*sum(results.OPT_SIZE_THAO'.*systemParams.TLlength)+(capexTL_DC)*sum(results.OPT_SIZE_TRTL'.*systemParams.TLlength);
    convCapex = (capexCS)*sum(results.OPT_SIZE_TRTL);
    opexCS_fix=1.8;
    linesOpex = (opexTL_AC_fix)*sum(results.OPT_SIZE_THAO'.*systemParams.TLlength)+(opexTL_DC_fix)*sum(results.OPT_SIZE_TRTL'.*systemParams.TLlength);
    convOpex = (opexCS_fix)*sum(results.OPT_SIZE_TRTL);
    transmOpex = linesOpex + convOpex;
    % exported/imported energy: first index (rows) represent sources, second index (cols) represent sinks
    transMatrix = LCOEremoteGeneration(results);

    % adding self supply to transMatrix
    transMatrix = transMatrix +  diag((prod_electricity_sys - sum(transMatrix,2)));

    for k=1:size(transMatrix,1)
        supplyRatio(:,k) = transMatrix(:,k) / sum(transMatrix(:,k),1);
    end
    supplyRatio(isnan(supplyRatio)) = 1; % set to 1 for single region systems
    maxTransmPerRegion = max(abs(results.GRID_DC+results.GRID_AC));
    % Call of funnction LCOT: Returns LCOT per region and for the whole
    % system; 1- split by exports, 0 - by imports, 0.5 - by half
    [lcot_sys,lcot_sys_average] = LCOT(results,systemParams,activeElements,1,prod_electricity_sys,maxTransmPerRegion,allDem_sys,costYear);
    lcot_sys(isnan(lcot_sys)) = 0; % set to zero for single region systems
    lcot_sys(isinf(lcot_sys)) = 0; % set to zero for zero demand regions - to be FIXED!

    lcot_tot = lcot_sys.* (allDem_sys)'./ (allDem_tot)';
    lcot_tot_average = lcot_sys_average.*(sum(sum(systemParams.ValueLoad(1,:,2,:))))./(sum(sum(systemParams.ValueLoad(1,:,2,:)))+sum(sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand,1)'));



catch
    capexTL_DC_CRF = 0;
    capexTL_AC_CRF = 0;
    capexCS_CRF = 0;
    capexTL_DC = 0;
    capexTL_AC = 0;
    capexCS = 0;
    opexTL_DC_fix = 0;
    opexTL_AC_fix = 0;
    opexCS_fix = 0;

    transmAnCost=zeros(regNumb,1);
    transmCapex=zeros(1,1);
    transmCapexNew=zeros(1,1);
    transmCapexLinesNew=zeros(1,1);
    transmOpex=0;
    linesCapex = 0;
    convCapex = 0;
    linesOpex = 0;
    convOpex = 0;
    results.OPT_SIZE_TRTL=0;
    results.OPT_SIZE_THAO=0;
    lcot_sys=zeros(regNumb,1)';
    lcot_sys_average = 0;

    lcot_tot=zeros(size(regNumb,1))';
    lcot_tot_average = 0;

    allDem_sys = demand_sys+demand_des+gasProdElCons+heatProdElCons+demand_transp+demand_industry+demand_export;
    allDem_sys(allDem_sys<=0) = 0;
    allDem_tot = demand_tot+demand_des+gasProdElCons+heatProdElCons+demand_transp+demand_industry+demand_export;
    allDem_tot(allDem_tot<=0) = 0;

    results.LINEpos_AC=0;
    results.LINEneg_AC=0;
    results.LINEpos_DC=0;
    results.LINEneg_DC=0;
    Share=zeros(regNumb,1);
    transmLosses_sys_DC = zeros(regNumb,1);
    transmLosses_sys_AC = zeros(regNumb,1);

    transmLosses_sys = zeros(regNumb,1);
end

%% Electricity

%%
prod_electricity_sys = sum(results.RPVO_EL,1)'+sum(results.RPVA_EL,1)'+...
    sum(results.RPBO_EL,1)'+sum(results.RPBA_EL,1)'+sum(results.RPBV_EL,1)'+sum(results.RPVF_EL,1)'+...
    sum(results.RWIN_EL,1)'+sum(results.RWIO_EL,1)'+sum(results.ROWI_EL,1)'+...
    sum(results.RWAV_EL,1)'+...
    sum(results.RRRI_EL,1)'+sum(results.HDAM_EL,1)'+...
    sum(results.TSTU_EL,1)'+sum(results.TGEO_EL,1)'+...
    sum(results.TBPP_EL,1)'+sum(results.TCHP_EL,1)'+sum(results.TMSW_EL,1)'+...
    sum(results.TBCS_EL,1)'+sum(results.TCBC_EL,1)'+sum(results.TCWC_EL,1)'+...
    sum(results.THPP_EL,1)'+sum(results.THCS_EL,1)'+sum(results.TICG_EL,1)'+sum(results.TICM_EL,1)'.*(shareOfICM_NGA+shareOfICM_BME+shareOfICM_SNGimp+shareOfICM_HYblack+shareOfICM_HYimpor+shareOfICM_Oil)+...
    sum(results.TCCG_EL+results.TOCG_EL+results.TGCS_EL+results.TCNG_EL,1)'.*(shareOfGT_NGA+shareOfGT_BME + shareOfGT_SNGimp+shareOfGT_HYimpor+shareOfGT_HYblack+shareOfGT_Oil)+...
    sum(results.TNUC_EL,1)'+...
    sum(results.TCCO_EL,1)'+sum(results.TCOI_EL,1)'+sum(results.TCBP_EL,1)';


prod_electricity_tot_RE = sum(resultsSC.RPVR_EL,1)'+ sum(resultsSC.RPVC_EL,1)'+ sum(resultsSC.RPVI_EL,1)'+...
    sum(results.RPVO_EL,1)'+sum(results.RPVA_EL,1)'+...
    sum(results.RPBO_EL,1)'+sum(results.RPBA_EL,1)'+sum(results.RPBV_EL,1)'+sum(results.RPVF_EL,1)'+...
    sum(results.RWIN_EL,1)'+sum(results.RWIO_EL,1)'+sum(results.ROWI_EL,1)'+...
    sum(results.RWAV_EL,1)'+...
    sum(results.RRRI_EL,1)'+sum(results.HDAM_EL,1)'+...
    sum(results.TSTU_EL,1)'+sum(results.TGEO_EL,1)'+...
    sum(results.TBPP_EL,1)'+sum(results.TCBP_EL,1)'+sum(results.TCHP_EL,1)'+sum(results.TMSW_EL,1)'+...
    sum(results.TBCS_EL,1)'+sum(results.TCBC_EL,1)'+sum(results.TCWC_EL,1)'+...
    sum(results.TICM_EL,1)'.*(shareOfICM_BME+ shareOfICM_SNGimp+shareOfICM_HYimpor)+...
    sum(results.TCCG_EL+results.TOCG_EL+results.TGCS_EL+results.TCNG_EL,1)'.*(shareOfGT_BME+ shareOfGT_SNGimp+shareOfGT_HYimpor);

prod_electricity_tot_Fos = sum(results.THPP_EL,1)'+sum(results.THCS_EL,1)'+sum(results.TICG_EL,1)'+sum(results.TICM_EL,1)'.*(shareOfICM_NGA+shareOfICM_Oil+shareOfICM_HYblack)+...
    sum(results.TCCG_EL+results.TOCG_EL+results.TGCS_EL+results.TCNG_EL,1)'.*(shareOfGT_NGA+shareOfGT_Oil+shareOfGT_HYblack)+...
    sum(results.TCCO_EL,1)'+sum(results.TCOI_EL,1)';

prod_electricity_tot_ocean = sum(results.RWAV_EL,1)';
prod_electricity_tot_hydro = sum(results.RRRI_EL,1)'+sum(results.HDAM_EL,1)';
prod_electricity_tot_wind = sum(results.RWIN_EL,1)'+sum(results.RWIO_EL,1)'+sum(results.ROWI_EL,1)';
prod_electricity_tot_PV =  sum(results.RPVO_EL,1)'+sum(results.RPVA_EL,1)'+sum(resultsSC.RPVR_EL,1)'+ sum(resultsSC.RPVC_EL,1)'+ sum(resultsSC.RPVI_EL,1)'+...
    sum(results.RPBO_EL,1)'+sum(results.RPBA_EL,1)'+sum(results.RPBV_EL,1)'+sum(results.RPVF_EL,1)';
prod_electricity_totcoll_biomass = sum(results.TBPP_EL,1)'+sum(results.TCBP_EL,1)'+...
    sum(results.TMSW_EL,1)'+sum(results.TCHP_EL,1)'+...
    sum(results.TBCS_EL,1)'+sum(results.TCBC_EL,1)'+sum(results.TCWC_EL,1)'+...
    sum(results.TICM_EL,1)'.*(shareOfICM_BME)+...
    sum(results.TCCG_EL+results.TOCG_EL+results.TCNG_EL,1)'.*shareOfGT_BME;
prod_electricity_tot_other = sum(results.TSTU_EL,1)'+sum(results.TGEO_EL,1)'+...
    sum(results.TICM_EL,1)'.*(shareOfICM_SNGimp+shareOfICM_HYimpor)+...
    sum(results.TCCG_EL+results.TOCG_EL+results.TCNG_EL,1)'.*(shareOfGT_SNGimp+shareOfGT_HYimpor);
prod_electricity_tot_Fos_gas = sum(results.TCCG_EL+results.TOCG_EL+results.TCNG_EL+results.TGCS_EL,1)'.*shareOfGT_NGA+sum(results.TICM_EL,1)'.*(shareOfICM_NGA+shareOfICM_HYblack);
prod_electricity_tot_Fos_oil = sum(results.TICG_EL,1)'+sum(results.TICM_EL,1)'.*shareOfICM_Oil+sum(results.TCOI_EL,1)'+sum(results.TCCG_EL+results.TOCG_EL+results.TGCS_EL+results.TCNG_EL,1)'.*(shareOfGT_Oil);
prod_electricity_tot_Fos_coal = sum(results.THPP_EL,1)'+sum(results.THCS_EL,1)'+sum(results.TCCO_EL,1)';
prod_electricity_tot_nuclear = sum(results.TNUC_EL,1)';

prod_electricity_tot = prod_electricity_sys+ sum(resultsSC.RPVR_EL,1)'+ sum(resultsSC.RPVC_EL,1)'+ sum(resultsSC.RPVI_EL,1)';
%%


primaryCost_EL = ...
    sum(opex_capex(:,:,ismember(allActiveExLoadLabels,'RPVO')).*capacity(:,:,ismember(allActiveExLoadLabels,'RPVO'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RPVO')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RPVO')).*permute(repmat(sum(results.RPVO_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RPVA')).*capacity(:,:,ismember(allActiveExLoadLabels,'RPVA'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RPVA')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RPVA')).*permute(repmat(sum(results.RPVA_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RPBO')).*capacity(:,:,ismember(allActiveExLoadLabels,'RPBO'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RPBO')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RPBO')).*permute(repmat(sum(results.RPBO_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RPBA')).*capacity(:,:,ismember(allActiveExLoadLabels,'RPBA'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RPBA')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RPBA')).*permute(repmat(sum(results.RPBA_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RPBV')).*capacity(:,:,ismember(allActiveExLoadLabels,'RPBV'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RPBV')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RPBV')).*permute(repmat(sum(results.RPBV_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RPVF')).*capacity(:,:,ismember(allActiveExLoadLabels,'RPVF'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RPVF')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RPVF')).*permute(repmat(sum(results.RPVF_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RWIN')).*capacity(:,:,ismember(allActiveExLoadLabels,'RWIN'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RWIN')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RWIN')).*permute(repmat(sum(results.RWIN_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.RWIN_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'ROWI')).*capacity(:,:,ismember(allActiveExLoadLabels,'ROWI'))+opex_var(:,:,ismember(allActiveExLoadLabels,'ROWI')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'ROWI')).*permute(repmat(sum(results.ROWI_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.ROWI_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RWIO')).*capacity(:,:,ismember(allActiveExLoadLabels,'RWIO'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RWIO')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RWIO')).*permute(repmat(sum(results.RWIO_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.RWIN_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RWAV')).*capacity(:,:,ismember(allActiveExLoadLabels,'RWAV'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RWAV')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RWAV')).*permute(repmat(sum(results.RWAV_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.RWAV_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RRRI')).*capacity(:,:,ismember(allActiveExLoadLabels,'RRRI'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RRRI')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RRRI')).*permute(repmat(sum(results.RRRI_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.RRRI_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'HDAM')).*capacity(:,:,ismember(allActiveExLoadLabels,'HDAM'))+opex_var(:,:,ismember(allActiveExLoadLabels,'HDAM')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'HDAM')).*permute(repmat(sum(results.HDAM_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.HDAM_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TSTU')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSTU'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSTU')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSTU')).*permute(repmat(sum(results.TSTU_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'RCSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'RCSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RCSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RCSP')).*permute(repmat(sum(results.RCSP_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CSPtoEl),1,length(systemParams.IndexYears)) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TGEO')).*capacity(:,:,ismember(allActiveExLoadLabels,'TGEO'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TGEO')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TGEO')).*permute(repmat(sum(results.TGEO_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TGEO_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TBPP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBPP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBPP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBPP')).*permute(repmat(sum(results.TBPP_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TBPP_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TBCS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBCS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBCS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBCS')).*permute(repmat(sum(results.TBCS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TBPP_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THPP')).*capacity(:,:,ismember(allActiveExLoadLabels,'THPP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'THPP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'THPP')).*permute(repmat(sum(results.THPP_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.THPP_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THCS')).*capacity(:,:,ismember(allActiveExLoadLabels,'THCS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'THCS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'THCS')).*permute(repmat(sum(results.THCS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.THPP_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TICG')).*capacity(:,:,ismember(allActiveExLoadLabels,'TICG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TICG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TICG')).*permute(repmat(sum(results.TICG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TICG_EL,1)' + ...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TICM')).*capacity(:,:,ismember(allActiveExLoadLabels,'TICM'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TICM')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TICM')).*permute(repmat(sum(results.TICM_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(shareOfICM_NGA+shareOfICM_BME+shareOfICM_SNGimp+shareOfICM_HYblack+shareOfICM_HYimpor+shareOfICM_Oil,1,length(systemParams.IndexYears)) +....*sum(results.TICG_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TNUC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TNUC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TNUC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TNUC')).*permute(repmat(sum(results.TNUC_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TNUC_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TMSW')).*capacity(:,:,ismember(allActiveExLoadLabels,'TMSW'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TMSW')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TMSW')).*permute(repmat(sum(results.TMSW_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TMSW_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCWC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCWC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCWC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCWC')).*permute(repmat(sum(results.TCWC_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TMSW_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCHP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCHP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCHP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCHP')).*permute(repmat(sum(results.TCHP_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TCHP_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCOI')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOI'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCOI')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCOI')).*permute(repmat(sum(results.TCOI_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCCO')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCCO'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCCO')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCCO')).*permute(repmat(sum(results.TCCO_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCBP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCBP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCBP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCBP')).*permute(repmat(sum(results.TCBP_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCBC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCBC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCBC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCBC')).*permute(repmat(sum(results.TCBC_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGD'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGD')).*permute(repmat(sum(results.RBGA_FU,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((shareOfBGA_El),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGU'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGU')).*permute(repmat(sum(results.TBGU_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((shareOfGasEl),1,length(systemParams.IndexYears)) +....*sum(results.RBGA_TBGD,1)').*(1-shareOfgasBGA) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TCNG')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCNG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCNG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCNG')).*permute(repmat(sum(results.TCNG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TGCS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TGCS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TGCS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TGCS')).*permute(repmat(sum(results.TGCS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TCCG_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCCG')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCCG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCCG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCCG')).*permute(repmat(sum(results.TCCG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TCCG_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TOCG')).*capacity(:,:,ismember(allActiveExLoadLabels,'TOCG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TOCG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TOCG')).*permute(repmat(sum(results.TOCG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(shareOfGT_NGA+shareOfGT_BME+shareOfGT_SNGimp+shareOfGT_HYimpor+shareOfGT_HYblack+shareOfGT_Oil,1,length(systemParams.IndexYears)),2)+....*sum(results.TOCG_EL,1)').*shareOfGT_NGA);
    opex_capexRehab_RoR(:,yearNumb).*(RoR_LL)'+...
    opex_capexRehab_Dam(:,yearNumb).*(Dam_LL)'+...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWOO')).*sum(results.RWOO_TBPP+results.RWOO_TCBP+results.RWOO_TBCS+results.RWOO_TCBC,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWWO')).*sum(results.RWWO_TBPP+results.RWWO_TCBP+results.RWWO_TBCS+results.RWWO_TCBC,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RHAR')).*sum(results.RHAR_THPP+results.RHAR_THCS+results.RHAR_TCCO,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RPET')).*sum(results.RPET_TICG+results.RPET_TICM+results.RPET_TCOI+results.RPET_TCCG+results.RPET_TOCG+results.RPET_TCNG+results.RPET_TGCS,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RURA')).*sum(results.RURA_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBMW')).*sum(results.RBMW_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.TCCG_GAS_FOS+results.TOCG_GAS_FOS+results.TCNG_GAS_FOS+results.TGCS_GAS_FOS+results.TICM_GAS_FOS,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBGA')).*sum(results.RBGA_FU,1)'.* shareOfBGA_El+ ...
    setup.importLNGCost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importLNG_regas,1)'.* shareOfGasEl+...
    setup.importH2Cost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importH2_regas,1)'.* share.H2toEL+...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RHAR_CO,1)'/1000.*shareOfElCoalGHG + ...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RPET_CO,1)'/1000.*shareOfElOilGHG + ...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RNGA_CO,1)'/1000.*shareOfElGasGHG+...
    setup.importELCost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importEL,1)';

primaryCost_EL_fuel = (opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWOO')).*sum(results.RWOO_TBPP+results.RWOO_TCBP+results.RWOO_TBCS+results.RWOO_TCBC,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWWO')).*sum(results.RWWO_TBPP+results.RWWO_TCBP+results.RWWO_TBCS+results.RWWO_TCBC,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RHAR')).*sum(results.RHAR_THPP+results.RHAR_THCS+results.RHAR_TCCO,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RPET')).*sum(results.RPET_TICG+results.RPET_TICM+results.RPET_TCOI+results.RPET_TCCG+results.RPET_TOCG+results.RPET_TCNG+results.RPET_TGCS,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RURA')).*sum(results.RURA_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBMW')).*sum(results.RBMW_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.TCCG_GAS_FOS+results.TOCG_GAS_FOS+results.TCNG_GAS_FOS+results.TGCS_GAS_FOS+results.TICM_GAS_FOS,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBGA')).*sum(results.RBGA_FU,1)'.* shareOfBGA_El+...
    setup.importLNGCost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importLNG_regas,1)'.* shareOfGasEl+...
    setup.importH2Cost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importH2_regas,1)'.* share.H2toEL+...
    0)./primaryCost_EL;
primaryCost_EL_fuelav = sum(opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWOO')).*sum(results.RWOO_TBPP+results.RWOO_TCBP+results.RWOO_TBCS+results.RWOO_TCBC,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWWO')).*sum(results.RWWO_TBPP+results.RWWO_TCBP+results.RWWO_TBCS+results.RWWO_TCBC,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RHAR')).*sum(results.RHAR_THPP+results.RHAR_THCS+results.RHAR_TCCO,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RPET')).*sum(results.RPET_TICG+results.RPET_TICM+results.RPET_TCOI+results.RPET_TCCG+results.RPET_TOCG+results.RPET_TCNG+results.RPET_TGCS,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RURA')).*sum(results.RURA_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBMW')).*sum(results.RBMW_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.TCCG_GAS_FOS+results.TOCG_GAS_FOS+results.TCNG_GAS_FOS+results.TGCS_GAS_FOS+results.TICM_GAS_FOS,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBGA')).*sum(results.RBGA_FU,1)'.* shareOfBGA_El+...
    setup.importLNGCost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importLNG_regas,1)'.* shareOfGasEl+...
    setup.importH2Cost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importH2_regas,1)'.* share.H2toEL+...
    0)./sum(primaryCost_EL);

primaryCost_EL_CO2 = (systemParams.fossilCO2Cost(yearNumb).*sum(results.RHAR_CO/10^3,1)'.*shareOfElCoalGHG + ...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RPET_CO/10^3,1)'.*shareOfElOilGHG + ...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RNGA_CO/10^3,1)'.*shareOfElGasGHG)./primaryCost_EL;
primaryCost_EL_CO2av = sum(systemParams.fossilCO2Cost(yearNumb).*sum(results.RHAR_CO/10^3,1)'.*shareOfElCoalGHG + ...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RPET_CO/10^3,1)'.*shareOfElOilGHG + ...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RNGA_CO/10^3,1)'.*shareOfElGasGHG)./sum(primaryCost_EL);

primaryCost_EL_main = 1-primaryCost_EL_fuel-primaryCost_EL_CO2;
primaryCost_EL_mainav = 1-primaryCost_EL_fuelav-primaryCost_EL_CO2av;



RE_capacity_tot = (results.OPT_SIZE_RPVO+results.OPT_SIZE_RPVA +...
    results.OPT_SIZE_RPBO+results.OPT_SIZE_RPBA + results.OPT_SIZE_RPBV+results.OPT_SIZE_RPVF +...
    resultsSC.OPT_SIZE_RPVR + resultsSC.OPT_SIZE_RPVC + resultsSC.OPT_SIZE_RPVI +...
    results.OPT_SIZE_RWIN + results.OPT_SIZE_RWIO + results.OPT_SIZE_ROWI +...
    results.OPT_SIZE_RRRI + results.OPT_SIZE_HDAM +...
    results.OPT_SIZE_TBPP + results.OPT_SIZE_TCBP + results.OPT_SIZE_TCHP + results.OPT_SIZE_TMSW + ...
    results.OPT_SIZE_TBCS + results.OPT_SIZE_TCBC + results.OPT_SIZE_TCWC + ...
    results.OPT_SIZE_TSTU + results.OPT_SIZE_TGEO)' +...
    (results.OPT_SIZE_TCCG+results.OPT_SIZE_TOCG+results.OPT_SIZE_TCNG)'.*(shareOfGT_BME+shareOfGT_SNGimp+shareOfGT_HYimpor) ;

primaryLCOEreal = primaryCost_EL./prod_electricity_sys;
primaryLCOEreal(isinf(primaryLCOEreal))=0;
primaryLCOEreal_average = sum(primaryCost_EL)./sum(prod_electricity_sys);

primaryLCOEsys = (primaryLCOEreal.*(prod_electricity_sys-sum(results.EL_EXCESS,1)'-storLosses_sys-transmLosses_sys))./((allDem_sys));
if sum(isinf(primaryLCOEsys))
    primaryLCOEsys(isinf(primaryLCOEsys))=0;
    warning('allDem_sys is negative')
end

primaryLCOEtot = (primaryLCOEsys.*(allDem_sys) + primaryLCOE_SC.*(sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand+resultsSC.EL_THHP+resultsSC.EL_THHR,1)'-el_cons_localpowerFromGrid-el_cons_localheatFromGrid))./(allDem_tot);
if sum(isinf(primaryLCOEtot))
    primaryLCOEtot(isinf(primaryLCOEtot))=0;
    warning('allDem_sys is negative')
end

Transfer_sys = - 0.02*sum(resultsSC.EL_EXCESS,1)'./((allDem_sys));

primaryLCOEsys_average = sum(primaryLCOEsys.*(allDem_sys))./sum((allDem_sys));
primaryLCOEtot_average = sum(primaryLCOEsys.*(allDem_sys) + primaryLCOE_SC.*(sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand+resultsSC.EL_THHP+resultsSC.EL_THHR,1)'-el_cons_localpowerFromGrid-el_cons_localheatFromGrid))./sum((allDem_tot));



LCOCsys = (primaryLCOEreal.*(sum(results.EL_EXCESS,1)'))./((allDem_sys));
if sum(isinf(LCOCsys))
    LCOCsys(isinf(LCOCsys))=0;
    warning('allDem_sys is negative')
end
LCOCtot = (LCOCsys.*(allDem_sys))./((allDem_tot));
LCOCtot(isnan(LCOCtot))=0;

LCOCsys_average = sum(LCOCsys.*(allDem_sys))./sum((allDem_sys));
LCOCtot_average = sum(LCOCtot.*(allDem_tot))./sum((allDem_tot));


allDem_sys = demand_sys+demand_des+gasProdElCons+heatProdElCons+demand_transp+demand_industry+demand_export+demand_CDR;
if sum(allDem_sys<0)
    warning('allDem_sys is negative')
end
allDem_sys(allDem_sys<=0) = 0;
allDem_tot = demand_tot+demand_des+gasProdElCons+heatProdElCons+demand_transp+demand_industry+demand_export+demand_CDR;
if sum(allDem_tot<0)
    warning('allDem_tot is negative')
end
allDem_tot(allDem_tot<=0) = 0;

allDemCons_sys = demand_El_Power_sys+demand_des+gasProdElCons+heatProdElCons+demand_transp+demand_industry+demand_export+demand_CDR;
allDemCons_tot = demand_El_Power_tot+demand_des+gasProdElCons+heatProdElCons+demand_transp+demand_industry+demand_export+demand_CDR;

q1 = (prod_electricity_tot);
q2 = (allDem_tot+sum(imports_AC')'+sum(imports_DC')'+sum(exports_AC')'+sum(exports_DC')')+2*(transmLosses_sys_AC)+2*(transmLosses_sys_DC)+(storLosses_sys)+(sum(results.EL_EXCESS)');
q3 = (allDemCons_tot+transmLosses_sys_AC+transmLosses_sys_DC+sum(resultsSC.EL_THHP+resultsSC.EL_THHR)') + (sum(resultsSC.EL_EXCESS,1)')+(storLosses_sys)+(sum(results.EL_EXCESS)');

lcot_sys = lcot_sys' + (primaryLCOEreal.*(transmLosses_sys))./((allDem_sys));

lcot_tot = lcot_sys.* (allDem_sys)./ (allDem_tot);

lcot_sys_average = sum(lcot_sys.* (allDem_sys))./sum((allDem_sys));
lcot_tot_average = sum(lcot_sys.* (allDem_sys))./ sum((allDem_tot));
lcot_sys(isnan(lcot_sys)) = 0;
lcot_tot(isnan(lcot_tot)) = 0;
lcot_sys_average(isnan(lcot_sys_average)) = 0;
lcot_tot_average(isnan(lcot_tot_average)) = 0;
lcot_sys(isinf(lcot_sys)) = 0;
lcot_tot(isinf(lcot_tot)) = 0;
lcot_sys_average(isinf(lcot_sys_average)) = 0;
lcot_tot_average(isinf(lcot_tot_average)) = 0;

try results.TCOS_GA;
catch
    results.TCOS_GA = results.TCOS_CO;
end

if setup.Heat.Flag
    shareOfHeatInElStorage = 0.8*(shareOfGasDH.*(sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'+sum(results.EL_TWEL,1)')+sum(results.EL_TDHR,1)'+sum(results.EL_TDHP,1)'+sum(results.EL_TDGE,1)'+sum(resultsSC.EL_THHR,1)'+sum(resultsSC.EL_THHP,1)')./prod_electricity_tot;

else
    shareOfHeatInElStorage = zeros(size(prod_electricity_tot));

end
storageCost = sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'SBAT')).*capacity(:,:,ismember(allActiveExLoadLabels,'SBAT'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SBAT')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SBAT')).*permute(repmat(sum(results.SBAT_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SBAT_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IBAT')).*capacity(:,:,ismember(allActiveExLoadLabels,'IBAT'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SACA')).*capacity(:,:,ismember(allActiveExLoadLabels,'SACA'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SACA')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SACA')).*permute(repmat(sum(results.SACA_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SACA_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IACA')).*capacity(:,:,ismember(allActiveExLoadLabels,'IACA'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SPHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SPHS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SPHS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SPHS')).*permute(repmat(sum(results.SPHS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SPHS_EL,1)' ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IPHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'IPHS'))+...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'SCO2'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SCO2')).*permute(repmat(sum(results.SCO2_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toCH4.*shareOfGasEl),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'ICO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'ICO2'))).*repmat((share.CO2toCH4.*shareOfGasEl),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHYG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHYG')).*permute(repmat(sum(results.SHYG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toEL),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'IHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHYG'))).*repmat((share.H2toEL),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SGAS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SGAS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SGAS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SGAS')).*permute(repmat(sum(results.SGAS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SGAS_TCCG,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IGAS')).*capacity(:,:,ismember(allActiveExLoadLabels,'IGAS'))).*repmat((1-shareOfHeatInElStorage),1,length(systemParams.IndexYears))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SBGA')).*capacity(:,:,ismember(allActiveExLoadLabels,'SBGA'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SBGA')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SBGA')).*permute(repmat(sum(results.SBGA_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SBGA_TCHP,1)' +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TDHR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDHR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDHR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDHR')).*permute(repmat(sum(results.TDHR_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDHP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDHP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDHP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDHP')).*permute(repmat(sum(results.TDHP_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IHOT')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHOT'))+...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SHOT')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHOT'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHOT')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHOT')).*permute(repmat(sum(results.SHOT_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SHOT_HE,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SDHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SDHS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SDHS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SDHS')).*permute(repmat(sum(results.SDHS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...sum(results.SHOT_HE,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IDHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'IDHS'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THHR')).*(resultsSC.Cap_THHR)+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THHP')).*(resultsSC.Cap_THHP))*(~setup.Heat.Flag)+...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TICM')).*capacity(:,:,ismember(allActiveExLoadLabels,'TICM'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TICM')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TICM')).*permute(repmat(sum(results.TICM_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(shareOfICM_SNG+shareOfICM_HY,1,length(systemParams.IndexYears)) +....*sum(results.TICG_EL,1)' + ...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TCNG')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCNG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCNG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCNG')).*permute(repmat(sum(results.TCNG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TGCS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TGCS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TGCS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TGCS')).*permute(repmat(sum(results.TGCS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TCCG_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCCG')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCCG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCCG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCCG')).*permute(repmat(sum(results.TCCG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.TCCG_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TOCG')).*capacity(:,:,ismember(allActiveExLoadLabels,'TOCG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TOCG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TOCG')).*permute(repmat(sum(results.TOCG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((shareOfGT_SNG+shareOfGT_HYgreen),1,length(systemParams.IndexYears)) +...sum(results.TOCG_EL,1)').*(1-shareOfGT_NGA) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.H2toEL,1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*shareOfGasEl +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMR')).*permute(repmat(sum(results.TSMR_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.H2toEL,1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*shareOfGasEl +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCOS')).*permute(repmat(sum(results.TCOS_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*shareOfGasEl,1,length(systemParams.IndexYears)) +...sum(results.TCOS_GA,1)').*shareOfGasEl +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSC')).*permute(repmat(sum(results.TPSC_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*shareOfGasEl,1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSP')).*permute(repmat(sum(results.TPSP_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*shareOfGasEl,1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TMET')).*capacity(:,:,ismember(allActiveExLoadLabels,'TMET'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TMET')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TMET')).*permute(repmat(sum(results.TMET_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(shareOfGasEl,1,length(systemParams.IndexYears)),2) +...sum(results.TMET_GA,1)').*shareOfGasEl);
    opex_capexRehab_PHS(:,yearNumb).*(PHS_LL)';



gasCost = sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGU'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGU')).*permute(repmat(sum(results.TBGU_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((shareOfGasInd),1,length(systemParams.IndexYears)) +...sum(results.RBME,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGD'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGD')).*permute(repmat(sum(results.RBGA_FU,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(shareOfBGA_Ind,1,length(systemParams.IndexYears)) +...sum(results.RBGA_TBGD,1)').*(shareOfgasBGA).*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toCH4.*shareOfGasInd),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMR')).*permute(repmat(sum(results.TSMR_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toCH4.*shareOfGasInd),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*shareOfGasEl +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCOS')).*permute(repmat(sum(results.TCOS_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toCH4.*shareOfGasInd),1,length(systemParams.IndexYears)) +....*sum(results.TCOS_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSC')).*permute(repmat(sum(results.TPSC_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toCH4.*shareOfGasInd),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSP')).*permute(repmat(sum(results.TPSP_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toCH4.*shareOfGasInd),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TMET')).*capacity(:,:,ismember(allActiveExLoadLabels,'TMET'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TMET')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TMET')).*permute(repmat(sum(results.TMET_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((shareOfGasInd),1,length(systemParams.IndexYears)),2) +...sum(results.TMET_GA,1)').*(1-shareOfGasEl))+...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBGA')).*sum(results.RBGA_FU,1)'.* shareOfBGA_Ind+ ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.LIGA_GAS_FOS,1)' +...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RNGA_CO,1)'/1000.*shareOfNGAInd;

hyCost = sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2]))+....*repmat((share.H2toM),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMR')).*permute(repmat(sum(results.TSMR_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2]))+....*repmat((share.H2toM),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*shareOfGasEl +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHYG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHYG')).*permute(repmat(sum(results.SHYG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]))+....*repmat((share.H2toM),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'IHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHYG'))),2)+....*repmat((share.H2toM),1,length(systemParams.IndexYears)),2) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMC')).*permute(repmat(sum(results.TSMC_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])),2)+....*repmat((share.H2toM),1,length(systemParams.IndexYears)),2) +...sum(results.TMET_GA,1)').*(1-shareOfGasEl))+...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.TSMR_GAS_FOS,1)'+....*share.H2toM +...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.TSMC_GAS_FOS,1)'+....*share.H2toM +...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RNGA_CO,1)'/1000.*shareOfNGAHy+...%CO2 flows in kg
    setup.importH2Cost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importH2_regas,1)';

hyCost_ren = sum(((opex_capex(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2]))+....*repmat((share.H2toM),1,length(systemParams.IndexYears))+...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHYG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHYG')).*permute(repmat(sum(results.SHYG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]))+....*repmat((share.H2toM),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'IHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHYG')))),2)+....*repmat((share.H2toM),1,length(systemParams.IndexYears))),2)+...
    setup.importH2Cost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importH2_regas,1)';%.* share.H2toM;

hyCost_fos = sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMR')).*permute(repmat(sum(results.TSMR_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])),2)+....*repmat((share.H2toM),1,length(systemParams.IndexYears)),2) +...sum(results.TMET_GA,1)').*(1-shareOfGasEl))+...
    sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMC')).*permute(repmat(sum(results.TSMC_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])),2)+....*repmat((share.H2toM),1,length(systemParams.IndexYears)),2) +...sum(results.TMET_GA,1)').*(1-shareOfGasEl))+...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.TSMR_GAS_FOS,1)'+....*share.H2toM +...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.TSMC_GAS_FOS,1)'+....*share.H2toM +...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RNGA_CO,1)'/1000.*shareOfNGAHy;%.*share.H2toM;

FTCost = sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'TFTU')).*capacity(:,:,ismember(allActiveExLoadLabels,'TFTU'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TFTU')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TFTU')).*permute(repmat(sum(results.TFTU_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCOS')).*permute(repmat(sum(results.TCOS_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toFTU),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSC')).*permute(repmat(sum(results.TPSC_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toFTU),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSP')).*permute(repmat(sum(results.TPSP_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toFTU),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toFTU),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMR')).*permute(repmat(sum(results.TSMR_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toFTU),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*shareOfGasEl +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'SCO2'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SCO2')).*permute(repmat(sum(results.SCO2_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toFTU),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'ICO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'ICO2'))).*repmat((share.CO2toFTU),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHYG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHYG')).*permute(repmat(sum(results.SHYG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toFTU),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'IHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHYG'))).*repmat((share.H2toFTU),1,length(systemParams.IndexYears)),2) +...sum(results.TMET_GA,1)').*(1-shareOfGasEl))+...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.TSMR_GAS_FOS,1)'.*share.H2toFTU +...
    setup.importH2Cost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importH2_regas,1)'.* share.H2toFTU+...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RNGA_CO,1)'/1000.*shareOfNGAHy.*share.H2toFTU;

FTCost_running = sum((opex_fix(:,:,ismember(allActiveExLoadLabels,'TFTU')).*capacity(:,:,ismember(allActiveExLoadLabels,'TFTU'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TFTU')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TFTU')).*permute(repmat(sum(results.TFTU_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_fix(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCOS')).*permute(repmat(sum(results.TCOS_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toFTU),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSC')).*permute(repmat(sum(results.TPSC_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toFTU),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSP')).*permute(repmat(sum(results.TPSP_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toFTU),1,length(systemParams.IndexYears)) +...
    (opex_fix(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toFTU),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_fix(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMR')).*permute(repmat(sum(results.TSMR_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toFTU),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*shareOfGasEl +...
    (opex_fix(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'SCO2'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SCO2')).*permute(repmat(sum(results.SCO2_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toFTU),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_fix(:,:,ismember(allActiveExLoadLabels,'ICO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'ICO2'))).*repmat((share.CO2toFTU),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_fix(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHYG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHYG')).*permute(repmat(sum(results.SHYG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toFTU),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_fix(:,:,ismember(allActiveExLoadLabels,'IHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHYG'))).*repmat((share.H2toFTU),1,length(systemParams.IndexYears)),2) +...sum(results.TMET_GA,1)').*(1-shareOfGasEl))+...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.TSMR_GAS_FOS,1)'.*share.H2toFTU +...
    setup.importH2Cost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importH2_regas,1)'.* share.H2toFTU+...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RNGA_CO,1)'/1000.*shareOfNGAHy.*share.H2toFTU;

LH2Cost = sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'TLH2')).*capacity(:,:,ismember(allActiveExLoadLabels,'TLH2'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TLH2')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TLH2')).*permute(repmat(sum(results.TLH2_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHYG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHYG')).*permute(repmat(sum(results.SHYG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toLH2),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'IHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHYG'))).*repmat((share.H2toLH2),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toLH2),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMR')).*permute(repmat(sum(results.TSMR_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toLH2),1,length(systemParams.IndexYears)),2) +...sum(results.TMET_GA,1)').*(1-shareOfGasEl))+...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.TSMR_GAS_FOS,1)'.*share.H2toLH2 +...
    setup.importH2Cost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importH2_regas,1)'.* share.H2toLH2+...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RNGA_CO,1)'/1000.*shareOfNGAHy.*share.H2toLH2;

PetrFuelCost = opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWOO')).*sum(results.RWOO_RBDS,1)' +...opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBDS')).*sum(results.RBDS_DI,1)' +...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWWO')).*sum(results.RWWO_RBDS,1)' +...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RPET')).*sum(results.RPET_DI+results.RPET_KE,1)' +...
    setup.importFTCost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importFT,1)'+...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RPET_CO,1)'/1000.*shareOfPETTrans;

LNGCost = sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'TLNG')).*capacity(:,:,ismember(allActiveExLoadLabels,'TLNG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TLNG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TLNG')).*permute(repmat(sum(results.TLNG_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...sum(results.RBME,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGU'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGU')).*permute(repmat(sum(results.TBGU_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((shareOfGasTrans),1,length(systemParams.IndexYears)) +...sum(results.RBME,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGD'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGD')).*permute(repmat(sum(results.RBGA_FU,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(shareOfBGA_Trans,1,length(systemParams.IndexYears)) +...sum(results.RBGA_TBGD,1)').*(shareOfgasBGA).*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toLNG),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMR')).*permute(repmat(sum(results.TSMR_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toLNG),1,length(systemParams.IndexYears)) +...sum(results.TMET_GA,1)').*(1-shareOfGasEl))+...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCOS')).*permute(repmat(sum(results.TCOS_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toLNG),1,length(systemParams.IndexYears)) +....*sum(results.TCOS_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSC')).*permute(repmat(sum(results.TPSC_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toLNG),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSP')).*permute(repmat(sum(results.TPSP_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toLNG),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'SCO2'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SCO2')).*permute(repmat(sum(results.SCO2_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toLNG),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'ICO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'ICO2'))).*repmat((share.CO2toLNG),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHYG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHYG')).*permute(repmat(sum(results.SHYG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toLNG),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'IHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHYG'))).*repmat((share.H2toLNG),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TMET')).*capacity(:,:,ismember(allActiveExLoadLabels,'TMET'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TMET')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TMET')).*permute(repmat(sum(results.TMET_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((shareOfGasTrans),1,length(systemParams.IndexYears)),2) +...sum(results.TMET_GA,1)').*(1-shareOfGasEl))+...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBGA')).*sum(results.RBGA_FU,1)'.* shareOfBGA_Trans +...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.TLNG_GAS_FOS,1)' +...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RNGA_CO,1)'/1000.*shareOfNGATrans +...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.TSMR_GAS_FOS,1)'.*share.H2toLNG +...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RNGA_CO,1)'/1000.*shareOfNGAHy.*share.H2toLNG;

SNGCost = sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toCH4),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMR')).*permute(repmat(sum(results.TSMR_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toCH4),1,length(systemParams.IndexYears)) +...sum(results.TMET_GA,1)').*(1-shareOfGasEl))+...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'SCO2'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SCO2')).*permute(repmat(sum(results.SCO2_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toCH4),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'ICO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'ICO2'))).*repmat((share.CO2toCH4),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHYG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHYG')).*permute(repmat(sum(results.SHYG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toCH4),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'IHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHYG'))).*repmat((share.H2toCH4),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCOS')).*permute(repmat(sum(results.TCOS_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toCH4),1,length(systemParams.IndexYears)) +....*sum(results.TCOS_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSC')).*permute(repmat(sum(results.TPSC_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toCH4),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSP')).*permute(repmat(sum(results.TPSP_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toCH4),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TMET')).*capacity(:,:,ismember(allActiveExLoadLabels,'TMET'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TMET')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TMET')).*permute(repmat(sum(results.TMET_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])),2) +...sum(results.TMET_GA,1)').*(1-shareOfGasEl))+...
    0;

CDRCost = sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toCDR),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMR')).*permute(repmat(sum(results.TSMR_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toCDR),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMR')).*permute(repmat(sum(results.TSMR_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toCDR),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'SCO2'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SCO2')).*permute(repmat(sum(results.SCO2_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toCDR),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'ICO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'ICO2'))).*repmat((share.CO2toCDR),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHYG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHYG')).*permute(repmat(sum(results.SHYG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toCDR),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'IHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHYG'))).*repmat((share.H2toCDR),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCOS')).*permute(repmat(sum(results.TCOS_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toCDR),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSC')).*permute(repmat(sum(results.TPSC_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toCDR),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSP')).*permute(repmat(sum(results.TPSP_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toCDR),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRSL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRSL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRSL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRSL')).*permute(repmat(sum(results.TRSL_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRSS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRSS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRSS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRSS')).*permute(repmat(sum(results.TRSS_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRMI')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRMI'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRMI')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRMI')).*permute(repmat(sum(results.TRMI_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRMO')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRMO'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRMO')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRMO')).*permute(repmat(sum(results.TRMO_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRME')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRME'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRME')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRME')).*permute(repmat(sum(results.TRME_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRSi')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRSi'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRSi')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRSi')).*permute(repmat(sum(results.TRSi_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRAP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRAP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRAP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRAP')).*permute(repmat(sum(results.TRAP_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRAD')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRAD'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRAD')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRAD')).*permute(repmat(sum(results.TRAD_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TREW')).*capacity(:,:,ismember(allActiveExLoadLabels,'TREW'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TREW')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TREW')).*permute(repmat(sum(results.TREW_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRBC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRBC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRBC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRBC')).*permute(repmat(sum(results.TRBC_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRGE')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRGE'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRGE')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRGE')).*permute(repmat(sum(results.TRGE_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])),2) +...
    0;

LCOSsys = (storageCost + primaryLCOEreal.*storLosses_sys)./((allDem_sys));

if sum(isinf(LCOSsys))
    LCOSsys(isinf(LCOSsys))=0;
    warning('allDem_sys is negative')
end
LCOStot = (LCOSsys.*(allDem_sys)+LCOS_SC.*(sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand+resultsSC.EL_THHP+resultsSC.EL_THHR,1)'-el_cons_localpowerFromGrid-el_cons_localheatFromGrid))./((allDem_tot));
LCOStot(isinf(LCOStot))=0;

%% to be checked and fixed, used for 'Electricity direct'; '[TWh]' calculations
stored_electricity_sys = sum(results.EL_SBAT,1)'+sum(results.EL_SACA,1)'+sum(results.EL_SPHS,1)'+sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'+sum(results.EL_TWEL,1)'+sum(results.EL_SHOT,1)';
stored_electricity_tot = stored_electricity_sys + sum(resultsSC.EL_SBAR,1)'+ sum(resultsSC.EL_SBAC,1)'+ sum(resultsSC.EL_SBAI,1)';
%%
LCOSsys_average = sum(LCOSsys.*(allDem_sys))./sum((allDem_sys));
LCOStot_average = sum(LCOStot.*(allDem_tot))./sum((allDem_tot));

LCOSBsys = (0.02*sum(resultsSC.EL_EXCESS,1)')./((allDem_sys));
if sum(isinf(LCOSBsys))
    LCOSBsys(isinf(LCOSBsys))=0;
    warning('allDem_sys is negative')
end

LCOSBsys_average = sum(0.02*sum(resultsSC.EL_EXCESS,1)')./sum((allDem_sys));

primaryLCOEsysR = (primaryLCOEsys.*(allDem_sys) + primaryLCOEsys_average * (sum(importPower)'-sum(-exportPower)'))./((allDem_sys+(sum(importPower)'-sum(-exportPower)')));
primaryLCOEsysR_average = sum(primaryLCOEsysR.*((allDem_sys+(sum(importPower)'-sum(-exportPower)'))))/sum((allDem_sys+(sum(importPower)'-sum(-exportPower)')));

subs_cost = 0.002*(sum(resultsSC.EL_EXCESS,1)');

fullLCOEtot = primaryLCOEtot+LCOCtot+LCOStot+lcot_tot;
fullLCOEtot_av = primaryLCOEtot_average+LCOCtot_average+LCOStot_average+lcot_tot_average;
fullLCOEsys = primaryLCOEsys+LCOCsys+LCOSsys+lcot_sys+LCOSBsys;
fullLCOEsys_av = primaryLCOEsys_average+LCOCsys_average+LCOSsys_average+lcot_sys_average+LCOSBsys_average;


total_Cons_ElCost = allDem_tot.*fullLCOEtot + sum((imports_AC+imports_DC).*repmat(fullLCOEtot',length(systemParams.IndexNodes),1),2)+sum(exports_AC+exports_DC,2).*fullLCOEtot+(transmLosses_sys_AC+transmLosses_sys_DC).*fullLCOEtot;
full_Cons_LCOEtot = total_Cons_ElCost./(allDem_tot+sum(imports_AC')'+sum(imports_DC')'+sum(exports_AC')'+sum(exports_DC')'+transmLosses_sys_AC+transmLosses_sys_DC);
full_Cons_LCOEtot_av = sum(total_Cons_ElCost)/sum(allDem_tot+sum(imports_AC')'+sum(imports_DC')'+sum(exports_AC')'+sum(exports_DC')'+transmLosses_sys_AC+transmLosses_sys_DC);

system_Cons_ElCost = allDem_sys.*fullLCOEsys + sum((imports_AC+imports_DC).*repmat(fullLCOEsys',length(systemParams.IndexNodes),1),2)+sum(exports_AC+exports_DC,2).*fullLCOEsys+(transmLosses_sys_AC+transmLosses_sys_DC).*fullLCOEsys;
full_Cons_LCOEsys = system_Cons_ElCost./(allDem_sys+sum(imports_AC')'+sum(imports_DC')'+sum(exports_AC')'+sum(exports_DC')'+transmLosses_sys_AC+transmLosses_sys_DC);
full_Cons_LCOEsys_av = sum(system_Cons_ElCost)/sum(allDem_sys+sum(imports_AC')'+sum(imports_DC')'+sum(exports_AC')'+sum(exports_DC')'+transmLosses_sys_AC+transmLosses_sys_DC);

%% Heat
%Total
prod_heat_tot = sum(results.TCNG_HE,1)'+sum(results.TCOI_HE,1)'+sum(results.TCCO_HE,1)'+sum(results.TCBP_HE,1)'+sum(results.TCHP_HE,1)'+sum(results.TMSW_HE,1)'+...
    sum(results.TCBC_HE,1)'+sum(results.TCWC_HE,1)'+...
    sum(results.TDHR_HE,1)'+sum(results.TDHP_HE,1)'+sum(results.TDNG_HE,1)'+sum(results.TDOI_HE,1)'+sum(results.TDCO_HE,1)'+sum(results.TDBP_HE,1)'+...
    sum(resultsSC.RRSH_HE,1)'+sum(resultsSC.THHR_HE,1)'+sum(resultsSC.THHP_HE,1)'+sum(resultsSC.THNG_HE,1)'+sum(resultsSC.THOI_HE,1)'+sum(resultsSC.THBP_HE,1)'+sum(resultsSC.THBG_HE,1)'+...
    sum(results.Loss_HE,1)';

el_cons_heat = sum(results.EL_TDHR,1)'+sum(results.EL_TDHP,1)'+sum(results.EL_TDGE,1)';

%District

shareOfHeCoalGHG = sum(results.RHAR_TDCO+results.RHAR_LHIN,1)'./sum(results.RHAR_FU+0.01,1)';
shareOfHeCoalGHG(isnan(shareOfHeCoalGHG))=0;
shareOfHeOilGHG = sum(results.RPET_TDOI+resultsSC.RPET_FU+results.RPET_LHIN,1)'./sum(results.RPET_FU+0.01,1)' ;
shareOfHeOilGHG(isnan(shareOfHeOilGHG))=0;
shareOfHeGasGHG = sum(results.TDNG_GAS_FOS+results.LHIN_GAS_FOS+results.Pros_GAS_FOS,1)'./sum(results.RNGA_FU+0.01,1)';
shareOfHeGasGHG(isnan(shareOfHeGasGHG))=0;

primaryCost_HEdistr = sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.H2toCH4.*(shareOfGasDH+shareOfGasINH)+share.H2toINH,1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*shareOfGasEl +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMR')).*permute(repmat(sum(results.TSMR_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.H2toCH4.*(shareOfGasDH+shareOfGasINH)+share.H2toINH,1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCOS')).*permute(repmat(sum(results.TCOS_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*(shareOfGasDH+shareOfGasINH),1,length(systemParams.IndexYears)) +...sum(results.TCOS_GA,1)').*shareOfGasEl +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSC')).*permute(repmat(sum(results.TPSC_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*(shareOfGasDH+shareOfGasINH),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSP')).*permute(repmat(sum(results.TPSP_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*(shareOfGasDH+shareOfGasINH),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TMET')).*capacity(:,:,ismember(allActiveExLoadLabels,'TMET'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TMET')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TMET')).*permute(repmat(sum(results.TMET_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((shareOfGasDH+shareOfGasINH),1,length(systemParams.IndexYears)) + ...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGU'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGU')).*permute(repmat(sum(results.TBGU_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(((shareOfGasDH+shareOfGasINH)),1,length(systemParams.IndexYears)) +...sum(results.RBME,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGD'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGD')).*permute(repmat(sum(results.RBGA_FU,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((shareOfBGA_DH+shareOfBGA_INH),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'SCO2'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SCO2')).*permute(repmat(sum(results.SCO2_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toCH4.*(shareOfGasDH+shareOfGasINH)),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'ICO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'ICO2'))).*repmat((share.CO2toCH4.*(shareOfGasDH+shareOfGasINH)),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHYG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHYG')).*permute(repmat(sum(results.SHYG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toCH4.*(shareOfGasDH+shareOfGasINH)+share.H2toINH),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'IHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHYG'))).*repmat((share.H2toCH4.*(shareOfGasDH+shareOfGasINH)+share.H2toINH),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'RCSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'RCSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RCSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RCSP')).*permute(repmat(sum(results.RCSP_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CSPtoHe),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'RDSH')).*capacity(:,:,ismember(allActiveExLoadLabels,'RDSH'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RDSH')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RDSH')).*permute(repmat(sum(results.RDSH_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDNG')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDNG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDNG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDNG')).*permute(repmat(sum(results.TDNG_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDOI')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDOI'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDOI')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDOI')).*permute(repmat(sum(results.TDOI_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDCO')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDCO'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDCO')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDCO')).*permute(repmat(sum(results.TDCO_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDBP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDBP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDBP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDBP')).*permute(repmat(sum(results.TDBP_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TDHR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDHR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDHR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDHR')).*permute(repmat(sum(results.TDHR_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDHP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDHP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDHP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDHP')).*permute(repmat(sum(results.TDHP_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IHOT')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHOT'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SHOT')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHOT'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHOT')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHOT')).*permute(repmat(sum(results.SHOT_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SHOT_HE,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SDHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SDHS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SDHS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SDHS')).*permute(repmat(sum(results.SDHS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SHOT_HE,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IDHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'IDHS'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SSHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SSHS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SSHS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SSHS')).*permute(repmat(sum(results.SSHS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SHOT_HE,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'ISHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'ISHS')))*(setup.Heat.Flag)+...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SBAT')).*capacity(:,:,ismember(allActiveExLoadLabels,'SBAT'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SBAT')).*permute(repmat(sum(results.SBAT_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SBAT_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IBAT')).*capacity(:,:,ismember(allActiveExLoadLabels,'IBAT'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SACA')).*capacity(:,:,ismember(allActiveExLoadLabels,'SACA'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SACA')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SACA')).*permute(repmat(sum(results.SACA_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SACA_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IACA')).*capacity(:,:,ismember(allActiveExLoadLabels,'IACA'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SPHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SPHS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SPHS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SPHS')).*permute(repmat(sum(results.SPHS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SPHS_EL,1)' ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IPHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'IPHS'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SGAS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SGAS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SGAS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SGAS')).*permute(repmat(sum(results.SGAS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SGAS_TCCG,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IGAS')).*capacity(:,:,ismember(allActiveExLoadLabels,'IGAS'))).*repmat(shareOfHeatInElStorage,1,length(systemParams.IndexYears))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDGE')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDGE'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDGE')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDGE')).*permute(repmat(sum(results.TDGE_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) ,2) +...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWOO')).*sum(results.RWOO_TDBP+results.RWOO_LHIN,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWWO')).*sum(results.RWWO_TDBP+results.RWWO_LHIN,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RHAR')).*sum(results.RHAR_TDCO+results.RHAR_LHIN,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RPET')).*sum(results.RPET_TDOI+results.RPET_LHIN,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBGA')).*sum(results.RBGA_FU,1)'.*shareOfBGA_DH + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.TDNG_GAS_FOS+results.LHIN_GAS_FOS,1)' +...
    setup.importLNGCost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importLNG_regas,1)'.* (shareOfGasDH+shareOfGasINH)+...
    setup.importH2Cost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importH2_regas,1)'.* share.H2toINH+...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RHAR_CO/10^3,1)'.*shareOfHeCoalGHG + ...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RPET_CO/10^3,1)'.*shareOfHeOilGHG + ...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RNGA_CO/10^3,1)'.*shareOfHeGasGHG;



prod_distheat_tot = sum(results.TCNG_HE,1)'+sum(results.TCOI_HE,1)'+sum(results.TCCO_HE,1)'+...
    sum(results.TCBP_HE,1)'+sum(results.TCHP_HE,1)'+sum(results.TMSW_HE,1)'+...
    sum(results.TCBC_HE,1)'+sum(results.TCWC_HE,1)'+...
    sum(results.TDHR_HE,1)'+sum(results.TDHP_HE,1)'+...
    sum(results.TDNG_HE,1)'+sum(results.TDOI_HE,1)'+sum(results.TDCO_HE,1)'+...
    sum(results.TDBP_HE,1)'+sum(results.TDGE_HE,1)'+...
    sum(results.RCSP_HE,1)'+...
    sum(results.RDSH_HE,1)'+...
    +sum(results.Loss_HE,1)';

el_cons_distrheat = sum(results.EL_TDHR,1)'+sum(results.EL_TDHP,1)'+sum(results.EL_TDGE,1)';

distrheatDemand = shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHIN'),:,1,:)),3)+...
    (shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHSP'),:,1,:)),3)+...
    shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHDW'),:,1,:)),3)).*systemParams.Heat.shareOfDistrHeat(:,find(systemParams.IndexYears == costYear))+...
    sum(results.HE_TICB+results.HE_TICI+results.HE_TISB+results.HE_TISH+results.HE_TISR+results.HE_TISE+results.HE_TIAA+results.HE_TIAR+results.HE_TIPP,1)'+...
    sum(results.HE_TRSi+results.HE_TRME+results.HE_TRGE,1)'+...
    sum(results.HE_TCOS	+results.HE_TPSC+results.HE_TPSP,1)';

distrHeatProductionCost = primaryCost_HEdistr + ...
    +(el_cons_distrheat+sum(results.EL_TWEL,1)'.*(share.H2toCH4.*(shareOfGasDH+shareOfGasINH)+share.H2toINH)).*(primaryLCOEsys+LCOCsys+LCOSsys+lcot_sys+LCOSBsys);
distrHeatProductionCost = primaryCost_HEdistr + ...
    +(el_cons_distrheat+sum(results.EL_TWEL,1)'.*(share.H2toCH4.*(shareOfGasDH+shareOfGasINH)+share.H2toINH)).*(full_Cons_LCOEsys);


distrHeatCost = distrHeatProductionCost./distrheatDemand;
distrHeatCost(isnan(distrHeatCost))=0;
distrHeatCost(isinf(distrHeatCost))=0;

%Local
% if setup.SC.Flag

genShares_SHHS = (resultsSC.Cap_SHHS)./repmat((resultsSC.OPT_SIZE_SHHS)',1,size(systemParams.IndexYears,2));
genShares_SHHS(isinf(genShares_SHHS))=0;
genShares_SHHS(isnan(genShares_SHHS))=0;
genShares_RRSH = (resultsSC.Cap_RRSH)./repmat((resultsSC.OPT_SIZE_RRSH)',1,size(systemParams.IndexYears,2));
genShares_RRSH(isinf(genShares_RRSH))=0;
genShares_RRSH(isnan(genShares_RRSH))=0;
genShares_THHR = (resultsSC.Cap_THHR)./repmat((resultsSC.OPT_SIZE_THHR)',1,size(systemParams.IndexYears,2));
genShares_THHR(isinf(genShares_THHR))=0;
genShares_THHR(isnan(genShares_THHR))=0;
genShares_THHP = (resultsSC.Cap_THHP)./repmat((resultsSC.OPT_SIZE_THHP)',1,size(systemParams.IndexYears,2));
genShares_THHP(isinf(genShares_THHP))=0;
genShares_THHP(isnan(genShares_THHP))=0;
genShares_THNG = (resultsSC.Cap_THNG)./repmat((resultsSC.OPT_SIZE_THNG)',1,size(systemParams.IndexYears,2));
genShares_THNG(isinf(genShares_THNG))=0;
genShares_THNG(isnan(genShares_THNG))=0;
genShares_THOI = (resultsSC.Cap_THOI)./repmat((resultsSC.OPT_SIZE_THOI)',1,size(systemParams.IndexYears,2));
genShares_THOI(isinf(genShares_THOI))=0;
genShares_THOI(isnan(genShares_THOI))=0;
genShares_THBP = (resultsSC.Cap_THBP)./repmat((resultsSC.OPT_SIZE_THBP)',1,size(systemParams.IndexYears,2));
genShares_THBP(isinf(genShares_THBP))=0;
genShares_THBP(isnan(genShares_THBP))=0;
genShares_THBG = (resultsSC.Cap_THBG)./repmat((resultsSC.OPT_SIZE_THBG)',1,size(systemParams.IndexYears,2));
genShares_THBG(isinf(genShares_THBG))=0;
genShares_THBG(isnan(genShares_THBG))=0;


primaryCost_HElocal = sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.H2toCH4.*shareOfGasIH,1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*shareOfGasEl +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMR')).*permute(repmat(sum(results.TSMR_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.H2toCH4.*(shareOfGasIH),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCOS')).*permute(repmat(sum(results.TCOS_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*shareOfGasIH,1,length(systemParams.IndexYears)) +...sum(results.TCOS_GA,1)').*shareOfGasEl +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSC')).*permute(repmat(sum(results.TPSC_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*shareOfGasIH,1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSP')).*permute(repmat(sum(results.TPSP_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*shareOfGasIH,1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TMET')).*capacity(:,:,ismember(allActiveExLoadLabels,'TMET'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TMET')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TMET')).*permute(repmat(sum(results.TMET_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(shareOfGasIH,1,length(systemParams.IndexYears)) + ...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGU'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGU')).*permute(repmat(sum(results.TBGU_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((shareOfGasIH),1,length(systemParams.IndexYears)) +...sum(results.RBME,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGD'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGD')).*permute(repmat(sum(results.RBGA_FU,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(shareOfBGA_IH,1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'SCO2'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SCO2')).*permute(repmat(sum(results.SCO2_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toCH4.*(shareOfGasIH)),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'ICO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'ICO2'))).*repmat((share.CO2toCH4.*(shareOfGasIH)),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHYG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHYG')).*permute(repmat(sum(results.SHYG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toCH4.*(shareOfGasIH)),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'IHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHYG'))).*repmat((share.H2toCH4.*(shareOfGasIH)),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SHHS')).*(resultsSC.Cap_SHHS) +opex_var(:,:,ismember(allActiveExLoadLabels,'SHHS')).*repmat(sum(resultsSC.SHHS_EL)',1,size(systemParams.IndexYears,2)).* genShares_SHHS +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IHHS')).*(resultsSC.Cap_IHHS) )*(setup.Heat.Flag)+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RRSH')).*(resultsSC.Cap_RRSH) +opex_var(:,:,ismember(allActiveExLoadLabels,'RRSH')).*repmat(sum(resultsSC.RRSH_HE)',1,size(systemParams.IndexYears,2)).* genShares_RRSH +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THHR')).*(resultsSC.Cap_THHR) +opex_var(:,:,ismember(allActiveExLoadLabels,'THHR')).*repmat(sum(resultsSC.THHR_HE)',1,size(systemParams.IndexYears,2)).* genShares_THHR +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THHP')).*(resultsSC.Cap_THHP) +opex_var(:,:,ismember(allActiveExLoadLabels,'THHP')).*repmat(sum(resultsSC.THHP_HE)',1,size(systemParams.IndexYears,2)).* genShares_THHP +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THNG')).*(resultsSC.Cap_THNG) +opex_var(:,:,ismember(allActiveExLoadLabels,'THNG')).*repmat(sum(resultsSC.THNG_HE)',1,size(systemParams.IndexYears,2)).* genShares_THNG +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THOI')).*(resultsSC.Cap_THOI) +opex_var(:,:,ismember(allActiveExLoadLabels,'THOI')).*repmat(sum(resultsSC.THOI_HE)',1,size(systemParams.IndexYears,2)).* genShares_THOI +...%    opex_capex(:,:,ismember(allActiveExLoadLabels,'THCO')).*resultsSC.Cap_THCO +opex_var(:,:,ismember(allActiveExLoadLabels,'THCO')).*sum(resultsSC.THCO_HE).*(resultsSC.Cap_THCO/resultsSC.OPT_SIZE_THCO) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THBP')).*(resultsSC.Cap_THBP) +opex_var(:,:,ismember(allActiveExLoadLabels,'THBP')).*repmat(sum(resultsSC.THBP_HE)',1,size(systemParams.IndexYears,2)).* genShares_THBP +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THBG')).*(resultsSC.Cap_THBG) +opex_var(:,:,ismember(allActiveExLoadLabels,'THBG')).*repmat(sum(resultsSC.THBG_HE)',1,size(systemParams.IndexYears,2)).* genShares_THBG,2) +...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWOO')).*sum(resultsSC.RWOO_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWWO')).*sum(resultsSC.RWWO_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RHAR')).*sum(resultsSC.RHAR_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RPET')).*sum(resultsSC.RPET_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBGA')).*sum(results.RBGA_FU,1)'.*shareOfBGA_IH + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(resultsSC.RNGA_FU,1)';


prod_localheat_tot = sum(resultsSC.RRSH_HE,1)'+...
    sum(resultsSC.THHR_HE,1)'+sum(resultsSC.THHP_HE,1)'+sum(resultsSC.THNG_HE,1)'+sum(resultsSC.THOI_HE,1)'+sum(resultsSC.THBP_HE,1)'+sum(resultsSC.THBG_HE,1)';

prod_heat_tot = prod_distheat_tot+prod_localheat_tot;

try results.HE_EXCESS_Local;
catch
    results.HE_EXCESS_Local = 0*resultsSC.HE_EXCESS_Local;
end
excess_localheat_tot = sum(resultsSC.HE_EXCESS_Local,1)'+sum(results.HE_EXCESS_Local,1)';

el_cons_localheat = sum(resultsSC.EL_THHR,1)'+sum(resultsSC.EL_THHP,1)';

localheatDemand = (shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHSP'),:,1,:)),3)+...
    shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHDW'),:,1,:)),3)).*(1-systemParams.Heat.shareOfDistrHeat(:,find(systemParams.IndexYears == costYear)))+...
    shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHBC'),:,1,:)),3);

localHeatProductionCost = primaryCost_HElocal + ...
    +(el_cons_localheat+sum(results.EL_TWEL,1)'.*(share.H2toCH4.*shareOfGasIH)).*(primaryLCOEsys+LCOCsys+LCOSsys+lcot_sys+LCOSBsys); % must be redone
localHeatProductionCost = primaryCost_HElocal + ...
    +(el_cons_localheat+sum(results.EL_TWEL,1)'.*(share.H2toCH4.*shareOfGasIH)).*(full_Cons_LCOEsys); % has been redone


localHeatCost = localHeatProductionCost./localheatDemand;
localHeatCost(isnan(localHeatCost))=0;
localHeatCost(isinf(localHeatCost))=0;

heatDemand_IND = shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHIN'),:,1,:)),3);
heatDemand_HSP = shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHSP'),:,1,:)),3);
heatDemand_HDW = shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHDW'),:,1,:)),3);
heatDemand_HBC = shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHBC'),:,1,:)),3);

HeatCost = (localHeatProductionCost+distrHeatProductionCost)./(localheatDemand+distrheatDemand);
HeatCost_av = sum(localHeatProductionCost+distrHeatProductionCost)./sum(localheatDemand+distrheatDemand);


HeatCost_local = (localHeatProductionCost)./(localheatDemand);
HeatCost_district = (distrHeatProductionCost)./(distrheatDemand);
heatDemand =(localheatDemand+distrheatDemand);



primaryCost_HEtot = sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.H2toCH4.*(shareOfGasDH+shareOfGasINH)+share.H2toINH,1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*shareOfGasEl +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMR')).*permute(repmat(sum(results.TSMR_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.H2toCH4.*(shareOfGasDH+shareOfGasINH)+share.H2toINH,1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCOS')).*permute(repmat(sum(results.TCOS_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*(shareOfGasDH+shareOfGasINH),1,length(systemParams.IndexYears)) +...sum(results.TCOS_GA,1)').*shareOfGasEl +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSC')).*permute(repmat(sum(results.TPSC_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*(shareOfGasDH+shareOfGasINH),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSP')).*permute(repmat(sum(results.TPSP_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*(shareOfGasDH+shareOfGasINH),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TMET')).*capacity(:,:,ismember(allActiveExLoadLabels,'TMET'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TMET')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TMET')).*permute(repmat(sum(results.TMET_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((shareOfGasDH+shareOfGasINH),1,length(systemParams.IndexYears)) + ...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGU'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGU')).*permute(repmat(sum(results.TBGU_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((shareOfGasDH+shareOfGasINH),1,length(systemParams.IndexYears)) +...sum(results.RBME,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGD'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGD')).*permute(repmat(sum(results.RBGA_FU,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((shareOfBGA_DH+shareOfBGA_INH),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'RCSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'RCSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RCSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RCSP')).*permute(repmat(sum(results.RCSP_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CSPtoHe),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'RDSH')).*capacity(:,:,ismember(allActiveExLoadLabels,'RDSH'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RDSH')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RDSH')).*permute(repmat(sum(results.RDSH_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDNG')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDNG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDNG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDNG')).*permute(repmat(sum(results.TDNG_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDOI')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDOI'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDOI')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDOI')).*permute(repmat(sum(results.TDOI_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDCO')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDCO'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDCO')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDCO')).*permute(repmat(sum(results.TDCO_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDBP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDBP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDBP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDBP')).*permute(repmat(sum(results.TDBP_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TDHR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDHR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDHR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDHR')).*permute(repmat(sum(results.TDHR_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDHP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDHP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDHP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDHP')).*permute(repmat(sum(results.TDHP_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IHOT')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHOT'))+...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SHOT')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHOT'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHOT')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHOT')).*permute(repmat(sum(results.SHOT_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SHOT_HE,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SDHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SDHS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SDHS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SDHS')).*permute(repmat(sum(results.SDHS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...sum(results.SHOT_HE,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IDHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'IDHS'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SSHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SSHS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SSHS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SSHS')).*permute(repmat(sum(results.SSHS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SHOT_HE,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'ISHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'ISHS')))*(setup.Heat.Flag)+...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SHHS')).*(resultsSC.Cap_SHHS) +opex_var(:,:,ismember(allActiveExLoadLabels,'SHHS')).*repmat(sum(resultsSC.SHHS_EL)',1,size(systemParams.IndexYears,2)).* genShares_SHHS +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IHHS')).*(resultsSC.Cap_IHHS) )*(setup.Heat.Flag)+...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'SCO2'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SCO2')).*permute(repmat(sum(results.SCO2_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.CO2toCH4.*(shareOfGasDH+shareOfGasIH+shareOfGasINH)),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'ICO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'ICO2'))).*repmat((share.CO2toCH4.*(shareOfGasDH+shareOfGasIH+shareOfGasINH)),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHYG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHYG')).*permute(repmat(sum(results.SHYG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((share.H2toCH4.*(shareOfGasDH+shareOfGasIH+shareOfGasINH)+share.H2toINH),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'IHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHYG'))).*repmat((share.H2toCH4.*(shareOfGasDH+shareOfGasIH+shareOfGasINH))+share.H2toINH,1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SBAT')).*capacity(:,:,ismember(allActiveExLoadLabels,'SBAT'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SBAT')).*permute(repmat(sum(results.SBAT_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SBAT_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IBAT')).*capacity(:,:,ismember(allActiveExLoadLabels,'IBAT'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SACA')).*capacity(:,:,ismember(allActiveExLoadLabels,'SACA'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SACA')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SACA')).*permute(repmat(sum(results.SACA_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SACA_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IACA')).*capacity(:,:,ismember(allActiveExLoadLabels,'IACA'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SPHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SPHS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SPHS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SPHS')).*permute(repmat(sum(results.SPHS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SPHS_EL,1)' ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IPHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'IPHS'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SGAS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SGAS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SGAS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SGAS')).*permute(repmat(sum(results.SGAS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SGAS_TCCG,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IGAS')).*capacity(:,:,ismember(allActiveExLoadLabels,'IGAS'))).*repmat(shareOfHeatInElStorage,1,length(systemParams.IndexYears))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDGE')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDGE'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDGE')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDGE')).*permute(repmat(sum(results.TDGE_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) ,2) +...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWOO')).*sum(results.RWOO_TDBP+results.RWOO_LHIN,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWWO')).*sum(results.RWWO_TDBP+results.RWWO_LHIN,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RHAR')).*sum(results.RHAR_TDCO+results.RHAR_LHIN,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RPET')).*sum(results.RPET_TDOI+results.RPET_LHIN,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBGA')).*sum(results.RBGA_FU,1)'.*shareOfBGA_DH + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.TDNG_GAS_FOS+results.LHIN_GAS_FOS,1)' +...ma
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RHAR_CO/10^3,1)'.*shareOfHeCoalGHG + ...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RPET_CO/10^3,1)'.*shareOfHeOilGHG + ...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RNGA_CO/10^3,1)'.*shareOfHeGasGHG +...
    sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.H2toCH4.*shareOfGasIH,1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*shareOfGasEl +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMR')).*permute(repmat(sum(results.TSMR_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.H2toCH4.*(shareOfGasIH),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCOS')).*permute(repmat(sum(results.TCOS_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*shareOfGasIH,1,length(systemParams.IndexYears)) +...sum(results.TCOS_GA,1)').*shareOfGasEl +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSC')).*permute(repmat(sum(results.TPSC_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*shareOfGasIH,1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSP')).*permute(repmat(sum(results.TPSP_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*shareOfGasIH,1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TMET')).*capacity(:,:,ismember(allActiveExLoadLabels,'TMET'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TMET')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TMET')).*permute(repmat(sum(results.TMET_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(shareOfGasIH,1,length(systemParams.IndexYears)) + ...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGU'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGU')).*permute(repmat(sum(results.TBGU_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((shareOfGasIH),1,length(systemParams.IndexYears)) +...sum(results.RBME,1)').*(1-shareOfGasEl) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGD'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGD')).*permute(repmat(sum(results.RBGA_FU,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(shareOfBGA_IH,1,length(systemParams.IndexYears)) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RRSH')).*(resultsSC.Cap_RRSH) +opex_var(:,:,ismember(allActiveExLoadLabels,'RRSH')).*repmat(sum(resultsSC.RRSH_HE)',1,size(systemParams.IndexYears,2)).* genShares_RRSH +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THHR')).*(resultsSC.Cap_THHR) +opex_var(:,:,ismember(allActiveExLoadLabels,'THHR')).*repmat(sum(resultsSC.THHR_HE)',1,size(systemParams.IndexYears,2)).* genShares_THHR +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THHP')).*(resultsSC.Cap_THHP) +opex_var(:,:,ismember(allActiveExLoadLabels,'THHP')).*repmat(sum(resultsSC.THHP_HE)',1,size(systemParams.IndexYears,2)).* genShares_THHP +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THNG')).*(resultsSC.Cap_THNG) +opex_var(:,:,ismember(allActiveExLoadLabels,'THNG')).*repmat(sum(resultsSC.THNG_HE)',1,size(systemParams.IndexYears,2)).* genShares_THNG +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THOI')).*(resultsSC.Cap_THOI) +opex_var(:,:,ismember(allActiveExLoadLabels,'THOI')).*repmat(sum(resultsSC.THOI_HE)',1,size(systemParams.IndexYears,2)).* genShares_THOI +...%    opex_capex(:,:,ismember(allActiveExLoadLabels,'THCO')).*resultsSC.Cap_THCO +opex_var(:,:,ismember(allActiveExLoadLabels,'THCO')).*sum(resultsSC.THCO_HE).*(resultsSC.Cap_THCO/resultsSC.OPT_SIZE_THCO) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THBP')).*(resultsSC.Cap_THBP) +opex_var(:,:,ismember(allActiveExLoadLabels,'THBP')).*repmat(sum(resultsSC.THBP_HE)',1,size(systemParams.IndexYears,2)).* genShares_THBP +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THBG')).*(resultsSC.Cap_THBG) +opex_var(:,:,ismember(allActiveExLoadLabels,'THBG')).*repmat(sum(resultsSC.THBG_HE)',1,size(systemParams.IndexYears,2)).* genShares_THBG,2) +...
    setup.importLNGCost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importLNG_regas,1)'.* (shareOfGasDH+shareOfGasINH)+...
    setup.importH2Cost(1:length(systemParams.IndexNodes),yearNumb)/1000 .*sum(results.importH2_regas,1)'.* share.H2toINH+...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWOO')).*sum(resultsSC.RWOO_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWWO')).*sum(resultsSC.RWWO_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RHAR')).*sum(resultsSC.RHAR_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RPET')).*sum(resultsSC.RPET_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBGA')).*sum(results.RBGA_FU,1)'.*shareOfBGA_IH + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(resultsSC.RNGA_FU,1)';


primaryCost_HE_stor = sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.H2toCH4.*shareOfGasDH,1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*shareOfGasEl +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCOS')).*permute(repmat(sum(results.TCOS_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*shareOfGasDH,1,length(systemParams.IndexYears)) +...sum(results.TCOS_GA,1)').*shareOfGasEl +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSC')).*permute(repmat(sum(results.TPSC_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*shareOfGasDH,1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSP')).*permute(repmat(sum(results.TPSP_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(share.CO2toCH4.*shareOfGasDH,1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TMET')).*capacity(:,:,ismember(allActiveExLoadLabels,'TMET'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TMET')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TMET')).*permute(repmat(sum(results.TMET_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(shareOfGasDH,1,length(systemParams.IndexYears)) + ...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGU'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGU')).*permute(repmat(sum(results.TBGU_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((shareOfGasDH),1,length(systemParams.IndexYears)) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'IHOT')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHOT'))+...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SHOT')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHOT'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHOT')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHOT')).*permute(repmat(sum(results.SHOT_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SHOT_HE,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SDHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SDHS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SDHS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SDHS')).*permute(repmat(sum(results.SDHS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...sum(results.SHOT_HE,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IDHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'IDHS')))*(setup.Heat.Flag)+...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'SBAT')).*capacity(:,:,ismember(allActiveExLoadLabels,'SBAT'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SBAT')).*permute(repmat(sum(results.SBAT_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SBAT_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IBAT')).*capacity(:,:,ismember(allActiveExLoadLabels,'IBAT'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SACA')).*capacity(:,:,ismember(allActiveExLoadLabels,'SACA'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SACA')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SACA')).*permute(repmat(sum(results.SACA_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SACA_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IACA')).*capacity(:,:,ismember(allActiveExLoadLabels,'IACA'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SPHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SPHS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SPHS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SPHS')).*permute(repmat(sum(results.SPHS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SPHS_EL,1)' ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IPHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'IPHS'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SGAS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SGAS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SGAS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SGAS')).*permute(repmat(sum(results.SGAS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SGAS_TCCG,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IGAS')).*capacity(:,:,ismember(allActiveExLoadLabels,'IGAS'))).*repmat(shareOfHeatInElStorage,1,length(systemParams.IndexYears)),2);

primaryCost_HE_main = (primaryCost_HEtot - primaryCost_HE_stor);


HeatCostStor = primaryCost_HE_stor./(localheatDemand+distrheatDemand);
HeatCostMain = (primaryCost_HE_main+el_cons_localheat.*(primaryLCOEsys+LCOCsys+LCOSsys+lcot_sys+LCOSBsys)+el_cons_distrheat.*(primaryLCOEsys+LCOCsys+LCOSsys+lcot_sys+LCOSBsys))./(localheatDemand+distrheatDemand);
HeatCostMain = (primaryCost_HE_main+el_cons_localheat.*(primaryLCOEtot+LCOCtot+LCOStot+lcot_tot)+el_cons_distrheat.*(primaryLCOEtot+LCOCtot+LCOStot+lcot_tot))./(localheatDemand+distrheatDemand);
HeatCostMain = (primaryCost_HE_main+el_cons_localheat.*(full_Cons_LCOEsys)+el_cons_distrheat.*(full_Cons_LCOEsys))./(localheatDemand+distrheatDemand);
HeatCostMain = (primaryCost_HE_main+el_cons_localheat.*(full_Cons_LCOEtot)+el_cons_distrheat.*(full_Cons_LCOEtot))./(localheatDemand+distrheatDemand);

%% Total ann cost

totAnCost = (sum(opex_capex(:,:,ismember(allActiveExLoadLabels,'RPVO')).*capacity(:,:,ismember(allActiveExLoadLabels,'RPVO'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RPVO')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RPVO')).*permute(repmat(sum(results.RPVO_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RPVA')).*capacity(:,:,ismember(allActiveExLoadLabels,'RPVA'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RPVA')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RPVA')).*permute(repmat(sum(results.RPVA_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RPBO')).*capacity(:,:,ismember(allActiveExLoadLabels,'RPBO'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RPBO')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RPBO')).*permute(repmat(sum(results.RPBO_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RPBA')).*capacity(:,:,ismember(allActiveExLoadLabels,'RPBA'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RPBA')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RPBA')).*permute(repmat(sum(results.RPBA_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RPVF')).*capacity(:,:,ismember(allActiveExLoadLabels,'RPVF'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RPVF')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RPVF')).*permute(repmat(sum(results.RPVF_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RPBV')).*capacity(:,:,ismember(allActiveExLoadLabels,'RPBV'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RPBV')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RPBV')).*permute(repmat(sum(results.RPBV_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RWIN')).*capacity(:,:,ismember(allActiveExLoadLabels,'RWIN'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RWIN')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RWIN')).*permute(repmat(sum(results.RWIN_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.RWIN_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'ROWI')).*capacity(:,:,ismember(allActiveExLoadLabels,'ROWI'))+opex_var(:,:,ismember(allActiveExLoadLabels,'ROWI')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'ROWI')).*permute(repmat(sum(results.ROWI_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.ROWI_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RWIO')).*capacity(:,:,ismember(allActiveExLoadLabels,'RWIO'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RWIO')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RWIO')).*permute(repmat(sum(results.RWIO_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.RWIN_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RWAV')).*capacity(:,:,ismember(allActiveExLoadLabels,'RWAV'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RWAV')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RWAV')).*permute(repmat(sum(results.RWAV_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.RWAV_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RRRI')).*capacity(:,:,ismember(allActiveExLoadLabels,'RRRI'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RRRI')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RRRI')).*permute(repmat(sum(results.RRRI_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.RRRI_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'HDAM')).*capacity(:,:,ismember(allActiveExLoadLabels,'HDAM'))+opex_var(:,:,ismember(allActiveExLoadLabels,'HDAM')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'HDAM')).*permute(repmat(sum(results.HDAM_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.HDAM_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'RCSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'RCSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RCSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RCSP')).*permute(repmat(sum(results.RCSP_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TCSP_HE,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TSTU')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSTU'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSTU')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSTU')).*permute(repmat(sum(results.TSTU_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TSTU_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TGEO')).*capacity(:,:,ismember(allActiveExLoadLabels,'TGEO'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TGEO')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TGEO')).*permute(repmat(sum(results.TGEO_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TGEO_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TBPP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBPP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBPP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBPP')).*permute(repmat(sum(results.TBPP_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TBPP_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TBCS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBCS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBCS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBCS')).*permute(repmat(sum(results.TBCS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TBPP_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THPP')).*capacity(:,:,ismember(allActiveExLoadLabels,'THPP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'THPP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'THPP')).*permute(repmat(sum(results.THPP_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.THPP_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THCS')).*capacity(:,:,ismember(allActiveExLoadLabels,'THCS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'THCS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'THCS')).*permute(repmat(sum(results.THCS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.THPP_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TICG')).*capacity(:,:,ismember(allActiveExLoadLabels,'TICG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TICG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TICG')).*permute(repmat(sum(results.TICG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TICG_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TICM')).*capacity(:,:,ismember(allActiveExLoadLabels,'TICM'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TICM')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TICM')).*permute(repmat(sum(results.TICM_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TICG_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TNUC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TNUC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TNUC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TNUC')).*permute(repmat(sum(results.TNUC_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TNUC_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TMSW')).*capacity(:,:,ismember(allActiveExLoadLabels,'TMSW'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TMSW')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TMSW')).*permute(repmat(sum(results.TMSW_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TMSW_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCWC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCWC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCWC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCWC')).*permute(repmat(sum(results.TCWC_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TMSW_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCHP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCHP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCHP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCHP')).*permute(repmat(sum(results.TCHP_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TCHP_EL,1)' + ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCNG')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCNG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCNG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCNG')).*permute(repmat(sum(results.TCNG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCOI')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOI'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCOI')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCOI')).*permute(repmat(sum(results.TCOI_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCCO')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCCO'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCCO')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCCO')).*permute(repmat(sum(results.TCCO_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCBP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCBP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCBP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCBP')).*permute(repmat(sum(results.TCBP_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCBC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCBC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCBC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCBC')).*permute(repmat(sum(results.TCBC_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGD'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGD')).*permute(repmat(sum(results.RBGA_FU,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGU'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGU')).*permute(repmat(sum(results.TBGU_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.RBGA_TBGD,1)').*(1-shareOfgasBGA) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TGCS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TGCS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TGCS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TGCS')).*permute(repmat(sum(results.TGCS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TCCG_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCCG')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCCG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCCG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCCG')).*permute(repmat(sum(results.TCCG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +....*sum(results.TCCG_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TOCG')).*capacity(:,:,ismember(allActiveExLoadLabels,'TOCG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TOCG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TOCG')).*permute(repmat(sum(results.TOCG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]),2)+....*sum(results.TOCG_EL,1)').*shareOfGT_NGA);
    opex_capexRehab_RoR(:,yearNumb).*(RoR_LL)'+...
    opex_capexRehab_Dam(:,yearNumb).*(Dam_LL)'+...
    sum(opex_capex(:,:,ismember(allActiveExLoadLabels,'SBAT')).*capacity(:,:,ismember(allActiveExLoadLabels,'SBAT'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SBAT')).*permute(repmat(sum(results.SBAT_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SBAT_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SACA')).*capacity(:,:,ismember(allActiveExLoadLabels,'SACA'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SACA')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SACA')).*permute(repmat(sum(results.SACA_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SACA_EL,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SPHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SPHS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SPHS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SPHS')).*permute(repmat(sum(results.SPHS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SPHS_EL,1)' ...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SHOT')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHOT'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHOT')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHOT')).*permute(repmat(sum(results.SHOT_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SHOT_HE,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SBGA')).*capacity(:,:,ismember(allActiveExLoadLabels,'SBGA'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SBGA')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SBGA')).*permute(repmat(sum(results.SBGA_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SBGA_TCHP,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SGAS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SGAS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SGAS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SGAS')).*permute(repmat(sum(results.SGAS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SGAS_TCCG,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'SHYG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SHYG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SHYG')).*permute(repmat(sum(results.SHYG_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SHOT_HE,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'SCO2'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SCO2')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SCO2')).*permute(repmat(sum(results.SCO2_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SBGA_TCHP,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'SDHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'SDHS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'SDHS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'SDHS')).*permute(repmat(sum(results.SDHS_EL,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.SBGA_TCHP,1)' +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IHYG')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHYG'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'ICO2')).*capacity(:,:,ismember(allActiveExLoadLabels,'ICO2'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IBAT')).*capacity(:,:,ismember(allActiveExLoadLabels,'IBAT'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IACA')).*capacity(:,:,ismember(allActiveExLoadLabels,'IACA'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IPHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'IPHS'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IGAS')).*capacity(:,:,ismember(allActiveExLoadLabels,'IGAS'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IHOT')).*capacity(:,:,ismember(allActiveExLoadLabels,'IHOT'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'IDHS')).*capacity(:,:,ismember(allActiveExLoadLabels,'IDHS'))+...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.TWEL_GA,1)').*shareOfGasEl +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCOS')).*permute(repmat(sum(results.TCOS_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.TCOS_GA,1)').*shareOfGasEl +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TMET')).*capacity(:,:,ismember(allActiveExLoadLabels,'TMET'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TMET')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TMET')).*permute(repmat(sum(results.TMET_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.TMET_GA,1)').*shareOfGasEl);
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TFTU')).*capacity(:,:,ismember(allActiveExLoadLabels,'TFTU'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TFTU')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TFTU')).*permute(repmat(sum(results.TFTU_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.TCOS_GA,1)').*shareOfGasEl +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TLNG')).*capacity(:,:,ismember(allActiveExLoadLabels,'TLNG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TLNG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TLNG')).*permute(repmat(sum(results.TLNG_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.TMET_GA,1)').*shareOfGasEl);
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TLH2')).*capacity(:,:,ismember(allActiveExLoadLabels,'TLH2'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TLH2')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TLH2')).*permute(repmat(sum(results.TLH2_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMR')).*permute(repmat(sum(results.TSMR_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...sum(results.TMET_GA,1)').*shareOfGasEl);
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TSMC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TSMC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TSMC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TSMC')).*permute(repmat(sum(results.TSMC_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2]),2) +...sum(results.TMET_GA,1)').*shareOfGasEl);
    opex_capexRehab_PHS(:,yearNumb).*(PHS_LL)'+...
    sum(opex_capex(:,:,ismember(allActiveExLoadLabels,'RDSH')).*capacity(:,:,ismember(allActiveExLoadLabels,'RDSH'))+opex_var(:,:,ismember(allActiveExLoadLabels,'RDSH')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'RDSH')).*permute(repmat(sum(results.RDSH_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDNG')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDNG'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDNG')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDNG')).*permute(repmat(sum(results.TDNG_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDOI')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDOI'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDOI')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDOI')).*permute(repmat(sum(results.TDOI_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDCO')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDCO'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDCO')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDCO')).*permute(repmat(sum(results.TDCO_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDBP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDBP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDBP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDBP')).*permute(repmat(sum(results.TDBP_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDHR')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDHR'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDHR')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDHR')).*permute(repmat(sum(results.TDHR_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDHP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDHP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDHP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDHP')).*permute(repmat(sum(results.TDHP_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'TDGE')).*capacity(:,:,ismember(allActiveExLoadLabels,'TDGE'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TDGE')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TDGE')).*permute(repmat(sum(results.TDGE_HE,1)',1,1,length(systemParams.IndexYears)),[1 3 2]) ,2) +...
    sum(opex_capex(:,:,ismember(allActiveExLoadLabels,'RRSH')).*(resultsSC.Cap_RRSH) +opex_var(:,:,ismember(allActiveExLoadLabels,'RRSH')).*repmat(sum(resultsSC.RRSH_HE)',1,size(systemParams.IndexYears,2)).* genShares_RRSH +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THHR')).*(resultsSC.Cap_THHR) +opex_var(:,:,ismember(allActiveExLoadLabels,'THHR')).*repmat(sum(resultsSC.THHR_HE)',1,size(systemParams.IndexYears,2)).* genShares_THHR +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THHP')).*(resultsSC.Cap_THHP) +opex_var(:,:,ismember(allActiveExLoadLabels,'THHP')).*repmat(sum(resultsSC.THHP_HE)',1,size(systemParams.IndexYears,2)).* genShares_THHP +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THNG')).*(resultsSC.Cap_THNG) +opex_var(:,:,ismember(allActiveExLoadLabels,'THNG')).*repmat(sum(resultsSC.THNG_HE)',1,size(systemParams.IndexYears,2)).* genShares_THNG +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THOI')).*(resultsSC.Cap_THOI) +opex_var(:,:,ismember(allActiveExLoadLabels,'THOI')).*repmat(sum(resultsSC.THOI_HE)',1,size(systemParams.IndexYears,2)).* genShares_THOI +...%    opex_capex(:,:,ismember(allActiveExLoadLabels,'THCO')).*resultsSC.Cap_THCO +opex_var(:,:,ismember(allActiveExLoadLabels,'THCO')).*sum(resultsSC.THCO_HE).*(resultsSC.Cap_THCO/resultsSC.OPT_SIZE_THCO) +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THBP')).*(resultsSC.Cap_THBP) +opex_var(:,:,ismember(allActiveExLoadLabels,'THBP')).*repmat(sum(resultsSC.THBP_HE)',1,size(systemParams.IndexYears,2)).* genShares_THBP +...
    opex_capex(:,:,ismember(allActiveExLoadLabels,'THBG')).*(resultsSC.Cap_THBG) +opex_var(:,:,ismember(allActiveExLoadLabels,'THBG')).*repmat(sum(resultsSC.THBG_HE)',1,size(systemParams.IndexYears,2)).* genShares_THBG,2)+...
    sum((opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSC')).*permute(repmat(sum(results.TPSC_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TPSP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TPSP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TPSP')).*permute(repmat(sum(results.TPSP_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRSL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRSL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRSL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRSL')).*permute(repmat(sum(results.TRSL_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRSS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRSS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRSS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRSS')).*permute(repmat(sum(results.TRSS_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRMI')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRMI'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRMI')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRMI')).*permute(repmat(sum(results.TRMI_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRMO')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRMO'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRMO')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRMO')).*permute(repmat(sum(results.TRMO_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRME')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRME'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRME')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRME')).*permute(repmat(sum(results.TRME_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRSi')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRSi'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRSi')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRSi')).*permute(repmat(sum(results.TRSi_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRAP')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRAP'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRAP')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRAP')).*permute(repmat(sum(results.TRAP_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRAD')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRAD'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRAD')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRAD')).*permute(repmat(sum(results.TRAD_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TREW')).*capacity(:,:,ismember(allActiveExLoadLabels,'TREW'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TREW')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TREW')).*permute(repmat(sum(results.TREW_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRBC')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRBC'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRBC')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRBC')).*permute(repmat(sum(results.TRBC_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])) +...
    (opex_capex(:,:,ismember(allActiveExLoadLabels,'TRGE')).*capacity(:,:,ismember(allActiveExLoadLabels,'TRGE'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TRGE')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TRGE')).*permute(repmat(sum(results.TRGE_GA/10^3,1)',1,1,length(systemParams.IndexYears)),[1 3 2])),2) +...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWOO')).*sum(results.RWOO_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWWO')).*sum(results.RWWO_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBGA')).*sum(results.RBGA_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBMW')).*sum(results.RBMW_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.RNGA_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RHAR')).*sum(results.RHAR_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RPET')).*sum(results.RPET_FU,1)' + ... %opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBDS')).*sum(results.RBDS_FU,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RURA')).*sum(results.RURA_FU,1)' + ...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RHAR_CO/10^3,1)' + ...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RPET_CO/10^3,1)' + ...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RNGA_CO/10^3,1)')+...
    sum(results.importFT,1)'.*setup.importFTCost(Regs,yearNumb)/1000+...
    sum(results.importH2,1)'.*setup.importH2Cost(Regs,yearNumb)/1000+...
    sum(results.importLNG,1)'.*setup.importLNGCost(Regs,yearNumb)/1000+...
    sum(results.importEL,1)'.*setup.importELCost(Regs,yearNumb)/1000+...
    sum(results.importMeOH,1)'.*setup.importMeOHCost(Regs,yearNumb)/1000+...
    sum(results.importNH3,1)'.*setup.importNH3Cost(Regs,yearNumb)/1000;

totAnCost_PEDAddition = opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.RNGA_OTHER,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RHAR')).*sum(results.RHAR_OTHER,1)' + ...
    opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RPET')).*sum(results.RPET_OTHER,1)';

sum(primaryCost_EL+storageCost+primaryCost_HEtot+gasCost+hyCost+FTCost+LH2Cost+LNGCost+opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RPET')).*sum(results.RPET_DI+results.RPET_KE,1)')
sum(primaryCost_EL+storageCost+primaryCost_HEtot+gasCost+hyCost+FTCost+LH2Cost+LNGCost+PetrFuelCost)
demand_all = demand+(sum(resultsSC.RES.SC_demand+resultsSC.COM.SC_demand+resultsSC.IND.SC_demand,1)'+sum(resultsSC.EL_EXCESS,1)')+demand_des+gasProdElCons+heatProdElCons+demand_transp+demand_industry+demand_export;

desLabels = systemParams.IndexID;
demand_water = sum(shiftdim(shiftdim(systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LPOT')),:,1,:),3),2),1)';
emand_water(demand_water<0) = 0;
pumpCapDemand = max(results.WROD_Desalination+results.WMSS_Desalination+results.WMSC_Desalination);


elprofit_WMSC = (sum(results.WMSC_EL,1)'- sum(results.EL_WMSC,1)').*((results.OPT_SIZE_WMSC'.*opex_capex(:,find(ismember(allActiveExLoadLabels,'WMSC'))))+ sum(results.GA_WMSC,1)'.*opex_var(:,find(ismember(allActiveExLoadLabels,'RNGA'))))./(sum(results.WMSC_EL,1)');
elprofit_WMSC(isnan(elprofit_WMSC)) = 0 ;
elprofit_WMDC = (sum(results.WMDC_EL,1)'- sum(results.EL_WMDC,1)').*((results.OPT_SIZE_WMDC'.*opex_capex(:,find(ismember(allActiveExLoadLabels,'WMDC'))))+ sum(results.GA_WMDC,1)'.*opex_var(:,find(ismember(allActiveExLoadLabels,'RNGA'))))./sum(results.WMDC_EL,1)';
elprofit_WMDC(isnan(elprofit_WMDC)) = 0 ;

desalinationCost = sum(capacity(:,:,ismember(allActiveExLoadLabels,'WROD')).*opex_capex(:,:,find(ismember(allActiveExLoadLabels,'WROD')))+...
    capacity(:,:,ismember(allActiveExLoadLabels,'WMSS')).*opex_capex(:,:,find(ismember(allActiveExLoadLabels,'WMSS')))+...
    capacity(:,:,ismember(allActiveExLoadLabels,'WMSC')).*opex_capex(:,:,find(ismember(allActiveExLoadLabels,'WMSC')))+...
    capacity(:,:,ismember(allActiveExLoadLabels,'WMDS')).*opex_capex(:,:,find(ismember(allActiveExLoadLabels,'WMDS')))+...
    capacity(:,:,ismember(allActiveExLoadLabels,'WMDC')).*opex_capex(:,:,find(ismember(allActiveExLoadLabels,'WMDC')))+...
    capacity(:,:,ismember(allActiveExLoadLabels,'SWAT')).*opex_capex(:,:,find(ismember(allActiveExLoadLabels,'SWAT'))),2)+...
    pumpCapDemand'.*systemParams.DesalinationParams.regionPumpLong(:,yearNumb).*opex_capex(:,yearNumb,find(ismember(allActiveExLoadLabels,'WHPU')))+...
    pumpCapDemand'.*systemParams.DesalinationParams.regionPumpUp(:,yearNumb).*opex_capex(:,yearNumb,find(ismember(allActiveExLoadLabels,'WVPU')))+...
    sum(results.GA_WDES,1)'.*opex_var(:,yearNumb,find(ismember(allActiveExLoadLabels,'RNGA')))+...
    sum(results.EL_WROD,1)'.*(full_Cons_LCOEsys)+...
    sum(results.EL_WMSS,1)'.*(full_Cons_LCOEsys)+...
    sum(results.EL_WMDS,1)'.*(full_Cons_LCOEsys)+...
    sum(results.EL_WMSC,1)'.*(full_Cons_LCOEsys)+...
    sum(results.EL_WMDC,1)'.*(full_Cons_LCOEsys)-...
    sum(results.WMSC_EL,1)'.*(full_Cons_LCOEsys)-...
    sum(results.WMDC_EL,1)'.*(full_Cons_LCOEsys);

desCAPEX = sum(capacity(:,:,ismember(allActiveExLoadLabels,'WROD')).*capex(:,:,find(ismember(allActiveExLoadLabels,'WROD')))+...
    capacity(:,:,ismember(allActiveExLoadLabels,'WMSS')).*capex(:,:,find(ismember(allActiveExLoadLabels,'WMSS')))+...
    capacity(:,:,ismember(allActiveExLoadLabels,'WMSC')).*capex(:,:,find(ismember(allActiveExLoadLabels,'WMSC')))+...
    capacity(:,:,ismember(allActiveExLoadLabels,'WMDS')).*capex(:,:,find(ismember(allActiveExLoadLabels,'WMDS')))+...
    capacity(:,:,ismember(allActiveExLoadLabels,'WMDC')).*capex(:,:,find(ismember(allActiveExLoadLabels,'WMDC')))+...
    capacity(:,:,ismember(allActiveExLoadLabels,'SWAT')).*capex(:,:,find(ismember(allActiveExLoadLabels,'SWAT'))),2)+...
    pumpCapDemand'.*systemParams.DesalinationParams.regionPumpLong(:,yearNumb).*capex(:,yearNumb,find(ismember(allActiveExLoadLabels,'WHPU')))+...
    pumpCapDemand'.*systemParams.DesalinationParams.regionPumpUp(:,yearNumb).*capex(:,yearNumb,find(ismember(allActiveExLoadLabels,'WVPU')));

desOPEX = sum(capacity(:,:,ismember(allActiveExLoadLabels,'WROD')).*opex_fix(:,:,find(ismember(allActiveExLoadLabels,'WROD')))+...
    capacity(:,:,ismember(allActiveExLoadLabels,'WMSS')).*opex_fix(:,:,find(ismember(allActiveExLoadLabels,'WMSS')))+...
    capacity(:,:,ismember(allActiveExLoadLabels,'WMSC')).*opex_fix(:,:,find(ismember(allActiveExLoadLabels,'WMSC')))+...
    capacity(:,:,ismember(allActiveExLoadLabels,'WMDS')).*opex_fix(:,:,find(ismember(allActiveExLoadLabels,'WMDS')))+...
    capacity(:,:,ismember(allActiveExLoadLabels,'WMDC')).*opex_fix(:,:,find(ismember(allActiveExLoadLabels,'WMDC')))+...
    capacity(:,:,ismember(allActiveExLoadLabels,'SWAT')).*opex_fix(:,:,find(ismember(allActiveExLoadLabels,'SWAT'))),2)+...
    pumpCapDemand'.*systemParams.DesalinationParams.regionPumpLong(:,yearNumb).*opex_fix(:,yearNumb,find(ismember(allActiveExLoadLabels,'WHPU')))+...
    pumpCapDemand'.*systemParams.DesalinationParams.regionPumpUp(:,yearNumb).*opex_fix(:,yearNumb,find(ismember(allActiveExLoadLabels,'WVPU')))+...
    sum(results.GA_WDES,1)'.*opex_var(:,yearNumb,find(ismember(allActiveExLoadLabels,'RNGA')))+...
    sum(results.EL_WROD,1)'.*(full_Cons_LCOEsys)+...
    sum(results.EL_WMSS,1)'.*(full_Cons_LCOEsys)+...
    sum(results.EL_WMDS,1)'.*(full_Cons_LCOEsys)+...
    sum(results.EL_WMSC,1)'.*(full_Cons_LCOEsys)+...
    sum(results.EL_WMDC,1)'.*(full_Cons_LCOEsys)-...
    sum(results.WMSC_EL,1)'.*(full_Cons_LCOEsys)-...
    sum(results.WMDC_EL,1)'.*(full_Cons_LCOEsys);


desalinationElCost = demand_des.*(primaryLCOEsys+LCOCsys+LCOSsys+lcot_sys);
desalinationElCost = demand_des.*(full_Cons_LCOEsys);
desalinationGasCost = sum(results.GA_WMSC,1)'.*opex_var(:,find(ismember(allActiveExLoadLabels,'RNGA')));
WaterCost = desalinationCost./(demand_water+sum(results.Wa_to_TRAD+results.Wa_to_TRBC+results.Wa_to_TRMI+results.Wa_to_TRMO)');
WaterCost(isnan(WaterCost))=0;
WaterCost(isinf(WaterCost))=0;

% !! to have appropriate water cost, excluding el production cost

demand_gas = shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LIGA'),:,1,:)),3);%.*GasFlag;

gasCAPEX = sum((capex(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGU'))).*repmat((1-shareOfGasEl),1,length(systemParams.IndexYears)) +...sum(results.RBME,1)').*(1-shareOfGasEl) +...
    (capex(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGD'))).*repmat(shareOfgasBGA,1,length(systemParams.IndexYears)) +...sum(results.RBGA_TBGD,1)').*(shareOfgasBGA).*(1-shareOfGasEl) +...
    (capex(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))).*repmat((1-shareOfGasEl),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (capex(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOS'))).*repmat((1-shareOfGasEl),1,length(systemParams.IndexYears)) +....*sum(results.TCOS_GA,1)').*(1-shareOfGasEl) +...
    (capex(:,:,ismember(allActiveExLoadLabels,'TMET')).*capacity(:,:,ismember(allActiveExLoadLabels,'TMET'))).*repmat((1-shareOfGasEl),1,length(systemParams.IndexYears)),2);


gasOPEX = sum((opex_fix(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGU'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGU')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGU')).*permute(repmat(sum(results.TBGU_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((1-shareOfGasEl),1,length(systemParams.IndexYears)) +...sum(results.RBME,1)').*(1-shareOfGasEl) +...
    (opex_fix(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capacity(:,:,ismember(allActiveExLoadLabels,'TBGD'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TBGD')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TBGD')).*permute(repmat(sum(results.RBGA_FU,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat(shareOfgasBGA,1,length(systemParams.IndexYears)) +...sum(results.RBGA_TBGD,1)').*(shareOfgasBGA).*(1-shareOfGasEl) +...
    (opex_fix(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capacity(:,:,ismember(allActiveExLoadLabels,'TWEL'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TWEL')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TWEL')).*permute(repmat(sum(results.TWEL_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((1-shareOfGasEl),1,length(systemParams.IndexYears)) +...sum(results.TWEL_GA,1)').*(1-shareOfGasEl) +...
    (opex_fix(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capacity(:,:,ismember(allActiveExLoadLabels,'TCOS'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TCOS')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TCOS')).*permute(repmat(sum(results.TCOS_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((1-shareOfGasEl),1,length(systemParams.IndexYears)) +....*sum(results.TCOS_GA,1)').*(1-shareOfGasEl) +...
    (opex_fix(:,:,ismember(allActiveExLoadLabels,'TMET')).*capacity(:,:,ismember(allActiveExLoadLabels,'TMET'))+opex_var(:,:,ismember(allActiveExLoadLabels,'TMET')).*capSharesByYears(:,:,ismember(allActiveExLoadLabels,'TMET')).*permute(repmat(sum(results.TMET_GA,1)',1,1,length(systemParams.IndexYears)),[1 3 2])).*repmat((1-shareOfGasEl),1,length(systemParams.IndexYears)),2) +...sum(results.TMET_GA,1)').*(1-shareOfGasEl))+...
    (opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBGA')).*sum(results.TBGU_GA,1)').*(shareOfgasBGA).*(1-shareOfGasEl) +...
    (opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA')).*sum(results.LIGA_GAS_FOS,1)') +...
    systemParams.fossilCO2Cost(yearNumb).*sum(results.RNGA_CO,1)'/1000.*shareOfNGAInd;


gasElCost = gasProdElCons.*(full_Cons_LCOEsys);

GasProductionCost = gasCost + ...
    +gasProdElCons.*(full_Cons_LCOEsys);

GasCost = GasProductionCost./demand_gas;%.*GasFlag;
GasCost(isnan(GasCost))=0;
GasCost(isinf(GasCost))=0;

demMat = shiftdim((systemParams.ValueLoad(1,:,1,:)),3);
demMat = shiftdim(demMat,2);

gridDirect = (results.GRID_AC+results.GRID_DC).*((results.GRID_AC+results.GRID_DC)>0).*demMat./(demMat+results.EL_SBAT+results.EL_SACA+results.EL_SPHS+results.EL_TCOS+results.EL_TWEL+results.EL_TDHP+results.EL_TDHR+results.EL_TDGE);

storIndirOut = sum(results.SBAT_EL+results.SACA_EL+resultsSC.SBAR_EL+resultsSC.SBAC_EL+resultsSC.SBAI_EL,1)'+sum(results.SPHS_EL,1)' +sum(results.TCCG_EL,1)'.*shareOfGT_SNG+sum(results.TOCG_EL,1)'.*shareOfGT_SNG+sum(results.TCNG_EL,1)'.*shareOfGT_SNG++sum(results.TSTU_EL,1)';
if setup.Heat.Flag
    storIndirOut = sum(results.SBAT_EL+results.SACA_EL+resultsSC.SBAR_EL+resultsSC.SBAC_EL+resultsSC.SBAI_EL,1)'+sum(results.SPHS_EL,1)' +sum(results.TCCG_EL,1)'.*shareOfGT_SNG+sum(results.TOCG_EL,1)'.*shareOfGT_SNG+sum(results.TCNG_EL,1)'.*shareOfGT_SNG;
end

storHeIndirOut = sum(results.SHOT_EL,1)'+sum(results.SDHS_EL,1)' +sum(results.TDNG_HE+results.TCNG_HE,1)'.*shareOfGT_SNG+sum(resultsSC.THNG_HE,1)'.*shareOfGT_SNG;



biomassCons = (sum(results.RWOO_FU)+sum(results.RWWO_FU)+...
    sum(results.RBMW_FU)+sum(results.RBGA_FU))';


fuel_cost = sum(results.RBGA_FU,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBGA'))+...
    sum(results.RWOO_FU,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWOO'))+...
    sum(results.RWWO_FU,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWWO'))+...
    sum(results.RBMW_FU,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBMW'))+...;
    sum(results.RNGA_FU,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA'))+...;
    sum(results.RHAR_FU,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RHAR'))+...;
    sum(results.RPET_FU,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RPET'))+...sum(results.RBDS_DI,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBDS'))+...
    sum(results.RURA_FU,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RURA'))+...
    sum(results.importFT,1)'.*setup.importFTCost(Regs,yearNumb)/1000+...
    sum(results.importH2,1)'.*setup.importH2Cost(Regs,yearNumb)/1000+...
    sum(results.importMeOH,1)'.*setup.importMeOHCost(Regs,yearNumb)/1000+...
    sum(results.importNH3,1)'.*setup.importNH3Cost(Regs,yearNumb)/1000+...
    sum(results.importLNG,1)'.*setup.importLNGCost(Regs,yearNumb)/1000;


foss_fuel_cost = sum(results.RNGA_FU,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA'))+...;
    sum(results.RHAR_FU,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RHAR'))+...;
    sum(results.RPET_FU,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RPET'));

nuc_fuel_cost = opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RURA')).*sum(results.RURA_FU,1)';


ren_fuel_cost = sum(results.RBGA_FU,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBGA'))+...
    sum(results.RWOO_FU,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWOO'))+...
    sum(results.RWWO_FU,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWWO'))+...
    sum(results.RBMW_FU,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RBMW'));

import_fuel_cost = sum(results.importFT,1)'.*setup.importFTCost(Regs,yearNumb)/1000+...
    sum(results.importH2,1)'.*setup.importH2Cost(Regs,yearNumb)/1000+...
    sum(results.importLNG,1)'.*setup.importLNGCost(Regs,yearNumb)/1000+...
    sum(results.importMeOH,1)'.*setup.importMeOHCost(Regs,yearNumb)/1000+...
    sum(results.importNH3,1)'.*setup.importNH3Cost(Regs,yearNumb)/1000;

%% El demand



demand_tot; % power direct
demand_des; % power to desalination
gasProdElCons; % power to industrial gas
heatProdElCons; % power to district heat (individual el heating is in demand_tot)
demand_transp;
el_demand_transp_dir = sum(results.EltoMobility,1)';
el_demand_transp_fuels = (sum(results.EL_TWEL,1)'.*share.H2toM)+(sum(results.EL_TWEL,1)'.*share.H2toLH2 + sum(results.EL_TLH2,1)')+(sum(results.EL_TWEL,1)'.*share.H2toLNG + sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*share.CO2toLNG + sum(results.EL_TFTU,1)'+ sum(results.EL_TLNG,1)')+(sum(results.EL_TWEL,1)'.*share.H2toFTU + sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*share.CO2toFTU);

el_toFuels_H2_p = sum(results.EL_TWEL,1)';
el_toFuels_H2_M = (sum(results.EL_TWEL,1)'.*share.H2toM);
el_toFuels_LH2 = (sum(results.EL_TWEL,1)'.*share.H2toLH2 + sum(results.EL_TLH2,1)');
el_toFuels_LNG = (sum(results.EL_TWEL,1)'.*share.H2toLNG + sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*share.CO2toLNG + sum(results.EL_TLNG,1)');
el_toFuels_FTU = (sum(results.EL_TWEL,1)'.*share.H2toFTU + sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*share.CO2toFTU+ sum(results.EL_TFTU,1)');

all_el_demand = demand_tot+demand_des+gasProdElCons+heatProdElCons+el_demand_transp_dir+el_demand_transp_fuels;

%% He demand
all_he_demand = localheatDemand+distrheatDemand + sum(results.HE_TCOS,1)'.*share.CO2toCH4.*shareOfGasInd + sum(results.HE_TCOS,1)'.*share.CO2toLNG + sum(results.HE_TCOS,1)'.*share.CO2toFTU;
%% transportation

%% TTW

baseAllGHGEmissionsGas = (sum(results.RNGA_FU,1)')*(systemParams.GasEmissionsMain)/10^3;
baseAllGHGEmissionsOil = (sum(results.RPET_FU,1)')*(systemParams.OilEmissionsMain)/10^3;
baseAllGHGEmissionsCoa = (sum(results.RHAR_FU,1)')*(systemParams.CoalEmissionsMain)/10^3;

baseAllGHGEmissions_CapturedCCS = ((sum(results.TGCS_GAS_FOS,1)'.*systemParams.EffGCS_CO(:,systemParams.IndexYears==costYear)+sum(results.TSMC_GAS_FOS,1)'.*systemParams.EffSMC_CO(1:length(systemParams.IndexNumNodes),systemParams.IndexYears==costYear))*(systemParams.GasEmissionsMain)/10^3+(sum(results.RPET_TGCS,1)'.*systemParams.EffGCS_CO(:,systemParams.IndexYears==costYear))*(systemParams.OilEmissionsMain)/10^3+(sum(results.RHAR_THCS,1)'.*systemParams.EffHCS_CO(:,systemParams.IndexYears==costYear))*(systemParams.CoalEmissionsMain)/10^3);

baseElGHGEmissionsGas = (sum(results.TCCG_GAS_FOS+results.TOCG_GAS_FOS+results.TCNG_GAS_FOS+results.TICM_GAS_FOS+results.TGCS_GAS_FOS,1)')*(systemParams.GasEmissionsMain)/10^3;
baseElGHGEmissionsOil = (sum(results.RPET_TICG+results.RPET_TICM+results.RPET_TCOI+results.RPET_TCCG+results.RPET_TOCG+results.RPET_TCNG+results.RPET_TGCS,1)')*(systemParams.OilEmissionsMain)/10^3;
baseElGHGEmissionsCoa = (sum(results.RHAR_THPP+results.RHAR_THCS+results.RHAR_TCCO,1)')*(systemParams.CoalEmissionsMain)/10^3;
baseElGHGEmissions = baseElGHGEmissionsGas+baseElGHGEmissionsOil+baseElGHGEmissionsCoa;
baseHeGHGEmissionsGas = (sum(results.TDNG_GAS_FOS+results.Pros_GAS_FOS,1)')*(systemParams.GasEmissionsMain)/10^3;
baseHeGHGEmissionsOil = (sum(results.RPET_TDOI+resultsSC.RPET_FU+results.RPET_LHIN,1)')*(systemParams.OilEmissionsMain)/10^3;
baseHeGHGEmissionsCoa = (sum(results.RHAR_TDCO,1)')*(systemParams.CoalEmissionsMain)/10^3;
baseHeGHGEmissions = baseHeGHGEmissionsGas+baseHeGHGEmissionsOil+baseHeGHGEmissionsCoa;
baseHyGHGEmissions = (sum(results.TSMR_GAS_FOS,1)')*(systemParams.GasEmissionsMain)/10^3;
baseGaGHGEmissions = (sum(results.TLNG_GAS_FOS,1)')*(systemParams.GasEmissionsMain)/10^3;
baseTrGHGEmissions = (sum(results.RPET_KE+results.RPET_DI,1)')*(systemParams.OilEmissionsMain)/10^3;
baseIGGHGEmissions = (sum(results.LIGA_GAS_FOS,1)')*(systemParams.GasEmissionsMain)/10^3;
baseDeGHGEmissions = (sum(results.WDES_GAS_FOS,1)')*(systemParams.GasEmissionsMain)/10^3;
baseINDGHGEmissions = (sum(results.HC_TISB+results.RHAR_LHIN+results.LHIN_GAS_FOS+results.RPET_LHIN,1)')*(systemParams.OilEmissionsMain)/10^3;

GHGel_dir = baseElGHGEmissions./all_el_demand.*demand_tot;
GHGel_des = baseElGHGEmissions./all_el_demand.*demand_des;
GHGel_gas = baseElGHGEmissions./all_el_demand.*gasProdElCons;
GHGel_he = baseElGHGEmissions./all_el_demand.*heatProdElCons;
GHGel_trans_dir = baseElGHGEmissions./all_el_demand.*el_demand_transp_dir;
GHGel_trans_H2 = baseElGHGEmissions./all_el_demand.*el_toFuels_H2_M;
GHGel_trans_LH2 = baseElGHGEmissions./all_el_demand.*el_toFuels_LH2;
GHGel_trans_LNG = baseElGHGEmissions./all_el_demand.*el_toFuels_LNG;
GHGel_trans_FTU = baseElGHGEmissions./all_el_demand.*el_toFuels_FTU;

GHGhe_dir = baseHeGHGEmissions./all_he_demand.*(localheatDemand+distrheatDemand);GHGhe_dir(isnan(GHGhe_dir)) = 0;
GHGhe_gas = baseHeGHGEmissions./all_he_demand.*(sum(results.HE_TCOS,1)'.*share.CO2toCH4.*shareOfGasInd);GHGhe_gas(isnan(GHGhe_gas)) = 0;
GHGhe_trans_LNG = baseHeGHGEmissions./all_he_demand.*sum(results.HE_TCOS,1)'.*share.CO2toLNG;GHGhe_trans_LNG(isnan(GHGhe_trans_LNG)) = 0;
GHGhe_trans_FTU = baseHeGHGEmissions./all_he_demand.*sum(results.HE_TCOS,1)'.*share.CO2toFTU;GHGhe_trans_FTU(isnan(GHGhe_trans_FTU)) = 0;

GHGel_dir+GHGel_des+GHGel_gas+GHGel_he+GHGel_trans_dir+GHGel_trans_H2+GHGel_trans_LH2+GHGel_trans_LNG+GHGel_trans_FTU;
GHGhe_dir+GHGhe_gas+GHGhe_trans_LNG+GHGhe_trans_FTU;

GHG_trans_El_tot = GHGel_trans_dir; % emissions from cars on hydrogen
GHG_trans_El_rel = GHG_trans_El_tot./(sum(results.EltoMobility+0.01,1)');
GHG_trans_H2_tot = GHGel_trans_H2 + baseHyGHGEmissions.*share.H2toM; % emissions from cars on hydrogen
GHG_trans_H2_rel = GHG_trans_H2_tot./(sum(results.HY_TRSP+0.01,1)');
GHG_trans_LH2_tot= GHGel_trans_LH2+ baseHyGHGEmissions.*share.H2toLH2; % emissions from LH2 based transp
GHG_trans_LH2_rel= GHG_trans_LH2_tot./(sum(results.LH2toMobility+0.01,1)');
GHG_trans_LNG_tot= GHGel_trans_LNG+ baseGaGHGEmissions + baseHyGHGEmissions.*share.H2toLNG; % emissions from LNG based transp
GHG_trans_LNG_rel= GHG_trans_LNG_tot./(sum(results.LNGtoMobility+0.01,1)');
GHG_trans_OIL_tot= GHGel_trans_FTU+ baseTrGHGEmissions; % emissions from petrol based transp
GHG_trans_OIL_rel= GHG_trans_OIL_tot./(sum(results.diesToMobility+results.kersToMobility+0.01,1)');

%% WTW

baseElGHGEmissionsGasWTW = (sum(results.TCCG_GAS_FOS+results.TOCG_GAS_FOS+results.TCNG_GAS_FOS+results.TICM_GAS_FOS+results.TGCS_GAS_FOS,1)')*(systemParams.GasEmissionsMain+systemParams.GasEmissionsAdditional)/10^3;
baseElGHGEmissionsOilWTW = (sum(results.RPET_TICG+results.RPET_TICM+results.RPET_TCOI+results.RPET_TCCG+results.RPET_TOCG+results.RPET_TCNG+results.RPET_TGCS,1)')*(systemParams.OilEmissionsMain+systemParams.OilEmissionsAdditional)/10^3;
baseElGHGEmissionsCoaWTW = (sum(results.RHAR_THPP+results.RHAR_THCS+results.RHAR_TCCO,1)')*(systemParams.CoalEmissionsMain+systemParams.CoalEmissionsAdditional)/10^3;
baseElGHGEmissionsWTW = baseElGHGEmissionsGasWTW+baseElGHGEmissionsOilWTW+baseElGHGEmissionsCoaWTW;
baseHeGHGEmissionsGasWTW = (sum(results.TDNG_GAS_FOS+results.Pros_GAS_FOS+results.LHIN_GAS_FOS,1)')*(systemParams.GasEmissionsMain+systemParams.GasEmissionsAdditional)/10^3;
baseHeGHGEmissionsOilWTW = (sum(results.RPET_TDOI+resultsSC.RPET_FU+results.RPET_LHIN,1)')*(systemParams.OilEmissionsMain+systemParams.OilEmissionsAdditional)/10^3;
baseHeGHGEmissionsCoaWTW = (sum(results.RHAR_TDCO+results.RHAR_LHIN,1)')*(systemParams.CoalEmissionsMain+systemParams.CoalEmissionsAdditional)/10^3;
baseHeGHGEmissionsWTW = baseHeGHGEmissionsGasWTW+baseHeGHGEmissionsOilWTW+baseHeGHGEmissionsCoaWTW;
baseHyGHGEmissionsWTW = (sum(results.TSMR_GAS_FOS,1)')*(systemParams.GasEmissionsMain+systemParams.GasEmissionsAdditional)/10^3;
baseGaGHGEmissionsWTW = (sum(results.TLNG_GAS_FOS,1)')*(systemParams.GasEmissionsMain+systemParams.GasEmissionsAdditional)/10^3;
baseTrGHGEmissionsWTW = (sum(results.RPET_KE+results.RPET_DI,1)')*(systemParams.OilEmissionsMain+systemParams.OilEmissionsAdditional)/10^3;
baseIGGHGEmissionsWTW = (sum(results.LIGA_GAS_FOS,1)')*(systemParams.GasEmissionsMain+systemParams.GasEmissionsAdditional)/10^3;
baseDeGHGEmissionsWTW = (sum(results.WDES_GAS_FOS,1)')*(systemParams.GasEmissionsMain+systemParams.GasEmissionsAdditional)/10^3;

GHGel_dirWTW = baseElGHGEmissionsWTW./all_el_demand.*demand_tot;
GHGel_desWTW = baseElGHGEmissionsWTW./all_el_demand.*demand_des;
GHGel_gasWTW = baseElGHGEmissionsWTW./all_el_demand.*gasProdElCons;
GHGel_heWTW = baseElGHGEmissionsWTW./all_el_demand.*heatProdElCons;
GHGel_trans_dirWTW = baseElGHGEmissionsWTW./all_el_demand.*el_demand_transp_dir;
GHGel_trans_H2WTW = baseElGHGEmissionsWTW./all_el_demand.*el_toFuels_H2_M;
GHGel_trans_LH2WTW = baseElGHGEmissionsWTW./all_el_demand.*el_toFuels_LH2;
GHGel_trans_LNGWTW = baseElGHGEmissionsWTW./all_el_demand.*el_toFuels_LNG;
GHGel_trans_FTUWTW = baseElGHGEmissionsWTW./all_el_demand.*el_toFuels_FTU;

GHGhe_dirWTW = baseHeGHGEmissionsWTW./all_he_demand.*(localheatDemand+distrheatDemand);GHGhe_dirWTW(isnan(GHGhe_dirWTW)) = 0;
GHGhe_gasWTW = baseHeGHGEmissionsWTW./all_he_demand.*(sum(results.HE_TCOS,1)'.*share.CO2toCH4.*shareOfGasInd);GHGhe_gasWTW(isnan(GHGhe_gasWTW)) = 0;
GHGhe_trans_LNGWTW = baseHeGHGEmissionsWTW./all_he_demand.*sum(results.HE_TCOS,1)'.*share.CO2toLNG;GHGhe_trans_LNGWTW(isnan(GHGhe_trans_LNGWTW)) = 0;
GHGhe_trans_FTUWTW = baseHeGHGEmissionsWTW./all_he_demand.*sum(results.HE_TCOS,1)'.*share.CO2toFTU;GHGhe_trans_FTUWTW(isnan(GHGhe_trans_FTUWTW)) = 0;



GHG_trans_El_totWTW = GHGel_trans_dirWTW; % emissions from cars on hydrogen
GHG_trans_El_relWTW = GHG_trans_El_totWTW./(sum(results.EltoMobility+0.01,1)');
GHG_trans_H2_totWTW = GHGel_trans_H2WTW + baseHyGHGEmissionsWTW.*share.H2toM; % emissions from cars on hydrogen
GHG_trans_H2_relWTW = GHG_trans_H2_totWTW./(sum(results.HY_TRSP+0.01,1)');
GHG_trans_LH2_totWTW= GHGel_trans_LH2WTW+ baseHyGHGEmissionsWTW.*share.H2toLH2; % emissions from LH2 based transp
GHG_trans_LH2_relWTW= GHG_trans_LH2_totWTW./(sum(results.LH2toMobility+0.01,1)');
GHG_trans_LNG_totWTW= GHGel_trans_LNGWTW+ baseGaGHGEmissionsWTW + baseHyGHGEmissionsWTW.*share.H2toLNG; % emissions from LNG based transp
GHG_trans_LNG_relWTW= GHG_trans_LNG_totWTW./(sum(results.LNGtoMobility+0.01,1)');
GHG_trans_OIL_totWTW= GHGel_trans_FTUWTW+ baseTrGHGEmissionsWTW; % emissions from petrol based transp
GHG_trans_OIL_relWTW= GHG_trans_OIL_totWTW./(sum(results.diesToMobility+results.kersToMobility+0.01,1)');
%%
results.RNGA_CO_WTW = results.RNGA_CO/systemParams.GasEmissionsMain*(systemParams.GasEmissionsMain+systemParams.GasEmissionsAdditional);
results.RPET_CO_WTW = results.RPET_CO/systemParams.OilEmissionsMain*(systemParams.OilEmissionsMain+systemParams.OilEmissionsAdditional);
results.RHAR_CO_WTW = results.RHAR_CO/systemParams.CoalEmissionsMain*(systemParams.CoalEmissionsMain+systemParams.CoalEmissionsAdditional);

GasEmissions_WTW = (systemParams.GasEmissionsMain+systemParams.GasEmissionsAdditional)/10^3;
OilEmissions_WTW = (systemParams.OilEmissionsMain+systemParams.OilEmissionsAdditional)/10^3;
CoalEmissions_WTW= (systemParams.CoalEmissionsMain+systemParams.CoalEmissionsAdditional)/10^3;

GasEmissions_TTW = (systemParams.GasEmissionsMain)/10^3;
OilEmissions_TTW = (systemParams.OilEmissionsMain)/10^3;
CoalEmissions_TTW= (systemParams.CoalEmissionsMain)/10^3;


GHG_tot_WTW = sum(results.LHIN_GAS_FOS +results.LIGA_GAS_FOS+results.WDES_GAS_FOS+results.TCCG_GAS_FOS+results.TOCG_GAS_FOS+results.TCNG_GAS_FOS+results.TICM_GAS_FOS+results.TGCS_GAS_FOS+results.TDNG_GAS_FOS+results.Pros_GAS_FOS+results.TSMR_GAS_FOS+results.TLNG_GAS_FOS,1)'.*GasEmissions_WTW +...
    sum(results.RPET_LHIN+results.RPET_TICG+results.RPET_TICM+results.RPET_TCOI+results.RPET_TCCG+results.RPET_TOCG+results.RPET_TCNG+results.RPET_TGCS+results.RPET_TDOI+resultsSC.RPET_FU+results.RPET_KE+results.RPET_DI,1)'.*OilEmissions_WTW +...
    sum(results.RHAR_LHIN +results.RHAR_THPP+results.RHAR_TCCO+results.RHAR_TDCO,1)'.*CoalEmissions_WTW;
GHG_El_WTW = sum(results.TCCG_GAS_FOS+results.TOCG_GAS_FOS+results.TCNG_GAS_FOS+results.TICM_GAS_FOS+results.TGCS_GAS_FOS,1)'.*GasEmissions_WTW +...
    sum(results.RPET_TICG+results.RPET_TICM+results.RPET_TCOI+results.RPET_TCCG+results.RPET_TOCG+results.RPET_TCNG+results.RPET_TGCS,1)'.*OilEmissions_WTW +...
    sum(results.RHAR_THPP+results.RHAR_THCS+results.RHAR_TCCO,1)'.*CoalEmissions_WTW;
GHG_Ind_WTW = sum(results.LIGA_GAS_FOS,1)'.*GasEmissions_WTW;
GHG_De_WTW = sum(results.WDES_GAS_FOS,1)'.*GasEmissions_WTW;
GHG_He_WTW = sum(results.LHIN_GAS_FOS+results.TDNG_GAS_FOS+results.Pros_GAS_FOS,1)'.*GasEmissions_WTW +...
    sum(results.RPET_LHIN+results.RPET_TDOI+resultsSC.RPET_FU,1)'.*OilEmissions_WTW +...
    sum(results.RHAR_LHIN+results.RHAR_TDCO,1)'.*CoalEmissions_WTW;


GHG_tot_TTW = sum(results.LHIN_GAS_FOS +results.LIGA_GAS_FOS+results.WDES_GAS_FOS+results.TCCG_GAS_FOS+results.TOCG_GAS_FOS+results.TCNG_GAS_FOS+results.TICM_GAS_FOS+results.TGCS_GAS_FOS+results.TDNG_GAS_FOS+results.Pros_GAS_FOS+results.TSMR_GAS_FOS+results.TLNG_GAS_FOS,1)'.*GasEmissions_TTW +...
    sum(results.RPET_LHIN+results.RPET_TICG+results.RPET_TICM+results.RPET_TCOI+results.RPET_TCCG+results.RPET_TOCG+results.RPET_TCNG+results.RPET_TGCS+results.RPET_TDOI+resultsSC.RPET_FU+results.RPET_KE+results.RPET_DI,1)'.*OilEmissions_TTW +...
    sum(results.RHAR_LHIN +results.RHAR_THPP+results.RHAR_TCCO+results.RHAR_TDCO,1)'.*CoalEmissions_TTW;
GHG_El_TTW = sum(results.TCCG_GAS_FOS+results.TOCG_GAS_FOS+results.TCNG_GAS_FOS+results.TICM_GAS_FOS+results.TGCS_GAS_FOS,1)'.*GasEmissions_TTW +...
    sum(results.RPET_TICG+results.RPET_TICM+results.RPET_TCOI+results.RPET_TCCG+results.RPET_TOCG+results.RPET_TCNG+results.RPET_TGCS,1)'.*OilEmissions_TTW +...
    sum(results.RHAR_THPP+results.RHAR_THCS+results.RHAR_TCCO,1)'.*CoalEmissions_TTW;
GHG_Ind_TTW = sum(results.LIGA_GAS_FOS,1)'.*GasEmissions_TTW;
GHG_De_TTW = sum(results.WDES_GAS_FOS,1)'.*GasEmissions_TTW;
GHG_He_TTW = sum(results.LHIN_GAS_FOS+results.TDNG_GAS_FOS+results.Pros_GAS_FOS,1)'.*GasEmissions_TTW +...
    sum(results.RPET_LHIN+results.RPET_TDOI+resultsSC.RPET_FU,1)'.*OilEmissions_TTW +...
    sum(results.RHAR_LHIN+results.RHAR_TDCO,1)'.*CoalEmissions_TTW;

GHG_El_WTW_rel = GHG_El_WTW.*10^6./prod_electricity_tot; %[gCO2eq/kWh,el]
GHG_El_TTW_rel = GHG_El_TTW.*10^6./prod_electricity_tot; %[gCO2eq/kWh,el]

GHG_He_WTW_rel = GHG_He_WTW.*10^6./prod_heat_tot;%[gCO2eq/kWh,th]
GHG_He_TTW_rel = GHG_He_TTW.*10^6./prod_heat_tot;%[gCO2eq/kWh,th]

GHG_tot_WTW_cost = systemParams.fossilCO2Cost(yearNumb)*GHG_tot_WTW;
GHG_tot_TTW_cost = systemParams.fossilCO2Cost(yearNumb)*GHG_tot_TTW;

%% Chemical industry



AmmoniaFromRE = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LINH'));
AmmoniaFromFossil = setup.AmmoniaDemand(Regs,1:find(ismember(systemParams.IndexYears ,setup.endYear)))-AmmoniaFromRE;

FossilsToChemicals = systemParams.Instalations(:,yearNumb,ismember(systemParams.IndexID,'LICH'));

GasToAmmonia = min(FossilsToChemicals.*setup.GasToChem(Regs),AmmoniaFromFossil(:,yearNumb)*setup.GasToAmmonia_Factor);
OilToAmmonia = (AmmoniaFromFossil(:,yearNumb)-GasToAmmonia/setup.GasToAmmonia_Factor)*setup.OilToAmmonia_Factor;

OtherGasLoss = (FossilsToChemicals.*setup.GasToChem(Regs) - GasToAmmonia)*setup.Loss;
OtherOilLoss = (FossilsToChemicals.*setup.OilToChem(Regs) - OilToAmmonia)*setup.Loss;

CoalToChemAll = FossilsToChemicals.*setup.CoalToChem(Regs);

GasToAmmonia_CO2_TTW = GasToAmmonia * systemParams.GasEmissionsMain; %kgCO2
OilToAmmonia_CO2_TTW = OilToAmmonia * systemParams.OilEmissionsMain;

OtherGasLoss_CO2_TTW = OtherGasLoss * systemParams.GasEmissionsMain;
OtherOilLoss_CO2_TTW = OtherOilLoss * systemParams.OilEmissionsMain;

CoalToChem_CO2_TTW = CoalToChemAll * systemParams.CoalEmissionsMain;

TotalChem_CO2_TTW = GasToAmmonia_CO2_TTW+OilToAmmonia_CO2_TTW+OtherGasLoss_CO2_TTW+OtherOilLoss_CO2_TTW+CoalToChem_CO2_TTW;
TotalChem_CO2_TTW(TotalChem_CO2_TTW<0) = 0;

GasToAmmonia_CO2_WTW = GasToAmmonia * (systemParams.GasEmissionsMain+systemParams.GasEmissionsAdditional); %kgCO2
OilToAmmonia_CO2_WTW = OilToAmmonia * (systemParams.OilEmissionsMain+systemParams.OilEmissionsAdditional);

OtherGasLoss_CO2_WTW = OtherGasLoss * (systemParams.GasEmissionsMain+systemParams.GasEmissionsAdditional);
OtherOilLoss_CO2_WTW = OtherOilLoss * (systemParams.OilEmissionsMain+systemParams.OilEmissionsAdditional);

CoalToChem_CO2_WTW = CoalToChemAll * (systemParams.CoalEmissionsMain+systemParams.CoalEmissionsAdditional);

TotalChem_CO2_WTW = GasToAmmonia_CO2_WTW+OilToAmmonia_CO2_WTW+OtherGasLoss_CO2_WTW+OtherOilLoss_CO2_WTW+CoalToChem_CO2_WTW;
TotalChem_CO2_WTW(TotalChem_CO2_WTW<0) = 0;


GasToChemFeedstock = FossilsToChemicals.*setup.GasToChem(Regs);
OilToChemFeedstock = FossilsToChemicals.*setup.OilToChem(Regs);
CoalToChemFeedstock = FossilsToChemicals.*setup.CoalToChem(Regs);

TotalChemFeedstockCost = GasToChemFeedstock.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RNGA'))+...;
    CoalToChemFeedstock.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RHAR'))+...;
    OilToChemFeedstock.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RPET'));

%% for slide 9
ElDirect = sum(results.EltoMobility,1)';
ElIndirHydrogen = sum(results.EL_TWEL,1)'.*(share.H2toM+share.H2toLH2) + sum(results.EL_TLH2,1)';
ElIndirMethane = sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC)'.*share.CO2toCH4.*shareOfGasTrans+sum(results.EL_TWEL,1)'.*share.H2toCH4.*shareOfGasTrans + sum(results.EL_TLNG,1)';
ElIndirLiqFuel = sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC)'.*share.CO2toFTU+sum(results.EL_TWEL,1)'.*share.H2toFTU;

%% costs
LCOEfull = primaryLCOEsys+LCOCsys+LCOSsys+lcot_sys+LCOSBsys;
LCOEfull = primaryLCOEtot+LCOCtot+LCOStot+lcot_tot;
LCOEfull(LCOEfull==0) = primaryLCOEtot_average+LCOCtot_average+LCOStot_average+lcot_tot_average;
LCOEfull_av = mean(primaryLCOEsys);
LCOEfull_av = primaryLCOEsys_average+LCOCsys_average+LCOSsys_average+lcot_sys_average+LCOSBsys_average;
LCOEfull_av = primaryLCOEtot_average+LCOCtot_average+LCOStot_average+lcot_tot_average;

LCOHy = (hyCost + LCOEfull.*(sum(results.EL_TWEL,1)'))./sum(results.TWEL_GA+results.TSMR_GA+results.TSMC_GA+results.importH2,1)';
LCOHyav = sum((hyCost + LCOEfull.*(sum(results.EL_TWEL,1)')))./sum(sum(results.TWEL_GA+results.TSMR_GA+results.TSMC_GA+results.importH2,1)');

LCOHy_fos = (hyCost_fos)./(sum(results.TSMR_GA+results.TSMC_GA+0.01,1)');
LCOHy_fos(isnan(LCOHy_fos))=0;
LCOHyav_fos = sum(hyCost_fos)./sum(sum(results.TSMR_GA+results.TSMC_GA+0.01,1)');
LCOHyav_fos(isnan(LCOHyav_fos))=0;

fosH2_CO2cost = (sum(results.TSMR_GAS_FOS,1)'*systemParams.GasEmissionsMain/1000*systemParams.fossilCO2Cost(yearNumb))./(sum(results.TSMR_GA+results.TSMC_GA+0.01,1)'); %[EUR/kWh]
fosH2_CO2av = sum(sum(results.TSMR_GAS_FOS,1)'*systemParams.GasEmissionsMain/1000*systemParams.fossilCO2Cost(yearNumb))./sum(sum(results.TSMR_GA+results.TSMC_GA+0.01,1)); %[EUR/kWh]

LCOHy_ren = (hyCost_ren + LCOEfull.*(sum(results.EL_TWEL,1)'))./(sum(results.TWEL_GA+results.importH2+0.01,1)');
LCOHy_ren(isnan(LCOHy_ren))=0;
LCOHy_ren(sum(results.TWEL_GA+results.importH2+0.01,1)'<10^5)=0;
LCOHyav_ren = sum(hyCost_ren + LCOEfull.*(sum(results.EL_TWEL,1)'))./sum(sum(results.TWEL_GA+results.importH2+0.01,1)');
LCOHyav_ren(isnan(LCOHyav_ren))=0;
LCOHy_ren(sum(sum(results.TWEL_GA+results.importH2+0.01,1)')<10^5)=0;

if costYear == 2025
    qq=1;
end

LCOLH2 = (LH2Cost + LCOEfull.*(sum(results.EL_TWEL,1)'.*share.H2toLH2 + sum(results.EL_TLH2,1)'))./(sum(results.LH2toMobility+0.01,1)');
LCOLH2av = sum((LH2Cost + LCOEfull.*(sum(results.EL_TWEL,1)'.*share.H2toLH2 + sum(results.EL_TLH2,1)')))./sum(sum(results.LH2toMobility+0.01,1)');

LCOLNG = (LNGCost + LCOEfull.*(sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*share.CO2toLNG + sum(results.EL_TWEL,1)'.*share.H2toLNG + sum(results.EL_TLNG,1)'))./(sum(results.LNGtoMobility+0.01,1)');
LCOLNGav = sum(LCOLNG.*(sum(results.LNGtoMobility,1)'))./sum(sum(results.LNGtoMobility+0.01,1)');

LCOSNG = (SNGCost + LCOEfull.*(sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*share.CO2toCH4 + sum(results.EL_TWEL,1)'.*share.H2toCH4))./(sum(results.TMET_GA+0.01,1)');
LCOSNG(isnan(LCOSNG))=0;
LCOSNGav = sum(LCOSNG.*(sum(results.TMET_GA,1)'))./sum(sum(results.TMET_GA+0.01,1)');
LCOSNGav(isnan(LCOSNGav))=0;

LCOFT = (FTCost + LCOEfull.*(sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*share.CO2toFTU + sum(results.EL_TWEL,1)'.*share.H2toFTU + sum(results.EL_TFTU,1)') )./(sum(results.TFTU_DI+results.TFTU_KE+results.TFTU_NA+0.01,1)');
LCOFTav = sum(LCOFT.*(sum(results.TFTU_DI+results.TFTU_KE,1)'))./sum(sum(results.TFTU_DI+results.TFTU_KE+0.01,1)');

sum(sum(results.TFTU_DI+results.TFTU_KE)+sum(results.RBDS_DI+results.RBDS_KE)-sum(results.diesToMobility+results.kersToMobility));
sum(sum(results.TFTU_KE)+sum(results.RBDS_KE)-sum(results.kersToMobility));
sum(sum(results.TFTU_DI)+sum(results.RBDS_DI)-sum(results.diesToMobility));

FT_AllCost = (FTCost + LCOEfull.*(sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*share.CO2toFTU + sum(results.EL_TWEL,1)'.*share.H2toFTU) );
FT_ElCons = ((sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*share.CO2toFTU + sum(results.EL_TWEL,1)'.*share.H2toFTU) );
FT_excessTot = sum(sum(results.TFTU_DI+results.TFTU_KE)+sum(results.RBDS_DI+results.RBDS_KE)-sum(results.diesToMobility+results.kersToMobility)-sum(results.exportFT_diesel+results.exportFT_kerosene));
FT_excessTot(FT_excessTot<0)=0;

FT_excess = (sum(results.TFTU_DI+results.TFTU_KE)+sum(results.RBDS_DI+results.RBDS_KE)-sum(results.diesToMobility+results.kersToMobility)-sum(results.exportFT_diesel+results.exportFT_kerosene));
FT_excess = (sum(results.TFTU_DI+results.TFTU_KE)-sum(results.FTdiesel+results.FTkerosene)-sum(results.exportFT_diesel+results.exportFT_kerosene));


FT_excessShareAv = sum(sum(results.TFTU_DI+results.TFTU_KE)+sum(results.RBDS_DI+results.RBDS_KE)-sum(results.diesToMobility+results.kersToMobility)-sum(results.exportFT_diesel+results.exportFT_kerosene))./...
    sum(sum(results.TFTU_DI+results.TFTU_KE));
FT_excessShareAv = sum(sum(results.TFTU_DI+results.TFTU_KE)-sum(results.FTdiesel+results.FTkerosene)-sum(results.exportFT_diesel+results.exportFT_kerosene))./...
    sum(sum(results.TFTU_DI+results.TFTU_KE));
FT_excessShareAv(isnan(FT_excessShareAv))=0;
FT_excessShareAv((FT_excessShareAv<0))=0;

FT_excessShare = (sum(results.TFTU_DI+results.TFTU_KE)+sum(results.RBDS_DI+results.RBDS_KE)-sum(results.diesToMobility+results.kersToMobility)-sum(results.exportFT_diesel+results.exportFT_kerosene))./...
    (sum(results.TFTU_DI+results.TFTU_KE));
FT_excessShare = (sum(results.TFTU_DI+results.TFTU_KE)-sum(results.FTdiesel+results.FTkerosene)-sum(results.exportFT_diesel+results.exportFT_kerosene))./...
    (sum(results.TFTU_DI+results.TFTU_KE));
FT_excessShare(isnan(FT_excessShare))=0;
FT_excessShare((FT_excessShare<0))=0;
FT_excessShare = 0*FT_excessShare;

LCOLQ = (PetrFuelCost + FTCost + LCOEfull.*(sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*share.CO2toFTU + sum(results.EL_TWEL,1)'.*share.H2toFTU) )./(sum(results.diesToMobility+results.kersToMobility+results.TFTU_NA+0.01,1)');
LCOLQav = sum(LCOLQ.*(sum(results.diesToMobility+results.kersToMobility+results.TFTU_NA,1)'))./sum(sum(results.diesToMobility+results.kersToMobility+results.TFTU_NA+0.01,1)');


LCOEfull.* sum(results.EltoMobility)'+(hyCost.*share.H2toM + LCOEfull.*(sum(results.EL_TWEL,1)'.*share.H2toM))+ (LH2Cost + LCOEfull.*(sum(results.EL_TWEL,1)'.*share.H2toLH2 + sum(results.EL_TLH2,1)'))+(LNGCost + LCOEfull.*(sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*share.CO2toLNG + sum(results.EL_TWEL,1)'.*share.H2toLNG + sum(results.EL_TLNG,1)'))+...
    (PetrFuelCost + FTCost + LCOEfull.*(sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC,1)'.*share.CO2toFTU + sum(results.EL_TWEL,1)'.*share.H2toFTU+ sum(results.EL_TFTU,1)') );



LCO_CDR = (CDRCost + LCOEfull.*demand_CDR)./(sum(results.TRSL_GA+results.TRSS_GA+results.TRMI_GA+results.TRMO_GA+results.TRME_GA+results.TRSi_GA+results.TRAP_GA+results.TRAD_GA+results.TREW_GA+results.TRBC_GA+results.TRGE_GA+0.001,1)');
LCO_CDR(isnan(LCO_CDR))=0;
LCO_CDRav = sum(LCO_CDR.*(sum(results.TRSL_GA+results.TRSS_GA+results.TRMI_GA+results.TRMO_GA+results.TRME_GA+results.TRSi_GA+results.TRAP_GA+results.TRAD_GA+results.TREW_GA+results.TRBC_GA+results.TRGE_GA+0.001,1)'))./sum(sum(results.TRSL_GA+results.TRSS_GA+results.TRMI_GA+results.TRMO_GA+results.TRME_GA+results.TRSi_GA+results.TRAP_GA+results.TRAD_GA+results.TREW_GA+results.TRBC_GA+results.TRGE_GA+0.001,1)');
LCO_CDRav(isnan(LCO_CDR))=0;

transpDem = sum(results.FU_MRLI+results.FU_MRBI+results.FU_MRWI+results.FU_MRHI+results.FU_MRMI+...%ice cars: LDV, Busses, 2 and 3 W ( these 3 passengers), HDV and MDV (these 2 freight)
    results.FU_MRLP+results.FU_MRBP+results.FU_MRWP+results.FU_MRHP+results.FU_MRMP+...%PHEV cars fuel demand: LDV, Busses, 2 and 3 W ( these 3 passengers), HDV and MDV (these 2 freight)
    results.EL_MRLP+results.EL_MRBP+results.EL_MRWP+results.EL_MRHP+results.EL_MRMP+...%PHEV cars elec demand: LDV, Busses, 2 and 3 W ( these 3 passengers), HDV and MDV (these 2 freight)
    results.EL_MRLB+results.EL_MRBB+results.EL_MRWB+results.EL_MRHB+results.EL_MRMB+...%BEV cars: LDV, Busses, 2 and 3 W ( these 3 passengers), HDV and MDV (these 2 freight)
    results.HY_MRLF+results.HY_MRBF+results.HY_MRWF+results.HY_MRHF+results.HY_MRMF+...%FCEV cars: LDV, Busses, 2 and 3 W ( these 3 passengers), HDV and MDV (these 2 freight)
    results.MRPF_Dem + results.MRPE_Dem +... % Rail passenger, fuel based and el based
    results.MRFF_Dem + results.MRFE_Dem +... % Rail freight, fuel based and el based
    results.MAPF_Dem + results.MAPE_Dem + results.MAPH_Dem +... % Avia passenger, fuel based, elec based and hydrogen based
    results.MAFF_Dem + results.MAFE_Dem + results.MAFH_Dem +... % Avia freight, fuel based, elec based and hydrogen based
    results.MMPF_Dem + results.MMPE_Dem + results.MMPH_Dem + results.MMPG_Dem+ results.MMPA_Dem + results.MMPM_Dem+... % Marine passenger, fuel based, elec based, hydrogen, and LNG based
    results.MMFF_Dem + results.MMFE_Dem + results.MMFH_Dem + results.MMFG_Dem+ results.MMFA_Dem + results.MMFM_Dem); % Marine freight, fuel based, elec based, hydrogen, and LNG based


sum(results.EltoMobility+results.diesToMobility+results.kersToMobility+results.HY_TRSP+results.LNGtoMobility+results.LH2toMobility);


transpCost = sum(results.FU_MRLI+results.FU_MRBI+results.FU_MRWI+results.FU_MRHI+results.FU_MRMI)'.*LCOLQ+...%ice cars: LDV, Busses, 2 and 3 W ( these 3 passengers), HDV and MDV (these 2 freight)
    sum(results.FU_MRLP+results.FU_MRBP+results.FU_MRWP+results.FU_MRHP+results.FU_MRMP)'.*LCOLQ+...%PHEV cars fuel demand: LDV, Busses, 2 and 3 W ( these 3 passengers), HDV and MDV (these 2 freight)
    sum(results.EL_MRLP+results.EL_MRBP+results.EL_MRWP+results.EL_MRHP+results.EL_MRMP)'.*LCOEfull+...%PHEV cars elec demand: LDV, Busses, 2 and 3 W ( these 3 passengers), HDV and MDV (these 2 freight)
    sum(results.EL_MRLB+results.EL_MRBB+results.EL_MRWB+results.EL_MRHB+results.EL_MRMB)'.*LCOEfull+...%BEV cars: LDV, Busses, 2 and 3 W ( these 3 passengers), HDV and MDV (these 2 freight)
    sum(results.HY_MRLF+results.HY_MRBF+results.HY_MRWF+results.HY_MRHF+results.HY_MRMF)'.*LCOHy+...%FCEV cars: LDV, Busses, 2 and 3 W ( these 3 passengers), HDV and MDV (these 2 freight)
    sum(results.MRPF_Dem)'.*LCOLQ + sum(results.MRPE_Dem)'.*LCOEfull +... % Rail passenger, fuel based and el based
    sum(results.MRFF_Dem)'.*LCOLQ + sum(results.MRFE_Dem)'.*LCOEfull +... % Rail freight, fuel based and el based
    sum(results.MAPF_Dem)'.*LCOLQ + sum(results.MAPE_Dem)'.*LCOEfull + sum(results.MAPH_Dem)'.*LCOLH2 +... % Avia passenger, fuel based, elec based and hydrogen based
    sum(results.MAFF_Dem)'.*LCOLQ + sum(results.MAFE_Dem)'.*LCOEfull + sum(results.MAFH_Dem)'.*LCOLH2 +... % Avia freight, fuel based, elec based and hydrogen based
    sum(results.MMPF_Dem)'.*LCOLQ + sum(results.MMPE_Dem)'.*LCOEfull + sum(results.MMPH_Dem)'.*LCOLH2 + sum(results.MMPG_Dem)'.*LCOLNG+... % Marine passenger, fuel based, elec based, hydrogen, and LNG based
    sum(results.MMFF_Dem)'.*LCOLQ + sum(results.MMFE_Dem)'.*LCOEfull + sum(results.MMFH_Dem)'.*LCOLH2 + sum(results.MMFG_Dem)'.*LCOLNG;

transpCostTot = sum(sum(results.FU_MRLI+results.FU_MRBI+results.FU_MRWI+results.FU_MRHI+results.FU_MRMI)').*LCOLQav+...%ice cars: LDV, Busses, 2 and 3 W ( these 3 passengers), HDV and MDV (these 2 freight)
    sum(sum(results.FU_MRLP+results.FU_MRBP+results.FU_MRWP+results.FU_MRHP+results.FU_MRMP)').*LCOLQav+...%PHEV cars fuel demand: LDV, Busses, 2 and 3 W ( these 3 passengers), HDV and MDV (these 2 freight)
    sum(sum(results.EL_MRLP+results.EL_MRBP+results.EL_MRWP+results.EL_MRHP+results.EL_MRMP)').*LCOEfull_av+...%PHEV cars elec demand: LDV, Busses, 2 and 3 W ( these 3 passengers), HDV and MDV (these 2 freight)
    sum(sum(results.EL_MRLB+results.EL_MRBB+results.EL_MRWB+results.EL_MRHB+results.EL_MRMB)').*LCOEfull_av+...%BEV cars: LDV, Busses, 2 and 3 W ( these 3 passengers), HDV and MDV (these 2 freight)
    sum(sum(results.HY_MRLF+results.HY_MRBF+results.HY_MRWF+results.HY_MRHF+results.HY_MRMF)').*LCOHyav+...%FCEV cars: LDV, Busses, 2 and 3 W ( these 3 passengers), HDV and MDV (these 2 freight)
    sum(sum(results.MRPF_Dem)').*LCOLQav + sum(sum(results.MRPE_Dem)').*LCOEfull_av +... % Rail passenger, fuel based and el based
    sum(sum(results.MRFF_Dem)').*LCOLQav + sum(sum(results.MRFE_Dem)').*LCOEfull_av +... % Rail freight, fuel based and el based
    sum(sum(results.MAPF_Dem)').*LCOLQav + sum(sum(results.MAPE_Dem)').*LCOEfull_av + sum(sum(results.MAPH_Dem)').*LCOLH2av +... % Avia passenger, fuel based, elec based and hydrogen based
    sum(sum(results.MAFF_Dem)').*LCOLQav + sum(sum(results.MAFE_Dem)').*LCOEfull_av + sum(sum(results.MAFH_Dem)').*LCOLH2av +... % Avia freight, fuel based, elec based and hydrogen based
    sum(sum(results.MMPF_Dem)').*LCOLQav + sum(sum(results.MMPE_Dem)').*LCOEfull_av + sum(sum(results.MMPH_Dem)').*LCOLH2av + sum(sum(results.MMPG_Dem)').*LCOLNGav+... % Marine passenger, fuel based, elec based, hydrogen, and LNG based
    sum(sum(results.MMFF_Dem)').*LCOLQav + sum(sum(results.MMFE_Dem)').*LCOEfull_av + sum(sum(results.MMFH_Dem)').*LCOLH2av + sum(sum(results.MMFG_Dem)').*LCOLNGav;

try sum(results.importEL,1)';
catch
    results.importEL = zeros(size(results.RWIN_EL));
end

%% primary energy calculations
primary_energy_tot_RE = sum(resultsSC.RPVR_EL,1)'+ sum(resultsSC.RPVC_EL,1)'+ sum(resultsSC.RPVI_EL,1)'+ ...
    sum(results.RPVO_EL,1)'+sum(results.RPVA_EL,1)'+...
    sum(results.RPBO_EL,1)' + sum(results.RPBA_EL,1)' + sum(results.RPBV_EL,1)' + sum(results.RPVF_EL,1)' + ...
    sum(results.RWAV_EL,1)'+...
    sum(results.RWIN_EL,1)'+sum(results.RWIO_EL,1)'+sum(results.ROWI_EL,1)'+...
    sum(results.RRRI_EL,1)'+sum(results.HDAM_EL,1)'+sum(results.RGEO_HE,1)'+...
    sum(results.RCSP_HE,1)'+sum(results.RDSH_HE,1)'+sum(resultsSC.RRSH_HE,1)'+...
    sum(results.RWOO_FU,1)'+sum(results.RWWO_FU,1)'+sum(results.RBMW_FU,1)'+sum(results.RBGA_FU,1)'+sum(results.RBDS_FU,1)'+...
    sum(results.importFT,1)'+sum(results.importH2,1)'+sum(results.importLNG,1)'+sum(results.importNH3,1)'+sum(results.importMeOH,1)'+sum(results.importEL,1)';

primary_energy_tot_Fossil = sum(results.RPET_FU,1)'+sum(results.RHAR_FU,1)'+sum(results.RNGA_FU,1)';%+sum(results.RURA_FU,1)'

primary_energy_tot_Nuclear = sum(results.RURA_FU,1)';

results.Import_RWOO = (sum(results.RWOO_FU,1)'-systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWOO'),:)').*((sum(results.RWOO_FU,1)'-systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWOO'),:)')>0);
results.Import_RWWO = (sum(results.RWWO_FU,1)'-systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWWO'),:)').*((sum(results.RWWO_FU,1)'-systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RWWO'),:)')>0);
results.Import_RBMW = (sum(results.RBMW_FU,1)'-systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBMW'),:)');%.*((sum(results.RBMW_FU,1)'-systemParams.ValueResourceTotal(ismember(systemParams.IndexIDR,'RBMW'),:)')>0);


primary_energy_tot_Local = sum(resultsSC.RPVR_EL,1)'+ sum(resultsSC.RPVC_EL,1)'+ sum(resultsSC.RPVI_EL,1)'+...
    sum(results.RPVO_EL,1)'+sum(results.RPVA_EL,1)'+...
    sum(results.RPBO_EL,1)' + sum(results.RPBA_EL,1)' + sum(results.RPBV_EL,1)' + sum(results.RPVF_EL,1)' + ...
    sum(results.RWAV_EL,1)'+...
    sum(results.RWIN_EL,1)'+sum(results.RWIO_EL,1)'+sum(results.ROWI_EL,1)'+...
    sum(results.RRRI_EL,1)'+sum(results.HDAM_EL,1)'+sum(results.RGEO_HE,1)'+...
    sum(results.RCSP_HE,1)'+sum(results.RDSH_HE,1)'+sum(resultsSC.RRSH_HE,1)'+...
    sum(results.RWOO_FU,1)'+sum(results.RWWO_FU,1)'+sum(results.RBMW_FU,1)'+sum(results.RBGA_FU,1)'-...
    results.Import_RBMW-results.Import_RWOO-results.Import_RWWO;


primary_energy_tot_Import = sum(results.RPET_FU,1)'+sum(results.RHAR_FU,1)'+sum(results.RNGA_FU,1)'+sum(results.RURA_FU,1)'+...
    sum(results.importFT,1)'+sum(results.importH2,1)'+sum(results.importLNG,1)'+sum(results.importEL,1)'+...
    results.Import_RBMW+results.Import_RWOO+results.Import_RWWO;

primary_energy_tot_ImportRE = sum(results.importFT,1)'+sum(results.importH2,1)'+sum(results.importLNG,1)'+sum(results.importEL,1)'+...
    results.Import_RBMW+results.Import_RWOO+results.Import_RWWO;

primary_energy_tot_ImportFoss = sum(results.RPET_FU,1)'+sum(results.RHAR_FU,1)'+sum(results.RNGA_FU,1)'+sum(results.RURA_FU,1)';

primary_energy_tot_ImportFoss_EL = sum(results.RPET_TCCG,1)'+sum(results.RPET_TOCG,1)'+sum(results.RPET_TGCS,1)'+sum(results.RPET_TCNG,1)'+sum(results.RPET_TICG,1)'+sum(results.RPET_TICM,1)'+sum(results.RPET_TCOI,1)'+...
    sum(results.RHAR_THPP,1)'+sum(results.RHAR_THCS,1)'+sum(results.RHAR_TCCO,1)'+...
    sum(results.TCCG_GAS_FOS,1)'+sum(results.TOCG_GAS_FOS,1)'+sum(results.TGCS_GAS_FOS,1)'+sum(results.TCNG_GAS_FOS,1)'+sum(results.TICM_GAS_FOS,1)';

primary_energy_tot_ImportNuc_EL = sum(results.RURA_FU,1)';

primary_energy_tot_ImportFoss_HE = sum(results.RPET_TDOI,1)'+sum(results.RPET_LHIN,1)'+...
    sum(results.RHAR_TDCO,1)'+sum(results.RHAR_LHIN,1)'+...
    sum(results.TDNG_GAS_FOS,1)'+sum(results.LHIN_GAS_FOS,1)'++sum(results.Pros_GAS_FOS,1)';

primary_energy_tot_ImportFoss_TR = sum(results.RPET_KE,1)'+sum(results.RPET_DI,1)'+...
    sum(results.TLNG_GAS_FOS,1)'+sum(results.TSMR_GAS_FOS,1)';

%Final energy demand

FinalEnergyDem_Trans = transpDem';

FinalEnergyDem_Des = (sum(results.EL_WROD,1)+sum(results.EL_WMSS,1)+sum(results.EL_WMSC,1)+sum(results.EL_WMDS,1)+sum(results.EL_WMDC,1)+sum(results.HE_WMSS,1)+sum(results.WDES_GAS_FOS,1)+sum(results.HE_WMDS,1)+...
    sum(results.WDES_GAS_REN,1)+sum(results.EL_WHPU+results.EL_WVPU,1))';

FinalEnergyDem_HeatDistr = (shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHSP'),:,1,:)),3)+...
    shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHDW'),:,1,:)),3)).*systemParams.Heat.shareOfDistrHeat(:,find(systemParams.IndexYears == costYear));
FinalEnergyDem_HeatIndiv = (shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHSP'),:,1,:)),3)+...
    shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHDW'),:,1,:)),3)).*(1-systemParams.Heat.shareOfDistrHeat(:,find(systemParams.IndexYears == costYear)))+...
    shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHBC'),:,1,:)),3);


FinalEnergyDem_Exports = sum(results.exportFT + results.exportH2 + results.exportLNG + results.exportMeOH + results.exportNH3,1)';
FinalEnergyDem_El = demand_origFinal+...
    sum(results.EL_TWEL,1)'.*share.H2toMeO+sum(results.EL_TWEL,1)'.*share.H2toNH3 + sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC+results.HE_TCOS+results.HE_TPSP+results.HE_TPSC,1)'.*share.CO2toMeO+sum(results.EL_TMeO+results.EL_TNH3,1)'+...
    sum(results.EL_TIAM+results.EL_TIAR+results.EL_TICB+results.EL_TICI+results.EL_TIPP+results.EL_TISB+results.EL_TISE+results.EL_TISH+results.EL_TISR,1)';


FinalEnergyDem_Industry = shiftdim(sum(systemParams.ValueLoad(ismember(systemParams.IndexIDL,'LHIN'),:,1,:)),3)+...
    sum(results.EL_TIAM+results.EL_TIAR+results.EL_TICB+results.EL_TICI+results.EL_TIPP+results.EL_TISB+results.EL_TISE+results.EL_TISH+results.EL_TISR,1)'+...
    sum(results.HE_TICB+results.HE_TICI+results.HE_TISB+results.HE_TISH+results.HE_TISR+results.HE_TISE+results.HE_TIAA+results.HE_TIAR+results.HE_TIPP,1)'+...
    sum(results.HY_TISH+results.HC_TISB,1)'+...
    sum(results.ME_LIME+results.NH_LINH,1)'+...
    systemParams.Instalations(:,systemParams.IndexYears==costYear,ismember(systemParams.IndexID,'LICH'))+...
    sum(results.CC_TISE+results.CC_TISH+results.CC_TISR,1)';

FinalEnergyDem_CDR = sum(results.EL_TRSS+results.EL_TRSL+results.EL_TRMI+results.EL_TRMO+results.EL_TRME+results.EL_TRSi+results.EL_TRAD+results.EL_TREW+results.EL_TRBC+results.EL_TRGE,1)'+...
    sum(results.EL_TCOS+results.EL_TPSP+results.EL_TPSC+results.HE_TCOS+results.HE_TPSP+results.HE_TPSC,1)'.*(share.CO2toCDR)+sum(results.HY_TRSi,1)'+...
    sum(results.HE_TRME+results.HE_TRSi+results.HE_TRGE,1)';

FinalEnergyDem_El = demand_origFinal;

FeedstockToSteel = sum(results.CC_TISE+results.CC_TISH+results.CC_TISR,1)';
SteelCharcFeedstockCost = sum(results.CC_TISE+results.CC_TISH+results.CC_TISR,1)'.*opex_var(:,yearNumb,ismember(allActiveExLoadLabels,'RWOO'));

if strcmp(type, 'Hist')
    mkdir([rootDir filesep 'projects' filesep 'resultsCalculated'])
    save([rootDir filesep 'projects' filesep 'resultsCalculated' filesep 'resultsCalculated_' pName '_' name],'-v7.3');
end



calc.primaryLCOEtot = primaryLCOEtot;
calc.primaryLCOEtot_average = primaryLCOEtot_average;

calc.LCOCtot = LCOCtot;
calc.LCOCtot_average = LCOCtot_average;

calc.lcot_tot = lcot_tot;
calc.lcot_tot_average = lcot_tot_average;

calc.LCOStot = LCOStot;
calc.LCOStot_average = LCOStot_average;

calc.WaterCost = WaterCost;
calc.GasCost = GasCost;


%LCO energy
calc.results.RES.Annual = results.RES.Annual;
calc.results.COM.Annual = results.COM.Annual;
calc.results.IND.Annual = results.IND.Annual;
calc.totAnCost = totAnCost;
calc.transmAnCost = transmAnCost;
calc.Share = Share;
calc.subs_cost = subs_cost;


