function [resultsSC,instalationsSC] = readProsResTrans(setup,systemParams,Reg1,Reg2,costYear,varargin)
%
% FUNCTION readProsResTrans(setup, systemParams, Reg1, Reg2, costYear, varargin)
%
% Reads and processes results for a transition scenario with prosumers.
%
%
% INPUT:
%            setup:           Structure that contains all necessary settings and data for processing.
%            systemParams:    Structure with system parameters related to the scenario.
%            Reg1:            Name or index of the first region.
%            Reg2:            Name or index of the second region.
%            costYear:        Year to which all costs are adjusted.
%            varargin:        ADD!!!!!!!!!!!!!!!!!!
%
% OUTPUT:
%            resultsSC:       Structure containing the processed scenario results for prosumers.
%            instalationsSC:  Installation data for prosumer technologies.
%
%Dmitrii Bogdanov
%last change 24.07.2025


if nargin == 6
    name = varargin{1};
    if length(name)>0
        nameAdd = ['_' varargin{1}]
    else
        nameAdd = ['']
    end
else
    name = '';
    nameAdd = ['']
end

[instalationsSC,activeSC] = ReadProsInstalations(setup,systemParams,Reg1,Reg2,name);



%try
if ~setup.OvernightFlag
    if setup.Heat.Flag
        projName = 'SC_HE';
    else
        projName = 'SC';
    end

    resCOM = load([setup.rootDir filesep 'projects' filesep projName filesep 'output' filesep 'results_' projName '_COM_' num2str(costYear) nameAdd]);
    resRES = load([setup.rootDir filesep 'projects' filesep projName filesep 'output' filesep 'results_' projName '_RES_' num2str(costYear) nameAdd]);
    resIND = load([setup.rootDir filesep 'projects' filesep projName filesep 'output' filesep 'results_' projName '_IND_' num2str(costYear) nameAdd]);

else
    if setup.Heat.Flag
        projName = ['SC_HE_Overnight_' num2str(setup.startYear)];
    else
        projName = ['SC_Overnight_' num2str(setup.startYear)];
    end

    resCOM = load([setup.rootDir filesep 'projects' filesep projName filesep 'output' filesep 'results_' projName '_COM_' num2str(costYear)]);
    resRES = load([setup.rootDir filesep 'projects' filesep projName filesep 'output' filesep 'results_' projName '_RES_' num2str(costYear)]);
    resIND = load([setup.rootDir filesep 'projects' filesep projName filesep 'output' filesep 'results_' projName '_IND_' num2str(costYear)]);

end

%% Capacities
for ii = 1:length(systemParams.IndexID)

    if activeSC(ii)==1 & ~sum(ismember({'RPVO','SBAT','IBAT'},systemParams.IndexID{ii}))

        eval(['resultsSC.OPT_SIZE_' systemParams.IndexID{ii} ' = resRES.results.OPT_SIZE_' systemParams.IndexID{ii} '(Reg1:Reg2) + resCOM.results.OPT_SIZE_' systemParams.IndexID{ii} '(Reg1:Reg2) + resIND.results.OPT_SIZE_' systemParams.IndexID{ii} '(Reg1:Reg2);'])

    end
end

resultsSC.OPT_SIZE_RPVR = resRES.results.OPT_SIZE_RPVO(Reg1:Reg2);
resultsSC.OPT_SIZE_RPVC = resCOM.results.OPT_SIZE_RPVO(Reg1:Reg2);
resultsSC.OPT_SIZE_RPVI = resIND.results.OPT_SIZE_RPVO(Reg1:Reg2);

resultsSC.OPT_SIZE_SBAR = resRES.results.OPT_SIZE_SBAT(Reg1:Reg2);
resultsSC.OPT_SIZE_SBAC = resCOM.results.OPT_SIZE_SBAT(Reg1:Reg2);
resultsSC.OPT_SIZE_SBAI = resIND.results.OPT_SIZE_SBAT(Reg1:Reg2);

resultsSC.OPT_SIZE_IBAR = resRES.results.OPT_SIZE_IBAT(Reg1:Reg2);
resultsSC.OPT_SIZE_IBAC = resCOM.results.OPT_SIZE_IBAT(Reg1:Reg2);
resultsSC.OPT_SIZE_IBAI = resIND.results.OPT_SIZE_IBAT(Reg1:Reg2);

%% Historical capacities
if costYear==2015
    instalMask = (repmat(systemParams.IndexYears,size(systemParams.IndexID,1),1)<=costYear);
else
    instalMask = ((repmat(systemParams.IndexYears,size(systemParams.IndexID,1),1)<=costYear).*(repmat(systemParams.IndexYears,size(systemParams.IndexID,1),1)>(costYear-systemParams.Lifetime)));

end

actCount = 0;
for i=1:length(activeSC)
    if activeSC(i) | sum([ismember({'RPVR' 'RPVC' 'RPVI' 'SBAR' 'SBAC' 'SBAI' 'IBAR' 'IBAC' 'IBAI'},systemParams.IndexID(i))])
        actCount = actCount + 1;
        fString = ['OPT_SIZE_' systemParams.IndexID{i}];
        if isfield(resultsSC,fString)
            if setup.OvernightFlag
                for reg = 1:size(systemParams.IndexNodes,1)
                    eval(['resultsSC.Cap_' systemParams.IndexID{i} '(reg,systemParams.IndexYears==costYear) = resultsSC.OPT_SIZE_' systemParams.IndexID{i} '(reg);']) % (regions,years)
                end
            else
                for reg = 1:size(systemParams.IndexNodes,1)
                    eval(['resultsSC.Cap_' systemParams.IndexID{i} '(reg,:) = instalationsSC(reg,:,i).*instalMask(i,:);']) % (regions,years)
                end
            end

        end
    end
end

if ~setup.Heat.Flag
    resultsSC.Cap_THBG = 0*resultsSC.Cap_THBG;
    resultsSC.Cap_THBP = 0*resultsSC.Cap_THBP;
    resultsSC.Cap_THHP = 0*resultsSC.Cap_THHP;
    resultsSC.Cap_THHR = 0*resultsSC.Cap_THHR;
    resultsSC.Cap_THNG = 0*resultsSC.Cap_THNG;
    resultsSC.Cap_THOI = 0*resultsSC.Cap_THOI;
    resultsSC.Cap_RRSH = 0*resultsSC.Cap_RRSH;
end

%% El/Heat Generation/Consuption
for ii = 1:length(systemParams.IndexID)

    if activeSC(ii)==1 & sum(ismember(systemParams.IndexIDT,systemParams.IndexID{ii}))
        try
            eval(['resultsSC.' systemParams.IndexID{ii} '_HE = resRES.results.' systemParams.IndexID{ii} '_HE(:,Reg1:Reg2) + resCOM.results.' systemParams.IndexID{ii} '_HE(:,Reg1:Reg2) + resIND.results.' systemParams.IndexID{ii} '_HE(:,Reg1:Reg2);']);
        end
        try
            eval(['resultsSC.EL_' systemParams.IndexID{ii} ' = resRES.results.EL_' systemParams.IndexID{ii} '(:,Reg1:Reg2) + resCOM.results.EL_' systemParams.IndexID{ii} '(:,Reg1:Reg2) + resIND.results.EL_' systemParams.IndexID{ii} '(:,Reg1:Reg2);']);
        end
    end
end

resultsSC.RRSH_HE = resRES.results.RRSH_HE(:,Reg1:Reg2)+resCOM.results.RRSH_HE(:,Reg1:Reg2)+resIND.results.RRSH_HE(:,Reg1:Reg2);
resultsSC.SHHS_EL = resRES.results.SHHS_EL(:,Reg1:Reg2)+resCOM.results.SHHS_EL(:,Reg1:Reg2)+resIND.results.SHHS_EL(:,Reg1:Reg2);
resultsSC.EL_SHHS = resRES.results.EL_SHHS(:,Reg1:Reg2)+resCOM.results.EL_SHHS(:,Reg1:Reg2)+resIND.results.EL_SHHS(:,Reg1:Reg2);
%% Fuel Consuption


resultsSC.RWOO_FU = resRES.results.RWOO_FU(:,Reg1:Reg2) + resCOM.results.RWOO_FU(:,Reg1:Reg2) + resIND.results.RWOO_FU(:,Reg1:Reg2);
resultsSC.RWWO_FU = resRES.results.RWWO_FU(:,Reg1:Reg2) + resCOM.results.RWWO_FU(:,Reg1:Reg2) + resIND.results.RWWO_FU(:,Reg1:Reg2);
resultsSC.RHAR_FU = resRES.results.RHAR_FU(:,Reg1:Reg2) + resCOM.results.RHAR_FU(:,Reg1:Reg2) + resIND.results.RHAR_FU(:,Reg1:Reg2);
resultsSC.RPET_FU = resRES.results.RPET_FU(:,Reg1:Reg2) + resCOM.results.RPET_FU(:,Reg1:Reg2) + resIND.results.RPET_FU(:,Reg1:Reg2);
resultsSC.RBGA_FU = resRES.results.RBGA_FU(:,Reg1:Reg2) + resCOM.results.RBGA_FU(:,Reg1:Reg2) + resIND.results.RBGA_FU(:,Reg1:Reg2);
resultsSC.RNGA_FU = resRES.results.RNGA_FU(:,Reg1:Reg2) + resCOM.results.RNGA_FU(:,Reg1:Reg2) + resIND.results.RNGA_FU(:,Reg1:Reg2);

%% PV Bat Generation/Consumption

resultsSC.RPVR_EL = resRES.results.RPVO_EL(:,Reg1:Reg2);
resultsSC.RPVC_EL = resCOM.results.RPVO_EL(:,Reg1:Reg2);
resultsSC.RPVI_EL = resIND.results.RPVO_EL(:,Reg1:Reg2);

resultsSC.SBAR_EL = resRES.results.SBAT_EL(:,Reg1:Reg2);
resultsSC.SBAC_EL = resCOM.results.SBAT_EL(:,Reg1:Reg2);
resultsSC.SBAI_EL = resIND.results.SBAT_EL(:,Reg1:Reg2);

resultsSC.EL_SBAR = resRES.results.EL_SBAT(:,Reg1:Reg2);
resultsSC.EL_SBAC = resCOM.results.EL_SBAT(:,Reg1:Reg2);
resultsSC.EL_SBAI = resIND.results.EL_SBAT(:,Reg1:Reg2);

resultsSC.SoC_SBAR = resRES.results.SoC_SBAT(:,Reg1:Reg2);
resultsSC.SoC_SBAC = resCOM.results.SoC_SBAT(:,Reg1:Reg2);
resultsSC.SoC_SBAI = resIND.results.SoC_SBAT(:,Reg1:Reg2);



%% Other

resultsSC.EL_GRID = resRES.results.EL_GRID(:,Reg1:Reg2) + resCOM.results.EL_GRID(:,Reg1:Reg2) + resIND.results.EL_GRID(:,Reg1:Reg2);
resultsSC.EL_EXCESS = resRES.results.EL_EXCESS(:,Reg1:Reg2) + resCOM.results.EL_EXCESS(:,Reg1:Reg2) + resIND.results.EL_EXCESS(:,Reg1:Reg2);
resultsSC.HE_EXCESS_Local = resRES.results.HE_EXCESS_Local(:,Reg1:Reg2) + resCOM.results.HE_EXCESS_Local(:,Reg1:Reg2) + resIND.results.HE_EXCESS_Local(:,Reg1:Reg2);

resultsSC.RES.SC_demand = resRES.results.SC_demand(:,Reg1:Reg2);
resultsSC.COM.SC_demand = resCOM.results.SC_demand(:,Reg1:Reg2);
resultsSC.IND.SC_demand = resIND.results.SC_demand(:,Reg1:Reg2);

%% prosumers sector calculations
resultsSC.RES.OPT_SIZE_RPVO = resRES.results.OPT_SIZE_RPVO(Reg1:Reg2);
resultsSC.RES.OPT_SIZE_SBAT = resRES.results.OPT_SIZE_SBAT(Reg1:Reg2);
resultsSC.RES.OPT_SIZE_IBAT = resRES.results.OPT_SIZE_IBAT(Reg1:Reg2);
resultsSC.RES.RPVO_EL = resRES.results.RPVO_EL(:,Reg1:Reg2);
resultsSC.RES.EL_SBAT = resRES.results.EL_SBAT(:,Reg1:Reg2);
resultsSC.RES.SBAT_EL = resRES.results.SBAT_EL(:,Reg1:Reg2);
resultsSC.RES.EL_EXCESS = resRES.results.EL_EXCESS(:,Reg1:Reg2);
resultsSC.RES.EL_GRID = resRES.results.EL_GRID(:,Reg1:Reg2);
resultsSC.RES.Dem = sum((resRES.results.RPVO_EL(:,Reg1:Reg2) - resRES.results.EL_SBAT(:,Reg1:Reg2) + resRES.results.SBAT_EL(:,Reg1:Reg2) - resRES.results.EL_EXCESS(:,Reg1:Reg2)),1);
resultsSC.RES.DemP = resRES.results.RPVO_EL(:,Reg1:Reg2) - resRES.results.EL_SBAT(:,Reg1:Reg2) + resRES.results.SBAT_EL(:,Reg1:Reg2) - resRES.results.EL_EXCESS(:,Reg1:Reg2);


resultsSC.COM.OPT_SIZE_RPVO = resCOM.results.OPT_SIZE_RPVO(Reg1:Reg2);
resultsSC.COM.OPT_SIZE_SBAT = resCOM.results.OPT_SIZE_SBAT(Reg1:Reg2);
resultsSC.COM.OPT_SIZE_IBAT = resCOM.results.OPT_SIZE_IBAT(Reg1:Reg2);
resultsSC.COM.RPVO_EL = resCOM.results.RPVO_EL(:,Reg1:Reg2);
resultsSC.COM.EL_SBAT = resCOM.results.EL_SBAT(:,Reg1:Reg2);
resultsSC.COM.SBAT_EL = resCOM.results.SBAT_EL(:,Reg1:Reg2);
resultsSC.COM.EL_EXCESS = resCOM.results.EL_EXCESS(:,Reg1:Reg2);
resultsSC.COM.EL_GRID = resCOM.results.EL_GRID(:,Reg1:Reg2);
resultsSC.COM.Dem = sum((resCOM.results.RPVO_EL(:,Reg1:Reg2) - resCOM.results.EL_SBAT(:,Reg1:Reg2) + resCOM.results.SBAT_EL(:,Reg1:Reg2) - resCOM.results.EL_EXCESS(:,Reg1:Reg2)),1);
resultsSC.COM.DemP = resCOM.results.RPVO_EL(:,Reg1:Reg2) - resCOM.results.EL_SBAT(:,Reg1:Reg2) + resCOM.results.SBAT_EL(:,Reg1:Reg2) - resCOM.results.EL_EXCESS(:,Reg1:Reg2);


resultsSC.IND.OPT_SIZE_RPVO = resIND.results.OPT_SIZE_RPVO(Reg1:Reg2);
resultsSC.IND.OPT_SIZE_SBAT = resIND.results.OPT_SIZE_SBAT(Reg1:Reg2);
resultsSC.IND.OPT_SIZE_IBAT = resIND.results.OPT_SIZE_IBAT(Reg1:Reg2);
resultsSC.IND.RPVO_EL = resIND.results.RPVO_EL(:,Reg1:Reg2);
resultsSC.IND.EL_SBAT = resIND.results.EL_SBAT(:,Reg1:Reg2);
resultsSC.IND.SBAT_EL = resIND.results.SBAT_EL(:,Reg1:Reg2);
resultsSC.IND.EL_EXCESS = resIND.results.EL_EXCESS(:,Reg1:Reg2);
resultsSC.IND.EL_GRID = resIND.results.EL_GRID(:,Reg1:Reg2);
resultsSC.IND.Dem = sum((resIND.results.RPVO_EL(:,Reg1:Reg2) - resIND.results.EL_SBAT(:,Reg1:Reg2) + resIND.results.SBAT_EL(:,Reg1:Reg2) - resIND.results.EL_EXCESS(:,Reg1:Reg2)),1);
resultsSC.IND.DemP = resIND.results.RPVO_EL(:,Reg1:Reg2) - resIND.results.EL_SBAT(:,Reg1:Reg2) + resIND.results.SBAT_EL(:,Reg1:Reg2) - resIND.results.EL_EXCESS(:,Reg1:Reg2);

