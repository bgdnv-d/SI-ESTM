function [resultsSC,instalationsSC] = readProsResTrans(setup,systemParams,Reg1,Regs,costYear,varargin)
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
%            Regs:            Regions for data collection.
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

[instalationsSC,activeSC] = ReadProsInstalations(setup,systemParams,Reg1,Regs,name);
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

        eval(['resultsSC.OPT_SIZE_' systemParams.IndexID{ii} ' = resRES.results.OPT_SIZE_' systemParams.IndexID{ii} '(Regs) + resCOM.results.OPT_SIZE_' systemParams.IndexID{ii} '(Regs) + resIND.results.OPT_SIZE_' systemParams.IndexID{ii} '(Regs);'])

    end
end

resultsSC.OPT_SIZE_RPVR = resRES.results.OPT_SIZE_RPVO(Regs);
resultsSC.OPT_SIZE_RPVC = resCOM.results.OPT_SIZE_RPVO(Regs);
resultsSC.OPT_SIZE_RPVI = resIND.results.OPT_SIZE_RPVO(Regs);

resultsSC.OPT_SIZE_SBAR = resRES.results.OPT_SIZE_SBAT(Regs);
resultsSC.OPT_SIZE_SBAC = resCOM.results.OPT_SIZE_SBAT(Regs);
resultsSC.OPT_SIZE_SBAI = resIND.results.OPT_SIZE_SBAT(Regs);

resultsSC.OPT_SIZE_IBAR = resRES.results.OPT_SIZE_IBAT(Regs);
resultsSC.OPT_SIZE_IBAC = resCOM.results.OPT_SIZE_IBAT(Regs);
resultsSC.OPT_SIZE_IBAI = resIND.results.OPT_SIZE_IBAT(Regs);

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
            %setfield() = getfield(resultsSC,fString)';
            if setup.OvernightFlag
                for reg = 1:size(systemParams.IndexNodes,1)
                    capacity(reg,:,actCount) = zeros(1,length(systemParams.IndexYears)); % (regions,years,tech)
                    capacity(reg,systemParams.IndexYears==costYear,actCount) = capacityCumul(reg,actCount);
                end
            else
                for reg = 1:size(systemParams.IndexNodes,1)
                    eval(['resultsSC.Cap_' systemParams.IndexID{i} '(reg,:) = instalationsSC(reg,:,i).*instalMask(i,:);']) % (regions,years)
                    %capacity(reg,systemParams.IndexYears==costYear,actCount) = capacityCumul(reg,actCount) - sum(capacity(reg,systemParams.IndexYears<costYear,actCount));
                end
            end
        else
            %capacity(:,:,actCount) = zeros(size(systemParams.IndexNodes,1),size(systemParams.IndexYears,2),1);
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
            eval(['resultsSC.' systemParams.IndexID{ii} '_HE = resRES.results.' systemParams.IndexID{ii} '_HE(:,Regs) + resCOM.results.' systemParams.IndexID{ii} '_HE(:,Regs) + resIND.results.' systemParams.IndexID{ii} '_HE(:,Regs);']);
        end
        try
            eval(['resultsSC.EL_' systemParams.IndexID{ii} ' = resRES.results.EL_' systemParams.IndexID{ii} '(:,Regs) + resCOM.results.EL_' systemParams.IndexID{ii} '(:,Regs) + resIND.results.EL_' systemParams.IndexID{ii} '(:,Regs);']);
        end
    end
end

resultsSC.RRSH_HE = resRES.results.RRSH_HE(:,Regs)+resCOM.results.RRSH_HE(:,Regs)+resIND.results.RRSH_HE(:,Regs);
resultsSC.SHHS_EL = resRES.results.SHHS_EL(:,Regs)+resCOM.results.SHHS_EL(:,Regs)+resIND.results.SHHS_EL(:,Regs);
resultsSC.EL_SHHS = resRES.results.EL_SHHS(:,Regs)+resCOM.results.EL_SHHS(:,Regs)+resIND.results.EL_SHHS(:,Regs);
%% Fuel Consuption


resultsSC.RWOO_FU = resRES.results.RWOO_FU(:,Regs) + resCOM.results.RWOO_FU(:,Regs) + resIND.results.RWOO_FU(:,Regs);
resultsSC.RWWO_FU = resRES.results.RWWO_FU(:,Regs) + resCOM.results.RWWO_FU(:,Regs) + resIND.results.RWWO_FU(:,Regs);
resultsSC.RHAR_FU = resRES.results.RHAR_FU(:,Regs) + resCOM.results.RHAR_FU(:,Regs) + resIND.results.RHAR_FU(:,Regs);
resultsSC.RPET_FU = resRES.results.RPET_FU(:,Regs) + resCOM.results.RPET_FU(:,Regs) + resIND.results.RPET_FU(:,Regs);
resultsSC.RBGA_FU = resRES.results.RBGA_FU(:,Regs) + resCOM.results.RBGA_FU(:,Regs) + resIND.results.RBGA_FU(:,Regs);
resultsSC.RNGA_FU = resRES.results.RNGA_FU(:,Regs) + resCOM.results.RNGA_FU(:,Regs) + resIND.results.RNGA_FU(:,Regs);

%% PV Bat Generation/Consumption

resultsSC.RPVR_EL = resRES.results.RPVO_EL(:,Regs);
resultsSC.RPVC_EL = resCOM.results.RPVO_EL(:,Regs);
resultsSC.RPVI_EL = resIND.results.RPVO_EL(:,Regs);

resultsSC.SBAR_EL = resRES.results.SBAT_EL(:,Regs);
resultsSC.SBAC_EL = resCOM.results.SBAT_EL(:,Regs);
resultsSC.SBAI_EL = resIND.results.SBAT_EL(:,Regs);

resultsSC.EL_SBAR = resRES.results.EL_SBAT(:,Regs);
resultsSC.EL_SBAC = resCOM.results.EL_SBAT(:,Regs);
resultsSC.EL_SBAI = resIND.results.EL_SBAT(:,Regs);

resultsSC.SoC_SBAR = resRES.results.SoC_SBAT(:,Regs);
resultsSC.SoC_SBAC = resCOM.results.SoC_SBAT(:,Regs);
resultsSC.SoC_SBAI = resIND.results.SoC_SBAT(:,Regs);



%% Other

resultsSC.EL_GRID = resRES.results.EL_GRID(:,Regs) + resCOM.results.EL_GRID(:,Regs) + resIND.results.EL_GRID(:,Regs);
resultsSC.EL_EXCESS = resRES.results.EL_EXCESS(:,Regs) + resCOM.results.EL_EXCESS(:,Regs) + resIND.results.EL_EXCESS(:,Regs);
resultsSC.HE_EXCESS_Local = resRES.results.HE_EXCESS_Local(:,Regs) + resCOM.results.HE_EXCESS_Local(:,Regs) + resIND.results.HE_EXCESS_Local(:,Regs);
% 
% resultsSC.RPVR_EL = resRES.results.RPVO_EL(:,Regs);
% resultsSC.RPVC_EL = resCOM.results.EL_THHP(:,Regs);
% resultsSC.RPVI_EL = resIND.results.RPVO_EL(:,Regs);

resultsSC.RES.SC_demand = resRES.results.SC_demand(:,Regs);
resultsSC.COM.SC_demand = resCOM.results.SC_demand(:,Regs);
resultsSC.IND.SC_demand = resIND.results.SC_demand(:,Regs);
    
%% prosumers sector calculations
resultsSC.RES.OPT_SIZE_RPVO = resRES.results.OPT_SIZE_RPVO(Regs);
resultsSC.RES.OPT_SIZE_SBAT = resRES.results.OPT_SIZE_SBAT(Regs);
resultsSC.RES.OPT_SIZE_IBAT = resRES.results.OPT_SIZE_IBAT(Regs);
resultsSC.RES.RPVO_EL = resRES.results.RPVO_EL(:,Regs);
resultsSC.RES.EL_SBAT = resRES.results.EL_SBAT(:,Regs);
resultsSC.RES.SBAT_EL = resRES.results.SBAT_EL(:,Regs);
resultsSC.RES.EL_EXCESS = resRES.results.EL_EXCESS(:,Regs);
resultsSC.RES.EL_GRID = resRES.results.EL_GRID(:,Regs);
resultsSC.RES.Dem = sum((resRES.results.RPVO_EL(:,Regs) - resRES.results.EL_SBAT(:,Regs) + resRES.results.SBAT_EL(:,Regs) - resRES.results.EL_EXCESS(:,Regs)),1);
resultsSC.RES.DemP = resRES.results.RPVO_EL(:,Regs) - resRES.results.EL_SBAT(:,Regs) + resRES.results.SBAT_EL(:,Regs) - resRES.results.EL_EXCESS(:,Regs);


resultsSC.COM.OPT_SIZE_RPVO = resCOM.results.OPT_SIZE_RPVO(Regs);
resultsSC.COM.OPT_SIZE_SBAT = resCOM.results.OPT_SIZE_SBAT(Regs);
resultsSC.COM.OPT_SIZE_IBAT = resCOM.results.OPT_SIZE_IBAT(Regs);
resultsSC.COM.RPVO_EL = resCOM.results.RPVO_EL(:,Regs);
resultsSC.COM.EL_SBAT = resCOM.results.EL_SBAT(:,Regs);
resultsSC.COM.SBAT_EL = resCOM.results.SBAT_EL(:,Regs);
resultsSC.COM.EL_EXCESS = resCOM.results.EL_EXCESS(:,Regs);
resultsSC.COM.EL_GRID = resCOM.results.EL_GRID(:,Regs);
resultsSC.COM.Dem = sum((resCOM.results.RPVO_EL(:,Regs) - resCOM.results.EL_SBAT(:,Regs) + resCOM.results.SBAT_EL(:,Regs) - resCOM.results.EL_EXCESS(:,Regs)),1);
resultsSC.COM.DemP = resCOM.results.RPVO_EL(:,Regs) - resCOM.results.EL_SBAT(:,Regs) + resCOM.results.SBAT_EL(:,Regs) - resCOM.results.EL_EXCESS(:,Regs);


resultsSC.IND.OPT_SIZE_RPVO = resIND.results.OPT_SIZE_RPVO(Regs);
resultsSC.IND.OPT_SIZE_SBAT = resIND.results.OPT_SIZE_SBAT(Regs);
resultsSC.IND.OPT_SIZE_IBAT = resIND.results.OPT_SIZE_IBAT(Regs);
resultsSC.IND.RPVO_EL = resIND.results.RPVO_EL(:,Regs);
resultsSC.IND.EL_SBAT = resIND.results.EL_SBAT(:,Regs);
resultsSC.IND.SBAT_EL = resIND.results.SBAT_EL(:,Regs);
resultsSC.IND.EL_EXCESS = resIND.results.EL_EXCESS(:,Regs);
resultsSC.IND.EL_GRID = resIND.results.EL_GRID(:,Regs);
resultsSC.IND.Dem = sum((resIND.results.RPVO_EL(:,Regs) - resIND.results.EL_SBAT(:,Regs) + resIND.results.SBAT_EL(:,Regs) - resIND.results.EL_EXCESS(:,Regs)),1);
resultsSC.IND.DemP = resIND.results.RPVO_EL(:,Regs) - resIND.results.EL_SBAT(:,Regs) + resIND.results.SBAT_EL(:,Regs) - resIND.results.EL_EXCESS(:,Regs);


