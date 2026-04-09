function [instalations,active] = readProsInstalations(setup,systemParams,Reg1,Regs, varargin)
%
% FUNCTION readProsInstalations(setup, systemParams, Reg1, Reg2, varargin)
%
% Reads installation data for prosumers across two regions.
%
%
% INPUT:
%            setup:         Structure that contains all necessary settings and data for processing.
%            systemParams:  Structure with system parameters related to prosumers.
%            Reg1:          Name or index of the first region.
%            Regs:          Regions for data collection.
%            varargin:      ADD!!!!!!!!!!!!!!!!!!
%
% OUTPUT:
%            instalations:  Structure with installation values.
%            active:        ADD!!!!!!!!!!!!!!!!!!
%
%Dmitrii Bogdanov
%last change 19.02.2026

if nargin == 5
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
instalations = zeros(size(Regs,2),size(systemParams.Instalations,2),size(systemParams.Instalations,3));

if setup.OvernightFlag
    if setup.Heat.Flag
        projName = ['SC_HE_Overnight_' num2str(setup.startYear)];
    else
        projName = ['SC_Overnight_' num2str(setup.startYear)];
    end
else
    if setup.Heat.Flag
        projName = 'SC_HE';
    else
        projName = 'SC';
    end
end

for Reg = 1:size(Regs,2)
    
    RegS = Reg-Reg1+1;
    instRES = load([setup.rootDir filesep 'projects' filesep projName filesep 'output' filesep 'Instalations' '_RES_' num2str(Regs(Reg)) nameAdd]);
    instCOM = load([setup.rootDir filesep 'projects' filesep projName filesep 'output' filesep 'Instalations' '_COM_' num2str(Regs(Reg)) nameAdd]);
    instIND = load([setup.rootDir filesep 'projects' filesep projName filesep 'output' filesep 'Instalations' '_IND_' num2str(Regs(Reg)) nameAdd]);

    for ii = 1:length(systemParams.IndexID)

        if instRES.systemParams.Active(ii)==1 & ~sum(ismember({'RPVO','SBAT','IBAT'},systemParams.IndexID{ii}))

            instalations(Reg,:,ii) = instalations(Reg,:,ii) + instRES.systemParams.Instalations(1,:,ii) + instCOM.systemParams.Instalations(1,:,ii) + instIND.systemParams.Instalations(1,:,ii);

        end
    end
    
    instalations(Reg,:,ismember(systemParams.IndexID,'RPVR')) = instRES.systemParams.Instalations(1,:,ismember(systemParams.IndexID,'RPVR'));
    instalations(Reg,:,ismember(systemParams.IndexID,'RPVC')) = instCOM.systemParams.Instalations(1,:,ismember(systemParams.IndexID,'RPVC'));
    instalations(Reg,:,ismember(systemParams.IndexID,'RPVI')) = instIND.systemParams.Instalations(1,:,ismember(systemParams.IndexID,'RPVI'));
    
    instalations(Reg,:,ismember(systemParams.IndexID,'SBAR')) = instRES.systemParams.Instalations(1,:,ismember(systemParams.IndexID,'SBAR'));
    instalations(Reg,:,ismember(systemParams.IndexID,'SBAC')) = instCOM.systemParams.Instalations(1,:,ismember(systemParams.IndexID,'SBAC'));
    instalations(Reg,:,ismember(systemParams.IndexID,'SBAI')) = instIND.systemParams.Instalations(1,:,ismember(systemParams.IndexID,'SBAI'));
    
    instalations(Reg,:,ismember(systemParams.IndexID,'IBAR')) = instRES.systemParams.Instalations(1,:,ismember(systemParams.IndexID,'IBAR'));
    instalations(Reg,:,ismember(systemParams.IndexID,'IBAC')) = instCOM.systemParams.Instalations(1,:,ismember(systemParams.IndexID,'IBAC'));
    instalations(Reg,:,ismember(systemParams.IndexID,'IBAI')) = instIND.systemParams.Instalations(1,:,ismember(systemParams.IndexID,'IBAI'));
        
end

active = instRES.systemParams.Active;