function active = activeProsumers(systemParams)
%
% FUNCTION activeComponents(struct)
%
% defines which components are active in the prosumers model
%
% INPUT:
%            struct:    system parameters data structure.
%
% OUTPUT:
%            active:  vector defining if the tecnology is active in the prosumers model
%
%Dmitrii Bogdanov
%last change 22.07.2025


active = ismember(systemParams.Type,{'LHF','LH'}) | ismember(systemParams.IndexID,{'RPVO','SBAT','IBAT','SHHS','IHHS'});