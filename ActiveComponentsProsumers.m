function activeElements = activeComponentsProsumers(struct,activeElements)
%
% FUNCTION activeComponentsProsumers(struct,activeElements)
%
% Redefines activeElements structure with elements active in the prosumers sub-model
% Returns activeElements structure containing all components active in prosumers model
%
%
% INPUT:
%            struct:    system parameters data structure.
%            activeElements:  activeElements structure defined for the centralised system model.
%
% OUTPUT:
%            activeElements:  activeElements structure redefined for the prosumers sub-model.
%
%Dmitrii Bogdanov
%last change 22.07.2025


activeElements.labels.elTransformer = {};
activeElements.labels.chpTransformer = {};
activeElements.labels.heatTransformer = activeElements.labels.localHeatTransformer;
activeElements.labels.localHeatTransformer;
activeElements.labels.distrHeatTransformer = struct.IndexID(activeElements.activeDHeatTransformer);
activeElements.labels.distrFHeatTransformer = {};
activeElements.labels.distrEHeatTransformer = {};
activeElements.labels.gasTransformer = {};
activeElements.labels.elFeedIn = {'RPVO'};
activeElements.labels.heatFeedIn = activeElements.labels.localHeatFeedIn;
activeElements.labels.localHeatFeedIn;
activeElements.labels.districtHeatFeedIn = {};
activeElements.labels.storage = [{'SBAT'};{'SHHS'}];

activeElements.labels.storageInterface = [{'IBAT'};{'IHHS'}];

activeElements.labels.resource = struct.IndexID(activeElements.activeResource);
activeElements.labels.hydro = {};
activeElements.labels.desalination = {};
activeElements.labels.load = [{'LELE'};{'LHSP'};{'LHDW'}];

if isfield(struct,'IndexIDTr')
    activeElements.labels.transmission = struct.IndexID(activeElements.activeTransmission);
end



% active components
activeElements.activeElTransformer=ismember(struct.IndexID,activeElements.labels.elTransformer);
activeElements.activeCHPTransformer = ismember(struct.IndexID,activeElements.labels.chpTransformer);
activeElements.activeLHeatTransformer = ismember(struct.IndexID,activeElements.labels.localHeatTransformer);
activeElements.activeDFHeatTransformer = ismember(struct.IndexID,activeElements.labels.distrFHeatTransformer);
activeElements.activeDEHeatTransformer = ismember(struct.IndexID,activeElements.labels.distrEHeatTransformer);
activeElements.activeDHeatTransformer = ismember(struct.IndexID,activeElements.labels.distrHeatTransformer);
activeElements.activeHeatTransformer = ismember(struct.IndexID,activeElements.labels.heatTransformer);
activeElements.activeGasTransformer = ismember(struct.IndexID,activeElements.labels.gasTransformer);
activeElements.activeElFeedIn = ismember(struct.IndexID,activeElements.labels.elFeedIn);
activeElements.activeLHeatFeedIn = ismember(struct.IndexID,activeElements.labels.localHeatFeedIn);
activeElements.activeDHeatFeedIn = ismember(struct.IndexID,activeElements.labels.districtHeatFeedIn);
activeElements.activeHeatFeedIn= ismember(struct.IndexID,activeElements.labels.heatFeedIn);


activeElements.activeHydro = ismember(struct.IndexID,activeElements.labels.hydro);
activeElements.activeStorage = ismember(struct.IndexID,activeElements.labels.storage);
activeElements.activeStorageInterface = ismember(struct.IndexID,activeElements.labels.storageInterface);

activeElements.activeResource = ismember(struct.IndexID,activeElements.labels.resource);
activeElements.activeDesalination = ismember(struct.IndexID,activeElements.labels.desalination);
activeElements.activeLoad = ismember(struct.IndexID,activeElements.labels.load);
if isfield(struct,'IndexIDTr')
    activeElements.activeTransmission = ismember(struct.IndexID,activeElements.labels.transmission);
end
activeElements.activeAllTransformer = activeElements.activeElTransformer + activeElements.activeCHPTransformer + activeElements.activeHeatTransformer + activeElements.activeGasTransformer;

struct.Active = activeElements.activeAllTransformer + activeElements.activeElFeedIn + activeElements.activeHeatFeedIn + activeElements.activeHydro + activeElements.activeStorage + activeElements.activeStorageInterface + activeElements.activeResource + activeElements.activeDesalination + activeElements.activeLoad + activeElements.activeTransmission;

activeElements.labels.allTransformer = struct.IndexID(logical(activeElements.activeAllTransformer));
activeElements.labels.all = struct.IndexID(logical(struct.Active));

