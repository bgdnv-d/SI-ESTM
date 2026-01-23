function activeElements = activeComponentsMain(struct)
%
% FUNCTION activeComponentsMain(struct)
%
% defines which components are active in the centralised model based on the system parametes
% Returns activeElements structure containing all components active in the centralised system model
%
%
% INPUT:
%            struct:    system parameters data structure.
%
% OUTPUT:
%            activeElements:  activeElements structure defined for the centralised system model.
%
%Dmitrii Bogdanov
%last change 22.07.2025


% technology type
electTransf = 'E';
chpTransf = 'C';
heatLTransf	= 'LH';
heatDFransf = 'DFH';
heatDETransf = 'DEH';
gasTransf = 'G';
electrFeedIn = 'EF';
heatLFeedIn ='LHF';
heatDFeedIn='DHF';

Hydro='HD';
Storage='S';
StorageInterface='I';
Resource='FU';
Desal='D';
TLcomponent='GR';
Load = 'L';


% active components
activeElements.activeElTransformer = and(ismember(struct.Type,electTransf),struct.Active);
activeElements.activeCHPTransformer = and(ismember(struct.Type,chpTransf),struct.Active);
activeElements.activeLHeatTransformer = and(ismember(struct.Type,heatLTransf),struct.Active);
activeElements.activeDFHeatTransformer = and(ismember(struct.Type,heatDFransf),struct.Active);
activeElements.activeDEHeatTransformer = and(ismember(struct.Type,heatDETransf),struct.Active);
activeElements.activeDHeatTransformer = or(activeElements.activeDEHeatTransformer,activeElements.activeDFHeatTransformer);
activeElements.activeHeatTransformer = or(activeElements.activeDHeatTransformer,activeElements.activeLHeatTransformer);
activeElements.activeHeatTransformer = (activeElements.activeDHeatTransformer);
activeElements.activeGasTransformer = and(ismember(struct.Type,gasTransf),struct.Active);
activeElements.activeElFeedIn = and(ismember(struct.Type,electrFeedIn),struct.Active);
activeElements.activeLHeatFeedIn = and(ismember(struct.Type,heatLFeedIn),struct.Active);
activeElements.activeDHeatFeedIn = and(ismember(struct.Type,heatDFeedIn),struct.Active);
%activeElements.activeHeatFeedIn= or(activeElements.activeDHeatFeedIn,activeElements.activeLHeatFeedIn);
activeElements.activeHeatFeedIn= (activeElements.activeDHeatFeedIn);


activeElements.activeHydro = and(ismember(struct.Type,Hydro),struct.Active);
activeElements.activeStorage = and(ismember(struct.Type,Storage),struct.Active);
try
    activeElements.activeStorageInterface = and(ismember(struct.Type,StorageInterface),struct.Active);
catch
end
activeElements.activeResource = and(ismember(struct.Type,Resource),struct.Active);
activeElements.activeDesalination = and(ismember(struct.Type,Desal),struct.Active);
activeElements.activeLoad = and(ismember(struct.Type,Load),struct.Active);

activeElements.activeAllTransformer = logical(activeElements.activeElTransformer+activeElements.activeCHPTransformer+...
    activeElements.activeHeatTransformer+activeElements.activeGasTransformer);
if isfield(struct,'IndexIDTr')
    activeElements.activeTransmission = and(ismember(struct.Type,TLcomponent),struct.Active);
end

% main storage update
activeElements.activeStorage(ismember(struct.IndexID,{'SHHS'})) = 0;
activeElements.activeStorageInterface(ismember(struct.IndexID,{'IHHS'})) = 0;

activeElements.labels.all = struct.IndexID(struct.Active);
activeElements.labels.elTransformer = struct.IndexID(activeElements.activeElTransformer);
activeElements.labels.chpTransformer = struct.IndexID(activeElements.activeCHPTransformer);
activeElements.labels.heatTransformer = struct.IndexID(activeElements.activeHeatTransformer);
activeElements.labels.localHeatTransformer = struct.IndexID(activeElements.activeLHeatTransformer);
activeElements.labels.distrHeatTransformer = struct.IndexID(activeElements.activeDHeatTransformer);
activeElements.labels.distrFHeatTransformer = struct.IndexID(activeElements.activeDFHeatTransformer);
activeElements.labels.distrEHeatTransformer = struct.IndexID(activeElements.activeDEHeatTransformer);
activeElements.labels.gasTransformer = struct.IndexID(activeElements.activeGasTransformer);
activeElements.labels.elFeedIn = struct.IndexID(activeElements.activeElFeedIn);
activeElements.labels.heatFeedIn = struct.IndexID(activeElements.activeHeatFeedIn);
activeElements.labels.localHeatFeedIn = struct.IndexID(activeElements.activeLHeatFeedIn);
activeElements.labels.districtHeatFeedIn = struct.IndexID(activeElements.activeDHeatFeedIn);
activeElements.labels.storage = struct.IndexID(activeElements.activeStorage);
try
    activeElements.labels.storageInterface = struct.IndexID(activeElements.activeStorageInterface);
catch
end
activeElements.labels.resource = struct.IndexID(activeElements.activeResource);
activeElements.labels.hydro = struct.IndexID(activeElements.activeHydro);
activeElements.labels.desalination = struct.IndexID(activeElements.activeDesalination);
activeElements.labels.load = struct.IndexID(activeElements.activeLoad);
activeElements.labels.allTransformer = struct.IndexID(activeElements.activeAllTransformer);
if isfield(struct,'IndexIDTr')
    activeElements.labels.transmission = struct.IndexID(activeElements.activeTransmission);
end
