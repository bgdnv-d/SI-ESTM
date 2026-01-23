function writeMiscProsumers(fileName,setup,struct,activeElements,transCapa,fossilLimit,storePeriod,el_cost,feedin_cost,shareOfSector,shareIndHeatPros,costYear)
%
% FUNCTION writeMiscProsumers(fileName, setup, struct, activeElements, transCapa, fossilLimit, storePeriod, el_cost, feedin_cost, shareOfSector, shareIndHeatPros, costYear)
%
% Writes miscellaneous data for prosumer scenarios to a file.
%
% INPUT:
%            fileName: Desired filename without extension
%            setup: Structure that contains all necessary settings and data for processing
%            struct: Parameter structure obtained from Excel sheets
%            activeElements: Structure with active elements in the system
%            transCapa: Transmission capacity data
%            fossilLimit: Limit on total fossil fuel use
%            storePeriod: Storage period
%            el_cost: Electricity cost for prosumers
%            feedin_cost: Feedin tariff or revenue for prosumers
%            shareOfSector: ADD!!!!!!!!!!!!!
%            shareIndHeatPros: ADD!!!!!!!!!!!!!
%            costYear: Year to which costs are adjusted
%
% OUTPUT:
%            -
%
%Dmitrii Bogdanov
%last change 24.07.2025


% input data in kW, kWh, kg for CO2 and tonnes of product for industry
try setup.units;
catch
    setup.units = 'kW';
end
switch setup.units
    case 'kW'
        m_factor = 1;
    case 'MW'
        m_factor = 0.001;
    case 'GW'
        m_factor = 0.000001;
end

WriteParamFileBegin(fileName);

WriteParams(fileName,'el_cost',round(el_cost,3),[1:length(struct.IndexNodes)]')
WriteParams(fileName,'feedin_cost',round(feedin_cost,3))

WriteParams(fileName,'transCapa',round(transCapa*m_factor,3))
WriteParams(fileName,'fossilLimit',round(fossilLimit*m_factor,3))


WriteParams(fileName,'feedInEfficiency',round(struct.FeedInEfficiencies(ismember(struct.IndexIDE,activeElements.labels.elFeedIn))',3),activeElements.labels.elFeedIn)

WriteParams(fileName,'urbanPop',round(struct.Urbanisation'),[1:length(struct.IndexNodes)]')

WriteParams(fileName,'shareOfSector',round(shareOfSector,3),[1:length(struct.IndexNodes)]')
WriteParams(fileName,'shareIndHeatPros',round(shareIndHeatPros,3),[1:length(struct.IndexNodes)]')
WriteParams(fileName,'shareOfSectorEl',round(struct.shareOfSectorEl,3),[1:length(struct.IndexNodes)]')
WriteParams(fileName,'shareOfDistrHeat',round(setup.Heat.shareOfDistrHeat(:,struct.IndexYears==costYear),3),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'shareOfLowHeatInd',round(setup.Heat.shareOfLowHeatInd,3),[1:length(struct.IndexNodes)]');
WriteParams(fileName,'shareOfHighHeatInd',round(setup.Heat.shareOfHighHeatInd,3),[1:length(struct.IndexNodes)]');


WriteParams(fileName,'fossilCO2Cost',round(struct.fossilCO2Cost(struct.IndexYears==costYear),3));

WriteParamFileEnd(fileName);