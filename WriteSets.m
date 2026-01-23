function writeSets(fileName,setup,struct,activeElements)
% specifies the sets with active elements of the centralise energy system
%
% last change Dmitrii Bogdanov 07.11.2023


for k=1:length(struct.IndexNodes)
    nodeNumber{k} = num2str(k);
end
TLlength = [];
for k=1:length(struct.TLlength)
    TLlength{k} = num2str(k);
end

WriteParamFileBegin(fileName);

SingleSet(fileName,'N',nodeNumber)
SingleSet(fileName,'L',TLlength)

SingleSet(fileName,'TDir',{'pos','neg'})
SingleSet(fileName,'TLcomponent',activeElements.labels.transmission)


SingleSet(fileName,'Resource',activeElements.labels.resource)

SingleSet(fileName,'Storage',activeElements.labels.storage)

try
    SingleSet(fileName,'StorInt',activeElements.labels.storageInterface)
catch
end

SingleSet(fileName,'Hydro',activeElements.labels.hydro)

SingleSet(fileName,'electTransf',activeElements.labels.elTransformer)

SingleSet(fileName,'heatTransf',activeElements.labels.distrHeatTransformer)

SingleSet(fileName,'heatDTransf',activeElements.labels.distrHeatTransformer)
SingleSet(fileName,'heatLTransf',activeElements.labels.localHeatTransformer)
SingleSet(fileName,'heatDFTransf',activeElements.labels.distrFHeatTransformer)
SingleSet(fileName,'heatDETransf',activeElements.labels.distrEHeatTransformer)

SingleSet(fileName,'chpTransf',activeElements.labels.chpTransformer)

SingleSet(fileName,'electFeedIn',activeElements.labels.elFeedIn)
SingleSet(fileName,'heatFeedIn',activeElements.labels.heatFeedIn)
SingleSet(fileName,'heatLFeedIn',activeElements.labels.localHeatFeedIn)
SingleSet(fileName,'heatDFeedIn',activeElements.labels.districtHeatFeedIn)

SingleSet(fileName,'Desal',activeElements.labels.desalination)

SingleSet(fileName,'ptgTransf',activeElements.labels.gasTransformer)


SingleSet(fileName,'TranspPrim',struct.Mobility.Cons_Names_Primary)
SingleSet(fileName,'TranspSecond',struct.Mobility.Cons_Names_Secondary)

SingleSet(fileName,'iceTypes',{'MRLI';'MRWI';'MRBI';'MRMI';'MRHI'})
SingleSet(fileName,'bevTypes',{'MRLB';'MRWB';'MRBB';'MRMB';'MRHB'})
SingleSet(fileName,'fceTypes',{'MRLF';'MRWF';'MRBF';'MRMF';'MRHF'})
SingleSet(fileName,'pheTypes',{'MRLP';'MRWP';'MRBP';'MRMP';'MRHP'})


if setup.New
    SingleSet(fileName,'Load',activeElements.labels.load)
end

SingleSet(fileName,'Fuel',struct.IndexFuels)

WriteParamFileEnd(fileName);



