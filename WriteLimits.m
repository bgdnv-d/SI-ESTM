function limits = writeLimits(fileName,setup,struct,activeElements)

% write limits to gmpl conform data file
% For a clear arrangement of all components: Please use the same component
% abbreviation and the same order, used in the Modelparams~.xls files,
% when adding new components.
%
%Dmitrii Bogdanov
%last change 24.07.2025

% input data in kW, kWh, kg for CO2 and tonnes of product for industry, m3
% per day for desalination
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

technames = char(struct.IndexID(struct.Active & ~activeElements.activeResource & ~activeElements.activeTransmission & ~activeElements.activeLoad));



WriteParamFileBegin(fileName);

struct.SizeLimits = floor(struct.SizeLimits);

% UPPER LIMITATIONS (lower limitations below)
% ------------------------------------------------------------------------------------------------------------------------------
for k=1:length(struct.IndexNodes)

    struct.SizeLimits(ismember(struct.IndexID,'TWEL'),2,k);
    for i = 1:size(technames,1)
        eval([technames(i,:) '_UpperLim(k) = ceil(struct.SizeLimits(ismember(struct.IndexID,' '''' technames(i,:) '''' '),2,k))*m_factor;']);
    end

end

% In gmpl it is not allowed to set a parameter with inf, which is why the "~isinf"-command is applied below
% only for the limit values and the relevant nodes:
NodesVektor = [1:length(struct.IndexNodes)];
for i = 1:size(technames,1)
    eval(['WriteParams(fileName,' '''' technames(i,:) '_UpperLim' '''' ', min(' technames(i,:) '_UpperLim(~isinf(' technames(i,:) '_UpperLim))' '''' ', 10^11),[NodesVektor(~isinf(' technames(i,:) '_UpperLim))])']);
end


% SETS FOR NODES WITH UPPER LIMITATIONS (lower limitations below)
% ------------------------------------------------------------------------------------------------------------------------------
% Use SingleSet for defining sets for all nodes which have limits (see rationale above)
% Relevant nodes must be defined and converted to string values

for i = 1:size(technames,1)
    eval([technames(i,:) '_UpperLimNodesVektor = NodesVektor(~isinf(' technames(i,:) '_UpperLim));if length(' technames(i,:) '_UpperLimNodesVektor) ~= 0,    for k=1:length(' technames(i,:) '_UpperLimNodesVektor),    ' technames(i,:) '_UpperLimNodes{k} = num2str(' technames(i,:) '_UpperLimNodesVektor(k));    end,    SingleSet(fileName,' '''' technames(i,:) '_UpperLimNodes' '''' ',' technames(i,:) '_UpperLimNodes), else SingleSet(fileName,' '''' technames(i,:) '_UpperLimNodes' '''' ',' '''' '''' '), end']);
end


% ##############################################################################################################################

% LOWER LIMITATIONS
% ------------------------------------------------------------------------------------------------------------------------------
for k=1:length(struct.IndexNodes)
    for i = 1:size(technames,1)
        eval([technames(i,:) '_LowerLim(k) = floor(struct.SizeLimits(ismember(struct.IndexID,' '''' technames(i,:) '''' '),1,k))*m_factor;']);
    end

end

% Components have always a lower limit (either 0 or a given value), which is why the find-command (instead of the ~isinf-command with upper limits)
% is applied in the following. Assigning lowerLimNodes through the find-command allows for a better handling within the glpk model
% (e.g. uncommenting of inactive components in the lower limits constraints is not necessary)
NodesVektor = [1:length(struct.IndexNodes)];

for i = 1:size(technames,1)
    eval(['WriteParams(fileName,' '''' technames(i,:) '_LowerLim' '''' ',' technames(i,:) '_LowerLim(~isinf(' technames(i,:) '_LowerLim))' '''' ',[NodesVektor(~isinf(' technames(i,:) '_LowerLim))])']);
end



% SETS FOR NODES WITH LOWER LIMITATIONS
% ------------------------------------------------------------------------------------------------------------------------------
% Use SingleSet for defining sets for all nodes which have limits (see rationale above)
% Relevant nodes must be defined and converted to string values

for i = 1:size(technames,1)
    eval([technames(i,:) '_LowerLimNodesVektor = NodesVektor;if length(' technames(i,:) '_LowerLimNodesVektor) ~= 0,    for k=1:length(' technames(i,:) '_LowerLimNodesVektor),    ' technames(i,:) '_LowerLimNodes{k} = num2str(' technames(i,:) '_LowerLimNodesVektor(k));    end,    SingleSet(fileName,' '''' technames(i,:) '_LowerLimNodes' '''' ',' technames(i,:) '_LowerLimNodes), else SingleSet(fileName,' '''' technames(i,:) '_LowerLimNodes' '''' ',' '''' '''' '), end']);
end



% ##############################################################################################################################

WriteParams(fileName,'totalFuelLimit',round(struct.ValueResourceTotal(find(ismember(struct.IndexIDR,activeElements.labels.resource)),:)'*m_factor,3),[1:length(struct.IndexNodes)],activeElements.labels.resource);

WriteParamFileEnd(fileName);

limits.q=0;
