function [distanceMatrix,linesDistances,gridMat,LowLimits_DC,UpLimits_DC,LowLimits_AC,UpLimits_AC] = getGridFromExcel(filename)
%
% FUNCTION getGridFromExcel(filename)
%
% Reads xls or xlsx file with power grids structure definition.
% Generates data structure for the model system parameters
%
%
% INPUTS:
%            filename:  xls or xlsx file with power grids structure definition.
%
% OUTPUTS:
%            distanceMatrix:  	grid connection distances in matrix form
%			 linesDistances:  	length of grid connection between regions
%			 gridMat:  			connected regions index
%			 LowLimits_DC:  	existing capacity of HVDC connections
%			 UpLimits_DC:  		upper limits for HVAC connections
%			 LowLimits_AC:  	existing capacity of HVAC connections
%			 UpLimits_AC:       upper limits for HVAC connections
%
%Dmitrii Bogdanov
%last change 22.07.2025


%% HVAC power lines upper limits
[num,txt,raw] = xlsread(filename,'UpLimits_AC');

% extract distances matrix
distanceMatrix = num(2:end,2:end);

% get transmission lines form matrix
UpLimits_AC = nan(prod(size(distanceMatrix)),1);
for k=1:size(distanceMatrix,1)
    for m=1:size(distanceMatrix,2)
        linesTmp((k-1)*size(distanceMatrix,2) + m,1) = distanceMatrix(k,m);
    end
end

UpLimits_AC = linesTmp(~isnan(linesTmp));

%% HVDC power lines upper limits
[num,txt,raw] = xlsread(filename,'UpLimits_DC');

% extract distances matrix
distanceMatrix = num(2:end,2:end);

% get transmission lines form matrix
UpLimits_DC = nan(prod(size(distanceMatrix)),1);
for k=1:size(distanceMatrix,1)
    for m=1:size(distanceMatrix,2)
        linesTmp((k-1)*size(distanceMatrix,2) + m,1) = distanceMatrix(k,m);
    end
end

UpLimits_DC = linesTmp(~isnan(linesTmp));


%% HVDC power lines lower limit
[num,txt,raw] = xlsread(filename,'LowLimits_DC');

% extract distances matrix
distanceMatrix = num(2:end,2:end);

% get transmission lines form matrix
LowLimits_DC = nan(prod(size(distanceMatrix)),1);
for k=1:size(distanceMatrix,1)
    for m=1:size(distanceMatrix,2)
        linesTmp((k-1)*size(distanceMatrix,2) + m,1) = distanceMatrix(k,m);
    end
end

LowLimits_DC = linesTmp(~isnan(linesTmp));

%% HVAC power lines lower limit
[num,txt,raw] = xlsread(filename,'LowLimits_AC');

% extract distances matrix
distanceMatrix = num(2:end,2:end);

% get transmission lines form matrix
LowLimits_AC = nan(prod(size(distanceMatrix)),1);
for k=1:size(distanceMatrix,1)
    for m=1:size(distanceMatrix,2)
        linesTmp((k-1)*size(distanceMatrix,2) + m,1) = distanceMatrix(k,m);
    end
end

LowLimits_AC = linesTmp(~isnan(linesTmp));


%% power lines lengths
[num,txt,raw] = xlsread(filename,'Length');

% extract distances matrix
distanceMatrix = num(2:end,2:end);

% get transmission lines form matrix
linesDistances = nan(prod(size(distanceMatrix)),1);
for k=1:size(distanceMatrix,1)
    for m=1:size(distanceMatrix,2)
        linesTmp((k-1)*size(distanceMatrix,2) + m,1) = distanceMatrix(k,m);
    end
end

linesDistances = linesTmp(~isnan(linesTmp));


% determine pairs of interconnected regions

% empty matrix of grid connections
gridMat = [];

for l=1:size(distanceMatrix,1)
    tmp = find(~isnan(distanceMatrix(l,:)));
    for m=1:length(tmp)
        gridMat = [gridMat; l tmp(m)];
    end
    tmp = [];
end