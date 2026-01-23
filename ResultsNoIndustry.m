function results = results_WO_Industry(results, numNodes, endHour)
%
% FUNCTION results_WO_Industry(results, numNodes, endHour)
%
% Removes industry-related results from the dataset.
%
%
% INPUT:
%            results:     Structure containing the original model results.
%            numNodes:    Number of nodes in the model.
%            endHour:     Last hour of the simulation time period.
%
% OUTPUT:
%            results:     Updated results structure without industry related data.
%
%Dmitrii Bogdanov
%last change 24.07.2025


results.TICB_PR = zeros(endHour,numNodes);
results.TICI_PR = zeros(endHour,numNodes);
results.TISB_PR = zeros(endHour,numNodes);
results.TISH_PR = zeros(endHour,numNodes);
results.TISR_PR = zeros(endHour,numNodes);
results.TISE_PR = zeros(endHour,numNodes);
results.TIAA_PR = zeros(endHour,numNodes);
results.TIAM_PR = zeros(endHour,numNodes);
results.TIAR_PR = zeros(endHour,numNodes);
results.TIPP_PR = zeros(endHour,numNodes);

results.EL_TICB = zeros(endHour,numNodes);
results.EL_TICI = zeros(endHour,numNodes);
results.EL_TISB = zeros(endHour,numNodes);
results.EL_TISH = zeros(endHour,numNodes);
results.EL_TISR = zeros(endHour,numNodes);
results.EL_TISE = zeros(endHour,numNodes);
results.EL_TIAM = zeros(endHour,numNodes);
results.EL_TIAR = zeros(endHour,numNodes);
results.EL_TIPP = zeros(endHour,numNodes);

results.HE_TICB = zeros(endHour,numNodes);
results.HE_TICI = zeros(endHour,numNodes);
results.HE_TISB = zeros(endHour,numNodes);
results.HE_TISH = zeros(endHour,numNodes);
results.HE_TISR = zeros(endHour,numNodes);
results.HE_TISE = zeros(endHour,numNodes);
results.HE_TIAA = zeros(endHour,numNodes);
results.HE_TIAR = zeros(endHour,numNodes);
results.HE_TIPP = zeros(endHour,numNodes);

results.HC_TISB = zeros(endHour,numNodes);
results.CC_TISH = zeros(endHour,numNodes);
results.CC_TISR = zeros(endHour,numNodes);
results.CC_TISE = zeros(endHour,numNodes);

results.TIAM_HE = zeros(endHour,numNodes);
results.TIAR_HE = zeros(endHour,numNodes);

results.EltoInd = zeros(endHour,numNodes);
results.HighHetoInd = zeros(endHour,numNodes);
results.MidHetoInd = zeros(endHour,numNodes);
results.MidHefromInd = zeros(endHour,numNodes);
results.LowHefromInd = zeros(endHour,numNodes);

results.TIPP_WWtot = zeros(endHour,numNodes);
