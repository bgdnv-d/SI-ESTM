function [LCOT,LCOT_average] = LCOT(results,systemParams,activeElements,z,prod_electricity,maxTransmPerRegion,demand,costYear)
%
% [LCOT,LCOT_AVERAGE] = LCOT(RESULTS,SYSTEMPARAMS,ACTIVEELEMENTS,Z,PROD_ELECTRICITY)
% Returns Levelized Cost of Transmission (LCOT) for each node and for the whole
% system (LCOT_AVERAGE). The input parameter z assigns cost to importing
% nodes/ regions (z=0) or exporting nodes/ regions (z=1).
%
%last change 13th Nov 2016 by Dmitrii Bogdanov


demand(demand<0) = 0;

yearNumb = find(systemParams.IndexYears==costYear);

CRF_TL = capitalRecoveryFactor(systemParams.WACC,systemParams.Lifetime(ismember(systemParams.IndexID,'TRTL')));
CRF_CS = capitalRecoveryFactor(systemParams.WACC,systemParams.Lifetime(ismember(systemParams.IndexID,'TRCS')));
% power on each transmission line in one variable
transPower_AC = results.LINEpos_AC - results.LINEneg_AC;
transPower_DC = results.LINEpos_DC - results.LINEneg_DC;

% This formulation includes cost of converter station
if isfield(activeElements,'activeTransmission')
    capexTL_DC_CRF = systemParams.Capex(ismember(systemParams.IndexID,'TRTL'),yearNumb)  .* capitalRecoveryFactor(systemParams.WACC,systemParams.Lifetime(ismember(systemParams.IndexID,'TRTL')));
    capexTL_AC_CRF = systemParams.Capex(ismember(systemParams.IndexID,'THAO'),yearNumb)  .* capitalRecoveryFactor(systemParams.WACC,systemParams.Lifetime(ismember(systemParams.IndexID,'THAO')));
    capexCS_CRF = systemParams.Capex(ismember(systemParams.IndexID,'TRCS'),yearNumb)  .* capitalRecoveryFactor(systemParams.WACC,systemParams.Lifetime(ismember(systemParams.IndexID,'TRCS')));
    capexTL_DC = systemParams.Capex(ismember(systemParams.IndexID,'TRTL'),yearNumb);
    capexTL_AC = systemParams.Capex(ismember(systemParams.IndexID,'THAO'),yearNumb);
    capexCS = systemParams.Capex(ismember(systemParams.IndexID,'TRCS'),yearNumb);
    opexTL_DC_fix = systemParams.Opex_fix(ismember(systemParams.IndexID,'TRTL'),yearNumb);
    opexTL_AC_fix = systemParams.Opex_fix(ismember(systemParams.IndexID,'THAO'),yearNumb);
    opexCS_fix = systemParams.Opex_fix(ismember(systemParams.IndexID,'TRCS'),yearNumb);

    opexTL_DC_var = systemParams.Opex_var(ismember(systemParams.IndexID,'TRTL'),yearNumb);
    opexTL_AC_var = systemParams.Opex_var(ismember(systemParams.IndexID,'THAO'),yearNumb);
    opexCS_var = systemParams.Opex_var(ismember(systemParams.IndexID,'TRCS'),yearNumb);


end

LCOT_total = ((capexTL_DC_CRF + opexTL_DC_fix) * sum(systemParams.TLlength .* results.OPT_SIZE_TRTL')) + sum(sum(abs(transPower_DC))) * opexTL_DC_var+...
    ((capexTL_AC_CRF + opexTL_AC_fix) * sum(systemParams.TLlength .* results.OPT_SIZE_THAO')) + sum(sum(abs(transPower_AC))) * opexTL_AC_var+...
    ((capexCS_CRF + opexCS_fix) * sum(results.OPT_SIZE_TRTL));

LCOT_total(isnan(LCOT_total)) = 0;
LCOT_average = LCOT_total / sum(demand);

% calculating share of each node at total transmission power
importPower = zeros(size(results.GRID_AC));
exportPower = zeros(size(results.GRID_AC));

importPower((results.GRID_AC + results.GRID_DC) >0) = (results.GRID_AC((results.GRID_AC + results.GRID_DC)>0) + results.GRID_DC((results.GRID_AC + results.GRID_DC)>0));
exportPower((results.GRID_AC + results.GRID_DC) <0) = (results.GRID_AC((results.GRID_AC + results.GRID_DC)<0) + results.GRID_DC((results.GRID_AC + results.GRID_DC)<0));


shareImport = sum(importPower,1) / sum(sum(importPower));
shareImport(isnan(shareImport)) = 0;
shareExport = sum(exportPower,1) / sum(sum(exportPower));
shareExport(isnan(shareExport)) = 0;

LCOT_total_per_Node = LCOT_total * (z * shareExport + (1-z) * shareImport);

if (sum(LCOT_total_per_Node)/LCOT_total>1.01)|(sum(LCOT_total_per_Node)/LCOT_total<0.99)
    warning('Problem with grid shares')
    LCOT_total_per_Node = LCOT_total*(demand/sum(demand))';
end


LCOT = LCOT_total_per_Node ./ (demand)';