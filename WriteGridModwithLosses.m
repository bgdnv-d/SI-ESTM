function writeGridModwithLosses(filename,setup,gridMat,systemParams)
% FUNCTION writeGridModwithLosses(filename, setup, gridMat, systemParams)
%
% Writes grid structure in .dat format(GMPL-readable)
%
% INPUT:
%            filename: Desired filename without extension
%            setup: Structure that contains all necessary settings and data for processing
%            gridMat: ADD!!!!!!!!!!!!!!!!!
%            systemParams: Structure with system parameters for the scenario
%
% OUTPUT:
%            -
%
%Dmitrii Bogdanov
%last change 24.07.2025


fID = fopen(filename,'w');
effCS = systemParams.TransmissionEfficiencies(ismember(systemParams.IndexIDTr,'TRCS'));
%HVDC
for nodeNo=1:max(max(gridMat))
    neg  = find(gridMat(:,1) == nodeNo);
    pos  = find(gridMat(:,2) == nodeNo);

    fprintf(fID,'s.t. gridDef%d{t in T}: transPower_DC[t,%d] = ',nodeNo,nodeNo);

    for k=1:length(neg)
        effTL = 1-(1-systemParams.TransmissionEfficiencies(ismember(systemParams.IndexIDTr,'TRTL')))*systemParams.TLlength(neg(k))/1000;
      	 if k==1
      		fprintf(fID,'%.3f*line_DC[t,%d,''neg''] - line_DC[t,%d,''pos''] ',round(effCS*effTL,3),neg(k),neg(k));
         elseif k>1
      		fprintf(fID,'+ %.3f*line_DC[t,%d,''neg''] - line_DC[t,%d,''pos''] ',round(effCS*effTL,3),neg(k),neg(k));
         end
    end


    for k=1:length(pos)
        effTL = 1-(1-systemParams.TransmissionEfficiencies(ismember(systemParams.IndexIDTr,'TRTL')))*systemParams.TLlength(pos(k))/1000;
      	 if (k==1 && isempty(neg))
      		fprintf(fID,'%.3f*line_DC[t,%d,''pos''] - line_DC[t,%d,''neg''] ',round(effCS*effTL,3),pos(k),pos(k));
         else
             fprintf(fID,'+ %.3f*line_DC[t,%d,''pos''] - line_DC[t,%d,''neg''] ',round(effCS*effTL,3),pos(k),pos(k));
         end
    end

    fprintf(fID,';\n');
    clear neg pos k;
end

%HVAC

for nodeNo=1:max(max(gridMat))
    neg  = find(gridMat(:,1) == nodeNo);
    pos  = find(gridMat(:,2) == nodeNo);

    fprintf(fID,'s.t. gridDef_AC%d{t in T}: transPower_AC[t,%d] = ',nodeNo,nodeNo);

    for k=1:length(neg)
        effTL = 1-(1-systemParams.TransmissionEfficiencies(ismember(systemParams.IndexIDTr,'THAO')))*systemParams.TLlength(neg(k))/1000;
      	 if k==1
      		fprintf(fID,'%.3f*line_AC[t,%d,''neg''] - line_AC[t,%d,''pos''] ',round(effTL,3),neg(k),neg(k));
         elseif k>1
      		fprintf(fID,'+ %.3f*line_AC[t,%d,''neg''] - line_AC[t,%d,''pos''] ',round(effTL,3),neg(k),neg(k));
         end
    end


    for k=1:length(pos)
        effTL = 1-(1-systemParams.TransmissionEfficiencies(ismember(systemParams.IndexIDTr,'THAO')))*systemParams.TLlength(pos(k))/1000;
      	 if (k==1 && isempty(neg))
      		fprintf(fID,'%.3f*line_AC[t,%d,''pos''] - line_AC[t,%d,''neg''] ',round(effTL,3),pos(k),pos(k));
         else
      		fprintf(fID,'+ %.3f*line_AC[t,%d,''pos''] - line_AC[t,%d,''neg''] ',round(effTL,3),pos(k),pos(k));
         end
    end

    fprintf(fID,';\n');
    clear neg pos k;
end

fclose(fID);