function systemParams = MobilityEfficiency(systemParams)
%
% FUNCTION MobilityEfficiency(systemParams)
%
% Calculates efficiency values for mobility technologies.
%
%
% INPUT:
%            systemParams:   Structure with system parameters to be updated.
%
% OUTPUT:
%            systemParams:   Updated structure with mobility efficiency values.
%
%Dmitrii Bogdanov
%last change 23.07.2025


%% Efficiencies
systemParams.Mobility.Cons_Primary_N = systemParams.Mobility.Cons_Values(ismember(systemParams.Mobility.Cons_Type,'P'),:);
systemParams.Mobility.Cons_Names_Primary = systemParams.Mobility.Cons_Names(ismember(systemParams.Mobility.Cons_Type,'P'));
systemParams.Mobility.Lifetime_Primary = systemParams.Mobility.Lifetime(ismember(systemParams.Mobility.Cons_Type,'P'));
systemParams.Mobility.Cons_Secondary_N = systemParams.Mobility.Cons_Values(ismember(systemParams.Mobility.Cons_Type,'S'),:);
systemParams.Mobility.Cons_Names_Secondary = systemParams.Mobility.Cons_Names(ismember(systemParams.Mobility.Cons_Type,'S'));
systemParams.Mobility.Lifetime_Secondary = systemParams.Mobility.Lifetime(ismember(systemParams.Mobility.Cons_Type,'S'));

systemParams.Mobility.Cons_Primary_S = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears),length(systemParams.Mobility.Cons_Names_Primary));
systemParams.Mobility.Cons_Secondary_S = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears),length(systemParams.Mobility.Cons_Names_Secondary));

%% Charge profiles generation
% Connection profiles in  MobilityInstalation
CarsDailyProf = [1 1 1 1 1 .2 .2 .2 .2 1 1 1 1 .2 .2 .2 .2 .2 .2 1 1 1 1 1];
BusDailyProf =  [1 1 1 1 1 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 1 1 1];
BaseProf = ones(1,24);

LDV_chProf = repmat(CarsDailyProf,length(systemParams.IndexNodes),365);
W23_chProf = repmat(CarsDailyProf,length(systemParams.IndexNodes),365);
BUS_chProf = repmat(BusDailyProf,length(systemParams.IndexNodes),365);
MDV_chProf = repmat(BusDailyProf,length(systemParams.IndexNodes),365);
HDV_chProf = repmat(BaseProf,length(systemParams.IndexNodes),365);

MP_chProf = repmat(BaseProf,length(systemParams.IndexNodes),365);
MF_chProf = repmat(BaseProf,length(systemParams.IndexNodes),365);
AP_chProf = repmat(BaseProf,length(systemParams.IndexNodes),365);
AF_chProf = repmat(BaseProf,length(systemParams.IndexNodes),365);

for i = 1:length(systemParams.IndexNodes)
    timeZone = systemParams.timeZone(i);
    if timeZone>=0

        systemParams.Mobility.LDV_chProf(i,:) = [LDV_chProf(i,1+timeZone:end) LDV_chProf(i,1:timeZone)];
        systemParams.Mobility.W23_chProf(i,:) = [W23_chProf(i,1+timeZone:end) W23_chProf(i,1:timeZone)];
        systemParams.Mobility.BUS_chProf(i,:) = [BUS_chProf(i,1+timeZone:end) BUS_chProf(i,1:timeZone)];
        systemParams.Mobility.MDV_chProf(i,:) = [MDV_chProf(i,1+timeZone:end) MDV_chProf(i,1:timeZone)];
        systemParams.Mobility.HDV_chProf(i,:) = [HDV_chProf(i,1+timeZone:end) HDV_chProf(i,1:timeZone)];
        systemParams.Mobility.MP_chProf(i,:) = [MP_chProf(i,1+timeZone:end) MP_chProf(i,1:timeZone)];
        systemParams.Mobility.MF_chProf(i,:) = [MF_chProf(i,1+timeZone:end) MF_chProf(i,1:timeZone)];
        systemParams.Mobility.AP_chProf(i,:) = [AP_chProf(i,1+timeZone:end) AP_chProf(i,1:timeZone)];
        systemParams.Mobility.AF_chProf(i,:) = [AF_chProf(i,1+timeZone:end) AF_chProf(i,1:timeZone)];

    else

        systemParams.Mobility.LDV_chProf(i,:) = [LDV_chProf(i,end+timeZone+1:end) LDV_chProf(i,1:end+timeZone)];
        systemParams.Mobility.W23_chProf(i,:) = [W23_chProf(i,end+timeZone+1:end) W23_chProf(i,1:end+timeZone)];
        systemParams.Mobility.BUS_chProf(i,:) = [BUS_chProf(i,end+timeZone+1:end) BUS_chProf(i,1:end+timeZone)];
        systemParams.Mobility.MDV_chProf(i,:) = [MDV_chProf(i,end+timeZone+1:end) MDV_chProf(i,1:end+timeZone)];
        systemParams.Mobility.HDV_chProf(i,:) = [HDV_chProf(i,end+timeZone+1:end) HDV_chProf(i,1:end+timeZone)];
        systemParams.Mobility.MP_chProf(i,:) = [MP_chProf(i,end+timeZone+1:end) MP_chProf(i,1:end+timeZone)];
        systemParams.Mobility.MF_chProf(i,:) = [MF_chProf(i,end+timeZone+1:end) MF_chProf(i,1:end+timeZone)];
        systemParams.Mobility.AP_chProf(i,:) = [AP_chProf(i,end+timeZone+1:end) AP_chProf(i,1:end+timeZone)];
        systemParams.Mobility.AF_chProf(i,:) = [AF_chProf(i,end+timeZone+1:end) AF_chProf(i,1:end+timeZone)];

    end



end