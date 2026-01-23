function systemParams = HourlyCOP(setup,systemParams)
%
% FUNCTION HourlyCOP(setup, systemParams)
%
% Calculates the hourly COP (Coefficient of Performance) for heat pump systems 
% based on ambient temperature profiles and assumptions on output heat temperature 
% and shares of different HP types.
%
% INPUT:
%            setup:         Structure that contains all necessary settings and data for processing.
%            systemParams:  Structure with system parameters that will be updated.
%
% OUTPUT:
%            systemParams:  Updated system parameters structure including hourly COP values.
%
%Dmitrii Bogdanov
%last change 23.07.2025


systemParams.HP.airHP_COP_IH = 6.08 - 0.0941*(setup.refTemp_IH-systemParams.TemperatureAir)+ 0.000464*(setup.refTemp_IH-systemParams.TemperatureAir).^2;
systemParams.HP.waterHP_COP_IH = 10.29 - 0.2084*(setup.refTemp_IH-systemParams.TemperatureWater)+ 0.001322*(setup.refTemp_IH-systemParams.TemperatureWater).^2;
systemParams.HP.groundHP_COP_IH = 9.99 - 0.2049*(setup.refTemp_IH-systemParams.TemperatureGround-5)+ 0.001249*(setup.refTemp_IH-systemParams.TemperatureGround-5).^2;

systemParams.HP.airHP_COP_DH = 6.08 - 0.0941*(setup.refTemp_DH-systemParams.TemperatureAir)+ 0.000464*(setup.refTemp_DH-systemParams.TemperatureAir).^2;
systemParams.HP.waterHP_COP_DH = 10.29 - 0.2084*(setup.refTemp_DH-systemParams.TemperatureWater)+ 0.001322*(setup.refTemp_DH-systemParams.TemperatureWater).^2;
systemParams.HP.groundHP_COP_DH = 9.99 - 0.2049*(setup.refTemp_DH-systemParams.TemperatureGround-5)+ 0.001249*(setup.refTemp_DH-systemParams.TemperatureGround-5).^2;

systemParams.HP.waterHP_COP_Ind = 10.29 - 0.2084*repmat(80,size(systemParams.TemperatureWater,1),size(systemParams.TemperatureWater,2))+ 0.001322*repmat(80,size(systemParams.TemperatureWater,1),size(systemParams.TemperatureWater,2)).^2;

for costYear = setup.startYear:setup.stepYear:setup.endYear

    systemParams.HP.COP_IH(:,:,systemParams.IndexYears==costYear) = 1./(repmat(setup.airHPshare_IH(:,systemParams.IndexYears==costYear)',setup.endHour,1)./systemParams.HP.airHP_COP_IH+repmat(setup.waterHPshare_IH(:,systemParams.IndexYears==costYear)',setup.endHour,1)./systemParams.HP.waterHP_COP_IH+repmat(setup.groundHPshare_IH(:,systemParams.IndexYears==costYear)',setup.endHour,1)./systemParams.HP.groundHP_COP_IH)*systemParams.Efficiency(ismember(systemParams.IndexID,'THHP'),systemParams.IndexYears==costYear);
    systemParams.HP.COP_DH(:,:,systemParams.IndexYears==costYear) = 1./(repmat(setup.airHPshare_DH(:,systemParams.IndexYears==costYear)',setup.endHour,1)./systemParams.HP.airHP_COP_DH+repmat(setup.waterHPshare_DH(:,systemParams.IndexYears==costYear)',setup.endHour,1)./systemParams.HP.waterHP_COP_DH+repmat(setup.groundHPshare_DH(:,systemParams.IndexYears==costYear)',setup.endHour,1)./systemParams.HP.groundHP_COP_DH)*systemParams.Efficiency(ismember(systemParams.IndexID,'TDHP'),systemParams.IndexYears==costYear);

    systemParams.HP.COP_IH_mean(:,:,systemParams.IndexYears==costYear) = repmat(mean(systemParams.HP.COP_IH(:,:,systemParams.IndexYears==costYear),1),8760,1);
    systemParams.HP.COP_DH_mean(:,:,systemParams.IndexYears==costYear) = repmat(mean(systemParams.HP.COP_IH(:,:,systemParams.IndexYears==costYear),1),8760,1);


end




qq = 1;