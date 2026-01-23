function systemParams = MobilityInstalation(setup, systemParams)
%
% FUNCTION MobilityInstalation(setup, systemParams)
%
% Adds installation data for mobility technologies to the system parameters.
%
%
% INPUT:
%            setup:         Structure that contains all necessary settings and data for processing.
%            systemParams:  Structure with system parameters to be updated.
%
% OUTPUT:
%            systemParams:  Updated structure with mobility installation data.
%
%Dmitrii Bogdanov
%last change 23.07.2025


% bat capacity		1960	1965	1970	1975	1980	1985	1990	1995	2000	2005	2010	2015	2020	2025	2030	2035	2040	2045	2050
LDV_PHEV =      [repmat([ 0	0	0	0	0	0	0	10	10	10	10	10	10	9	8	8	7	7	6],length(systemParams.IndexNodes),1)];
LDV_BEV =       repmat([ 0	0	0	0	0	0	0	70	70	70	70	70	70	70	70	70	70	70	70],length(systemParams.IndexNodes),1);
W23_PHEV =      repmat([ 0	0	0	0	0	0	0	2	2	2	2	2	2	2	2	2	2	2	2],length(systemParams.IndexNodes),1);
W23_BEV =       repmat([ 0	0	0	0	0	0	0	10	10	10	10	10	10	10	10	10	10	10	10],length(systemParams.IndexNodes),1);
BUS_PHEV =      repmat([ 0	0	0	0	0	0	0	50	50	50	50	50	50	50	50	50	50	50	50],length(systemParams.IndexNodes),1);
BUS_BEV =       repmat([ 0	0	0	0	0	0	0	333	333	333	333	333	333	333	333	333	333	333	333],length(systemParams.IndexNodes),1);
MDV_PHEV =      repmat([ 0	0	0	0	0	0	0	17	17	17	17	17	17	17	17	17	17	17	17],length(systemParams.IndexNodes),1);
MDV_BEV =       repmat([ 0	0	0	0	0	0	0	120	120	120	120	120	120	120	120	120	120	120	120],length(systemParams.IndexNodes),1);
HDV_PHEV =      repmat([ 0	0	0	0	0	0	0	36	36	36	36	36	36	36	36	36	36	36	36],length(systemParams.IndexNodes),1);
HDV_BEV =       repmat([ 0	0	0	0	0	0	0	900	900	900	900	900	900	900	900	900	900	900	900],length(systemParams.IndexNodes),1);
if length(systemParams.IndexYears)>19
    LDV_PHEV = [LDV_PHEV repmat(LDV_PHEV(:,end),1,length(systemParams.IndexYears)-19)];
    LDV_BEV = [LDV_BEV repmat(LDV_BEV(:,end),1,length(systemParams.IndexYears)-19)];
    W23_PHEV = [W23_PHEV repmat(W23_PHEV(:,end),1,length(systemParams.IndexYears)-19)];
    W23_BEV = [W23_BEV repmat(W23_BEV(:,end),1,length(systemParams.IndexYears)-19)];
    BUS_PHEV = [BUS_PHEV repmat(BUS_PHEV(:,end),1,length(systemParams.IndexYears)-19)];
    BUS_BEV = [BUS_BEV repmat(BUS_BEV(:,end),1,length(systemParams.IndexYears)-19)];
    MDV_PHEV = [MDV_PHEV repmat(MDV_PHEV(:,end),1,length(systemParams.IndexYears)-19)];
    MDV_BEV = [MDV_BEV repmat(MDV_BEV(:,end),1,length(systemParams.IndexYears)-19)];
    HDV_PHEV = [HDV_PHEV repmat(HDV_PHEV(:,end),1,length(systemParams.IndexYears)-19)];
    HDV_BEV = [HDV_BEV repmat(HDV_BEV(:,end),1,length(systemParams.IndexYears)-19)];
end


%% Efficiencies
systemParams.Mobility.Cons_Primary_N = systemParams.Mobility.Cons_Values(ismember(systemParams.Mobility.Cons_Type,'P'),:);
systemParams.Mobility.Cons_Names_Primary = systemParams.Mobility.Cons_Names(ismember(systemParams.Mobility.Cons_Type,'P'));
systemParams.Mobility.Lifetime_Primary = systemParams.Mobility.Lifetime(ismember(systemParams.Mobility.Cons_Type,'P'));
systemParams.Mobility.Cons_Secondary_N = systemParams.Mobility.Cons_Values(ismember(systemParams.Mobility.Cons_Type,'S'),:);
systemParams.Mobility.Cons_Names_Secondary = systemParams.Mobility.Cons_Names(ismember(systemParams.Mobility.Cons_Type,'S'));
systemParams.Mobility.Lifetime_Secondary = systemParams.Mobility.Lifetime(ismember(systemParams.Mobility.Cons_Type,'S'));

systemParams.Mobility.Cons_Primary_S = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears),length(systemParams.Mobility.Cons_Names_Primary));
systemParams.Mobility.Cons_Secondary_S = zeros(length(systemParams.IndexNodes),length(systemParams.IndexYears),length(systemParams.Mobility.Cons_Names_Secondary));



%% Road

%LDV
LDVNumber = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMLL'))./systemParams.Mobility.LDV_Pass./systemParams.Mobility.LDV_km;
% LDVNumber = [0	0	0	0	0	0	0	0	300	435	654	867	1205	1500	1700	1900	2050	2170	2250]
NewLDV = NewCars15(setup,LDVNumber,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRLI')) = NewLDV .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRLI'));
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRLB')) = NewLDV .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRLB'));
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRLF')) = NewLDV .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRLF'));
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRLP')) = NewLDV .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRLP'));

systemParams.Instalations(:,:,ismember(systemParams.IndexID,'SLDV')) = (LDV_BEV.*systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRLB'))+LDV_PHEV.*systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRLP')));%setup.DumpChargeShare.*
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'ILDV')) = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'SLDV'))/systemParams.EnPoRatioStorage(ismember(systemParams.IndexIDS,'SLDV'),1);

%2,3W
W23Number = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMLW'))./systemParams.Mobility.W23_Pass./systemParams.Mobility.W23_km;

NewW23 = NewCars10(setup,W23Number,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRWI')) = NewW23 .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRWI'));
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRWB')) = NewW23 .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRWB'));
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRWF')) = NewW23 .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRWF'));
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRWP')) = NewW23 .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRWP'));

systemParams.Instalations(:,:,ismember(systemParams.IndexID,'SW23')) = (W23_BEV.*systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRWB'))+W23_PHEV.*systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRWP')));%setup.DumpChargeShare.*
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'IW23')) = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'SW23'))/systemParams.EnPoRatioStorage(ismember(systemParams.IndexIDS,'SW23'),1);

%BUS
BUSNumber = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMLB'))./systemParams.Mobility.BUS_Pass./systemParams.Mobility.BUS_km;

NewBUS = NewCars15(setup,BUSNumber,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRBI')) = NewBUS .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRBI'));
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRBB')) = NewBUS .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRBB'));
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRBF')) = NewBUS .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRBF'));
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRBP')) = NewBUS .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRBP'));

systemParams.Instalations(:,:,ismember(systemParams.IndexID,'SBUS')) = (BUS_BEV.*systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRBB'))+BUS_PHEV.*systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRBP')));%setup.DumpChargeShare.*
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'IBUS')) = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'SBUS'))/systemParams.EnPoRatioStorage(ismember(systemParams.IndexIDS,'SBUS'),1);

%MDV
MDVNumber = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMLM'))./systemParams.Mobility.MDV_Tonne./systemParams.Mobility.MDV_km;

NewMDV = NewCars15(setup,MDVNumber,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRMI')) = NewMDV .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRMI'));
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRMB')) = NewMDV .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRMB'));
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRMF')) = NewMDV .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRMF'));
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRMP')) = NewMDV .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRMP'));

systemParams.Instalations(:,:,ismember(systemParams.IndexID,'SMDV')) = (MDV_BEV.*systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRMB'))+MDV_PHEV.*systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRMP')));%setup.DumpChargeShare.*
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'IMDV')) = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'SMDV'))/systemParams.EnPoRatioStorage(ismember(systemParams.IndexIDS,'SMDV'),1);

%HDV
HDVNumber = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMLH'))./systemParams.Mobility.HDV_Tonne./systemParams.Mobility.HDV_km;

NewHDV = NewCars15(setup,HDVNumber,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRHI')) = NewHDV .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRHI'));
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRHB')) = NewHDV .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRHB'));
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRHF')) = NewHDV .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRHF'));
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRHP')) = NewHDV .*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRHP'));

systemParams.Instalations(:,:,ismember(systemParams.IndexID,'SHDV')) = (HDV_BEV.*systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRHB'))+HDV_PHEV.*systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRHP')));%setup.DumpChargeShare.*
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'IHDV')) = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'SHDV'))/systemParams.EnPoRatioStorage(ismember(systemParams.IndexIDS,'SHDV'),1);

%% Others 30years lifetime
%% Rail
%Rail Pass
MRPF = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMRP')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRPF'));
NewMRPF = NewCars30(setup,MRPF,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRPF')) = NewMRPF;

MRPE = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMRP')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRPE'));
NewMRPE = NewCars30(setup,MRPE,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRPE')) = NewMRPE;
%Rail Freight
MRFF = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMRF')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRFF'));
NewMRFF = NewCars30(setup,MRFF,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRFF')) = NewMRFF;

MRFE = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMRF')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MRFE'));
NewMRFE = NewCars30(setup,MRFE,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MRFE')) = NewMRFE;

%% Marine
%Marine Pass
MMPF = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMMP')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MMPF'));
NewMMPF = NewCars30(setup,MMPF,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MMPF')) = NewMMPF;

MMPE = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMMP')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MMPE'));
NewMMPE = NewCars30(setup,MMPE,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MMPE')) = NewMMPE;

MMPH = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMMP')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MMPH'));
NewMMPH = NewCars30(setup,MMPH,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MMPH')) = NewMMPH;

MMPG = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMMP')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MMPG'));
NewMMPG = NewCars30(setup,MMPG,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MMPG')) = NewMMPG;

MMPA = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMMP')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MMPA'));
NewMMPA = NewCars30(setup,MMPA,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MMPA')) = NewMMPA;

MMPM = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMMP')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MMPM'));
NewMMPM = NewCars30(setup,MMPM,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MMPM')) = NewMMPM;

%Marine Freight
MMFF = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMMF')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MMFF'));
NewMMFF = NewCars30(setup,MMFF,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MMFF')) = NewMMFF;

MMFE = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMMF')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MMFE'));
NewMMFE = NewCars30(setup,MMFE,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MMFE')) = NewMMFE;

MMFH = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMMF')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MMFH'));
NewMMFH = NewCars30(setup,MMFH,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MMFH')) = NewMMFH;

MMFG = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMMF')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MMFG'));
NewMMFG = NewCars30(setup,MMFG,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MMFG')) = NewMMFG;

MMFA = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMMF')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MMFA'));
NewMMFA = NewCars30(setup,MMFA,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MMFA')) = NewMMFA;

MMFM = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMMF')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MMFM'));
NewMMFM = NewCars30(setup,MMFM,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MMFM')) = NewMMFM;

%% Avia
%Avia Pass
MAPF = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMAP')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MAPF'));
NewMAPF = NewCars30(setup,MAPF,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MAPF')) = NewMAPF;

MAPE = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMAP')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MAPE'));
NewMAPE = NewCars30(setup,MAPE,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MAPE')) = NewMAPE;

MAPH = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMAP')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MAPH'));
NewMAPH = NewCars30(setup,MAPH,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MAPH')) = NewMAPH;

%Avia Freight
MAFF = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMAF')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MAFF'));
NewMAFF = NewCars30(setup,MAFF,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MAFF')) = NewMAFF;

MAFE = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMAF')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MAFE'));
NewMAFE = NewCars30(setup,MAFE,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MAFE')) = NewMAFE;

MAFH = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMAF')).*systemParams.Mobility.Shares(:,:,ismember(systemParams.IndexID,'MAFH'));
NewMAFH = NewCars30(setup,MAFH,systemParams.IndexYears);
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'MAFH')) = NewMAFH;



for costYear = setup.startYear:setup.stepYear:setup.endYear
    Y_ind = find(systemParams.IndexYears == costYear);
    % LDV
    techM = systemParams.IndexIDM(1:4);
    totMob = [];
    for i = 1:length(techM)
        totMob(:,i) = sum(systemParams.Instalations(:,Y_ind-systemParams.Mobility.Lifetime_Primary(ismember(systemParams.IndexIDM,techM(i)))/setup.stepYear+1:Y_ind,ismember(systemParams.IndexID,techM(i))),2);
    end
    for i = 1:length(techM)
        systemParams.Mobility.SharesTot(:,Y_ind,ismember(systemParams.IndexID,techM(i))) = totMob(:,i)./sum(totMob,2);
    end

    % 23W
    techM = systemParams.IndexIDM(5:8);
    totMob = [];
    for i = 1:length(techM)
        totMob(:,i) = sum(systemParams.Instalations(:,Y_ind-systemParams.Mobility.Lifetime_Primary(ismember(systemParams.IndexIDM,techM(i)))/setup.stepYear+1:Y_ind,ismember(systemParams.IndexID,techM(i))),2);
    end
    for i = 1:length(techM)
        systemParams.Mobility.SharesTot(:,Y_ind,ismember(systemParams.IndexID,techM(i))) = totMob(:,i)./sum(totMob,2);
    end

    %Bus
    techM = systemParams.IndexIDM(9:12);
    totMob = [];
    for i = 1:length(techM)
        totMob(:,i) = sum(systemParams.Instalations(:,Y_ind-systemParams.Mobility.Lifetime_Primary(ismember(systemParams.IndexIDM,techM(i)))/setup.stepYear+1:Y_ind,ismember(systemParams.IndexID,techM(i))),2);
    end
    for i = 1:length(techM)
        systemParams.Mobility.SharesTot(:,Y_ind,ismember(systemParams.IndexID,techM(i))) = totMob(:,i)./sum(totMob,2);
    end

    %MDV
    techM = systemParams.IndexIDM(13:16);
    totMob = [];
    for i = 1:length(techM)
        totMob(:,i) = sum(systemParams.Instalations(:,Y_ind-systemParams.Mobility.Lifetime_Primary(ismember(systemParams.IndexIDM,techM(i)))/setup.stepYear+1:Y_ind,ismember(systemParams.IndexID,techM(i))),2);
    end
    for i = 1:length(techM)
        systemParams.Mobility.SharesTot(:,Y_ind,ismember(systemParams.IndexID,techM(i))) = totMob(:,i)./sum(totMob,2);
    end
    %HDV
    techM = systemParams.IndexIDM(17:20);
    totMob = [];
    for i = 1:length(techM)
        totMob(:,i) = sum(systemParams.Instalations(:,Y_ind-systemParams.Mobility.Lifetime_Primary(ismember(systemParams.IndexIDM,techM(i)))/setup.stepYear+1:Y_ind,ismember(systemParams.IndexID,techM(i))),2);
    end
    for i = 1:length(techM)
        systemParams.Mobility.SharesTot(:,Y_ind,ismember(systemParams.IndexID,techM(i))) = totMob(:,i)./sum(totMob,2);
    end

    %train pass
    techM = systemParams.IndexIDM(21:22);
    totMob = [];
    for i = 1:length(techM)
        totMob(:,i) = sum(systemParams.Instalations(:,Y_ind-systemParams.Mobility.Lifetime_Primary(ismember(systemParams.IndexIDM,techM(i)))/setup.stepYear+1:Y_ind,ismember(systemParams.IndexID,techM(i))),2);
    end
    for i = 1:length(techM)
        systemParams.Mobility.SharesTot(:,Y_ind,ismember(systemParams.IndexID,techM(i))) = totMob(:,i)./sum(totMob,2);
    end
    %train freight
    techM = systemParams.IndexIDM(23:24);
    totMob = [];
    for i = 1:length(techM)
        totMob(:,i) = sum(systemParams.Instalations(:,Y_ind-systemParams.Mobility.Lifetime_Primary(ismember(systemParams.IndexIDM,techM(i)))/setup.stepYear+1:Y_ind,ismember(systemParams.IndexID,techM(i))),2);
    end
    for i = 1:length(techM)
        systemParams.Mobility.SharesTot(:,Y_ind,ismember(systemParams.IndexID,techM(i))) = totMob(:,i)./sum(totMob,2);
    end

    %marine pass
    techM = systemParams.IndexIDM(25:30);
    totMob = [];
    for i = 1:length(techM)
        totMob(:,i) = sum(systemParams.Instalations(:,Y_ind-systemParams.Mobility.Lifetime_Primary(ismember(systemParams.IndexIDM,techM(i)))/setup.stepYear+1:Y_ind,ismember(systemParams.IndexID,techM(i))),2);
    end
    for i = 1:length(techM)
        systemParams.Mobility.SharesTot(:,Y_ind,ismember(systemParams.IndexID,techM(i))) = totMob(:,i)./sum(totMob,2);
    end
    %marine freight
    techM = systemParams.IndexIDM(31:36);
    totMob = [];
    for i = 1:length(techM)
        totMob(:,i) = sum(systemParams.Instalations(:,Y_ind-systemParams.Mobility.Lifetime_Primary(ismember(systemParams.IndexIDM,techM(i)))/setup.stepYear+1:Y_ind,ismember(systemParams.IndexID,techM(i))),2);
    end
    for i = 1:length(techM)
        systemParams.Mobility.SharesTot(:,Y_ind,ismember(systemParams.IndexID,techM(i))) = totMob(:,i)./sum(totMob,2);
    end
    %avia pass
    techM = systemParams.IndexIDM(37:39);
    totMob = [];
    for i = 1:length(techM)
        totMob(:,i) = sum(systemParams.Instalations(:,Y_ind-systemParams.Mobility.Lifetime_Primary(ismember(systemParams.IndexIDM,techM(i)))/setup.stepYear+1:Y_ind,ismember(systemParams.IndexID,techM(i))),2);
    end
    for i = 1:length(techM)
        systemParams.Mobility.SharesTot(:,Y_ind,ismember(systemParams.IndexID,techM(i))) = totMob(:,i)./sum(totMob,2);
    end
    %avia freight
    techM = systemParams.IndexIDM(40:42);
    totMob = [];
    for i = 1:length(techM)
        totMob(:,i) = sum(systemParams.Instalations(:,Y_ind-systemParams.Mobility.Lifetime_Primary(ismember(systemParams.IndexIDM,techM(i)))/setup.stepYear+1:Y_ind,ismember(systemParams.IndexID,techM(i))),2);
    end
    for i = 1:length(techM)
        systemParams.Mobility.SharesTot(:,Y_ind,ismember(systemParams.IndexID,techM(i))) = totMob(:,i)./sum(totMob,2);
    end
end


systemParams.Mobility.SharesTot(isnan(systemParams.Mobility.SharesTot))=0;


systemParams.Instalations(isnan(systemParams.Instalations))=0;

for costYear = setup.startYear:setup.stepYear:setup.endYear

    for i = 1:length(systemParams.Mobility.Cons_Names_Primary)
        for reg = 1:length(systemParams.IndexNodes)

            active = (systemParams.IndexYears<=costYear)&(systemParams.IndexYears>costYear-systemParams.Mobility.Lifetime_Primary(i));

            systemParams.Mobility.Cons_Primary_S(reg,systemParams.IndexYears == costYear,i)=...
                sum(max(systemParams.Instalations(reg,active,ismember(systemParams.IndexID,systemParams.Mobility.Cons_Names_Primary(i))),0).*repmat(systemParams.Mobility.Cons_Primary_N(i,active),1,1))./...
                sum(max(systemParams.Instalations(reg,active,ismember(systemParams.IndexID,systemParams.Mobility.Cons_Names_Primary(i))),0),2);
        end
    end

    for i = 1:length(systemParams.Mobility.Cons_Names_Secondary)
        for reg = 1:length(systemParams.IndexNodes)
            active = (systemParams.IndexYears<=costYear)&(systemParams.IndexYears>costYear-systemParams.Mobility.Lifetime_Secondary(i));
            systemParams.Mobility.Cons_Secondary_S(reg,systemParams.IndexYears == costYear,i)=...
                sum(systemParams.Instalations(reg,active,ismember(systemParams.IndexID,systemParams.Mobility.Cons_Names_Secondary(i))).*repmat(systemParams.Mobility.Cons_Secondary_N(i,active),1,1))./...
                sum(systemParams.Instalations(reg,active,ismember(systemParams.IndexID,systemParams.Mobility.Cons_Names_Secondary(i))),2);
        end

    end
end
systemParams.Mobility.Cons_Primary_S(isnan(systemParams.Mobility.Cons_Primary_S))=0;
systemParams.Mobility.Cons_Secondary_S(isnan(systemParams.Mobility.Cons_Secondary_S))=0;
%% Update road instalations (p*km to km)

systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMLL')) = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMLL'))./systemParams.Mobility.LDV_Pass; % in km
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMLW')) = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMLW'))./systemParams.Mobility.W23_Pass; % in km
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMLB')) = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMLB'))./systemParams.Mobility.BUS_Pass; % in km
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMLM')) = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMLM'))./systemParams.Mobility.MDV_Tonne; % in km
systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMLH')) = systemParams.Instalations(:,:,ismember(systemParams.IndexID,'LMLH'))./systemParams.Mobility.HDV_Tonne; % in km
systemParams.Instalations(isnan(systemParams.Instalations)) = 0;

%% Charge profiles generation

CarsDailyProf = [1 1 1 1 1 .2 .2 .2 .2 .8 .85 .85 .8 .2 .2 .2 .2 .2 .2 1 1 1 1 1]-0.05; %CarsDailyProf = CarsDailyProf/sum(CarsDailyProf);
BusDailyProf =  [1 1 1 1 1 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 .2 1 1 1]-0.05;

CarsDailyProf_Charge = [0.12	0.06	0.03	0.02	0.02	0.09	0.28	0.57	0.78	0.62	0.56	0.62	0.68	0.68	0.67	0.72	0.83	0.97	1.00	0.86	0.67	0.52	0.42	0.24];
BusDailyProf_Charge =  [0.12	0.06	0.03	0.02	0.02	0.09	0.28	0.57	0.78	0.62	0.56	0.62	0.68	0.68	0.67	0.72	0.83	0.97	1.00	0.86	0.67	0.52	0.42	0.24];

BaseProf = ones(1,24)*0.5;

LDV_connProf = repmat(CarsDailyProf,length(systemParams.IndexNodes),365);
W23_connProf = repmat(CarsDailyProf,length(systemParams.IndexNodes),365);
BUS_connProf = repmat(BusDailyProf,length(systemParams.IndexNodes),365);
MDV_connProf = repmat(BusDailyProf,length(systemParams.IndexNodes),365);
HDV_connProf = repmat(BaseProf,length(systemParams.IndexNodes),365);

LDV_chProf = repmat(CarsDailyProf_Charge,length(systemParams.IndexNodes),365);
W23_chProf = repmat(CarsDailyProf_Charge,length(systemParams.IndexNodes),365);
BUS_chProf = repmat(BusDailyProf_Charge,length(systemParams.IndexNodes),365);
MDV_chProf = repmat(BusDailyProf_Charge,length(systemParams.IndexNodes),365);
HDV_chProf = repmat(BaseProf,length(systemParams.IndexNodes),365);

MP_chProf = repmat(BaseProf,length(systemParams.IndexNodes),365);
MF_chProf = repmat(BaseProf,length(systemParams.IndexNodes),365);
AP_chProf = repmat(BaseProf,length(systemParams.IndexNodes),365);
AF_chProf = repmat(BaseProf,length(systemParams.IndexNodes),365);

for i = 1:length(systemParams.IndexNodes)
    timeZone = ceil(systemParams.timeZone(i));

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

        systemParams.Mobility.LDV_connProf(i,:) = [LDV_connProf(i,1+timeZone:end) LDV_connProf(i,1:timeZone)];
        systemParams.Mobility.W23_connProf(i,:) = [W23_connProf(i,1+timeZone:end) W23_connProf(i,1:timeZone)];
        systemParams.Mobility.BUS_connProf(i,:) = [BUS_connProf(i,1+timeZone:end) BUS_connProf(i,1:timeZone)];
        systemParams.Mobility.MDV_connProf(i,:) = [MDV_connProf(i,1+timeZone:end) MDV_connProf(i,1:timeZone)];
        systemParams.Mobility.HDV_connProf(i,:) = [HDV_connProf(i,1+timeZone:end) HDV_connProf(i,1:timeZone)];


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

        systemParams.Mobility.LDV_connProf(i,:) = [LDV_connProf(i,end+timeZone+1:end) LDV_connProf(i,1:end+timeZone)];
        systemParams.Mobility.W23_connProf(i,:) = [W23_connProf(i,end+timeZone+1:end) W23_connProf(i,1:end+timeZone)];
        systemParams.Mobility.BUS_connProf(i,:) = [BUS_connProf(i,end+timeZone+1:end) BUS_connProf(i,1:end+timeZone)];
        systemParams.Mobility.MDV_connProf(i,:) = [MDV_connProf(i,end+timeZone+1:end) MDV_connProf(i,1:end+timeZone)];
        systemParams.Mobility.HDV_connProf(i,:) = [HDV_connProf(i,end+timeZone+1:end) HDV_connProf(i,1:end+timeZone)];


    end

    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMLL')),:,1,i)= 1-systemParams.Mobility.LDV_chProf(i,:)';
    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMLW')),:,1,i)= 1-systemParams.Mobility.W23_chProf(i,:)';
    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMLB')),:,1,i)= 1-systemParams.Mobility.BUS_chProf(i,:)';
    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMLM')),:,1,i)= 1-systemParams.Mobility.MDV_chProf(i,:)';
    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMLH')),:,1,i)= 1-systemParams.Mobility.HDV_chProf(i,:)';
    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMMP')),:,1,i)= 1-systemParams.Mobility.MP_chProf(i,:)';
    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMMF')),:,1,i)= 1-systemParams.Mobility.MF_chProf(i,:)';
    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMAP')),:,1,i)= 1-systemParams.Mobility.AP_chProf(i,:)';
    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMAF')),:,1,i)= 1-systemParams.Mobility.AF_chProf(i,:)';

    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMLL')),:,2,i)= 1-systemParams.Mobility.LDV_chProf(i,:)';
    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMLW')),:,2,i)= 1-systemParams.Mobility.W23_chProf(i,:)';
    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMLB')),:,2,i)= 1-systemParams.Mobility.BUS_chProf(i,:)';
    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMLM')),:,2,i)= 1-systemParams.Mobility.MDV_chProf(i,:)';
    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMLH')),:,2,i)= 1-systemParams.Mobility.HDV_chProf(i,:)';
    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMMP')),:,2,i)= 1-systemParams.Mobility.MP_chProf(i,:)';
    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMMF')),:,2,i)= 1-systemParams.Mobility.MF_chProf(i,:)';
    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMAP')),:,2,i)= 1-systemParams.Mobility.AP_chProf(i,:)';
    systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LMAF')),:,2,i)= 1-systemParams.Mobility.AF_chProf(i,:)';


end

