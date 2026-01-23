function NewCar = NewCars30(setup,CarNumber,IndexYears)
%
% FUNCTION NewCars30(setup, CarNumber, IndexYears)
%
% Estimates new car additions based on a 30-year vehicle lifetime.
%
%
% INPUT:
%            setup:       Structure that contains all necessary settings and data for processing.
%            CarNumber:   Total number of cars to consider in the calculation.
%            IndexYears:  List of years to calculate new car additions for.
%
% OUTPUT:
%            NewCar:      Vector with number of new cars for each year.
%
%Dmitrii Bogdanov
%last change 23.07.2025


NewCar = zeros(size(CarNumber));
CarNumber(:,IndexYears==1995) = 0.85*CarNumber(:,IndexYears==2000);
CarNumber(:,IndexYears==1990) = 0.70*CarNumber(:,IndexYears==2000);
CarNumber(:,IndexYears==1985) = 0.5*CarNumber(:,IndexYears==2000);
NewCar(:,IndexYears<2005) = 0;

for costYear = setup.startYear-25:setup.stepYear:setup.startYear
    NewCar(:,IndexYears==costYear) = (CarNumber(:,IndexYears==costYear)-CarNumber(:,IndexYears==costYear-5))+...
        (CarNumber(:,IndexYears==setup.startYear-30)*5/30);
end
for costYear = setup.startYear+5:setup.stepYear:setup.endYear
    NewCar(:,IndexYears==costYear) = NewCar(:,IndexYears==costYear-30)+(CarNumber(:,IndexYears==costYear)-CarNumber(:,IndexYears==(costYear-5)));
end

