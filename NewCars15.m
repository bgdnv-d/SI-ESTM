function NewCar = NewCars15(setup,CarNumber,IndexYears)
%
% FUNCTION NewCars15(setup, CarNumber, IndexYears)
%
% Estimates new car additions based on a 15-year vehicle lifetime.
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
for costYear = setup.startYear-10:setup.stepYear:setup.startYear
    NewCar(:,IndexYears==costYear) = (CarNumber(:,IndexYears==costYear)-CarNumber(:,IndexYears==costYear-5))+...
        (CarNumber(:,IndexYears==setup.startYear-15)*0.33);
end
for costYear = setup.startYear+5:setup.stepYear:setup.endYear
    NewCar(:,IndexYears==costYear) = NewCar(:,IndexYears==costYear-15)+(CarNumber(:,IndexYears==costYear)-CarNumber(:,IndexYears==(costYear-5)));
end

NewCar(NewCar<0)=0;