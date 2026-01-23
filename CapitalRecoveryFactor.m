function crf=capitalRecoveryFactor(WACC,lifetime)
%
% FUNCTION capitalRecoveryFactor(WACC, lifetime)
%
% Calculates the capital recovery factor (crf) for one or multiple asset lifetimes.
%
%
% INPUT:
%            WACC:      Interest rate used to calculate yearly payments.
%            lifetime:  How many years the investment lasts.
%
% OUTPUT:
%            crf:       Column of yearly payment factors for each lifetime.
%
%Dmitrii Bogdanov
%last change 22.07.2025


crf = 0;
for k=1:length(lifetime)
    crf(k)  = (WACC * (1 + WACC)^lifetime(k)) / ((1 + WACC)^lifetime(k) - 1);
end

crf = crf';





