function a2prime=ElectrifyHousing_a2primeFn(installpv, solarpv, energy_pct_cost)
% Because we use vfoptions.refine_d, the decision variables for aprimeFn must follow the ordering d2,d3
% Also, experience assets must be listed first in aprimeFn

a2prime=0;
if installpv
    % Average New Zealander uses 20kW/day for house and 10kW/day for EV
    % 
    a2prime=(solarpv / 30)*energy_pct_cost;
end

end