function solarpv_prime=ElectrifyHousing_a2primeFn(installpv, solarpv)
% Because we use vfoptions.refine_d, the decision variables for aprimeFn must follow the ordering d2,d3
% Also, experience assets must be listed first in aprimeFn

solarpv_prime=solarpv*0.6;
if installpv
    % Average New Zealander uses 20kW/day for house and 10kW/day for EV
    % 
    if solarpv < 40
        solarpv_prime = randi([solarpv+1, min(solarpv+20,40)]);
    end
end


end