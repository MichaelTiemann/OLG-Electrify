function solarpv_prime=ElectrifyHousing_aprimeFn(installpv,solarpv)
% Because we use vfoptions.refine_d, the decision variables for aprimeFn must follow the ordering d2,d3
% Also, experience assets must be listed first in aprimeFn

solarpv_prime=-Inf;

if installpv
    % If we can install more solar, do so now
    if solarpv<=30
        solarpv_prime = 10 * randi([1,floor((40-solarpv)/10)]);
    end
    % Else -Inf return value will cause the ReturnFn to return -Inf
else
    % The slow degradation of installed solar capacity
    solarpv_prime=solarpv * 0.99;
end


end