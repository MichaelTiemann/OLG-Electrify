function solarpv_prime=ElectrifyHousing_aprimeFn(buyhouse,solarpv)
% Because we use vfoptions.refine_d, the decision variables for aprimeFn must follow the ordering d2,d3
% Also, experience assets must be listed last in aprimeFn

solarpv_prime=-Inf;

switch buyhouse
    case 0
        % Cannot install in house we don't own
        return
    case 1
        % We start from scratch with a new house
        solarpv=0;
    case 2
        % We get a random amount of solar
        solarpv_prime = randi([1,5]);
    case 3
        % The slow degradation of installed solar capacity
        solarpv_prime=solarpv * 0.99;
    case 4
        % If we can install more solar, do so now
        if solarpv<=4
            solarpv_prime = solarpv+randi([1,floor(5-solarpv)]);
        % Else -Inf return value will cause the ReturnFn to return -Inf
        end        
end


end