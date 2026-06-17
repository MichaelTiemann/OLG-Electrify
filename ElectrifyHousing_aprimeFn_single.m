function solarpv_prime=ElectrifyHousing_aprimeFn_single(buyhouse,solarpv,ypp,pv_delta)
% Because we use vfoptions.refine_d, the decision variables for aprimeFn must follow the ordering d2,d3
% Also, experience assets must be listed last in aprimeFn

solarpv_prime=single(-Inf);
single_0=single(0); single_1=single(1); single_5=single(5);

% Must ensure that hprime>h when buyhouse>0...

switch buyhouse
    case 0
        % Cannot install in house we don't own
        solarpv_prime=single_0;
    case 1
        % We start from scratch with a new house
        solarpv_prime=single_0;
    case 2
        % Buy new house, get a random amount of solar
        solarpv_prime = single(randi([single_1,single_5]));
    case 3
        % Keep house, experience the slow degradation of solarpv capacity
        solarpv_prime=solarpv*(single_1-pv_delta)^ypp;
    case 4
        % Keep house, install more solarpv if we can
        if solarpv<=4
            solarpv_prime = solarpv+single(randi([single_1,floor(single_5-solarpv)]));
            % Else -Inf return value will cause the ReturnFn to return -Inf
        end
end


end