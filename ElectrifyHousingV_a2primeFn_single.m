function solarpv_prime = ElectrifyHousingV_a2primeFn_single(installpv, solarpv)

% Because we use vfoptions.refine_d, the decision variables for aprimeFn must follow the ordering d2,d3
% Also, experience assets must be listed first in aprimeFn

% 1. Create a master sizing tensor to ensure implicit expansion works safely
master_tensor = installpv + solarpv;
solarpv_prime = -inf(size(master_tensor), 'like', master_tensor);

% 2. Condition: If installing solar for the first time
install_mask = (installpv == 1) & (solarpv == 0);

% (Safely generate a GPU-compatible matrix of random values)
rand_sizes = cast(10 * randi([1, 4], size(master_tensor)), 'like', master_tensor);
solarpv_prime(install_mask) = rand_sizes(install_mask);

% 3. Condition: The slow degradation of installed solar capacity
no_install_mask = (installpv == 0);
solarpv_prime(no_install_mask) = solarpv(no_install_mask) * 0.99;


end