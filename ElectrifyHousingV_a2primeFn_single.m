function solarpv_prime = ElectrifyHousingV_a2primeFn_single(installpv, solarpv)

% 1. Wrap constants in single() to prevent phantom double promotion
single_0_99 = single(0.99);

% We replace the random generator with a deterministic install size (e.g., 30kW)
% so the backward induction can mathematically converge.
single_install_size = single(30);

% 2. Pure Arithmetic Masking (No pre-allocations!)
install_mask    = (installpv == 1) & (solarpv == 0);
no_install_mask = (installpv == 0);

% The JIT compiler implicitly expands this natively inside the GPU registers
solarpv_prime = install_mask .* single_install_size + ...
    no_install_mask .* (solarpv .* single_0_99);

% (Note: The ReturnFn already bans states where installpv == 1 & solarpv > 0 with -Inf,
% so we don't need to explicitly handle their transitions here; they will naturally drop out).


end