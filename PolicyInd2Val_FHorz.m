function PolicyValues = PolicyInd2Val_FHorz(Policy,n_d,n_a,n_z,N_j,d_grid,a_grid,simoptions,varargin)
% VFIToolkit Shadow Override for SemiZ Dimension Bugs
% Safely translates indices to values and reshapes to EXACTLY 4D
% [NumPolicies, N_a_flat, N_z_flat, N_j] to satisfy LifeCycleProfiles

Policy_cpu = gather(Policy);
sz = size(Policy_cpu);

% 1. Dynamically Flatten Trailing Dimensions for mapping
num_pol = sz(1);
Pol_flat = reshape(Policy_cpu, [num_pol, prod(sz(2:end))]);
Pol_flat = max(1, round(Pol_flat)); % Clamp indices safely

% 2. Map Grids
PolVal_flat = zeros(size(Pol_flat), 'like', a_grid);

l_d = length(n_d);
l_aprime = num_pol - l_d; % Dynamically detect how many standard assets were chosen

% Map Decisions
cum_d = 0;
for i = 1:l_d
    grid_chunk = d_grid(cum_d + 1 : cum_d + n_d(i));
    PolVal_flat(i, :) = grid_chunk(Pol_flat(i, :));
    cum_d = cum_d + n_d(i);
end

% Map Assets (ONLY standard assets stored in the Policy tensor)
cum_a = 0;
for i = 1:l_aprime
    row_idx = l_d + i;
    grid_chunk = a_grid(cum_a + 1 : cum_a + n_a(i));
    PolVal_flat(row_idx, :) = grid_chunk(Pol_flat(row_idx, :));
    cum_a = cum_a + n_a(i);
end

% 3. Reshape to EXACTLY 4 Dimensions for LifeCycleProfiles
N_a_flat = prod(n_a);
N_j_flat = N_j;

% Dynamically calculate N_z_flat based on the total elements to guarantee no mismatch
N_z_flat = prod(sz(2:end)) / (N_a_flat * N_j_flat);

PolicyValues = reshape(PolVal_flat, [num_pol, N_a_flat, N_z_flat, N_j_flat]);


end