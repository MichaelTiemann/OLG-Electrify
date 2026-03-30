function [d_grid,a_grid,z_grid,pi_z,jequaloneDist,share_grid,k_grid,pv_grid_firm,Params,vfoptions,simoptions]=Electrify_GridSetup(scenario, ypp, n_d, n_a, n_z, small_z_no_e, Params, vfoptions, simoptions)

%% Grids for household

% Grid for labour choice
labor_grid=linspace(0,1,n_d.household(1))'; % Notice that it is imposing the 0<=h<=1 condition implicitly

% Grid for share holdings, always > 0
% For later scenarios, shrink the grid for more accuracy
% s_grid_cubed=linspace(0,1,ceil(n_a.household(1)/3)).^3; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.
% s_grid_linear=linspace(1,10,floor(n_a.household(1)*2/3)+1);
% share_grid=[s_grid_cubed, s_grid_linear(2:end)]';

% Set up d for VFI Toolkit
if scenario<3
    % One decision variable: labor hours percentage
    d_grid.household=labor_grid;
    % Grid for share holdings, always > 0
    s_grid_cubed=linspace(0,1,ceil(n_a.household(1)/3)).^3; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.
    s_grid_linear=linspace(1,16,floor(n_a.household(1)*2/3)+1);
    share_grid=[s_grid_cubed, s_grid_linear(2:end)]';
    a_grid.household=share_grid;
    Params.minhouse=1;
    pv_grid_hh=NaN;

    % This is a default, but we set explicitly to make this reentrant
    vfoptions.experienceasset.household=0;
    simoptions.experienceasset.household=0;
else
    % Grid for share holdings, always > 0; Small max due to other assets
    share_grid=8*linspace(0,1,n_a.household(1))';

    % Grid for bank account; a negative balance implies a mortgage
    a_grid_cubed=linspace(-1,0,ceil(n_a.household(2)/2)-1).^3;
    a_grid_linear=linspace(0,16,floor(n_a.household(2)/2)+2);
    asset_grid=[a_grid_cubed, a_grid_linear(2:end)]';
    
    % Make it so that there is a zero assets
    % Find closest to zero assets
    [~,zeroassetindex]=min(abs(asset_grid));
    asset_grid(zeroassetindex)=0;

    % PV grid is 5 kW per grid element (approx 20kWh/day)
    if scenario<4
        car_grid=zeros(0);
        house_grid=(0:1:n_a.household(3)-1)';
        pv_grid_hh=(0:1:n_a.household(4)-1)';
    else
        car_grid=(0:1:n_a.household(3)-1)'; % car assets: no car; petrol car; EV car
        house_grid=(0:1:n_a.household(4)-1)';
        pv_grid_hh=(0:1:n_a.household(5)-1)';
    end

    Params.minhouse=house_grid(2); % first is zero (no house)
    
    % buyhouse decisions
    %  0=no house
    %  1=buy house w/o pv this period
    %  2=keep house; no pv upgrade
    %  3=buy house w/ pv this period
    %  4=keep house; pv upgrade (if possible)
    %  5=testing (not used)
    buyhouse_grid=(0:1:n_d.household(2)-1)';
    
    d_grid.household=[labor_grid; buyhouse_grid];
    a_grid.household=[share_grid; asset_grid; car_grid; house_grid; pv_grid_hh];
    
    %% Solar PV is an experience asset
    vfoptions.experienceasset.household=1;
    simoptions.experienceasset.household=1;
    
    %% Define aprime function used for the experience asset
    
    % experienceasset: aprime_val=aprimeFn(d,a)
    % vfoptions.refine_d: the decision variables input to aprimeFn are d3
    aprimeFn=@(buyhouse, solarpv, ypp) ElectrifyHousing_aprimeFn(buyhouse, solarpv, ypp); % Will return the value of aprime (solarpv)
    
    %% Put the experience asset into vfoptions and simoptions
    vfoptions.aprimeFn.household=aprimeFn;
    % vfoptions.n_u=n_u;
    % vfoptions.u_grid=u_grid;
    % vfoptions.pi_u=pi_u;
    simoptions.aprimeFn.household=aprimeFn;
    % simoptions.n_u=n_u;
    % simoptions.u_grid=u_grid;
    % simoptions.pi_u=pi_u;
    % Because a_grid and d_grid are involved in experience assets, but are not
    % normally needed for agent distriubiton simulation, we have to also
    % include these in simoptions
    simoptions.a_grid.household=a_grid.household;
    simoptions.d_grid.household=d_grid.household;
    
    %% To speed up the use of experienceasset we use 'refine_d', which requires us to set the decision variables in a specific order
    if scenario<4
        vfoptions.refine_d.household=[1,0,1]; % tell the code how many d1, d2, and d3 there are
    else
        vfoptions.refine_d.household=[2,0,1]; % tell the code how many d1, d2, and d3 there are
    end
    % Idea is to distinguish three categories of decision variable:
    %  d1: decision is in the ReturnFn but not in aprimeFn
    %  d2: decision is in the aprimeFn but not in ReturnFn
    %  d3: decision is in both ReturnFn and in aprimeFn
    % Note: ReturnFn must use inputs (d1,d3,..) 
    %       aprimeFn must use inputs (d2,d3,..)
    % n_d must be set up as n_d=[n_d1, n_d2, n_d3]
    % d_grid must be set up as d_grid=[d1_grid; d2_grid; d3_grid];
    simoptions.refine_d=vfoptions.refine_d;
end
if small_z_no_e
    z_grid_J=zeros(n_z.household,Params.J);
    pi_z_J=ones(n_z.household,n_z.household,Params.J);
else
    % First, z, the AR(1) with age-dependent parameters
    [z_grid_J, pi_z_J] = discretizeLifeCycleAR1_FellaGallipoliPan(Params.rho_z,Params.sigma_epsilon_z,n_z.household,Params.J);
    % z_grid_J is n_z-by-J, so z_grid_J(:,j) is the grid for age j
    % pi_z_J is n_z-by-n_z-by-J, so pi_z_J(:,:,j) is the transition matrix for age j

    % Second, e, the iid normal with age-dependent parameters
    [e_grid_J, pi_e_J] = discretizeLifeCycleAR1_FellaGallipoliPan(zeros(1,Params.J),Params.sigma_e,vfoptions.n_e.household,Params.J); % Note: AR(1) with rho=0 is iid normal
    % Because e is iid we actually just use
    pi_e_J=shiftdim(pi_e_J(1,:,:),1);

    % Any (iid) e variable always has to go into vfoptions and simoptions
    vfoptions.e_grid.household=e_grid_J;
    vfoptions.pi_e.household=pi_e_J;
    simoptions.n_e.household=vfoptions.n_e.household;
    simoptions.e_grid.household=e_grid_J;
    simoptions.pi_e.household=pi_e_J;
end

% z_grid and pi_z for household
z_grid.household=z_grid_J;
pi_z.household=pi_z_J;


%% Grids for firm
if scenario<4
    d_grid.firm=linspace(0,1+floor(log(ypp)),n_d.firm)'; % Notice that it is imposing the d>=0 condition implicitly
    % k_max=10 replicates OLGModel14; K>4=infeasible when ypp=1, but need more as ypp increases
    k_max=[10,6+ceil(log(ypp)),10+ceil(log(ypp)),10+ceil(log(ypp))];
    k_grid_cubed=linspace(0,1,ceil(n_a.firm/2)).^3; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.
    k_grid_linear=linspace(1,k_max(scenario),floor(n_a.firm/2)+1);
    k_grid=[k_grid_cubed, k_grid_linear(2:end)];
    a_grid.firm=k_grid';
    pv_grid_firm=NaN;
else
    d_grid.firm=linspace(0,1,n_d.firm(1))'; % Electrification investment
    % k_max=10 replicates OLGModel14; K>4=infeasible when ypp=1, but need more as ypp increases
    k_max=6+ceil(log(ypp));
    k_grid_cubed=linspace(0,1,ceil(n_a.firm(1)/2)).^3; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.
    k_grid_linear=linspace(1,k_max,floor(n_a.firm(1)/2)+1);
    k_grid=[k_grid_cubed, k_grid_linear(2:end)];
    % 300 * 200GWh PV = 60 TWh solar generation of 69 TWh current fossil sources
    pv_grid_firm=linspace(0,100*(1+1/(n_a.firm(2)-1)),n_a.firm(2))-100/(n_a.firm(2)-1);
    a_grid.firm=[k_grid'; pv_grid_firm'];
end

if n_z.firm==1
    z_grid.firm=zeros(n_z.firm,1);
    pi_z.firm=ones(n_z.firm,n_z.firm);
else
    [z_grid.firm,pi_z.firm] = discretizeAR1_FarmerToda(0,Params.rho_z_firm,Params.sigma_z_e_firm,n_z.firm);
end
z_grid.firm=exp(z_grid.firm);


%% Grids for energy
if scenario < 4
    d_grid.energy=0; % Notice that it is imposing the d>=0 condition implicitly
    a_grid.energy=linspace(0,1,n_a.energy)'; % Nothing in particular
    pv_grid_energy=NaN;
else
    d_grid.energy=linspace(0,1,n_d.energy)'; % Notice that it is imposing the d>=0 condition implicitly
    % 300 * 200GWh PV = 60 TWh solar generation of 69 TWh current fossil sources
    pv_grid_energy=linspace(0,200*(1+1/(n_a.energy-1)),n_a.energy)-100/(n_a.energy-1);
    a_grid.energy=pv_grid_energy'; % PV assets
end

if small_z_no_e
    z_grid.energy=zeros(n_z.energy,1);
    pi_z.energy=ones(n_z.energy,n_z.energy);
else
    [z_grid.energy,pi_z.energy] = discretizeAR1_FarmerToda(0,Params.rho_z_firm,Params.sigma_z_e_energy,n_z.energy);
end
z_grid.energy=exp(z_grid.energy);

%% Initial distribution of agents at birth (j=1)
% Before we plot the life-cycle profiles we have to define how agents are
% at age j=1. We will give them all zero shares (and possibly zero assets, no house, no solarpv).
if small_z_no_e
    jequaloneDist.household=zeros([n_a.household,n_z.household],'gpuArray'); % Put no households anywhere on grid
    if scenario<3
        % All agents start with zero shares, and the median shocks
        jequaloneDist.household(1,floor((n_z.household+1)/2))=1;
    elseif scenario<4
        % All agents start with zero shares, assets, houses, solarpv, and median shocks
        jequaloneDist.household(1,zeroassetindex,1,1,floor((n_z.household+1)/2))=1;
    else
        % All agents start with zero shares, assets, cars, houses, solarpv, and median shocks
        jequaloneDist.household(1,zeroassetindex,1,1,1,floor((n_z.household+1)/2))=1;
    end
else
    jequaloneDist.household=zeros([n_a.household,n_z.household,vfoptions.n_e.household],'gpuArray'); % Put no households anywhere on grid
    if scenario<3
        % All agents start with zero shares, and the median shocks
        jequaloneDist.household(1,floor((n_z.household+1)/2),floor((simoptions.n_e.household+1)/2))=1;
    elseif scenario<4
        % All agents start with zero shares, assets, houses, solarpv, and median shocks
        jequaloneDist.household(1,zeroassetindex,1,1,floor((n_z.household+1)/2),floor((simoptions.n_e.household+1)/2))=1;
    else
        % All agents start with zero shares, assets, cars, houses, solarpv, and median shocks
        jequaloneDist.household(1,zeroassetindex,1,1,1,floor((n_z.household+1)/2),floor((simoptions.n_e.household+1)/2))=1;
    end
end

% Note that because the firms are infinite horizon they do not have an age=1 distribution

% We cannot store these values in a structure because we cannot pass structures via Params to applyfun.
Params.pv_max_firm=pv_grid_firm(end);
Params.pv_max_household=pv_grid_hh(end);
Params.pv_max_energy=pv_grid_energy(end);

% Last aspect of grid: can we divide and conquer?
vfoptions.divideandconquer.household = logical(scenario<3);


end