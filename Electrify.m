%% OLG Electrification (based on OLGModels14: Heterogenous households and heterogeneous firms
%% and also Life-Cycle Model 35: Portfolio-Choice with Housing)
% See https://www.vfitoolkit.com/updates-blog/2021/an-introduction-to-life-cycle-models/
% OLGModel14.m in the repo https://github.com/vfitoolkit/IntroToOLGModels
% and LifeCycleModel35.m in the repo https://github.com/vfitoolkit/IntroToLifeCycleModels

solve_GE=true;

Names_i={'household','firm'};
PTypeDistParamNames={'ptypemass'};
Params.ptypemass=[1,1]; % Mass of households and firms are each equal to one

% A line some need for running on the Server
addpath(genpath('./MatlabToolkits/'))

%% Parameters for household (3 scenarios)
% Scenario 1: no housing, no assets, no inflation
% Scenario 2: add rental+energy costs, but no housing/assets/inflation
% Scenario 3: add housing/assets/inflation
Params.scenario=3;

% To be able to solve such a big problem, I switched to 5 year model period.
% Note that ypp (years-per-period) must be at most 15 (for kappa_j labor productivity evolutions).
% Discounting parameters (beta_pp and sj) defined in terms of ypp
Params.ypp=1; % model period, in years (just used this to modify some parameters from annual to model period)

%% Begin setting up to use VFI Toolkit to solve
% vfoptions.howardsgreedy=0;
% vfoptions.howards=80;
% vfoptions.maxhowards=200;
if Params.scenario<3
    vfoptions.tolerance=10^(-9);
else
    vfoptions.tolerance=10^(-4);
end
% Note that simoptions.tolerance is used very differently than vfoptions.tolerance

% The user can experiment with gridinterplayer=0 (pure discretization) or gridinterplayer=1 (linear interpolation b/w grid points).
% If gridinterplayer=1, then you must set vfoptions.divideandconquer=1 (required for transition).
vfoptions.gridinterplayer  = 0;
vfoptions.ngridinterp      = 10;
vfoptions.divideandconquer.household = (Params.scenario<3);
vfoptions.divideandconquer.firm = 0;
vfoptions.level1n=11;
simoptions.gridinterplayer = vfoptions.gridinterplayer;
simoptions.ngridinterp     = vfoptions.ngridinterp;

% Grid sizes to use for household

% Lets model agents from age 20 to age 100, so 81 periods (or 61 for scenario 3)
max_J=[100,100,80];
Params.agejshifter=19; % Age 20 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=ceil((max_J(Params.scenario)-Params.agejshifter)/Params.ypp); % =60/ypp, Number of period in life-cycle
if Params.scenario<3
    n_d.household=101;
    n_a.household=201;
    n_z.household=3+2*floor(1.7*log(min(Params.J,60))); % AR(1) with age-dependent params = 15 with 60 periods
    vfoptions.lowmemory.household=0;
else
    n_d.household=[31,5]; % Decisions: labor, buyhouse (5)
    n_a.household=[13,15,5,4]; % Endogenous shares, assets (>=6), housing (>=2), and solarpv (4) assets (0-45 kW generation)
    n_z.household=1+2*floor(1.2*log(min(Params.J,60))); % AR(1) with age-dependent params = 7 with 60 periods
    vfoptions.lowmemory.household=2;
end
% Exogenous labor productivity units shocks (next two lines)
vfoptions.n_e.household=3; % iid
N_j.household=Params.J; % Number of periods in finite horizon

% Grids to use for firm
n_d.firm=101; % Dividend payment
n_a.firm=201; % Capital holdings
n_z.firm=3+2*floor(log(min(Params.J,60))); % Productivity shock; scaled to model, not firm horizon
N_j.firm=Inf; % Infinite horizon
vfoptions.lowmemory.firm=0;

%% Global parameters (applies to household and firm)
% Annual risk-free rate of return
r=0.05; Params.r_pp=(1+r)^Params.ypp-1;
r_wedge=0.05; Params.r_r_wedge_pp=(1+r+r_wedge)^Params.ypp-1;

% Housing
% Params.minhouse % set below, is the minimum value of house that can be purchased
rentprice=[0,0.3,0.3]; % I figured setting rent a decent fraction of income is sensible
Params.rentprice=rentprice(Params.scenario);
houseservices=[1,1,1]; % housing services as a fraction of house value
Params.houseservices=houseservices(Params.scenario);
energy_pct_cost=[0,0.05,0.07]; % Electricity: 3%; Gas: 1-2%; Petrol: 1-2%
Params.energy_pct_cost=energy_pct_cost(Params.scenario);
if Params.scenario==3
    Params.f_htc=0.03; % transaction cost of buying/selling house (is a percent of h+hprime)
    f_coll=[0,0,0.5]; % collateral contraint (fraction of house value that can be borrowed)
    Params.f_coll=f_coll(Params.scenario);
    pv_pct_cost=[0,0,0.05]; % modeling a $30K install for a $600K house
    Params.pv_pct_cost=pv_pct_cost(Params.scenario);
end

% Discount rate; Changed to get S to increase nearer to 1 given r=0.05
% (ran it with beta=0.99, got S=0.3, so increased this; note that it interacts with sj to give the actual discount factor)
beta=[0.95,0.95,0.98];
Params.beta_pp = beta(Params.scenario)^Params.ypp;
% Preferences
Params.sigma = 2; % Coeff of relative risk aversion (curvature of consumption)
sigma_h=[0,0,0.5];
Params.sigma_h=sigma_h(Params.scenario); % Relative importance of housing services (vs consumption) in utility
Params.eta = 1.5; % Curvature of leisure (This will end up being 1/Frisch elasty)
psi = [2, 2, 2]; % Weight on leisure
Params.psi=psi(Params.scenario);
% Labor productivity at start, peak, and end of working life
k_j1 = [0.5, 0.5, 0.5];
k_j2 = [2, 2, 2];
k_j2_length = [0,0,5];
k_j3 = [1, 1, 1];

% Demographics
Params.agej=1:1:Params.J; % Is a vector of all the periods: 1,2,3,...,J
Params.Jr=round((65-Params.agejshifter)/Params.ypp); % Age 65 (period 10 is ages 65-69 in the 5 year case)
% Population growth rate
n=0.02; Params.n_pp=(1+n)^Params.ypp-1; % percentage rate (expressed as fraction) of population growth per period

% Age-dependent labor productivity units
% Stage 1: starting out (typ. first 25-30 years)
% Stage 2: peak earnings (typ. years 25-30 (meaning ages 45-50))
% Stage 3: winding down (typ. last 14 years before retirement (ages 50-64))
% Stage r: retirement
if Params.Jr>5
    kappa_j12=linspace(k_j1(Params.scenario),k_j2(Params.scenario),Params.Jr-round((15+k_j2_length(Params.scenario))/Params.ypp));
    kappa_j2s=k_j2(Params.scenario)*ones(1,ceil(k_j2_length(Params.scenario)/Params.ypp));
    kappa_j23=linspace(k_j2(Params.scenario),k_j3(Params.scenario),ceil(14/Params.ypp));
else
    kappa_j12=linspace(k_j1(Params.scenario),k_j2(Params.scenario),Params.Jr-1-min(k_j2_length(Params.scenario),1));
    kappa_j2s=k_j2(Params.scenario)*ones(1,min(k_j2_length(Params.scenario),1)); % At most one period of max wage
    kappa_j23=k_j3(Params.scenario)*ones(1,1); % One period of "pre-retirement" work
end
kappa_jr=zeros(1,Params.J-Params.Jr+1);
kappa_j=[kappa_j12, kappa_j2s, kappa_j23, kappa_jr];

% If Params.J is rounded up, don't add extra zeros
Params.kappa_j=kappa_j(1:Params.J);

% Life-cycle AR(1) process z, on (log) labor productivity units
% Chosen following Karahan & Ozkan (2013) [as used by Fella, Gallipoli & Pan (2019)]
% Note that 37 covers 24 to 60 inclusive (as in the original)
% Now repeat the first and last values to fill in working age, and put zeros for retirement (where it is anyway irrelevant)
ones_pp4y=ones(1,ceil(4/Params.ypp));
rho_z=0.7596+0.2039*((1:Params.ypp:37)/10)-0.0535*((1:Params.ypp:37)/10).^2+0.0028*((1:Params.ypp:37)/10).^3; % Chosen following Karahan & Ozkan (2013) [as used by Fella, Gallipoli & Pan (2019)]
sigma_epsilon_z=0.0518-0.0405*((1:Params.ypp:37)/10)+0.0105*((1:Params.ypp:37)/10).^2-0.0002*((1:Params.ypp:37)/10).^3; % Chosen following Karahan & Ozkan (2013) [as used by Fella, Gallipoli & Pan (2019)]

% Here we allow one period each at the start and end of working age, followed by retirement
Params.rho_z=[rho_z(1)*ones_pp4y, ...
    rho_z, ...
    rho_z(end)*ones_pp4y, ...
    zeros(1,Params.J-Params.Jr+1)];
Params.sigma_epsilon_z=[sigma_epsilon_z(1)*ones_pp4y, ...
    sigma_epsilon_z, ...
    sigma_epsilon_z(end)*ones_pp4y, ...
    sigma_epsilon_z(end)*ones(1,Params.J-Params.Jr+1)];

% Transitory iid shock
sigma_e=0.0410+0.0221*((24:Params.ypp:60)/10)-0.0069*((24:Params.ypp:60)/10).^2+0.0008*((24:Params.ypp:60)/10).^3;
Params.sigma_e=[sigma_e(1)*ones_pp4y, ...
    sigma_e, ...
    sigma_e(end)*ones_pp4y, ...
    sigma_e(end)*ones(1,Params.J-Params.Jr+1)];

% Note: These iid shocks will interact with the endogenous labor so the final labor
% earnings process will not equal that of Karahan & Ozkan (2013)
% Note: Karahan & Ozkan (2013) also have a fixed effect (which they call alpha) and which I ignore here.

% Age-dependent universal, permanent supply-side shocks to housing, energy, and pv installation costs
shock_period=ceil(10/Params.ypp);
shock_years=(1+shock_period:shock_period:ceil(101/Params.ypp));
% Exponentially increasing every shock_period years from from ~2% to ~7% after initial shock-free period
shock_pct=[zeros(1,length(shock_years));zeros(1,length(shock_years));exp(shock_years/(Params.Jr-1))/50];
Params.agej_pct_cost=[zeros(1,1+shock_period), repelem(shock_pct(Params.scenario,:),shock_period)];
% Only worry about the shocks within the model horizon
Params.agej_pct_cost=Params.agej_pct_cost(1:Params.J);

% Conditional survival probabilities: sj is the probability of surviving to be age j+1, given alive at age j
% Most countries have calculations of these (as they are used by the government departments that oversee pensions)
% In fact I will here get data on the conditional death probabilities, and then survival is just 1-death.
% Here I just use them for the US, taken from "National Vital Statistics Report, volume 58, number 10, March 2010."
% I took them from first column (qx) of Table 1 (Total Population)
% Conditional death probabilities
dj=[0.006879, 0.000463, 0.000307, 0.000220, 0.000184, 0.000172, 0.000160, 0.000149, 0.000133, 0.000114, 0.000100, 0.000105, 0.000143, 0.000221, 0.000329, 0.000449, 0.000563, 0.000667, 0.000753, 0.000823,...
    0.000894, 0.000962, 0.001005, 0.001016, 0.001003, 0.000983, 0.000967, 0.000960, 0.000970, 0.000994, 0.001027, 0.001065, 0.001115, 0.001154, 0.001209, 0.001271, 0.001351, 0.001460, 0.001603, 0.001769, 0.001943, 0.002120, 0.002311, 0.002520, 0.002747, 0.002989, 0.003242, 0.003512, 0.003803, 0.004118, 0.004464, 0.004837, 0.005217, 0.005591, 0.005963, 0.006346, 0.006768, 0.007261, 0.007866, 0.008596, 0.009473, 0.010450, 0.011456, 0.012407, 0.013320, 0.014299, 0.015323,...
    0.016558, 0.018029, 0.019723, 0.021607, 0.023723, 0.026143, 0.028892, 0.031988, 0.035476, 0.039238, 0.043382, 0.047941, 0.052953, 0.058457, 0.064494,...
    0.071107, 0.078342, 0.086244, 0.094861, 0.104242, 0.114432, 0.125479, 0.137427, 0.150317, 0.164187, 0.179066, 0.194979, 0.211941, 0.229957, 0.249020, 0.269112, 0.290198, 0.312231, 1.000000]; 
dj=resize(dj,101+Params.ypp,FillValue=1);
% dj covers Ages 0-100, plus extras at end to make it period-friendly
Params.sj=prod(1-reshape(dj(1:Params.ypp*ceil(101/Params.ypp)),[Params.ypp,ceil(101/Params.ypp)]),1); % p5-year survival rates
Params.sj=Params.sj(1+ceil(20/Params.ypp):ceil(20/Params.ypp)+N_j.household); % Just the ages we are using (20yo and up)
Params.sj(end)=0; % In the present model the last period (j=J) value of sj is actually irrelevant

% Warm glow of bequest
Params.warmglow1=0.3; % (relative) importance of bequests
Params.warmglow2=3; % bliss point of bequests (essentially, the target amount)
Params.warmglow3=Params.sigma; % By using the same curvature as the utility of consumption it makes it much easier to guess appropraite parameter values for the warm glow

% Taxes
Params.tau_l = 0.2; % Tax rate on labour income

%% Parameters for firm
% Production
Params.alpha_k=0.311; % diminishing returns to capital input
Params.alpha_l=0.650; % diminishing returns to labor input
Params.delta=0.05; % Annual depreciation of physical capital
% Capital adjustment costs
Params.capadjconstant=1.21; % term in the capital adjustment cost
% Tax
Params.tau_corp=0.34; % Tax rate on corporate earnings
Params.phi=0.5; % Fraction of capital adjustment costs that can be deducted from corporate earnings
Params.tau_d=0.2; % Tax rate on dividends
Params.tau_cg=0.2; % Tax rate on capital gains
% Idiosyncatic productivity shocks
Params.rho_z_firm=0.767;
Params.sigma_z_e_firm=0.211;

% Set the firm discount factor below (as it is determined in general eqm)
% Params.firmbeta=1/(1+Params.r_pp)/(1-Params.tau_cg)); % 1/(1+r_pp) but returns net of capital gains tax

%% Also, we want to target a capital-output ratio
% This is relevant to the general equilibrium conditions
Params.TargetKdivL=2.03;

%% Remaining parameters

% Some initial values/guesses for variables that will be determined in general eqm
Params.pension=0.4; % Initial guess (this will be determined in general eqm)
Params.w=1;
AccidentBeqS= [0.02,0.02,0.02]; % Accidental bequests (this is the lump sum transfer of shares)
AccidentBeqAH= [0,0,0.02]; % Accidental bequests (this is the lump sum transfer of assets+house value)
Params.AccidentBeqS_pp=AccidentBeqS(Params.scenario)*Params.ypp;
Params.AccidentBeqAH_pp=AccidentBeqAH(Params.scenario)*Params.ypp;
Params.G_pp=0.1*Params.ypp; % Government expenditure
Params.firmbeta=1/(1+Params.r_pp/(1-Params.tau_cg)); % 1/(1+r_pp) but returns net of capital gains tax
Params.D_pp=(1+0.2)^Params.ypp-1; % Per-period dividends rate expected/received by households
Params.P0=1;
Lhscale=[0.21,0.21,0.29]; % Scaling the household labor supply
Params.Lhscale=Lhscale(Params.scenario);

%% Grids for household

% Grid for labour choice
labor_grid=linspace(0,1,n_d.household(1))'; % Notice that it is imposing the 0<=h<=1 condition implicitly

% Grid for share holdings, always > 0
% For later scenarios, shrink the grid for more accuracy
s_grid_cubed=linspace(0,1,ceil(n_a.household(1)/2)).^3; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.
s_grid_linear=linspace(1,10-4*(Params.scenario>1),floor(n_a.household(1)/2)+1);
share_grid=[s_grid_cubed, s_grid_linear(2:end)]';

% Set up d for VFI Toolkit (is the two decision variables)
if Params.scenario<3
    d_grid.household=labor_grid;
    a_grid.household=share_grid;
    Params.minhouse=1;
else
    % Grid for bank account; a negative balance implies a mortgage
    a_grid_cubed=linspace(-1,0,ceil(n_a.household(2)/2)-1).^3;
    a_grid_linear=linspace(0,3,floor(n_a.household(2)/2)+2);
    asset_grid=[0.5*a_grid_cubed, a_grid_linear(2:end)]';
    
    % Make it so that there is a zero assets
    % Find closest to zero assets
    [~,zeroassetindex]=min(abs(asset_grid));
    asset_grid(zeroassetindex)=0;
    
    % age20avgincome=Params.w*Params.kappa_j(1);
    % house_grid=[0; logspace(2*age20avgincome, 12*age20avgincome, 5)'];
    house_grid=(0:1:n_a.household(3)-1)';
    % Note, we can see from w*kappa_j*z and the values of these, that average
    % income is going to be around one, so will just use this simpler house grid
    % [We can think about the values of the house_grid as being relative the average income (or specifically average at a given age)]
    Params.minhouse=house_grid(2); % first is zero (no house)
    
    % buyhouse decisions
    %  0=no house
    %  1=buy house w/o pv this period
    %  2=buy house w/ pv this period
    %  3=keep house; no pv upgrade
    %  4=keep house; pv upgrade (if possible)
    %  5=testing (not used)
    
    buyhouse_grid=(0:1:n_d.household(2)-1)';
    
    % kW of solar generation installed, 10 kW per grid element
    solarpv_grid=(0:1:n_a.household(4)-1)';
    
    d_grid.household=[labor_grid; buyhouse_grid];
    
    a_grid.household=[share_grid; asset_grid; house_grid; solarpv_grid];
    
    %% Solar PV is an experience asset
    vfoptions.experienceasset.household=1;
    simoptions.experienceasset.household=1;
    
    %% Define aprime function used for the experience asset
    
    % experienceasset: aprime_val=aprimeFn(d,a)
    % vfoptions.refine_d: the decision variables input to aprimeFn are d3
    aprimeFn=@(buyhouse, solarpv) ElectrifyHousing_aprimeFn(buyhouse, solarpv); % Will return the value of aprime (solarpv)
    
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
    vfoptions.refine_d.household=[1,0,1]; % tell the code how many d1, d2, and d3 there are
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
% First, z, the AR(1) with age-dependent parameters
[z_grid_J, pi_z_J] = discretizeLifeCycleAR1_FellaGallipoliPan(Params.rho_z,Params.sigma_epsilon_z,n_z.household,Params.J);
% z_grid_J is n_z-by-J, so z_grid_J(:,j) is the grid for age j
% pi_z_J is n_z-by-n_z-by-J, so pi_z_J(:,:,j) is the transition matrix for age j

% Second, e, the iid normal with age-dependent parameters
[e_grid_J, pi_e_J] = discretizeLifeCycleAR1_FellaGallipoliPan(zeros(1,Params.J),Params.sigma_e,vfoptions.n_e.household,Params.J); % Note: AR(1) with rho=0 is iid normal
% Because e is iid we actually just use
pi_e_J=shiftdim(pi_e_J(1,:,:),1);

% z_grid and pi_z for household
z_grid.household=z_grid_J;
pi_z.household=pi_z_J;

% Any (iid) e variable always has to go into vfoptions and simoptions
vfoptions.e_grid.household=e_grid_J;
vfoptions.pi_e.household=pi_e_J;
simoptions.n_e.household=vfoptions.n_e.household;
simoptions.e_grid.household=e_grid_J;
simoptions.pi_e.household=pi_e_J;


%% Grids for firm
d_grid.firm=linspace(0,1,n_d.firm)'; % Notice that it is imposing the d>=0 condition implicitly
% k_max=10 replicates OLGModel14; K>4=infeasible when ypp=1, but need more as ypp increases
k_max=[10,6,6];
k_grid_cubed=linspace(0,1,ceil(n_a.firm/2)).^3; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.
k_grid_linear=linspace(1,k_max(Params.scenario),floor(n_a.firm/2)+1);
a_grid.firm=[k_grid_cubed, k_grid_linear(2:end)]';

[z_grid.firm,pi_z.firm] = discretizeAR1_FarmerToda(0,Params.rho_z_firm,Params.sigma_z_e_firm,n_z.firm);
z_grid.firm=exp(z_grid.firm);


%% Now, create the return function

% For households
DiscountFactorParamNames.household={'beta_pp','sj'};

if Params.scenario<3
    ReturnFn.household=@( ...
        labor,sprime,s,z,e, ...
        pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
        sigma,psi,eta,sigma_h,kappa_j,tau_l,tau_d,tau_cg,warmglow1,warmglow2,ypp,agej,Jr,J,...
        scenario,r_pp,r_r_wedge_pp,minhouse,rentprice,houseservices,energy_pct_cost ...
    ) Electrify_HouseholdReturnFn( ...
        labor,0,sprime,0,0,s,0,0,0,z,e, ...
        pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
        sigma,psi,eta,sigma_h,kappa_j,tau_l,tau_d,tau_cg,warmglow1,warmglow2,ypp,agej,Jr,J,...
        scenario,r_pp,r_r_wedge_pp,0,minhouse,rentprice,0,houseservices,0,0,energy_pct_cost ...
    );
else
    % Notice we use 'Electrify_HouseholdReturnFn'
    ReturnFn.household=@( ...
            labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
            pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
            sigma,psi,eta,sigma_h,kappa_j,tau_l,tau_d,tau_cg,warmglow1,warmglow2,ypp,agej,Jr,J,...
            scenario,r_pp,r_r_wedge_pp,f_htc,minhouse,rentprice,f_coll,houseservices,agej_pct_cost,pv_pct_cost,energy_pct_cost ...
        ) Electrify_HouseholdReturnFn( ...
            labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
            pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
            sigma,psi,eta,sigma_h,kappa_j,tau_l,tau_d,tau_cg,warmglow1,warmglow2,ypp,agej,Jr,J,...
            scenario,r_pp,r_r_wedge_pp,f_htc,minhouse,rentprice,f_coll,houseservices,agej_pct_cost,pv_pct_cost,energy_pct_cost ...
        );
end

% For firms
DiscountFactorParamNames.firm={'firmbeta'};
% Notice we use 'Electrify_FirmReturnFn'
ReturnFn.firm=@( ...
        d,kprime,k,z, ...
        w,D_pp, ...
        ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg ...
    ) Electrify_FirmReturnFn( ...
        d,kprime,k,z, ...
        w,D_pp, ...
        ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg ...
    );

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium
disp('Test ValueFnIter')
tic;
% Note: z_grid and pi_z, this will be ignored due to presence of vfoptions.z_grid_J and vfoptions.pi_z_J
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
toc

if Params.scenario<3
    checkfeasible_household_case12(V,Policy,a_grid.household,1,n_d,n_a,n_z,N_j,d_grid,a_grid,Params,vfoptions);
else
    checkfeasible_household_case3(V,Policy,asset_grid,zeroassetindex,n_d,n_a,n_z,N_j,d_grid,a_grid,Params,vfoptions);
end

%% Initial distribution of agents at birth (j=1)
% Before we plot the life-cycle profiles we have to define how agents are
% at age j=1. We will give them all zero shares (and possibly zero assets, no house, no solarpv).
jequaloneDist.household=zeros([n_a.household,n_z.household,vfoptions.n_e.household],'gpuArray'); % Put no households anywhere on grid
if Params.scenario<3
    % All agents start with zero shares, and the median shocks
    jequaloneDist.household(1,floor((n_z.household+1)/2),floor((simoptions.n_e.household+1)/2))=1;
else
    % All agents start with zero shares, assets, houses, solarpv, and median shocks
    jequaloneDist.household(1,zeroassetindex,1,1,floor((n_z.household+1)/2),floor((simoptions.n_e.household+1)/2))=1;
end

% Note that because the firms are infinite horizon they do not have an age=1 distribution

%% Agents age distribution
% Many OLG models include some kind of population growth, and perhaps some other things that create a weighting of different ages that needs to
% be used to calculate the stationary distribution and aggregate variable.
mewj=ones(1,Params.J); % Marginal distribution of households over age
for jj=2:length(mewj)
    mewj(jj)=Params.sj(jj-1)*mewj(jj-1)/(1+Params.n_pp);
end
Params.mewj=mewj./sum(mewj); % Normalize to one

AgeWeightsParamNames=struct('household',{{'mewj'}}); % So VFI Toolkit knows which parameter is the mass of agents of each age

%% Test
disp('Test StationaryDist')
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);

%% General eqm variables
if Params.scenario<3
    GEPriceParamNames={'pension','AccidentBeqS_pp','G_pp','w','firmbeta','D_pp','P0'};
else
    GEPriceParamNames={'pension','AccidentBeqS_pp','AccidentBeqAH_pp','G_pp','w','firmbeta','D_pp','P0'};
end
heteroagentoptions.constrainpositive={'w','D_pp','P0'};

% We don't need P
% We can get P from the equation that defines r as the return to the mutual fund
% 1+r_pp = (P0 +(1-tau_d)D_pp - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have
% P=((1-tau_cg)*P0 + (1-tau_d)*D_pp)/(1+r_pp-tau_cg);

%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)
% Note: we need to add z & e to FnsToEvaluate inputs for households,
% whereas firm only has z (it is just coincidence/lazy that I call them
% both z).

% Stationary Distribution Aggregates from households (important that ordering of Names and Functions is the same)
if Params.scenario<3
    FnsToEvaluate.L_h.household = @(labor,sprime,s,z,e,kappa_j,Lhscale) ...
        labor*kappa_j*exp(z+e)*Lhscale;  % Aggregate labour supply in efficiency units, not scaled by ypp
    FnsToEvaluate.S.household = @(labor,sprime,s,z,e) s; % Aggregate share holdings
    FnsToEvaluate.PensionSpending.household = @(labor,sprime,s,z,e,pension,ypp,agej,Jr) ...
        (agej>=Jr)*pension*ypp; % Total spending on pensions
    FnsToEvaluate.PayrollTaxRevenue.household = @(labor,sprime,s,z,e,ypp,agej,Jr,tau_l,w,kappa_j,Lhscale) ...
        (agej<Jr)*tau_l*labor*w*kappa_j*exp(z+e)*ypp*Lhscale; % Total spending on payroll taxes
    FnsToEvaluate.AccidentalBeqSLeft_pp.household = @(labor,sprime,s,z,e,sj) ...
        sprime*(1-sj); % Accidental share bequests left by people who die
    FnsToEvaluate.CapitalGainsTaxRevenue.household = @(labor,sprime,s,z,e,tau_cg,P0,D_pp,tau_d,r_pp) ...
        tau_cg*(P0-(((1-tau_cg)*P0 + (1-tau_d)*D_pp)/(1+r_pp-tau_cg)))*s; % tau_cg*(P0-Plag)*s, but substitute P=Plag, and then substitute for P
else
    FnsToEvaluate.L_h.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,kappa_j,Lhscale) ...
        labor*kappa_j*exp(z+e)*Lhscale;  % Aggregate labour supply in efficiency units, not scaled by ypp
    FnsToEvaluate.S.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) s; % Aggregate share holdings
    FnsToEvaluate.A.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) a; % Aggregate asset/mortgage holdings
    FnsToEvaluate.H.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) h; % Aggregate house holdings
    FnsToEvaluate.PV.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) solarpv; % Aggregate solarpv holdings
    FnsToEvaluate.PensionSpending.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,pension,ypp,agej,Jr) ...
        (agej>=Jr)*pension*ypp; % Total spending on pensions
    FnsToEvaluate.PayrollTaxRevenue.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,ypp,agej,Jr,tau_l,w,kappa_j,Lhscale) ...
        (agej<Jr)*tau_l*labor*w*kappa_j*exp(z+e)*Lhscale*ypp; % Total spending on payroll taxes
    FnsToEvaluate.AccidentalBeqSLeft_pp.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,sj) ...
        sprime*(1-sj); % Accidental share bequests left by people who die
    % AccidentalBeqAHLeft is zero (if in debt) or accidental asset+house bequests left by people who die
    FnsToEvaluate.AccidentalBeqAHLeft_pp.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,scenario,sj,agej_pct_cost) ...
        max(0,(aprime+(1+agej_pct_cost)*hprime)*(1-sj));
    % BadDebt is the debt somebody accidentally leaves behind, or zero if net worth is positive
    FnsToEvaluate.BadDebt_pp.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,scenario,sj,agej_pct_cost) ...
        min(0,(aprime+(1+agej_pct_cost)*hprime)*(1-sj));
    FnsToEvaluate.CapitalGainsTaxRevenue.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,tau_cg,P0,D_pp,tau_d,r_pp) ...
        tau_cg*(P0-(((1-tau_cg)*P0 + (1-tau_d)*D_pp)/(1+r_pp-tau_cg)))*s+(1-tau_d)*r_pp*max(a,0); % tau_cg*(P0-Plag)*s + deposit interest, but substitute P=Plag, and then substitute for P
end

% From firms
FnsToEvaluate.Output.firm = @(d,kprime,k,z,w,ypp,alpha_k,alpha_l) ...
    z*(k^alpha_k)*((w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)))^alpha_l*ypp; % Production function z*(k^alpha_k)*(l^alpha_l) (substituting for l)
FnsToEvaluate.L_f.firm = @(d,kprime,k,z,w,alpha_k,alpha_l) ...
    (w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % (effective units of) labor demanded by firm, not scaled by ypp
FnsToEvaluate.K.firm = @(d,kprime,k,z,w,alpha_k,alpha_l) k; % physical capital
FnsToEvaluate.DividendPaid_pp.firm = @(d,kprime,k,z,ypp) (1+d)^ypp-1; % dividend paid by firm
FnsToEvaluate.Sissued.firm = @(d,kprime,k,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi) ...
    Electrify_FirmShareIssuance(d,kprime,k,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi); % Share issuance
FnsToEvaluate.CorpTaxRevenue.firm = @(d,kprime,k,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi) ...
    Electrify_FirmCorporateTaxRevenue(d,kprime,k,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi); % revenue from the corporate profits tax

% General Equilibrium conditions (these should evaluate to zero in general equilbrium)
GeneralEqmEqns.sharemarket = @(S) S-1; % mass of all shares equals one
GeneralEqmEqns.labormarket = @(L_h,L_f) (Params.scenario+1)*(L_h-L_f)*Params.ypp; % labor supply of households equals labor demand of firms (scaled by ypp)
GeneralEqmEqns.pensions = @(PensionSpending,PayrollTaxRevenue) PensionSpending-PayrollTaxRevenue; % Retirement benefits equal Payroll tax revenue: pension*fractionretired-tau*w*H
GeneralEqmEqns.bequestsS_pp = @(AccidentalBeqSLeft_pp,AccidentBeqS_pp,n_pp) AccidentalBeqSLeft_pp/(1+n_pp)-AccidentBeqS_pp; % Accidental share bequests received equal accidental share bequests left
if Params.scenario==3
    GeneralEqmEqns.bequestsAH_pp = @(AccidentalBeqAHLeft_pp,AccidentBeqAH_pp,n_pp) AccidentalBeqAHLeft_pp/(1+n_pp)-AccidentBeqAH_pp; % Accidental asset+house bequests received equal accidental asset+house bequests left
end
GeneralEqmEqns.govbudget = @(G_pp,tau_d,D_pp,CapitalGainsTaxRevenue,CorpTaxRevenue) G_pp-tau_d*D_pp-CapitalGainsTaxRevenue-CorpTaxRevenue; % G is equal to the target, GdivYtarget*Y
GeneralEqmEqns.firmdiscounting = @(firmbeta,r_pp,tau_cg) firmbeta-1/(1+r_pp/(1-tau_cg)); % Firms discount rate is related to market return rate
GeneralEqmEqns.dividends = @(D_pp,DividendPaid_pp) D_pp-DividendPaid_pp; % That the dividend households receive equals that which firms give
GeneralEqmEqns.ShareIssuance = @(Sissued,P0,D_pp,tau_cg,tau_d,r_pp) ...
    P0-((((1-tau_cg)*P0 + (1-tau_d)*D_pp)/(1+r_pp-tau_cg))-Sissued); % P0=P-S, but substitute for P (see derivation inside the return fn)
GeneralEqmEqns.CapitalOutputRatio =@(K,L_f,TargetKdivL) K/L_f-TargetKdivL; % Ratio not based on ypp

%% Test
% Note: Because we used simoptions we must include this as an input
disp('Test AggVars')
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params, n_d, n_a, n_z,N_j,Names_i, d_grid, a_grid, z_grid,simoptions);

% Next few lines were used to try a few parameter values so as to get a
% decent initial guess before actually solving the general equilbrium
fprintf('Check: L_h, L_f, K \n')
[AggVars.L_h.Mean,AggVars.L_f.Mean,AggVars.K.Mean]
fprintf('Check: K/L_f (should be about 2.03) \n')
AggVars.K.Mean/AggVars.L_f.Mean
if Params.scenario<3
    fprintf('Check: S \n')
    [AggVars.S.Mean]
else
    fprintf('Check: S, A, H, PV \n')
    [AggVars.S.Mean,AggVars.A.Mean,AggVars.H.Mean,AggVars.PV.Mean]
end
fprintf('Check: ShareIssuance GE condition \n')
Params.P0-((((1-Params.tau_cg)*Params.P0 + (1-Params.tau_d)*Params.D_pp)/(1+Params.r_pp-Params.tau_cg))-AggVars.S.Mean)


%% Solve for the General Equilibrium
if solve_GE
    % heteroagentoptions.fminalgo=4 % CMA-ES algorithm 
    
    heteroagentoptions.verbose=1;
    if Params.scenario<3
        heteroagentoptions.toleranceGEprices=10^(-4);
        heteroagentoptions.toleranceGEcondns=10^(-4); % This is the hard one
    else
        heteroagentoptions.toleranceGEprices=10^(-3);
        heteroagentoptions.toleranceGEcondns=10^(-3); % This is the hard one
        heteroagentoptions.maxiter=600;
    end
        
    p_eqm=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, [], pi_z, d_grid, a_grid, z_grid,jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightsParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
    % p_eqm contains the general equilibrium parameter values
    % Put this into Params so we can calculate things about the initial equilibrium
    Params.pension=p_eqm.pension;
    Params.AccidentBeqS_pp=p_eqm.AccidentBeqS_pp;
    if Params.scenario==3
        Params.AccidentBeqAH_pp=p_eqm.AccidentBeqAH_pp;
    end
    Params.G_pp=p_eqm.G_pp;
    Params.w=p_eqm.w;
    Params.firmbeta=p_eqm.firmbeta;
    Params.D_pp=p_eqm.D_pp;
    Params.P0=p_eqm.P0;

    % Re-Calculate a few things related to the general equilibrium.
    [V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j, Names_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
    StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
end % Otherwise we can use initial guesses set for Params

% Can just use the same FnsToEvaluate as before.
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

%% Plot the life cycle profiles of capital and labour for the inital and final eqm.

figure_c=figure(10);
if Params.scenario<3
    subplot(2,1,1); plot(1:1:Params.J,AgeConditionalStats.L_h.Mean)
    title('Life Cycle Profile: Effective Labour Supply')
    subplot(2,1,2); plot(1:1:Params.J,AgeConditionalStats.S.Mean)
    title('Life Cycle Profile: Share holdings')
else
    subplot(3,2,1); plot(1:1:Params.J,AgeConditionalStats.L_h.Mean)
    title('Life Cycle Profile: Effective Labour Supply')
    subplot(3,2,3); plot(1:1:Params.J,AgeConditionalStats.S.Mean)
    title('Life Cycle Profile: Share holdings')
end
if Params.scenario==3
    subplot(3,2,2); plot(1:1:Params.J,AgeConditionalStats.A.Mean)
    title('Life Cycle Profile: Asset holdings')
    subplot(3,2,4); plot(1:1:Params.J,AgeConditionalStats.H.Mean)
    title('Life Cycle Profile: House holdings')
    subplot(3,2,6); plot(1:1:Params.J,AgeConditionalStats.PV.Mean)
    title('Life Cycle Profile: Solar PV installed')
end
saveas(figure_c,'./SavedOutput/Graphs/Electrify_LifeCycleProfiles','pdf')

%% Calculate some aggregates and print findings about them

% Add consumption to the FnsToEvaluate
if Params.scenario<3
    FnsToEvaluate.Consumption.household=@( ...
            labor,sprime,s,z,e, ...
            pension,AccidentBeqS_pp,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp, r_r_wedge_pp,rentprice,energy_pct_cost ...
        ) Electrify_HouseholdConsumptionFn( ...
            labor,0,sprime,0,0,s,0,0,0,z,e, ...
            pension,AccidentBeqS_pp,0,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,r_r_wedge_pp,0,rentprice,0,0,energy_pct_cost);
    FnsToEvaluate.Income.household=@( ...
            labor,sprime,s,z,e, ...
            pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,r_r_wedge_pp ...
        ) Electrify_HouseholdIncomeFn( ...
            labor,0,sprime,0,0,s,0,0,0,z,e, ...
            pension,AccidentBeqS_pp,0,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,r_r_wedge_pp,0);
else
    FnsToEvaluate.Consumption.household=@( ...
            labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
            pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,r_r_wedge_pp,f_htc,rentprice,agej_pct_cost,pv_pct_cost,energy_pct_cost ...
        ) Electrify_HouseholdConsumptionFn( ...
            labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
            pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,r_r_wedge_pp,f_htc,rentprice,agej_pct_cost,pv_pct_cost,energy_pct_cost);
    FnsToEvaluate.Income.household=@( ...
            labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
            pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,r_r_wedge_pp,agej_pct_cost ...
        ) Electrify_HouseholdIncomeFn( ...
            labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
            pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,r_r_wedge_pp,agej_pct_cost);
end

AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params, n_d, n_a, n_z,N_j, Names_i, d_grid, a_grid, z_grid,simoptions);

Y=AggVars.Output.Mean;

P=((1-Params.tau_cg)*Params.P0 + (1-Params.tau_d)*Params.D_pp)/(1+Params.r_pp-Params.tau_cg);

% Calculate the aggregate TFP as output/((capital^alpha_k)*(labor^alpha_l))
AggregateTFP=Y/((AggVars.K.Mean^Params.alpha_k)*(AggVars.L_f.Mean^Params.alpha_l));

% Total value of firms
temp=V.firm.*StationaryDist.firm;
temp(StationaryDist.firm==0)=0; % Get rid of points that have V=-inf but zero mass which would give nan
TotalValueOfFirms=sum(temp(isfinite(temp)));

fileID = fopen('SavedOutput\aggs.txt','w');
fprintf(fileID,'Following are some aggregates of the model economy (Scenario %d): \n', Params.scenario);
fprintf(fileID,'Output: Y=%8.2f \n',AggVars.Output.Mean);
fprintf(fileID,'Aggregate TFP: Y=%8.2f \n',AggregateTFP);
fprintf(fileID,'Capital-Output ratio (firm side): K/Y=%8.2f \n',AggVars.K.Mean/Y);
if Params.scenario<3
    fprintf(fileID,'Total share value (HH side): P*S (%.2f) = %8.2f\n',P*AggVars.S.Mean,P*AggVars.S.Mean);
else
    fprintf(fileID,'Total share+asset value (HH side): P*S (%.2f) + A (%.2f) = %8.2f\n',P*AggVars.S.Mean,AggVars.A.Mean,P*AggVars.S.Mean+AggVars.A.Mean);
    fprintf(fileID,'Total house value (HH side): H=%8.2f \n',AggVars.H.Mean);
    fprintf(fileID,'Total bad debt (HH side): P*S=%8.2f \n',AggVars.BadDebt_pp.Mean);
end
fprintf(fileID,'Total firm value (firm side): Value of firm=%8.2f \n',TotalValueOfFirms);
fprintf(fileID,'Consumption-Output ratio: C/Y=%8.2f \n',AggVars.Consumption.Mean/Y);
fprintf(fileID,'Government-to-Output ratio: G/Y=%8.2f \n', Params.G_pp/Y);
fprintf(fileID,'Wage: w=%8.2f \n',Params.w);
fclose(fileID);

type 'SavedOutput\aggs.txt'

function feasible=checkfeasible_household_case12(V,Policy,asset_grid,zeroassetindex,n_d,n_a,n_z,N_j,d_grid,a_grid,Params,vfoptions)

feasible=true;
fhh_options.tolerance=vfoptions.tolerance;
fhh_options.lowmemory=vfoptions.lowmemory.household;
fhh_options.n_e=vfoptions.n_e.household;
fhh_options.e_grid=vfoptions.e_grid.household;
fhh_options.pi_e=vfoptions.pi_e.household;
fhh=PolicyInd2Val_FHorz(Policy.household,n_d.household,n_a.household,n_z.household,N_j.household,d_grid.household,a_grid.household,fhh_options);

z_idx=ceil(n_z.household/2);
e_idx=ceil(fhh_options.n_e/2);
aprime_index_last=1;
for fhh_agej=1:5
    fhh_next=squeeze(fhh(:,:,z_idx,e_idx,fhh_agej));
    aprime_index_next=find(asset_grid==fhh_next(2,aprime_index_last));
    if fhh_next(1,aprime_index_next)==0
        feasible=false;
        error("labor strike")
    end

    income = Electrify_HouseholdIncomeFn( ...
        1,0,asset_grid(aprime_index_next),0,0,asset_grid(aprime_index_last),0,0,0,0,0, ...
        Params.pension,Params.AccidentBeqS_pp,0,Params.w,Params.P0,Params.D_pp, ...
        Params.kappa_j(fhh_agej),Params.tau_l,Params.tau_d,Params.tau_cg,Params.ypp,fhh_agej,Params.Jr, ...
        Params.r_pp,0,0);
    if income<0
        feasible=false;
        error("income negative")
    end
    consumption=Electrify_HouseholdConsumptionFn( ...
        1,0,asset_grid(aprime_index_next),0,0,asset_grid(aprime_index_last),0,0,0,0,0, ...
        Params.pension,Params.AccidentBeqS_pp,0,Params.w,Params.P0,Params.D_pp, ...
        Params.kappa_j(fhh_agej),Params.tau_l,Params.tau_d,Params.tau_cg,Params.ypp,fhh_agej,Params.Jr, ...
        Params.r_pp,0,0,Params.rentprice,0,0,Params.energy_pct_cost);
    if consumption<0
        feasible=false;
        error("consumption negative")
    end
    if aprime_index_next==1
        F=Electrify_HouseholdReturnFn( ...
            1,0,asset_grid(aprime_index_next),0,0,asset_grid(aprime_index_last),0,0,0,0,0, ...
            Params.pension,Params.AccidentBeqS_pp,Params.AccidentBeqAH_pp,Params.w,Params.P0,Params.D_pp, ...
            Params.sigma,Params.psi,Params.eta,Params.sigma_h,Params.kappa_j(fhh_agej),Params.tau_l,Params.tau_d,Params.tau_cg,Params.warmglow1,Params.warmglow2,Params.ypp,fhh_agej,Params.Jr,Params.J, ...
            Params.scenario,Params.r_pp,0,0,Params.minhouse,Params.rentprice,0,Params.houseservices,0,0,Params.energy_pct_cost);
        if isfinite(F)
            fprintf("F = %.2f \n", F);
        else
            feasible=false;
            error("infeasible last->next")
        end
    end
    aprime_index_last=aprime_index_next;
end

% all(squeeze(Policy.household(:,1,1:zeroassetindex-2,1,1,4:6,:,2:20))==1,[3 4 5])
% if all(squeeze(Policy.household(:,1,1:zeroassetindex-2,1,1,3:7,:,2:20))==1,'all')
%     warning("infeasible Stationary Distribution")
% end

end

function feasible=checkfeasible_household_case3(V,Policy,asset_grid,zeroassetindex,n_d,n_a,n_z,N_j,d_grid,a_grid,Params,vfoptions)

feasible=true;
fhh_options.tolerance=vfoptions.tolerance;
fhh_options.lowmemory=vfoptions.lowmemory.household;
fhh_options.n_e=vfoptions.n_e.household;
fhh_options.experienceasset=vfoptions.experienceasset.household;
fhh_options.aprimeFn=vfoptions.aprimeFn.household;
fhh_options.refine_d=vfoptions.refine_d.household;
fhh_options.e_grid=vfoptions.e_grid.household;
fhh_options.pi_e=vfoptions.pi_e.household;
fhh=PolicyInd2Val_FHorz(Policy.household,n_d.household,n_a.household,n_z.household,N_j.household,d_grid.household,a_grid.household,fhh_options);

z_idx=ceil(n_z.household/2);
e_idx=ceil(fhh_options.n_e/2);
aprime_index_last=zeroassetindex;
for fhh_agej=1:5
    fhh_next=squeeze(fhh(:,1,:,1,1,z_idx,e_idx,fhh_agej));
    aprime_index_next=find(asset_grid==fhh_next(4,aprime_index_last));
    if fhh_next(1,aprime_index_next)==0
        feasible=false;
        warning("labor strike")
    end
    income = Electrify_HouseholdIncomeFn( ...
        1,0,0,asset_grid(aprime_index_next),0,0,asset_grid(aprime_index_last),0,0,0,0, ...
        Params.pension,Params.AccidentBeqS_pp,Params.AccidentBeqAH_pp,Params.w,Params.P0,Params.D_pp, ...
        Params.kappa_j(fhh_agej),Params.tau_l,Params.tau_d,Params.tau_cg,Params.ypp,fhh_agej,Params.Jr, ...
        Params.r_pp,Params.r_r_wedge_pp,Params.agej_pct_cost(fhh_agej));
    if income<0
        feasible=false;
        error("income negative")
    end
    consumption=Electrify_HouseholdConsumptionFn( ...
        1,0,0,asset_grid(aprime_index_next),0,0,asset_grid(aprime_index_last),0,0,0,0, ...
        Params.pension,Params.AccidentBeqS_pp,Params.AccidentBeqAH_pp,Params.w,Params.P0,Params.D_pp, ...
        Params.kappa_j(fhh_agej),Params.tau_l,Params.tau_d,Params.tau_cg,Params.ypp,fhh_agej,Params.Jr, ...
        Params.r_pp,Params.r_r_wedge_pp,Params.f_htc,Params.rentprice,Params.agej_pct_cost(fhh_agej),Params.pv_pct_cost,Params.energy_pct_cost);
    if consumption<0
        feasible=false;
        error("consumption negative")
    end
    if aprime_index_next==1
        sprime=0;
        s=0;
        F=Electrify_HouseholdReturnFn( ...
            1,0,sprime,asset_grid(aprime_index_next),0,s,asset_grid(aprime_index_last),0,0,0,0, ...
            Params.pension,Params.AccidentBeqS_pp,Params.AccidentBeqAH_pp,Params.w,Params.P0,Params.D_pp, ...
            Params.sigma,Params.psi,Params.eta,Params.sigma_h,Params.kappa_j(fhh_agej),Params.tau_l,Params.tau_d,Params.tau_cg,Params.warmglow1,Params.warmglow2,Params.ypp,fhh_agej,Params.Jr,Params.J, ...
            Params.scenario,Params.r_pp,Params.r_r_wedge_pp,Params.f_htc,Params.minhouse,Params.rentprice,Params.f_coll,Params.houseservices,Params.agej_pct_cost(fhh_agej),Params.pv_pct_cost,Params.energy_pct_cost);
        if isfinite(F)
            fprintf("F = %.2f \n", F);
        else
            feasible=false;
            error("infeasible last->next")
        end
    end
    aprime_index_last=aprime_index_next;
end

% all(squeeze(Policy.household(:,1,1:zeroassetindex-2,1,1,4:6,:,2:20))==1,[3 4 5])
% if all(squeeze(Policy.household(:,1,1:zeroassetindex-2,1,1,3:7,:,2:20))==1,'all')
%     warning("infeasible Stationary Distribution")
% end

end
