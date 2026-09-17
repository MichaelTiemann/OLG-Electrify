%% Rooftop Solar Life-Cycle Model (based on Life-Cycle Model 35: Portfolio-Choice with Housing)
% Modify Life-Cycle Model 35semiz, adding semi-exogenous shocks.
% Four semi-exo with experienceasset:
%  first is house price as markov before purchase
%  second is house price as markov after purchase
%  third is years since purchase, which is used for mortgages
%  fourth is the downpayment when house was purchased
%  investment in PV Solar is an experienceasset (could become
%  experienceassetu with shocks being improvements to solar technologies)

% As always, semi-exo goes after endogenous states, and before and markov
% (z) or i.i.d. (e) exogenous states.

% Semi-exogenous states evolves based on a decision variable. In this model
% we want them to change when you buy a house, or hold a house. Since house 
% is an endogenous state, we will add a decision variable, that must be one
% when buying a house, zero if don't buy, and two if hold (we can easily enforce this
% in the return function) [actually it takes more values, we change it so hold house
% is four, and make one-to-three be buy as explained below].
% The decision variable that is relevant to the semi-exogenous state is 
% assumed to be the 'last' decision variable. But here we are using 
% vfoptions.refine_d and so it additionally is assumed to be the 'd4' 
% decision variable. [refine_d with riskyasset has d1,d2,d3, when also
% using semiz there is also d4]

% Recall that LifeCycleModel10 showed how we could make labor (hours
% worked) an exogenous variable (z) and thus not a decision variable.

% The decision variable that determines semi-exo state transitions is
% called 'buyhouse' and takes three values: 0=don't own house, 1-to-3=buying a
% house this period, 4=own house. The three different values for buying a
% house related to the downpayment size, which is 20, 40, 60% of
% the price of the house.

% To be able to solve such a big problem, I switched to 5 year model period.
% Note that p5 must be at least 3 (for Farmer-Toda) so years-owned >= 2.
% p5 must be at most 15 (for kappa_j labor productivity evolutions).
p5=1; % model period, in years (just used this to modify some parameters from annual to model period)

%% How does VFI Toolkit think about this?
%
% Two decision variables: installpv (experience asset), buyhouse (semi-exo)
% Two endogenous state variables: a and h (assets and housing)
% One experienceasset (solarpv)
% Four semi-exogenous state variables: pbefore,pafter,yearsowned,olddownpayment
% One stochastic exogenous state variable: z, an AR(1) process (in logs), idiosyncratic shock to labor productivity units
% Age: j (which is actually a period number, since periods span multiple years)

%% Begin setting up to use VFI Toolkit to solve
% Lets model agents from age 20 to age 79, in five year periods (so first
% period is ages 20-24, and last period is ages 75-79 when p5==5.

Params.agejshifter=19; % Age 20 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=ceil((79-Params.agejshifter)/p5); % =60/p5, Number of period in life-cycle

% Grid sizes to use
% --- The 9.2 Million State Grid ---
n_d = [2, 5];             % Decisions: PV (2), BuyHouse (5)
n_a = [15, 4, 5];         % Endogenous: Assets (15), Housing (4 sizes), SolarPV (5 sizes)
n_semiz = [7, 7, 30, 3];  % Semi-exog: PBefore (7), PAfter (7), Mortgage Years (30), Downpayment (3)
n_z = 7;                  % Exogenous: Labor productivity (7)
N_j = Params.J;

% LifeCycleModel35 had risky assets, but we delete that in this example
% vfoptions.riskyasset=1; % riskyasset aprime(d,u)
% simoptions.riskyasset=1;
% When there is more than one endogenous state, the riskyasset is the last one

%% 'refine_d' requires us to set the decision variables in a specific order
vfoptions.refine_d=[0,1,1]; % tell the code how many d1, d2, d3 and d4 there are
% Idea is to distinguish three categories of decision variable:
%  d1: decision is in the ReturnFn but not in aprimeFn
%  d2: decision is in the aprimeFn but not in ReturnFn
%  d3: decision is in both ReturnFn and in aprimeFn (installpv, an experienceasset)
% Note: ReturnFn must use inputs (d1,d3,..) 
%       aprimeFn must use inputs (d2,d3,..)
% n_d must be set up as n_d=[n_d1, n_d2, n_d3]
% d_grid must be set up as d_grid=[d1_grid; d2_grid; d3_grid];
% It is possible to solve models without any d1, as is the case here.
simoptions.refine_d=vfoptions.refine_d;

vfoptions.gridinterplayer=[0,0,0];
vfoptions.ngridinterp=5;
simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;
%% Parameters

% Housing
Params.f_htc=0.1; % transaction cost of buying/selling house (is a percent of h+prime)
% Params.minhouse % set below, is the minimum value of house that can be purchased
Params.rentprice=0.3; % I figured setting rent a decent fraction of income is sensible
Params.houseservices=0.3; % housing services as a fraction of house value
Params.pv_pct_cost=0.03; % modeling a $20K install for a $600K house
Params.energy_pct_cost=0.07; % Electricity: 3%; Gas: 2%; Petrol: 2%

% Discount rate
Params.beta = 0.96^p5;
Params.beta0 = 0.80^p5; % <-- NEW: Quasi-Hyperbolic present-bias parameter

% Preferences (Core)
Params.sigma = 2.0; 
Params.eta   = 0.5;  
Params.phi   = 10;   % HOW DOES THIS GET INTO EZ?  Or is it Params.psi?  Or what?

% Preferences (Epstein-Zin)
Params.ez_risk_aversion = 4.0; % High Risk Aversion 
Params.ez_eis           = 0.5; % Elasticity of Intertemporal Substitution

%% Exotic Preferences: Quasi-Hyperbolic Epstein-Zin (QH-EZ)
vfoptions.exoticpreferences    = 'QHEpsteinZin';
vfoptions.quasi_hyperbolic     = 'Sophisticated'; 
vfoptions.QHadditionaldiscount = 'beta0';         

% Epstein-Zin Aggregation Settings
vfoptions.EZriskaversion    = 'ez_risk_aversion'; % <-- Pass the string name!
vfoptions.EZeis             = 'ez_eis';           % <-- Pass the string name!
vfoptions.EZpositiveutility = 0;   
vfoptions.EZutils           = 1;

% Preferences
Params.sigma=10; % Coeff of relative risk aversion (curvature of consumption)
Params.sigma_h=0.5; % Relative importance of housing services (vs consumption) in utility


% Prices
Params.w=1; % Wage

% Asset returns
Params.r=(1.25^p5)-1; % Rate of return on risk free asset
% u is the stochastic component of the excess returns to the risky asset
% Params.rp=(1.03^p5)-1; % Mean excess returns to the risky asset (so the mean return of the risky asset will be r+rp)
% Params.sigma_u=0.025; % Standard deviation of innovations to the risky asset
% Params.rho_u=0; % Asset return risk component is modeled as iid (if you regressed, e.g., the percent change in S&P500 on it's one year lag you get a coefficient of essentially zero)
% [u_grid, pi_u]=discretizeAR1_FarmerToda(Params.rp,Params.rho_u,Params.sigma_u,n_u);
% pi_u=pi_u(1,:)'; % This is iid

% Demographics
Params.agej=1:1:Params.J; % Is a vector of all the agej periods: 1,2,3,...,J
Params.Jr=round((65-Params.agejshifter)/p5); % Age 65 (period 10 is ages 65-69 in the 5 year case)

% Pensions
Params.pension=0.4; % Increased to be greater than rental costs

% Age-dependent labor productivity units (Smoothed for annual 60-period life)
working_years = Params.Jr - 1; 
Params.kappa_j = [linspace(0.5, 2.0, working_years - 10), ...
                  linspace(2.0, 1.0, 10), ...
                  zeros(1, Params.J - working_years)];

% Annualized Exogenous shock process: AR1 on labor productivity units
Params.rho_z = 0.97;              % Increased persistence for annual wage shocks
Params.sigma_epsilon_z = 0.015;   % Lower annual variance

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
dj=resize(dj,101+p5,FillValue=1);
% dj covers Ages 0-100, plus extras at end to make it period-friendly
Params.sj=prod(1-reshape(dj(1:p5*ceil(101/p5)),[p5,ceil(101/p5)]),1); % p5-year survival rates
Params.sj=Params.sj(1+ceil(20/p5):ceil(20/p5)+N_j); % Just the ages we are using (20yo and up)
Params.sj(end)=0; % In the present model the last period (j=J) value of sj is actually irrelevant

%% Mortgages

% In periods 0-5 you make a mortgage repayment. Year '100' is an absorbing
% state, which indicates mortgage has been paid off. This is tracked by the
% 'yearsowned' semi-exo state. Note that the following param needs to align
% with the grid on yearsowned (here both are set for 30 year mortgages; 6 periods of 5 years per period).
Params.mortgageduration=n_semiz(3)-1;


%% House prices (Annualized Drift)
Params.probhousepricerise = 0.08; % 8% chance house prices rise a tier this year
Params.probhousepricefall = 0.08; % 8% chance house prices fall a tier this year
% remaining 84% probability that house price is unchanged from last year
% remaining 1-probhousepricerise-probhousepricefall probability that house
% price is unchanged from previous period

% The grids on house prices (pbefore_grid and pafter_grid are below).

%% Grids
vfoptions.precision='single'; simoptions.precision=vfoptions.precision;
zero=cast(0,vfoptions.precision);

% The ^3 means that there are more points near 0 than near 1. We know from
% theory that the value function will be more 'curved' near zero assets,
% and putting more points near curvature (where the derivative changes the most) increases accuracy of results.
asset_grid=10*(linspace(zero,1,n_a(1)))'; % Note, I use equal spacing (normally would put most points near zero)
% note: will go from 0 to 10
% assetprime_grid=10*(linspace(0,1,n_d(2)))'; % Want to let n_d(2) have different number of grid points from n_a(1).

% age20avgincome=Params.w*Params.kappa_j(1);
% house_grid=[0; logspace(2*age20avgincome, 12*age20avgincome, 5)'];
house_grid=(zero:1:n_a(2)-1)';
% Note, we can see from w*kappa_j*z and the values of these, that average
% income is going to be around one, so will just use this simpler house grid
% [We can think about the values of the house_grid as being relative the average income (or specifically average at a given age)]
Params.minhouse=house_grid(2); % first is zero (no house)

% kWh of solar generation installed, 10kW per grid element
solarpv_grid=10*(zero:1:n_a(3)-1)';

% First, the AR(1) process z
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,Params.rho_z,Params.sigma_epsilon_z,n_z);
z_grid=exp(z_grid); % Take exponential of the grid
[mean_z,~,~,~]=MarkovChainMoments(z_grid,pi_z); % Calculate the mean of the grid so as can normalise it
z_grid=z_grid./mean_z; % Normalise the grid on z (so that the mean of z is exactly 1)

% Share of assets invested in the risky asset
% riskyshare_grid=linspace(0,1,n_d(x))'; % Share of assets, from 0 to 1

% buyhouse
buyhouse_grid=(zero:1:n_d(2)-1)';

% installpv is a binary choice
installpv_grid=cast([0; 1], vfoptions.precision);

% Set up d for VFI Toolkit (is the two decision variables)
d_grid=[installpv_grid; buyhouse_grid];

a_grid=[asset_grid; house_grid; solarpv_grid];

% Expanded 7-point House Price Market (evenly spaced by 0.15)
% [0.70, 0.85, 1.00, 1.15, 1.30, 1.45, 1.60]
pbefore_grid = cast(0.70 : 0.15 : 1.60, vfoptions.precision)'; 
pafter_grid  = cast(0.70 : 0.15 : 1.60, vfoptions.precision)'; 

yearsowned_grid = [(zero : 1 : (n_semiz(3) - 2))'; 100]; 
downpayment_grid = cast([0.2, 0.4, 0.6], vfoptions.precision)'; 

semiz_grid = [pbefore_grid; pafter_grid; yearsowned_grid; downpayment_grid];

% (Keep the spacing checks that follow here...)
% Note, SemiExoStateFn hardcodes that the grid spacing for pbefore_grid
% must be evenly spaced, and same for pafter_grid.
Params.pbeforespacing=pbefore_grid(2)-pbefore_grid(1);
if any(abs(pbefore_grid(2:end)-pbefore_grid(1:end-1)-Params.pbeforespacing) > 1e-7)
    error('pbefore_grid must be evenly spaced (is hardcoded in SemiExoStateFn)')
end
Params.pafterspacing=pafter_grid(2)-pafter_grid(1);
if any(abs(pafter_grid(2:end)-pafter_grid(1:end-1)-Params.pafterspacing) > 1e-7)
    error('pafter_grid must be evenly spaced (is hardcoded in SemiExoStateFn)')
end
% need to store max/min of pbefore and pafter grids, so we can use them in
% SemiExoStateFn to avoid leaving the grid
Params.maxpbefore=max(pbefore_grid);
Params.minpbefore=min(pbefore_grid);
Params.maxpafter=max(pafter_grid);
Params.minpafter=min(pafter_grid);
% For initial agent distribution: 1.00 is the 3rd element in our new 0.70:0.15:1.60 grid
Params.pbefore1 = 3; 
Params.pafter1 = 3;

%% Solar PV is an experienceasset
vfoptions.experienceasset=1;
simoptions.experienceasset=1;

%% Define aprime function used for the riskyasset (value of next period assets, determined by this period decision, and u shock)
% This must all be adjusted if/when we add riskyassets back in
% riskyasset: aprime_val=aprimeFn(d,u)
% vfoptions.refine_d: the decision variables input to aprimeFn are d2,d3

% Experience assets must be listed first in aprime
if strcmp(vfoptions.precision, 'single')
    a2primeFn=@(installpv, solarpv, pbefore, pafter, yearsowned, olddownpayment) ElectrifyHousingV_a2primeFn_single(installpv, solarpv); % Will return the value of aprime
else
    a2primeFn=@(installpv, solarpv, pbefore, pafter, yearsowned, olddownpayment) ElectrifyHousingV_a2primeFn(installpv, solarpv); % Will return the value of aprime
end
% Note that u is risky asset excess return and effectively includes both the (excess) mean and standard deviation of risky assets

%% Put the risky asset/experienceasset into vfoptions and simoptions
vfoptions.aprimeFn=a2primeFn;
% vfoptions.n_u=n_u;
% vfoptions.u_grid=u_grid;
% vfoptions.pi_u=pi_u;
simoptions.aprimeFn=vfoptions.aprimeFn;
% simoptions.n_u=n_u;
% simoptions.u_grid=u_grid;
% simoptions.pi_u=pi_u;
% Because a_grid and d_grid are involved in risky assets and experienceassets, but are not
% normally needed for agent distribution simulation, we have to also
% include these in simoptions
% And we need to include z_grid to support later semiz bypass logic
simoptions.a_grid=a_grid;
simoptions.d_grid=d_grid;
simoptions.z_grid=z_grid;

%% Setup for how the semi-exogenous states evolve
vfoptions.l_dsemiz = 1; % or 2 depending on how many decision variables control semiz
vfoptions.n_semiz = n_semiz;
vfoptions.semiz_grid = semiz_grid;

% Define the transition probabilities function handle
vfoptions.SemiExoStateFn = @(pbefore,pafter,yearsowned,downpayment,pbeforeprime,pafterprime,yearsownedprime,downpaymentprime,buyhouse, ...
    probhousepricerise,probhousepricefall,pbeforespacing,pafterspacing,maxpbefore,minpbefore,maxpafter,minpafter,mortgageduration)...
    ElectrifyHousing_SemiExoStateFn(pbefore,pafter,yearsowned,downpayment,pbeforeprime,pafterprime,yearsownedprime,downpaymentprime,buyhouse, ...
    probhousepricerise,probhousepricefall,pbeforespacing,pafterspacing,maxpbefore,minpbefore,maxpafter,minpafter,mortgageduration);

%% Initialize time-invariant Semi-Exogenous transition tensor
disp('Initializing memory-optimized Semi-Exogenous transition tensor...');

% Unlock the setup temporarily
vfoptions.alreadygridvals_semiexo = 0;

% Restrict generation to 2 periods to allocate a single time-invariant transition 
% tensor, bypassing dense time-varying allocation limits for static housing grids.
vfoptions_temp = SemiExogShockSetup_FHorz(n_d, 2, d_grid, Params, vfoptions, 3);

% Extract the efficiently generated grids and transition matrix
vfoptions.semiz_gridvals_J = vfoptions_temp.semiz_gridvals_J;
vfoptions.pi_semiz_J = vfoptions_temp.pi_semiz_J;

% Lock grids to prevent time-varying expansion during ValueFnIter
vfoptions.alreadygridvals_semiexo = 1;

% We also need to tell simoptions about the semi-exogenous states
simoptions.SemiExoStateFn = vfoptions.SemiExoStateFn;
simoptions.n_semiz = vfoptions.n_semiz;
simoptions.semiz_grid = vfoptions.semiz_grid;
simoptions.l_dsemiz = vfoptions.l_dsemiz;
simoptions.semiz_gridvals_J = vfoptions.semiz_gridvals_J;
simoptions.pi_semiz_J = vfoptions.pi_semiz_J;
simoptions.alreadygridvals_semiexo = 1;

%% Now, create the return function 
% % There is not much agreement on how to handle mortality risk with Epstein-Zin preferences
% % We can treat them as a risk
% vfoptions.survivalprobability='sj';
% DiscountFactorParamNames={'beta'};
% % Or we can treat them as a discount factor
DiscountFactorParamNames={'beta','sj'};

% Use 'ElectrifyHousing_ReturnFn'
if strcmp(vfoptions.precision,'single')
    ReturnFn=@(installpv,buyhouse,aprime,hprime,a,h,solarpv, ...
        pbefore,pafter,yearsowned,olddownpayment, z, ...
        w,r,sigma,agej,Jr,pension,kappa_j,sigma_h,f_htc,minhouse,rentprice,houseservices,mortgageduration,pv_pct_cost,energy_pct_cost) ...
        ElectrifyHousingsemizV_ReturnFn_single(installpv,buyhouse,aprime,hprime,a,h,solarpv, ...
        pbefore,pafter,yearsowned,olddownpayment, z, ...
        w,r,sigma,agej,Jr,pension,kappa_j,sigma_h,f_htc,minhouse,rentprice,houseservices,mortgageduration,pv_pct_cost,energy_pct_cost);
else
    ReturnFn=@(installpv,buyhouse,aprime,hprime,a,h,solarpv, ...
            pbefore,pafter,yearsowned,olddownpayment, z, ...
            w,r,sigma,agej,Jr,pension,kappa_j,sigma_h,f_htc,minhouse,rentprice,houseservices,mortgageduration,pv_pct_cost,energy_pct_cost) ...
        ElectrifyHousingsemizV_ReturnFn(installpv,buyhouse,aprime,hprime,a,h,solarpv, ...
            pbefore,pafter,yearsowned,olddownpayment, z, ...
            w,r,sigma,agej,Jr,pension,kappa_j,sigma_h,f_htc,minhouse,rentprice,houseservices,mortgageduration,pv_pct_cost,energy_pct_cost);
end
% vfoptions.refine_d, with semiz: only (d1,d3,..) are input to ReturnFn [this model has no d1, so here just d3]

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium
disp('Solve ValueFnIter')
vfoptions.verbose=1;
vfoptions.lowmemory=1;
tic;
[V, Policy]=ValueFnIter_Case1_VFHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
toc

% V is now (a,z,j). This was already true, just that previously z was trivial (a single point) 
% Compare
size(V)
% with
[n_a,n_z,N_j]
% there are the same.
% Policy is
size(Policy)
% which is the same as
[length(n_d)+1,n_a,n_z,N_j]
% The n_a,n_z,N_j represent the state on which the decisions/policys
% depend, and there is one decision for each decision variable 'd' plus one
% more for the standard asset

%% Now, we want to graph Life-Cycle Profiles

%% Initial distribution of agents at birth (j=1)
% Before we plot the life-cycle profiles we have to define how agents are
% at age j=1. We will give them all zero assets.
jequaloneDist=zeros([n_a,n_semiz,n_z],vfoptions.precision,'gpuArray'); % Put no households anywhere on grid
jequaloneDist(1,1,1,Params.pbefore1,Params.pafter1,1,1,ceil(n_z/2))=1; 
% All agents start with zero assets, no house, zero solarpv
% note: yearsowned=0 and downpayment=0.2 initial values are anyway irrelevant
% and the median z shock

%% We now compute the 'stationary distribution' of households
% Start with a mass of one at initial age, use the conditional survival
% probabilities sj to calculate the mass of those who survive to next
% period, repeat. Once done for all ages, normalize to one
Params.mewj=ones(1,Params.J); % Marginal distribution of households over age
for jj=2:length(Params.mewj)
    Params.mewj(jj)=Params.sj(jj-1)*Params.mewj(jj-1);
end
Params.mewj=Params.mewj./sum(Params.mewj); % Normalize to one
AgeWeightsParamNames={'mewj'}; % So VFI Toolkit knows which parameter is the mass of agents of each age

StationaryDist=StationaryDist_VFHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
% riskyasset requires the grids when simulating the agent distribution to be able to handle aprime(d,u). The grids are passed in simoptions.


%% FnsToEvaluate are how we say what we want to graph the life-cycles of
% Takes all the d, then relevant aprime, then a, then semiz, then z
% FnsToEvaluate.riskyshare=@(savings,riskyshare,buyhouse,hprime,a,h,solarpv,pbefore,pafter,yearsowned,olddownpayment,z) riskyshare; % riskyshare, is the fraction of savings invested in the risky asset
FnsToEvaluate.earnings=@(installpv,buyhouse,aprime,hprime,a,h,solarpv,pbefore,pafter,yearsowned,olddownpayment,z,w,kappa_j) w*kappa_j*z + solarpv; % labor earnings + PV generation
FnsToEvaluate.assets=@(installpv,buyhouse,aprime,hprime,a,h,solarpv,pbefore,pafter,yearsowned,olddownpayment,z) a; % a is the current asset holdings
FnsToEvaluate.housing=@(installpv,buyhouse,aprime,hprime,a,h,solarpv,pbefore,pafter,yearsowned,olddownpayment,z) h; % h is housing holdings
FnsToEvaluate.solarpv=@(installpv,buyhouse,aprime,hprime,a,h,solarpv,pbefore,pafter,yearsowned,olddownpayment,z) solarpv; % solarpv is the current PV capacity


FnsToEvaluate.buyhouse=@(installpv,buyhouse,aprime,hprime,a,h,solarpv,pbefore,pafter,yearsowned,olddownpayment,z) buyhouse; % h is housing holdings
FnsToEvaluate.pbefore=@(installpv,buyhouse,aprime,hprime,a,h,solarpv,pbefore,pafter,yearsowned,olddownpayment,z) pbefore; % h is housing holdings
FnsToEvaluate.pafter=@(installpv,buyhouse,aprime,hprime,a,h,solarpv,pbefore,pafter,yearsowned,olddownpayment,z) pafter; % h is housing holdings
FnsToEvaluate.olddownpayment=@(installpv,buyhouse,aprime,hprime,a,h,solarpv,pbefore,pafter,yearsowned,olddownpayment,z) olddownpayment; % h is housing holdings
FnsToEvaluate.yearsowned=@(installpv,buyhouse,aprime,hprime,a,h,solarpv,pbefore,pafter,yearsowned,olddownpayment,z) yearsowned; % yearsowned, note, goes a bit silly due to the 100 being 5+ 
% FnsToEvaluate.yearsowned=@(installpv,buyhouse,aprime,hprime,a,h,solarpv,pbefore,pafter,yearsowned,olddownpayment,z) yearsowned*(yearsowned~=100); % yearsowned, note, goes a bit silly due to the 100 being 5+ 

% notice that we have called these earnings and assets

%% Calculate the life-cycle profiles
%% Calculate the life-cycle profiles
AgeConditionalStats = LifeCycleProfiles_VFHorz_ExpAssetsemiz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);

% For example
% AgeConditionalStats.earnings.Mean
% There are things other than Mean, but in our current deterministic model
% in which all agents are born identical the rest are meaningless.

%% Plot the life cycle profiles of fraction-of-time-worked, earnings, and assets

figure(5)
subplot(4,1,1); plot(1:1:Params.J,AgeConditionalStats.solarpv.Mean)
title('Life Cycle Profile: SolarPV installation (solarpv)')
subplot(4,1,2); plot(1:1:Params.J,AgeConditionalStats.earnings.Mean)
title('Life Cycle Profile: Labor Earnings (w kappa_j z)')
subplot(4,1,3); plot(1:1:Params.J,AgeConditionalStats.assets.Mean)
title('Life Cycle Profile: Assets (a)')
subplot(4,1,4); plot(1:1:Params.J,AgeConditionalStats.housing.Mean)
title('Life Cycle Profile: Housing (h)')

%%
figure(6)
subplot(5,1,1); plot(1:1:Params.J,AgeConditionalStats.buyhouse.Mean)
title('Life Cycle Profile: 0=no house, 4=hold house, 1/2/3 are all buying and reflect downpayment when buying (buyhouse)')
subplot(5,1,2); plot(1:1:Params.J,AgeConditionalStats.pbefore.Mean)
title('Life Cycle Profile: Price of house at purchase (pbefore)')
subplot(5,1,3); plot(1:1:Params.J,AgeConditionalStats.pafter.Mean)
title('Life Cycle Profile: House price now relative to purchase (pafter)')
subplot(5,1,4); plot(1:1:Params.J,AgeConditionalStats.olddownpayment.Mean)
title('Life Cycle Profile: Down-payment (lag) (olddownpayment)')
subplot(5,1,5); plot(1:1:Params.J,AgeConditionalStats.yearsowned.Mean)
title('Life Cycle Profile: Years owned current house (yearsowned)')
