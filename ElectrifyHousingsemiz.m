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
p5=6; % model period, in years (just used this to modify some parameters from annual to model period)

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
n_d=[2,5]; % Decisions: installpv, buyhouse (note, SemiExoStateFn hardcodes that buyhouse is 5 points, and it must be last)
n_a=[13,3,5]; % Endogenous asset, housing, and solarpv (0-40kW generation) 
n_semiz=[5,5,ceil(30/p5),3]; % Semi-exog: house prices before/after purchase, years since purchase (one minus this is the 30y duration of mortgages in model periods), and downpayment
n_z=7; % Exogenous labor productivity units shock (maybe later also solarpv productivity shock?)
% n_u=p5; % Between period i.i.d. shock (is this really p5 or is it 5 to match with SemiExoStateFn buyhouse?)
N_j=Params.J; % Number of periods in finite horizon

% LifeCycleModel35 had risky assets, but we delete that in this example
% vfoptions.riskyasset=1; % riskyasset aprime(d,u)
% simoptions.riskyasset=1;
% When there is more than one endogenous state, the riskyasset is the last one


% % Specify Epstein-Zin preferences
% vfoptions.exoticpreferences='EpsteinZin';
% vfoptions.EZpositiveutility=0; % Epstein-Zin preferences in utility-units have to be handled differently depending on whether the utility funciton is positive or negative valued (this is all done internally, you just need to use vfoptions to specify which)
% vfoptions.EZriskaversion='phi'; % additional risk-aversion
% % Params.phi is set below

%% To speed up the use of riskyasset we use 'refine_d', which requires us to set the decision variables in a specific order
% NOTE: semiz adds an n_d4, which are in ReturnFn but not in aprimeFn, and which determine semi-exogenous transitions
vfoptions.refine_d=[0,0,1,1]; % tell the code how many d1, d2, d3 and d4 there are
% Idea is to distinguish three categories of decision variable:
%  d1: decision is in the ReturnFn but not in aprimeFn
%  d2: decision is in the aprimeFn but not in ReturnFn
%  d3: decision is in both ReturnFn and in aprimeFn (installpv, an experienceasset)
%  d4: decision is in the ReturnFn but not in aprimeFn, and is in semiz (buyhouse)
% Note: ReturnFn must use inputs (d1,d3,d4,..) 
%       aprimeFn must use inputs (d2,d3,..)
% n_d must be set up as n_d=[n_d1, n_d2, n_d3, n_d4]
% d_grid must be set up as d_grid=[d1_grid; d2_grid; d3_grid; d4_grid];
% It is possible to solve models without any d1, as is the case here.
simoptions.refine_d=vfoptions.refine_d;

vfoptions.gridinterplayer=0;
vfoptions.ngridinterp=5;
simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;
%% Parameters

% Housing
Params.f_htc=0.1; % transaction cost of buying/selling house (is a percent of h+prime)
% Params.minhouse % set below, is the minimum value of house that can be purchased
Params.rentprice=0.3; % I figured setting rent a decent fraction of income is sensible
Params.houseservices=0.3; % housing services as a fraction of house value
Params.pv_pct_cost=0.05; % modeling a $30K install for a $600K house
Params.energy_pct_cost=0.07; % Electricity: 3%; Gas: 2%; Petrol: 2%

% Discount rate
Params.beta = 0.96^p5;
% Preferences
Params.sigma=10; % Coeff of relative risk aversion (curvature of consumption)
Params.phi=10; % Additional risk aversion (from Epstein-Zin preferences)
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

% Age-dependent labor productivity units
if Params.Jr>5
    Params.kappa_j=[linspace(0.5,2,Params.Jr-3),linspace(2,1,2),zeros(1,Params.J-Params.Jr+1)];
else
    Params.kappa_j=[linspace(0.5,2,Params.Jr-2),ones(1,1),zeros(1,Params.J-Params.Jr+1)];
end
% Exogenous shock process: AR1 on labor productivity units
Params.rho_z=0.9;
Params.sigma_epsilon_z=0.03;

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


%% House prices
Params.probhousepricerise=0.3; % increase one grid point
Params.probhousepricefall=0.2; % decrease one grid point
% remaining 1-probhousepricerise-probhousepricefall probability that house
% price is unchanged from previous period

% The grids on house prices (pbefore_grid and pafter_grid are below).

%% Grids
% The ^3 means that there are more points near 0 than near 1. We know from
% theory that the value function will be more 'curved' near zero assets,
% and putting more points near curvature (where the derivative changes the most) increases accuracy of results.
asset_grid=10*(linspace(0,1,n_a(1)))'; % Note, I use equal spacing (normally would put most points near zero)
% note: will go from 0 to 10
% assetprime_grid=10*(linspace(0,1,n_d(2)))'; % Want to let n_d(2) have different number of grid points from n_a(1).

% age20avgincome=Params.w*Params.kappa_j(1);
% house_grid=[0; logspace(2*age20avgincome, 12*age20avgincome, 5)'];
house_grid=(0:1:n_a(2)-1)';
% Note, we can see from w*kappa_j*z and the values of these, that average
% income is going to be around one, so will just use this simpler house grid
% [We can think about the values of the house_grid as being relative the average income (or specifically average at a given age)]
Params.minhouse=house_grid(2); % first is zero (no house)

% kWh of solar generation installed, 10kW per grid element
solarpv_grid=10*(0:1:n_a(3)-1)';

% First, the AR(1) process z
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,Params.rho_z,Params.sigma_epsilon_z,n_z);
z_grid=exp(z_grid); % Take exponential of the grid
[mean_z,~,~,~]=MarkovChainMoments(z_grid,pi_z); % Calculate the mean of the grid so as can normalise it
z_grid=z_grid./mean_z; % Normalise the grid on z (so that the mean of z is exactly 1)

% Share of assets invested in the risky asset
% riskyshare_grid=linspace(0,1,n_d(x))'; % Share of assets, from 0 to 1

% buyhouse
buyhouse_grid=(0:1:n_d(2)-1)';

% installpv is a binary choice
installpv_grid=[0; 1];

% Set up d for VFI Toolkit (is the two decision variables)
d_grid=[buyhouse_grid; installpv_grid];

a_grid=[asset_grid; house_grid; solarpv_grid];

% Now the semi-exogenous states, we define SemiExoStateFn later, for now just some grids
pbefore_grid=[0.8,1,1.2,1.4,1.6]'; % 1 represents price when agent is 'born'
pafter_grid=[0.8,1,1.2,1.4,1.6]'; % 1 represents price when house is purchased
% Note: is purely coincidence that pbefore and pafter use same grids (both
% must be equally spaced, but no need to be same values, nor same number of points)
yearsowned_grid=[(0:1:(n_semiz(3)-2))';100]; % note: 100 is an absorbing state representing 30+ years (so mortgage is fully repaid)
downpayment_grid=[0.2,0.4,0.6]'; % downpayment for new house must be 20%, 40%, 60%.
% Should the fact of solar PV installation affect pafter_grid?
semiz_grid=[pbefore_grid; pafter_grid; yearsowned_grid; downpayment_grid];
% Note, SemiExoStateFn hardcodes that the grid spacing for pbefore_grid
% must be evenly spaced, and same for pafter_grid.
Params.pbeforespacing=pbefore_grid(2)-pbefore_grid(1);
if any(abs(pbefore_grid(2:end)-pbefore_grid(1:end-1)-Params.pbeforespacing) > 1e-14)
    error('pbefore_grid must be evenly spaced (is hardcoded in SemiExoStateFn)')
end
Params.pafterspacing=pafter_grid(2)-pafter_grid(1);
if any(abs(pafter_grid(2:end)-pafter_grid(1:end-1)-Params.pafterspacing) > 1e-14)
    error('pafter_grid must be evenly spaced (is hardcoded in SemiExoStateFn)')
end
% need to store max/min of pbefore and pafter grids, so we can use them in
% SemiExoStateFn to avoid leaving the grid
Params.maxpbefore=max(pbefore_grid);
Params.minpbefore=min(pbefore_grid);
Params.maxpafter=max(pafter_grid);
Params.minpafter=min(pafter_grid);
% for initial agent distribution, which elements in pbefore and pafter_grid are the initial prices
Params.pbefore1=2; % second element is 1, which is where we want to start
Params.pafter1=2; % second element is 1, which is where we want to start

%% Solar PV is an experienceasset
vfoptions.experienceasset=1;
simoptions.experienceasset=1;

%% Define aprime function used for the riskyasset (value of next period assets, determined by this period decision, and u shock)
% This must all be adjusted if/when we add riskyassets back in
% riskyasset: aprime_val=aprimeFn(d,u)
% vfoptions.refine_d: the decision variables input to aprimeFn are d2,d3

% Experience assets must be listed first in aprime
a2primeFn=@(installpv, solarpv) ElectrifyHousing_a2primeFn(installpv, solarpv); % Will return the value of aprime
% Note that u is risky asset excess return and effectively includes both the (excess) mean and standard deviation of risky assets

%% Put the risky asset/experienceasset into vfoptions and simoptions
vfoptions.aprimeFn=a2primeFn;
% vfoptions.n_u=n_u;
% vfoptions.u_grid=u_grid;
% vfoptions.pi_u=pi_u;
simoptions.aprimeFn=a2primeFn;
% simoptions.n_u=n_u;
% simoptions.u_grid=u_grid;
% simoptions.pi_u=pi_u;
% Because a_grid and d_grid are involved in risky assets and experienceassets, but are not
% normally needed for agent distriubiton simulation, we have to also
% include these in simoptions
simoptions.a_grid=a_grid;
simoptions.d_grid=d_grid;

%% Setup for how the semi-exogneous states evolve

% Note: with riskyasset, the decision variables for the semi-exo states are determined by d4 in vftopions.refine_d
% Set up the semi-exogneous states
vfoptions.n_semiz=n_semiz;
vfoptions.semiz_grid=semiz_grid;
% Define the transition probabilities of the semi-exogenous states
vfoptions.SemiExoStateFn=@(pbefore,pafter,yearsowned,downpayment,pbeforeprime,pafterprime,yearsownedprime,downpaymentprime,buyhouse, ...
        probhousepricerise,probhousepricefall, pbeforespacing,pafterspacing, maxpbefore,minpbefore,maxpafter,minpafter, mortgageduration)...
    ElectrifyHousing_SemiExoStateFn(pbefore,pafter,yearsowned,downpayment,pbeforeprime,pafterprime,yearsownedprime,downpaymentprime,buyhouse, ...
        probhousepricerise,probhousepricefall, pbeforespacing,pafterspacing, maxpbefore,minpbefore,maxpafter,minpafter, mortgageduration);

% We also need to tell simoptions about the semi-exogenous states
simoptions.SemiExoStateFn=vfoptions.SemiExoStateFn;
simoptions.n_semiz=vfoptions.n_semiz;
simoptions.semiz_grid=vfoptions.semiz_grid;

%% Now, create the return function 
% % There is not much agreement on how to handle mortality risk with Epstein-Zin preferences
% % We can treat them as a risk
% vfoptions.survivalprobability='sj';
% DiscountFactorParamNames={'beta'};
% % Or we can treat them as a discount factor
DiscountFactorParamNames={'beta','sj'};

% Use 'ElectrifyHousing_ReturnFn'
ReturnFn=@(installpv,buyhouse,aprime,hprime,a,h,solarpv, ...
        pbefore,pafter,yearsowned,olddownpayment, z, ...
        w,r,sigma,agej,Jr,pension,kappa_j,sigma_h,f_htc,minhouse,rentprice,houseservices,mortgageduration,pv_pct_cost,energy_pct_cost) ...
    ElectrifyHousing_ReturnFn(installpv,buyhouse,aprime,hprime,a,h,solarpv, ...
        pbefore,pafter,yearsowned,olddownpayment, z, ...
        w,r,sigma,agej,Jr,pension,kappa_j,sigma_h,f_htc,minhouse,rentprice,houseservices,mortgageduration,pv_pct_cost,energy_pct_cost);
% vfoptions.refine_d, with semiz: only (d1,d3,d4,..) are input to ReturnFn [this model has no d1, so here just d3,d4]

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium
disp('Solve ValueFnIter')
vfoptions.verbose=1;
vfoptions.lowmemory=1;
tic;
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
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
jequaloneDist=zeros([n_a,n_semiz,n_z],'gpuArray'); % Put no households anywhere on grid
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
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
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
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);

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
