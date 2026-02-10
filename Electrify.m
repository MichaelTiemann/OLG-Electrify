%% OLG Electrification (based on OLGModels14: Heterogenous households and heterogeneous firms)
% See https://www.vfitoolkit.com/updates-blog/2021/an-introduction-to-life-cycle-models/
% and OLGModel14.m in the repo https://github.com/vfitoolkit/IntroToOLGModels

Names_i={'household','firm'};
PTypeDistParamNames={'ptypemass'};
Params.ptypemass=[1,1]; % Mass of households and firms are each equal to one

% A line some need for running on the Server
addpath(genpath('./MatlabToolkits/'))

%% Begin setting up to use VFI Toolkit to solve
vfoptions.lowmemory.household=2;
vfoptions.lowmemory.firm=0;
vfoptions.tolerance=10^(-1);
simoptions.tolerance=vfoptions.tolerance;

% The user can experiment with gridinterplayer=0 (pure discretization) or gridinterplayer=1 (linear interpolation b/w grid points).
% If gridinterplayer=1, then you must set vfoptions.divideandconquer=1 (required for transition).
vfoptions.gridinterplayer  = 0;
vfoptions.ngridinterp      = 10;
vfoptions.divideandconquer.household = 0;
vfoptions.divideandconquer.firm = 0;
vfoptions.level1n=11;
simoptions.gridinterplayer = vfoptions.gridinterplayer;
simoptions.ngridinterp     = vfoptions.ngridinterp;

% Grid sizes to use for household

% Lets model agents from age 20 to age 100, so 81 periods
Params.agejshifter=19; % Age 20 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=100-Params.agejshifter; % Number of period in life-cycle
n_d.household=[17,5]; % Decisions: labor, buyhouse
n_a.household=[17,8,5,4]; % Endogenous shares, assets, housing, and solarpv assets (0-45 kW generation)
% Exogenous labor productivity units shocks (next two lines)
n_z.household=7; % AR(1) with age-dependent params
vfoptions.n_e.household=3; % iid
N_j.household=Params.J; % Number of periods in finite horizon

% Grids to use for firm
n_d.firm=51; % Dividend payment
n_a.firm=101; % Capital holdings
n_z.firm=11; % Productivity shock
N_j.firm=Inf; % Infinite horizon

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

%% Parameters for household

% Housing
Params.f_htc=0.03; % transaction cost of buying/selling house (is a percent of h+hprime)
% Params.minhouse % set below, is the minimum value of house that can be purchased
Params.rentprice=0.3; % I figured setting rent a decent fraction of income is sensible
Params.f_coll=0.5; % collateral contraint (fraction of house value that can be borrowed)
Params.houseservices=0.5; % housing services as a fraction of house value
Params.pv_pct_cost=0.05; % modeling a $30K install for a $600K house
Params.energy_pct_cost=0.07; % Electricity: 3%; Gas: 2%; Petrol: 2%

% Discount rate
Params.beta = 1.016; % Changed to get S to increase nearer to 1 given r=0.05 (ran it with beta=0.99, got S=0.3, so increased this; note that it interacts with sj to give the actual discount factor)
% Preferences
Params.sigma = 2; % Coeff of relative risk aversion (curvature of consumption)
Params.sigma_h=0.5; % Relative importance of housing services (vs consumption) in utility
Params.eta = 1.5; % Curvature of leisure (This will end up being 1/Frisch elasty)
Params.psi = 0.5; % Weight on leisure

% Asset returns
Params.r_wedge=0.02; % Rate difference between deposits and loans

% Demographics
Params.agej=1:1:Params.J; % Is a vector of all the periods: 1,2,3,...,J
Params.Jr=65-Params.agejshifter; % Period number of retirement
% Population growth rate
Params.n=0.02; % percentage rate (expressed as fraction) of population growth per period

% Age-dependent labor productivity units
Params.kappa_j=[linspace(0.5,2,Params.Jr-15),linspace(2,1,14),zeros(1,Params.J-Params.Jr+1)];

% Life-cycle AR(1) process z, on (log) labor productivity units
% Chosen following Karahan & Ozkan (2013) [as used by Fella, Gallipoli & Pan (2019)]
% Note that 37 covers 24 to 60 inclusive (as in the original)
% Now repeat the first and last values to fill in working age, and put zeros for retirement (where it is anyway irrelevant)
rho_z=0.7596+0.2039*((1:1:37)/10)-0.0535*((1:1:37)/10).^2+0.0028*((1:1:37)/10).^3; % Chosen following Karahan & Ozkan (2013) [as used by Fella, Gallipoli & Pan (2019)]
Params.rho_z=[rho_z(1)*ones(1,4),rho_z,rho_z(end)*ones(1,4),zeros(1,101-65)];
sigma_epsilon_z=0.0518-0.0405*((1:1:37)/10)+0.0105*((1:1:37)/10).^2-0.0002*((1:1:37)/10).^3; % Chosen following Karahan & Ozkan (2013) [as used by Fella, Gallipoli & Pan (2019)]
Params.sigma_epsilon_z=[sigma_epsilon_z(1)*ones(1,4),sigma_epsilon_z,sigma_epsilon_z(end)*ones(1,4),sigma_epsilon_z(end)*ones(1,101-65)];

% Transitory iid shock
sigma_e=0.0410+0.0221*((24:1:60)/10)-0.0069*((24:1:60)/10).^2+0.0008*((24:1:60)/10).^3;
Params.sigma_e=[sigma_e(1)*ones(1,4),sigma_e,sigma_e(end)*ones(1,4),sigma_e(end)*ones(1,101-65)];

% Note: These iid shocks will interact with the endogenous labor so the final labor
% earnings process will not equal that of Karahan & Ozkan (2013)
% Note: Karahan & Ozkan (2013) also have a fixed effect (which they call alpha) and which I ignore here.

% Age-dependent universal, permanent supply-side shocks to housing, energy, and pv installation costs
shock_period=10;
shock_years=(1+shock_period:shock_period:80);
% Exponentially increasing every shock_period years from from ~2% to ~7% after initial shock-free period
shock_pct=exp(shock_years/60)/50;
Params.agej_pct_cost=[zeros(1,1+shock_period), repelem(shock_pct,shock_period)];

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
% dj covers Ages 0 to 100
Params.sj=1-dj(21:101); % Conditional survival probabilities
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
Params.delta=0.05; % Depreciation of physical capital
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
% Params.firmbeta=1/(1+Params.r/(1-Params.tau_cg)); % 1/(1+r) but returns net of capital gains tax

%% Also, we want to target a capital-output ratio
% This is relevant to the general equilibrium conditions
Params.TargetKdivL=2.03;

%% Remaining parameters

% Some initial values/guesses for variables that will be determined in general eqm
Params.pension=0.4; % Initial guess (this will be determined in general eqm)
Params.r=0.05;
Params.w=1;
Params.AccidentBeqS= 0.03; % Accidental bequests (this is the lump sum transfer of shares)
Params.AccidentBeqAH= 0.03; % Accidental bequests (this is the lump sum transfer of assets+house value)
Params.G=0.1; % Government expenditure
Params.firmbeta=1/(1+Params.r/(1-Params.tau_cg)); % 1/(1+r) but returns net of capital gains tax
Params.D=0.2; % Dividends
Params.P0=1;
Params.Lhscale=0.62; % Scaling the household labor supply

%% Grids for household

% Grid for labour choice
labor_grid=linspace(0,1,n_d.household(1))'; % Notice that it is imposing the 0<=h<=1 condition implicitly

% Grid for share holdings, always > 0
share_grid=10*(linspace(0,1,n_a.household(1)).^3)';

% Grid for bank account; a negative balance implies a mortgage
asset_grid=(-1+3*(linspace(0,1,n_a.household(2))).^3)';

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

buyhouse_grid=(0:1:n_d.household(2)-1)';

% kW of solar generation installed, 10 kW per grid element
solarpv_grid=(0:1:n_a.household(4)-1)';

% Set up d for VFI Toolkit (is the two decision variables)
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
a_grid.firm=10*(linspace(0,1,n_a.firm).^3)'; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.

[z_grid.firm,pi_z.firm] = discretizeAR1_FarmerToda(0,Params.rho_z_firm,Params.sigma_z_e_firm,n_z.firm);
z_grid.firm=exp(z_grid.firm);


%% Now, create the return function

% For households
DiscountFactorParamNames.household={'beta','sj'};
% Notice we use 'Electrify_HouseholdReturnFn'
ReturnFn.household=@( ...
        labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
        r,pension,AccidentBeqS,AccidentBeqAH,w,P0,D,Lhscale, ...
        sigma,psi,eta,sigma_h,kappa_j,tau_l,tau_d,tau_cg,warmglow1,warmglow2,agej,Jr,J,...
        r_wedge,f_htc,minhouse,rentprice,f_coll,houseservices,agej_pct_cost,pv_pct_cost,energy_pct_cost ...
    ) Electrify_HouseholdReturnFn( ...
        labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
        r,pension,AccidentBeqS,AccidentBeqAH,w,P0,D,Lhscale, ...
        sigma,psi,eta,sigma_h,kappa_j,tau_l,tau_d,tau_cg,warmglow1,warmglow2,agej,Jr,J,...
        r_wedge,f_htc,minhouse,rentprice,f_coll,houseservices,agej_pct_cost,pv_pct_cost,energy_pct_cost ...
    );

% For firms
DiscountFactorParamNames.firm={'firmbeta'};
% Notice we use 'Electrify_FirmReturnFn'
ReturnFn.firm=@( ...
        d,kprime,k,z, ...
        w, ...
        delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg ...
    ) Electrify_FirmReturnFn( ...
        d,kprime,k,z, ...
        w, ...
        delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg ...
    );

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium
disp('Test ValueFnIter')
tic;
% Note: z_grid and pi_z, this will be ignored due to presence of vfoptions.z_grid_J and vfoptions.pi_z_J
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
toc

%% Initial distribution of agents at birth (j=1)
% Before we plot the life-cycle profiles we have to define how agents are
% at age j=1. We will give them all zero assets, no house, no solarpv.
jequaloneDist.household=zeros([n_a.household,n_z.household,vfoptions.n_e.household],'gpuArray'); % Put no households anywhere on grid
jequaloneDist.household(1,zeroassetindex,1,1,floor((n_z.household+1)/2),floor((simoptions.n_e.household+1)/2))=1; % All agents start with zero shares, assets, houses, solarpv, and the median shock

% Note that because the firms are infinite horizon they do not have an age=1 distribution

%% Agents age distribution
% Many OLG models include some kind of population growth, and perhaps some other things that create a weighting of different ages that needs to
% be used to calculate the stationary distribution and aggregate variable.
mewj=ones(1,Params.J); % Marginal distribution of households over age
for jj=2:length(mewj)
    mewj(jj)=Params.sj(jj-1)*mewj(jj-1)/(1+Params.n);
end
Params.mewj.household=mewj./sum(mewj); % Normalize to one

AgeWeightsParamNames=struct('household',{{'mewj'}}); % So VFI Toolkit knows which parameter is the mass of agents of each age

%% Test
disp('Test StationaryDist')
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);

%% General eqm variables
GEPriceParamNames={'r','pension','AccidentBeqS','AccidentBeqAH','G','w','firmbeta','D','P0','Lhscale'}; 
% We don't need P
% We can get P from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have
% P=((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg);

%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)
% Note: we need to add z & e to FnsToEvaluate inputs for households,
% whereas firm only has z (it is just coincidence/lazy that I call them
% both z).

% Stationary Distribution Aggregates from households (important that ordering of Names and Functions is the same)
FnsToEvaluate.L_h.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,kappa_j,Lhscale) labor*kappa_j*exp(z+e)*Lhscale;  % Aggregate labour supply in efficiency units 
FnsToEvaluate.S.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) s; % Aggregate share holdings
FnsToEvaluate.A.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) a; % Aggregate asset/mortgage holdings
FnsToEvaluate.H.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) h; % Aggregate house holdings
FnsToEvaluate.PV.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) solarpv; % Aggregate solarpv holdings
FnsToEvaluate.PensionSpending.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,pension,agej,Jr) (agej>=Jr)*pension; % Total spending on pensions
FnsToEvaluate.PayrollTaxRevenue.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,agej,Jr,tau_l,w,kappa_j,Lhscale) (agej<Jr)*tau_l*labor*w*kappa_j*exp(z+e)*Lhscale; % Total spending on pensions
FnsToEvaluate.AccidentalBeqSLeft.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,sj) sprime*(1-sj); % Accidental share bequests left by people who die
FnsToEvaluate.AccidentalBeqAHLeft.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,sj) (aprime+hprime)*(1-sj); % Accidental asset+house bequests left by people who die
FnsToEvaluate.CapitalGainsTaxRevenue.household = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,tau_cg,P0,D,tau_d,r) tau_cg*(P0-(((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg)))*(s+min(a,0)); % tau_cg*(P0-Plag)*(s,min(a,0)), but substitute P=Plag, and then substitute for P
% From firms
FnsToEvaluate.Output.firm = @(d,kprime,k,z,w,alpha_k,alpha_l) z*(k^alpha_k)*((w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)))^alpha_l; % Production function z*(k^alpha_k)*(l^alpha_l) (substituting for l)
FnsToEvaluate.L_f.firm = @(d,kprime,k,z,w,alpha_k,alpha_l) (w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % (effective units of) labor demanded by firm
FnsToEvaluate.K.firm = @(d,kprime,k,z,w,alpha_k,alpha_l) k; % physical capital
FnsToEvaluate.DividendPaid.firm = @(d,kprime,k,z,w) d; % dividend paid by firm
FnsToEvaluate.Sissued.firm = @(d,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi) Electrify_FirmShareIssuance(d,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi); % Share issuance
FnsToEvaluate.CorpTaxRevenue.firm = @(d,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi) Electrify_FirmCorporateTaxRevenue(d,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi); % revenue from the corporate profits tax

% General Equilibrium conditions (these should evaluate to zero in general equilbrium)
GeneralEqmEqns.sharemarket = @(S) S-1; % mass of all shares equals one
GeneralEqmEqns.labormarket = @(L_h,L_f) L_h-L_f; % labor supply of households equals labor demand of firms
GeneralEqmEqns.pensions = @(PensionSpending,PayrollTaxRevenue) PensionSpending-PayrollTaxRevenue; % Retirement benefits equal Payroll tax revenue: pension*fractionretired-tau*w*H
GeneralEqmEqns.bequestsS = @(AccidentalBeqSLeft,AccidentBeqS,n) AccidentalBeqSLeft/(1+n)-AccidentBeqS; % Accidental share bequests received equal accidental share bequests left
GeneralEqmEqns.bequestsAH = @(AccidentalBeqAHLeft,AccidentBeqAH,n) AccidentalBeqAHLeft/(1+n)-AccidentBeqAH; % Accidental asset+house bequests received equal accidental asset+house bequests left
GeneralEqmEqns.govbudget = @(G,tau_d,D,CapitalGainsTaxRevenue,CorpTaxRevenue) G-tau_d*D-CapitalGainsTaxRevenue-CorpTaxRevenue; % G is equal to the target, GdivYtarget*Y
GeneralEqmEqns.firmdiscounting = @(firmbeta,r,tau_cg) firmbeta-1/(1+r/(1-tau_cg)); % Firms discount rate is related to market return rate
GeneralEqmEqns.dividends = @(D,DividendPaid) D-DividendPaid; % That the dividend households receive equals that which firms give
GeneralEqmEqns.ShareIssuance = @(Sissued,P0,D,tau_cg,tau_d,r) P0-((((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg))-Sissued); % P0=P-S, but substitute for P (see derivation inside the return fn)
GeneralEqmEqns.CapitalOutputRatio =@(K,L_f,TargetKdivL) K/L_f-TargetKdivL;

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
fprintf('Check: S, A, H, PV \n')
[AggVars.S.Mean,AggVars.A.Mean,AggVars.H.Mean,AggVars.PV.Mean]
fprintf('Check: ShareIssuance GE condition \n')
Params.P0-((((1-Params.tau_cg)*Params.P0 + (1-Params.tau_d)*Params.D)/(1+Params.r-Params.tau_cg))-AggVars.S.Mean)


%% Solve for the General Equilibrium
% heteroagentoptions.fminalgo=4 % CMA-ES algorithm 

heteroagentoptions.verbose=1;
heteroagentoptions.toleranceGEprices=10^(-2);
heteroagentoptions.toleranceGEcondns=10^(-1); % This is the hard one
heteroagentoptions.maxiter=300;

p_eqm=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, [], pi_z, d_grid, a_grid, z_grid,jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightsParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
% p_eqm contains the general equilibrium parameter values
% Put this into Params so we can calculate things about the initial equilibrium
Params.r=p_eqm.r;
Params.pension=p_eqm.pension;
Params.AccidentBeqS=p_eqm.AccidentBeqS;
Params.AccidentBeqAH=p_eqm.AccidentBeqAH;
Params.G=p_eqm.G;
Params.w=p_eqm.w;
Params.firmbeta=p_eqm.firmbeta;
Params.D=p_eqm.D;
Params.P0=p_eqm.P0;
Params.Lhscale=p_eqm.Lhscale;

% Calculate a few things related to the general equilibrium.
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j, Names_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
% Can just use the same FnsToEvaluate as before.
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

%% Plot the life cycle profiles of capital and labour for the inital and final eqm.

figure(1)
subplot(4,1,1); plot(1:1:Params.J,AgeConditionalStats.L_h.Mean)
title('Life Cycle Profile: Effective Labour Supply')
subplot(4,1,2); plot(1:1:Params.J,AgeConditionalStats.S.Mean)
title('Life Cycle Profile: Share holdings')
subplot(4,1,3); plot(1:1:Params.J,AgeConditionalStats.A.Mean)
title('Life Cycle Profile: Asset holdings')
subplot(4,1,4); plot(1:1:Params.J,AgeConditionalStats.H.Mean)
title('Life Cycle Profile: House holdings')
% saveas(figure_c,'./SavedOutput/Graphs/OLGModel6_LifeCycleProfiles','pdf')

%% Calculate some aggregates and print findings about them

% Add consumption to the FnsToEvaluate
FnsToEvaluate.Consumption.household=@( ...
        labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
        r,pension,AccidentBeqS,AccidentBeqAH,w,P0,D,Lhscale, ...
        kappa_j,tau_l,tau_d,tau_cg,agej,Jr,...
        r_wedge,f_htc,rentprice,agej_pct_cost,pv_pct_cost,energy_pct_cost ...
    ) Electrify_HouseholdConsumptionFn( ...
        labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
        r,pension,AccidentBeqS,AccidentBeqAH,w,P0,D,Lhscale, ...
        kappa_j,tau_l,tau_d,tau_cg,agej,Jr,...
        r_wedge,f_htc,rentprice,agej_pct_cost,pv_pct_cost,energy_pct_cost);
FnsToEvaluate.Income.household=@( ...
        labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
        r,pension,AccidentBeqS,AccidentBeqAH,w,P0,D,Lhscale, ...
        kappa_j,tau_l,tau_d,tau_cg,agej,Jr ...
    ) Electrify_HouseholdIncomeFn( ...
        labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
        r,pension,AccidentBeqS,AccidentBeqAH,w,P0,D,Lhscale, ...
        kappa_j,tau_l,tau_d,tau_cg,agej,Jr);

AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params, n_d, n_a, n_z,N_j, Names_i, d_grid, a_grid, z_grid,simoptions);

Y=AggVars.Output.Mean;

P=((1-Params.tau_cg)*Params.P0 + (1-Params.tau_d)*Params.D)/(1+Params.r-Params.tau_cg);

% Calculate the aggregate TFP as output/((capital^alpha_k)*(labor^alpha_l))
AggregateTFP=Y/((AggVars.K.Mean^Params.alpha_k)*(AggVars.L_f.Mean^Params.alpha_l));

% Total value of firms
temp=V.firm.*StationaryDist.firm;
temp(StationaryDist.firm==0)=0; % Get rid of points that have V=-inf but zero mass which would give nan
TotalValueOfFirms=sum(temp(isfinite(temp)));

fprintf('Following are some aggregates of the model economy: \n')
fprintf('Output: Y=%8.2f \n',AggVars.Output.Mean)
fprintf('Aggregate TFP: Y=%8.2f \n',AggregateTFP)
fprintf('Capital-Output ratio (firm side): K/Y=%8.2f \n',AggVars.K.Mean/Y)
fprintf('Total share value (HH side): P*S=%8.2f \n',P*AggVars.S.Mean)
fprintf('Total asset value (HH side): A=%8.2f \n',AggVars.A.Mean)
fprintf('Total house value (HH side): H=%8.2f \n',AggVars.H.Mean)
fprintf('Total net worth (HH side): P*S=%8.2f \n',P*AggVars.S.Mean+AggVars.A.Mean+AggVars.H.Mean)
fprintf('Total firm value (firm side): Value of firm=%8.2f \n',TotalValueOfFirms)
fprintf('Consumption-Output ratio: C/Y=%8.2f \n',AggVars.Consumption.Mean/Y)
fprintf('Government-to-Output ratio: G/Y=%8.2f \n', Params.G/Y)
fprintf('Wage: w=%8.2f \n',Params.w)










