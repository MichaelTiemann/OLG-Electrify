%% OLG Electrification (based on OLGModels14: Heterogenous households and heterogeneous firms
%% and also Life-Cycle Model 35: Portfolio-Choice with Housing)
% See https://www.vfitoolkit.com/updates-blog/2021/an-introduction-to-life-cycle-models/
% OLGModel14.m in the repo https://github.com/vfitoolkit/IntroToOLGModels
% and LifeCycleModel35.m in the repo https://github.com/vfitoolkit/IntroToLifeCycleModels
% See https://github.com/MichaelTiemann/OLG-Electrify/blob/main/README.md for more info

% A line some need for running on the Server
addpath(genpath('./MatlabToolkits/'))

solve_setup=true;
test_Lhscale=false;
solve_GE=3; % 0: skip GE; 1: solve initial, 2: solve final, 3: solve both
solve_TPath=false;
small_z_no_e=true; % n_z=1; n_e=0
small_model=false; % Minimal vs. maximal grid sizes
small_T=0; % small_T==1 means just do T=1, T=2 (or smallest not-to-be-confused-with-dimension); small_T==2 means use jpT

if solve_setup

Names_i={'firm','household','energy'};
PTypeDistParamNames={'ptypemass'};
Params.ptypemass=[1,1,1]; % Mass of households and firms are each equal to one

%% Parameters for household (4 scenarios)
% Scenario 1: no housing, no assets, no inflation
% Scenario 2: add rental+energy costs, but no housing/assets/inflation
% Scenario 3: add housing/assets/pv/inflation
% Scenario 4: add cars/detailed energy
Params.scenario=4;

% To be able to solve such a big problem, I switched to 5 year model period.
% Note that ypp (years-per-period) must be at most 15 (for kappa_j labor productivity evolutions).
% Discounting parameters (beta_pp and sj) defined in terms of ypp
Params.ypp=5; % model period, in years (just used this to modify some parameters from annual to model period)

% Lets model agents from age 20 to age 100, so 81 periods (or 61 for scenario 3)
max_age=100;
agejshifter=19; % Age 20 minus one. Makes keeping track of actual age easy in terms of model age

%% Global parameters (applies to household and firm)
% Note: with w=1, labor tax=20%, kappa_j(1)=0.5, agents have 0.4 budget to start
% If rent=0.3 energy=0.1, they live, but cannot save; they strike. Adjust rent to 0.3*sqrt(kappa_j)

% Annual risk-free rate of return
r=0.05; % We will discover risk-free rate of return per period in GE
r_wedge=0.05;

%% Parameters for households
% Discount rate; Changed to get S to increase nearer to 1 given r=0.05
% (ran it with beta=0.99, got S=0.3, so increased this; note that it interacts with sj to give the actual discount factor)
beta=[0.95,0.95,0.99,0.99];

Params.sigma = 2; % Coeff of relative risk aversion (curvature of consumption)

energy_pct_cost=[0,0.07,0.07,0.05]; % Electricity: 3%; Gas: 1-2%; Petrol: 1-2%; Scenario 4 disaggregates petrol from this cost

% Demographics
% Population growth rate
n=0.02; % percentage rate (expressed as fraction) of population growth per period

% Age-dependent labor productivity units
% Stage 1: starting out (typ. first 25-30 years)
% Stage 2: peak earnings (typ. years 25-30 (meaning ages 45-50))
% Stage 3: winding down (typ. last 14 years before retirement (ages 50-64))
% Stage r: retirement
% Labor productivity at start, peak, and end of working life
k_j1 = [0.5, 0.5, 0.5, 0.5];
k_j2 = [2, 2, 2, 2];
k_j2_length = [0,0,5,5];
k_j3 = [1, 1, 1, 1];

% Note: These iid shocks will interact with the endogenous labor so the final labor
% earnings process will not equal that of Karahan & Ozkan (2013)
% Note: Karahan & Ozkan (2013) also have a fixed effect (which they call alpha) and which I ignore here.

% Warm glow of bequest
Params.warmglow1=0.3; % (relative) importance of bequests
Params.warmglow2=3; % bliss point of bequests (essentially, the target amount)
Params.warmglow3=Params.sigma; % By using the same curvature as the utility of consumption it makes it much easier to guess appropraite parameter values for the warm glow

% The warmglow parameters will help us find the GE solution to actual bequest rates/values
AccidentBeqS=[0.02,0.02,0.02,0.02]; % Accidental bequests (this is the lump sum transfer of shares)
AccidentBeqAH=[0,0,0.02,0.02]; % Accidental bequests (this is the lump sum transfer of assets+house value)

% Preferences
% Relative importance of housing services (vs consumption) in utility
sigma_h=[0,0,0.5,0.2];
% Relative importance of car services (vs consumption) in utility
sigma_c=[0,0,0.5,0.3];
Params.eta=1.5; % Curvature of leisure (This will end up being 1/Frisch elasty)
psi = [2, 1, 1, 1]; % Weight on leisure

%% Energy parameters
if Params.scenario==4
    Params.Ek=1; Params.ek=1;
    Params.carbon_tax=35;
    Params.energy_pct_brown=0.8;
end

%% Government spending (to be found in GE)
G=0.1;
% Dividend target (to be found in GE)
D=0.2;

%% Taxes
% Household labor tax
Params.tau_l = 0.2; % Tax rate on labour income
% Firm Taxes
Params.tau_corp=0.34; % Tax rate on corporate earnings
Params.phi=0.5; % Fraction of capital adjustment costs that can be deducted from corporate earnings
Params.tau_d=0.2; % Tax rate on dividends
Params.tau_cg=0.2; % Tax rate on capital gains

Params=Electrify_Scenario_YPP_Setup(Params,Params.scenario,Params.ypp,small_z_no_e,max_age,agejshifter,r,r_wedge,beta,n,k_j1,k_j2,k_j2_length,k_j3,sigma_h,sigma_c,psi,Params.tau_cg,energy_pct_cost,G,D,AccidentBeqS,AccidentBeqAH);

% Housing (ignored/overwritten if no housing in scenario)
% Params.minhouse % set below, is the minimum value of house that can be purchased
Params.rentprice=0.3; % To make real fraction of income, must be multiplied by kappa_j in scenarios 3 & 4
% housing services as a fraction of house value (ignored in scenarios w/o housing)
Params.houseservices=0.5;
Params.f_htc=0.05; % transaction cost of buying/selling house (is a percent of h+hprime)
Params.f_coll=0.5; % collateral contraint (fraction of house value that can be borrowed)
Params.pv_pct_cost=0.05; % modeling a $30K install for a $600K house

%% Parameters for firm
% Production
Params.alpha_k=0.311; % diminishing returns to capital and energy inputs
Params.alpha_l=0.650; % diminishing returns to labor input
Params.delta=0.054; % Annual depreciation of physical capital
% Capital adjustment costs
Params.capadjconstant=1.21; % term in the capital adjustment cost

% Idiosyncatic productivity shocks
Params.rho_z_firm=0.767;
Params.sigma_z_e_firm=0.211;

%% Parameters for energy
% Idiosyncatic productivity shocks
Params.rho_z_energy=0.767;
Params.sigma_z_e_energy=0.211;

% Set the firm discount factor below (as it is determined in general eqm)
% Params.firmbeta=1/(1+Params.r_pp/(1-Params.tau_cg)); % 1/(1+r_pp) but returns net of capital gains tax

%% Create our Grids from Scenario and Parameters
vfoptions=struct(); simoptions=struct();
[n_d,n_a,n_z,N_j,vfoptions]=Electrify_GridSizeSetup(Params.scenario, Params.J, small_z_no_e, small_model, vfoptions);
[d_grid,a_grid,z_grid,pi_z,jequaloneDist,share_grid,k_grid,pv_grid_firm,Params,vfoptions,simoptions]=Electrify_GridSetup(Params.scenario, Params.ypp, n_d, n_a, n_z, small_z_no_e, Params, vfoptions, simoptions);

% Set up Transition Path control parameters
if small_T==2
    jpT=3;
else
    jpT=1; % Default: one transition period=1 time period; Could have multiple j's per T
end

T=ceil(Params.J*1.4/jpT);
if T==length(Names_i)
    T=T+1;
end
if T==Params.J
    % The toolkit thinks that T and J must be different (T larger to reach equilibrium post J)
    T=T+1;
end

%% Remaining Parameters will be set in GE below

%% Now, create the return function

% For households
DiscountFactorParamNames.household={'beta_pp','sj'};

% Hardwire buyhouse, hprime, aprime, h, a, and solarpv to zero
ReturnFn_12.household=@( ...
    labor,sprime,s,z,e, ...
    pension,AccidentBeqS_pp,w,P0,D_pp, ...
    sigma,psi,eta,sigma_h,kappa_j,tau_l,tau_d,tau_cg,warmglow1,warmglow2,ypp,agej,Jr,J,...
    scenario,r_pp,r_wedge_pp,minhouse,rentprice,houseservices,energy_pct_cost ...
) Electrify_HouseholdReturnFn( ...
    labor,0,sprime,0,0,s,0,0,0,z,e, ...
    pension,AccidentBeqS_pp,0,w,P0,D_pp, ...
    sigma,psi,eta,sigma_h,kappa_j,tau_l,tau_d,tau_cg,warmglow1,warmglow2,ypp,agej,Jr,J,...
    scenario,r_pp,r_wedge_pp,0,minhouse,rentprice,0,houseservices,0,0,energy_pct_cost ...
);
% Notice we use 'Electrify_HouseholdReturnFn'
ReturnFn_3.household=@( ...
        labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
        pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
        sigma,psi,eta,sigma_h,kappa_j,tau_l,tau_d,tau_cg,warmglow1,warmglow2,ypp,agej,Jr,J,...
        scenario,r_pp,r_wedge_pp,f_htc,minhouse,rentprice,f_coll,houseservices,cpi,pv_pct_cost,energy_pct_cost ...
    ) Electrify_HouseholdReturnFn( ...
        labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
        pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
        sigma,psi,eta,sigma_h,kappa_j,tau_l,tau_d,tau_cg,warmglow1,warmglow2,ypp,agej,Jr,J,...
        scenario,r_pp,r_wedge_pp,f_htc,minhouse,rentprice,f_coll,houseservices,cpi,pv_pct_cost,energy_pct_cost ...
    );
ReturnFn_4.household=@( ...
        labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e, ...
        pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
        sigma,psi,eta,sigma_h,sigma_c,kappa_j,tau_l,tau_d,tau_cg,warmglow1,warmglow2,ypp,agej,Jr,J,...
        r_pp,r_wedge_pp,f_htc,minhouse,rentprice,f_coll,houseservices,carservices_j,cpi_energy,pv_pct_cost,energy_pct_cost,energy_pct_brown,carbon_tax ...
    ) Electrify_4HouseholdReturnFn( ...
        labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e, ...
        pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
        sigma,psi,eta,sigma_h,sigma_c,kappa_j,tau_l,tau_d,tau_cg,warmglow1,warmglow2,ypp,agej,Jr,J,...
        r_pp,r_wedge_pp,f_htc,minhouse,rentprice,f_coll,houseservices,carservices_j,cpi_energy,pv_pct_cost,energy_pct_cost,energy_pct_brown,carbon_tax ...
    );

if Params.scenario<3
    ReturnFn.household=ReturnFn_12.household;
elseif Params.scenario<4
    ReturnFn.household=ReturnFn_3.household;
else
    ReturnFn.household=ReturnFn_4.household;
end

% For firms
DiscountFactorParamNames.firm={'firmbeta'};

% Notice we use 'Electrify_FirmReturnFn'
ReturnFn_123.firm=@( ...
        d,kprime,k,z, ...
        w, D_pp, ...
        ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg ...
    ) Electrify_FirmReturnFn( ...
        d,kprime,k,z, ...
        w, D_pp, ...
        ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg ...
    );
% Notice we use 'Electrify_4FirmReturnFn'
ReturnFn_4.firm=@( ...
        kprime,pvprime,k,pv,z, ...
        w, ...
        ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg,Ek,ek,pv_max_firm,carbon_tax ...
    ) Electrify_4FirmReturnFn( ...
        0,kprime,pvprime,k,pv,z, ...
        w, ...
        ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg,Ek,ek,pv_max_firm,carbon_tax ...
    );

if Params.scenario<4
    ReturnFn.firm=ReturnFn_123.firm;
else
    ReturnFn.firm=ReturnFn_4.firm;
end

% For energy
DiscountFactorParamNames.energy={};

% Notice we use 'Electrify_EnergyReturnFn'
ReturnFn_123.energy=@( ...
        aprime,a,z ...
    ) Electrify_EnergyReturnFn( ...
        aprime,a,z ...
    );
% Notice we use 'Electrify_EnergyReturnFn'
ReturnFn_4.energy=@( ...
        d,aprime,a,z ...
    ) Electrify_4EnergyReturnFn( ...
        d,aprime,a,z ...
    );

if Params.scenario<4
    ReturnFn.energy=ReturnFn_123.energy;
else
    ReturnFn.energy=ReturnFn_4.energy;
end

%% Begin setting up to use VFI Toolkit to solve
% vfoptions.howardsgreedy=0;
% vfoptions.howards=80;
% vfoptions.maxhowards=200;
if Params.scenario<3 && small_model==false
    vfoptions.tolerance=10^(-9);
else
    vfoptions.tolerance=10^(-4);
end
% Note that simoptions.tolerance is used very differently than vfoptions.tolerance

% The user can experiment with gridinterplayer=0 (pure discretization) or gridinterplayer=1 (linear interpolation b/w grid points).
% If gridinterplayer=1, then you must set vfoptions.divideandconquer=1 (required for transition).
vfoptions.gridinterplayer.household  = 0;
vfoptions.level1n.household          = 5;
vfoptions.divideandconquer.household = logical(Params.scenario<3);
vfoptions.gridinterplayer.firm       = 0;
vfoptions.divideandconquer.firm      = 0;
simoptions.gridinterplayer = vfoptions.gridinterplayer;
% simoptions.ngridinterp     = vfoptions.ngridinterp;

%% Remaining parameters
% Steering GE
Params.TargetKdivL=2.03;

% Parameters set/evolved by Transition Paths
Params.cpi=0; % Initial condition
Params.cpi_energy=0; % Initial condition

% Scaling the household labor supply; we scale model and GE finds its own equilibrium
% This is vaguely scenario-by-ypp
Lhscale=[[0.25,0.21,0.20,116];
    [0.37,0.27,0.24,2];    % 2
    [0.54,0.40,0.30,2.3];
    [0.57,0.44,0.31,2.3];  % 4
    [0.65,1.1,0.32,2.3];
    [3.3,2.8,1.6,2.3];     % 6
    [5.0,3.2,2.1,2.6];
    [6.2,4.2,2.9,3.0];     % 8
    [7.5,4.7,3.5,3.2];
    [7.9,5.0,3.5,3.2];    % 10
    [8.7,5.8,3.6,3.4];
    [9.5,6.3,3.7,3.6];    % 12
    ]; 
if Params.ypp<=size(Lhscale,1)
    Lhscale_final=Lhscale(Params.ypp,Params.scenario)*1.1;
    ParamPath.Lhscale=linspace(Lhscale(Params.ypp,Params.scenario),Lhscale_final,T);
else
    Lhscale_final=Lhscale(end,Params.scenario)*1.1;
    ParamPath.Lhscale=linspace(Lhscale(end,Params.scenario),Lhscale_final,T);
end
Params.Lhscale=ParamPath.Lhscale(1);

% Solved by GE

% Some initial values/guesses for variables that will be determined in general eqm
Params.P0=1;
Params.w=1;
Params.pension=0.4; % Initial guess (this will be determined in general eqm)
Params.G_pp=0.1*Params.ypp; % Government expenditure

% And some initial values/guesses for AggVar values that will be calculated while calculating the general eqm
Params.D_pp=(1+0.20)^Params.ypp-1; % The dividends paid by the firm per period
Params.EnergyCosts_h=0.3; % Energy used by households
Params.EnergyCosts_f=0.7; % Energy used by firms
Params.CarbonCosts_h=0.05; % Carbon tax paid by households
Params.CarbonCosts_f=0.3; % Cabron tax by firms

%% General eqm variables
if Params.scenario<3
    GEPriceParamNames={'w','D_pp','P0','pension','G_pp','AccidentBeqS_pp'};
elseif Params.scenario<4
    GEPriceParamNames={'w','D_pp','P0','pension','G_pp','AccidentBeqS_pp','AccidentBeqAH_pp'};
else
    GEPriceParamNames={'w','P0','pension','G_pp','AccidentBeqS_pp','AccidentBeqAH_pp'};
end
heteroagentoptions.constrainpositive=GEPriceParamNames;

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
% Note also we must differentiate based on Scenarios...

% Stationary Distribution Aggregates from households (important that ordering of Names and Functions is the same)
FnsToEvaluate_12.L_h.household=@(labor,sprime,s,z,e,kappa_j,Lhscale) ...
    labor*kappa_j*exp(z+e)*Lhscale;  % Aggregate labour supply in efficiency units, not scaled by ypp
FnsToEvaluate_12.S.household=@(labor,sprime,s,z,e) s; % Aggregate share holdings
FnsToEvaluate_12.PensionSpending.household=@(labor,sprime,s,z,e,pension,ypp,agej,Jr) ...
    (agej>=Jr)*pension*ypp; % Total spending on pensions
FnsToEvaluate_12.PayrollTaxRevenue.household=@(labor,sprime,s,z,e,ypp,agej,Jr,tau_l,w,kappa_j,Lhscale) ...
    (agej<Jr)*tau_l*labor*w*kappa_j*exp(z+e)*ypp*Lhscale; % Total spending on payroll taxes
FnsToEvaluate_12.CapitalGainsTaxRevenue.household=@(labor,sprime,s,z,e,tau_cg,P0,D_pp,tau_d,r_pp) ...
    tau_cg*(P0-(((1-tau_cg)*P0 + (1-tau_d)*D_pp)/(1+r_pp-tau_cg)))*s; % tau_cg*(P0-Plag)*s, but substitute P=Plag, and then substitute for P
FnsToEvaluate_12.BeqleftS_pp.household=@(labor,sprime,s,z,e,sj) ...
    sprime*(1-sj); % Accidental share bequests left by people who die
FnsToEvaluate_3.L_h.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,kappa_j,Lhscale) ...
    labor*kappa_j*exp(z+e)*Lhscale;  % Aggregate labour supply in efficiency units, not scaled by ypp
FnsToEvaluate_3.S.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) s; % Aggregate share holdings
FnsToEvaluate_3.A.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) a; % Aggregate share holdings
FnsToEvaluate_3.PensionSpending.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,pension,ypp,agej,Jr) ...
    (agej>=Jr)*pension*ypp; % Total spending on pensions
FnsToEvaluate_3.PayrollTaxRevenue.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,ypp,agej,Jr,tau_l,w,kappa_j,Lhscale) ...
    (agej<Jr)*tau_l*labor*w*kappa_j*exp(z+e)*Lhscale*ypp; % Total spending on payroll taxes
FnsToEvaluate_3.CapitalGainsTaxRevenue.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,tau_cg,P0,D_pp,tau_d,r_pp) ...
    tau_cg*(P0-(((1-tau_cg)*P0 + (1-tau_d)*D_pp)/(1+r_pp-tau_cg)))*s+(1-tau_d)*r_pp*max(a,0); % tau_cg*(P0-Plag)*s + deposit interest, but substitute P=Plag, and then substitute for P
FnsToEvaluate_3.BeqleftS_pp.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,sj) ...
    sprime*(1-sj); % Accidental share bequests left by people who die
% AccidentalBeqAHLeft is zero (if in debt) or accidental asset+house bequests left by people who die
FnsToEvaluate_3.BeqleftAH_pp.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,scenario,sj,cpi) ...
    max(0,(aprime+(1+cpi)*hprime)*(1-sj));
% BadDebt is the debt somebody accidentally leaves behind, or zero if net worth is positive
FnsToEvaluate_4.L_h.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,kappa_j,Lhscale) ...
    labor*kappa_j*exp(z+e)*Lhscale;  % Aggregate labour supply in efficiency units, not scaled by ypp
FnsToEvaluate_4.S.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e) s; % Aggregate share holdings
FnsToEvaluate_4.PensionSpending.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,pension,ypp,agej,Jr) ...
    (agej>=Jr)*pension*ypp; % Total spending on pensions
FnsToEvaluate_4.PayrollTaxRevenue.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,ypp,agej,Jr,tau_l,w,kappa_j,Lhscale) ...
    (agej<Jr)*tau_l*labor*w*kappa_j*exp(z+e)*Lhscale*ypp; % Total spending on payroll taxes
FnsToEvaluate_4.CapitalGainsTaxRevenue.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,tau_cg,P0,D_pp,tau_d,r_pp) ...
    tau_cg*(P0-(((1-tau_cg)*P0 + (1-tau_d)*D_pp)/(1+r_pp-tau_cg)))*s+(1-tau_d)*r_pp*max(a,0); % tau_cg*(P0-Plag)*s + deposit interest, but substitute P=Plag, and then substitute for P
FnsToEvaluate_4.BeqleftS_pp.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,sj) ...
    sprime*(1-sj); % Accidental share bequests left by people who die
% AccidentalBeqAHLeft is zero (if in debt) or accidental asset+house bequests left by people who die
FnsToEvaluate_4.BeqleftAH_pp.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,scenario,sj,cpi) ...
    max(0,(aprime+(1+cpi)*hprime)*(1-sj));
% BadDebt is the debt somebody accidentally leaves behind, or zero if net worth is positive
FnsToEvaluate_4.EnergyCosts_h.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,w,ypp,cpi_energy,energy_pct_cost,energy_pct_brown,carbon_tax) ...
    Electrify_4HouseholdEnergyCosts(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,w,ypp,cpi_energy,energy_pct_cost,energy_pct_brown,carbon_tax);
FnsToEvaluate_4.CarbonCosts_h.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,w,ypp,cpi_energy,energy_pct_cost,energy_pct_brown,carbon_tax) ...
    Electrify_4HouseholdCarbonCosts(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,w,ypp,cpi_energy,energy_pct_cost,energy_pct_brown,carbon_tax);

% From firms
FnsToEvaluate_12.L_f.firm=@(d,kprime,k,z,w,alpha_k,alpha_l) ...
    (w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % (effective units of) labor demanded by firm, not scaled by ypp
FnsToEvaluate_12.K.firm=@(d,kprime,k,z,w,alpha_k,alpha_l) k; % physical capital
FnsToEvaluate_12.dividend_pp.firm=@(d,kprime,k,z,ypp) (1+d)^ypp-1; % dividend paid by firm
FnsToEvaluate_12.Sissued.firm=@(d,kprime,k,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi) ...
    Electrify_FirmShareIssuance(d,kprime,k,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi); % Share issuance
FnsToEvaluate_12.CorpTaxRevenue.firm=@(d,kprime,k,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi) ...
    Electrify_FirmCorporateTaxRevenue(d,kprime,k,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi); % revenue from the corporate profits tax
fnnames=fieldnames(FnsToEvaluate_12);
for ff=1:length(fnnames)
    if isfield(FnsToEvaluate_12.(fnnames{ff}), 'firm')
        FnsToEvaluate_3.(fnnames{ff}).firm=FnsToEvaluate_12.(fnnames{ff}).firm;
    end
end
FnsToEvaluate_4.L_f.firm=@(kprime,pvprime,k,pv,z,w,alpha_k,alpha_l) ...
    (w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % (effective units of) labor demanded by firm, not scaled by ypp
FnsToEvaluate_4.K.firm=@(kprime,pvprime,k,pv,z,w,alpha_k,alpha_l) k; % physical capital
FnsToEvaluate_4.PV_f.firm=@(kprime,pvprime,k,pv,z,w,alpha_k,alpha_l) pv; % firm's solarPV generation capacity
FnsToEvaluate_4.D_pp.firm=@(kprime,pvprime,k,pv,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,Ek,ek,pv_max_firm,carbon_tax) ...
    Electrify_4FirmDividend(0,kprime,pvprime,k,pv,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,Ek,ek,pv_max_firm,carbon_tax); % dividend paid by firm
FnsToEvaluate_4.Sissued.firm=@(kprime,pvprime,k,pv,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,Ek,ek,pv_max_firm,carbon_tax) ...
    Electrify_4FirmShareIssuance(0,kprime,pvprime,k,pv,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,Ek,ek,pv_max_firm,carbon_tax); % Share issuance
FnsToEvaluate_4.CorpTaxRevenue.firm=@(kprime,pvprime,k,pv,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,Ek,ek,pv_max_firm,carbon_tax) ...
    Electrify_4FirmCorporateTaxRevenue(0,kprime,pvprime,k,pv,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,Ek,ek,pv_max_firm,carbon_tax); % revenue from the corporate profits tax
FnsToEvaluate_4.EnergyCosts_f.firm=@(kprime,pvprime,k,pv,z,w,ypp,alpha_k,alpha_l,Ek,ek,pv_max_firm,carbon_tax) ...
    Electrify_4FirmEnergyCosts(0,kprime,pvprime,k,pv,z,w,ypp,alpha_k,alpha_l,Ek,ek,pv_max_firm,carbon_tax);
FnsToEvaluate_4.CarbonCosts_f.firm=@(kprime,pvprime,k,pv,z,w,ypp,alpha_k,alpha_l,Ek,ek,pv_max_firm,carbon_tax) ...
    Electrify_4FirmCarbonCosts(0,kprime,pvprime,k,pv,z,w,ypp,alpha_k,alpha_l,Ek,ek,pv_max_firm,carbon_tax);

% From energy -- there must be at least one
FnsToEvaluate_12.EnergyRevenue.energy=@(aprime,a,z) 0;
FnsToEvaluate_3.EnergyRevenue.energy=@(aprime,a,z) 0;
FnsToEvaluate_4.EnergyRevenue.energy=@(invest,aprime,a,z,EnergyCosts_h,EnergyCosts_f) EnergyCosts_h+EnergyCosts_f;
FnsToEvaluate_4.TransitionInvestment.energy=@(invest,aprime,a,z,CarbonCosts_h,CarbonCosts_f) CarbonCosts_h+CarbonCosts_f;

% General Equilibrium conditions (these should evaluate to zero in general equilbrium)
GeneralEqmEqns.sharemarket=@(S) S-1; % mass of all shares equals one
GeneralEqmEqns.labormarket=@(L_h,L_f) (L_h-L_f)*max(2,Params.ypp); % labor supply of households equals labor demand of firms (scaled by ypp)
GeneralEqmEqns.pensions=@(PensionSpending,PayrollTaxRevenue) PensionSpending-PayrollTaxRevenue; % Retirement benefits equal Payroll tax revenue: pension*fractionretired-tau*w*H
GeneralEqmEqns.govbudget=@(G_pp,tau_d,D_pp,CapitalGainsTaxRevenue,CorpTaxRevenue) G_pp-tau_d*D_pp-CapitalGainsTaxRevenue-CorpTaxRevenue; % G is equal to the target, GdivYtarget*Y
% GeneralEqmEqns.firmdiscounting=@(firmbeta,r_pp,tau_cg) firmbeta-1/(1+r_pp/(1-tau_cg)); % Firms discount rate is related to market return rate
if Params.scenario<4
    GeneralEqmEqns.dividends=@(dividend_pp,D_pp) (dividend_pp-D_pp); % That the dividend households receive equals that which firms give
end
GeneralEqmEqns.ShareIssuance=@(Sissued,P0,D_pp,tau_cg,tau_d,r_pp) ...
    P0-((((1-tau_cg)*P0 + (1-tau_d)*D_pp)/(1+r_pp-tau_cg))-Sissued); % P0=P-S, but substitute for P (see derivation inside the return fn)
GeneralEqmEqns.CapitalOutputRatio=@(K,L_f,TargetKdivL) (K/L_f-TargetKdivL)/100; % Ratio not based on ypp
GeneralEqmEqns.bequestsS_pp=@(BeqleftS_pp,AccidentBeqS_pp,n_pp) BeqleftS_pp/(1+n_pp)-AccidentBeqS_pp; % Accidental share bequests received equal accidental share bequests left
if Params.scenario>2
    GeneralEqmEqns.bequestsAH_pp=@(BeqleftAH_pp,AccidentBeqAH_pp,n_pp) BeqleftAH_pp/(1+n_pp)-AccidentBeqAH_pp; % Accidental asset+house bequests received equal accidental asset+house bequests left
end

% For analysing the model
FnsToEvaluate2_12=FnsToEvaluate_12;
FnsToEvaluate2_12.earnings.household=@(labor,aprime,a,z,e,w,kappa_j,Lhscale) w*kappa_j*labor*exp(z+e)*Lhscale; % w*kappa_j is the labor earnings
FnsToEvaluate2_12.A.household=@(labor,aprime,a,z,e,w,kappa_j,Lhscale) a; % w*kappa_j is the labor earnings
FnsToEvaluate2_12.BeqleftS_pp.household=@(labor,aprime,a,z,e,sj) aprime*(1-sj); % Accidental asset bequests left by people who die
FnsToEvaluate2_3=FnsToEvaluate_3;
FnsToEvaluate2_3.earnings.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,w,kappa_j,Lhscale) w*kappa_j*labor*exp(z+e)*Lhscale; % w*kappa_j is the labor earnings
FnsToEvaluate2_3.A.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) a; % Aggregate asset/mortgage holdings
FnsToEvaluate2_3.S.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) s; % Aggregate share holdings
FnsToEvaluate2_3.H.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) h; % Aggregate house holdings
FnsToEvaluate2_3.PV_h.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) solarpv; % Aggregate solarpv holdings
FnsToEvaluate2_3.BeqleftS_pp.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,sj) sprime*(1-sj); % Accidental share bequests left by people who die
FnsToEvaluate2_3.BeqleftAH_pp.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,scenario,sj,cpi) max(0,(aprime+(1+cpi)*hprime)*(1-sj));
FnsToEvaluate2_3.BadDebt_pp.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,scenario,sj,cpi) ...
    min(0,(aprime+(1+cpi)*hprime)*(1-sj));
FnsToEvaluate2_4=FnsToEvaluate_4;
FnsToEvaluate2_4.earnings.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,w,kappa_j,Lhscale) w*kappa_j*labor*exp(z+e)*Lhscale; % w*kappa_j is the labor earnings
FnsToEvaluate2_4.A.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e) a; % Aggregate asset/mortgage holdings
FnsToEvaluate2_4.S.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e) s; % Aggregate share holdings
FnsToEvaluate2_4.Car.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e) car; % Aggregate house holdings
FnsToEvaluate2_4.H.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e) h; % Aggregate house holdings
FnsToEvaluate2_4.PV_h.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e) solarpv; % Aggregate solarpv holdings
FnsToEvaluate2_4.BeqleftS_pp.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,sj) sprime*(1-sj); % Accidental share bequests left by people who die
FnsToEvaluate2_4.BeqleftAH_pp.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,scenario,sj,cpi) max(0,(aprime+(1+cpi)*hprime)*(1-sj));
FnsToEvaluate2_4.BadDebt_pp.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,scenario,sj,cpi) ...
    min(0,(aprime+(1+cpi)*hprime)*(1-sj));
FnsToEvaluate2_12.Output.firm=@(d,kprime,k,z,w,ypp,alpha_k,alpha_l) ...
    z*(k^alpha_k)*((w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)))^alpha_l*ypp; % Production function z*(k^alpha_k)*(l^alpha_l) (substituting for l)
FnsToEvaluate2_3.Output.firm=FnsToEvaluate2_12.Output.firm;
FnsToEvaluate2_4.Output.firm=@(kprime,pvprime,k,pv,z,w,ypp,alpha_k,alpha_l,Ek,ek) ...
    (Ek*ek)*z*(k^alpha_k)*((w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)))^alpha_l*ypp; % Production function z*(k^alpha_k)*(l^alpha_l) (substituting for l)
FnsToEvaluate2_4.CarbonCosts_f.firm=@(kprime,pvprime,k,pv,z,w,ypp,alpha_k,alpha_l,Ek,ek,pv_max_firm,carbon_tax) ...
    Electrify_4FirmCarbonCosts(0,kprime,pvprime,k,pv,z,w,ypp,alpha_k,alpha_l,Ek,ek,pv_max_firm,carbon_tax);

% Note: I keep the FnsToEvaluate use in general eqm to a minimum (to reduce
% runtimes) and then use FnsToEvaluate2 to analyse model with more stats.
% Note: FnsToEvaluate may need 'e' grids, but AggVars and other stats use a
% a joint ze grid (which reads as z in their parameter lists).

if Params.scenario==3
    [V_f, Policy_f]=ValueFnIter_Case1_PType(n_d,n_a,n_z, {'firm'}, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
    % We can plot V as a 3d plot (surf is matlab command for 3d plot)
    figure(5)
    subplot(2,1,1);
    % Plot F as K vs D
    s_grid_firm=zeros(length(k_grid),length(d_grid.firm));
    for D=1:length(d_grid.firm)
        for K=1:length(k_grid)
            if k_grid(K)>d_grid.firm(D)
                % Pay dividend from capital stock
                s_grid_firm(K,D)=0;
            else
                % Exhaust capital stock and issue shares
                s_grid_firm(K,D)=Electrify_FirmShareIssuance(d_grid.firm(D),k_grid(K),k_grid(K),1,Params.w,Params.ypp,Params.delta,Params.alpha_k,Params.alpha_l,Params.capadjconstant,Params.tau_corp,Params.phi);
            end
        end
    end
    surf(d_grid.firm,k_grid,s_grid_firm)
    title('Value function: S as K vs D')
    xlabel('D')
    ylabel('K')
elseif Params.scenario==4
    [V_f, Policy_f]=ValueFnIter_Case1_PType(n_d,n_a,n_z, {'firm'}, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
    % We can plot V as a 3d plot (surf is matlab command for 3d plot)
    figure(5)
    subplot(2,1,1);
    % Plot F as K vs PV
    if small_z_no_e
        surf(pv_grid_firm,k_grid,V_f.firm)
    else
        surf(pv_grid_firm,k_grid,V_f.firm(:,:,4))
    end
    title('Value function: F as K vs PV')
    xlabel('PV')
    ylabel('K')
end

if Params.scenario>=3 % shares vs. assets not possible in Scenarios 1 and 2
    % household = [S, A, H, PV, z, agej]
    [V_h, Policy_h]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j, {'household'}, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
    share_grid=a_grid.household(1:n_a.household(1));
    asset_grid=a_grid.household(n_a.household(1)+1:n_a.household(1)+n_a.household(2));
    % We can plot V as a 3d plot (surf is matlab command for 3d plot)
    figure(6)
    for row=1:3
        for col=1:3
            subplot(3,3,(row-1)*3+col);
            agej=(row-1)*6+col*2-1;
            % Plot F as K vs PV
            % Plot F as S vs H
            if small_z_no_e
                if Params.scenario<4
                    surf(asset_grid,share_grid,V_h.household(:,:,2,2,1,agej), V_h.household(:,:,2,2,1,agej)-V_h.household(:,:,1,1,1,agej))
                else
                    surf(asset_grid,share_grid,V_h.household(:,:,2,2,2,1,agej), V_h.household(:,:,2,2,2,1,agej)-V_h.household(:,:,1,1,1,1,agej))
                end
                colorbar
            else
                if Params.scenario<4
                    surf(asset_grid,share_grid,V_h.household(:,:,1,1,4,Params.Jr))
                else
                    surf(asset_grid,share_grid,V_h.household(:,:,1,1,1,4,Params.Jr))
                end
            end
            title(sprintf('Value function: F as S vs A at age %d', (agej-1)*Params.ypp+agejshifter+1))
            xlabel('A')
            ylabel('S')
        end
    end
end


%% Agents age distribution
AgeWeightsParamNames=struct('household',{{'mewj'}}); % So VFI Toolkit knows which parameter is the mass of agents of each age

% Find parameters for Lhscale; note this is all small_z_no_e or all ~small_z_no_e
if test_Lhscale
    if small_z_no_e
        small_z_no_e_string="true";
    else
        small_z_no_e_string="false";
    end
    for scenario=1:4
        if scenario<3
            ReturnFn_Lh.household=ReturnFn_12.household;
        elseif scenario<4
            ReturnFn_Lh.household=ReturnFn_3.household;
        else
            ReturnFn_Lh.household=ReturnFn_4.household;
        end
        if scenario<4
            ReturnFn_Lh.firm=ReturnFn_123.firm;
            ReturnFn_Lh.energy=ReturnFn_123.energy;
        else
            ReturnFn_Lh.firm=ReturnFn_4.firm;
            ReturnFn_Lh.energy=ReturnFn_4.energy;
        end
        for ypp=[7,8,9,10,12]
            vfoptions_Lh=struct(); simoptions_Lh=struct();
            Params_Lh=Electrify_Scenario_YPP_Setup(Params,scenario,ypp,small_z_no_e,max_age,agejshifter,r,r_wedge,beta,n,k_j1,k_j2,k_j2_length,k_j3,sigma_h,sigma_c,psi,Params.tau_cg,energy_pct_cost,G,D,AccidentBeqS,AccidentBeqAH);
            [n_d_Lh,n_a_Lh,n_z_Lh,N_j_Lh,vfoptions_Lh]=Electrify_GridSizeSetup(scenario, Params_Lh.J, small_z_no_e, small_model, vfoptions_Lh);
            [d_grid_Lh,a_grid_Lh,z_grid_Lh,pi_z_Lh,jequaloneDist_Lh,share_grid_Lh,k_grid_Lh,pv_grid_firm_Lh,Params_Lh,vfoptions_Lh,simoptions_Lh]=Electrify_GridSetup(scenario, ypp, n_d_Lh, n_a_Lh, n_z_Lh, small_z_no_e, Params_Lh, vfoptions_Lh, simoptions_Lh);
            if ypp<=size(Lhscale,1)
                Params_Lh.Lhscale=Lhscale(ypp,scenario);
            else
                Params_Lh.Lhscale=Lhscale(end,scenario);
            end
            [V_Lh, Policy_Lh]=ValueFnIter_Case1_FHorz_PType(n_d_Lh,n_a_Lh,n_z_Lh,N_j_Lh,Names_i, d_grid_Lh, a_grid_Lh, z_grid_Lh, pi_z_Lh, ReturnFn_Lh, Params_Lh, DiscountFactorParamNames, vfoptions_Lh);
            StationaryDist_Lh=StationaryDist_Case1_FHorz_PType(jequaloneDist_Lh,AgeWeightsParamNames,PTypeDistParamNames,Policy_Lh,n_d_Lh,n_a_Lh,n_z_Lh,N_j_Lh,Names_i,pi_z_Lh,Params_Lh,simoptions_Lh);
            if scenario<3
                FnsToEvaluate2_Lh=FnsToEvaluate2_12;
            elseif scenario<4
                FnsToEvaluate2_Lh=FnsToEvaluate2_3;
            else
                FnsToEvaluate2_Lh=FnsToEvaluate2_4;
            end
            AggVars_Lh=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_Lh, Policy_Lh, FnsToEvaluate2_Lh, Params_Lh, n_d_Lh, n_a_Lh, n_z_Lh,N_j_Lh,Names_i, d_grid_Lh, a_grid_Lh, z_grid_Lh,simoptions_Lh);
            fprintf('Check Scenario %d; ypp= %d; small_z_no_e = %s: L_hscale = %.2f; L_h, L_f, K are ', scenario, ypp, small_z_no_e_string, Params_Lh.Lhscale)
            if AggVars_Lh.L_f.Mean~=0 && abs(1-AggVars_Lh.L_h.Mean/AggVars_Lh.L_f.Mean)<0.1
                cprintf('', "%10.4f %10.4f %10.4f \n", AggVars_Lh.L_h.Mean,AggVars_Lh.L_f.Mean,AggVars_Lh.K.Mean);
            else
                cprintf('err', "%10.4f %10.4f %10.4f \n", AggVars_Lh.L_h.Mean,AggVars_Lh.L_f.Mean,AggVars_Lh.K.Mean);
            end
            clear vfoptions_Lh simoptions_Lh Params_Lh FnsToEvaluate2_Lh V_Lh Policy_Lh StationaryDist_Lh AggVars_Lh
        end
        clear ReturnFn_Lh
    end
else

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium
disp('Test ValueFnIter')
tic;
% Note: z_grid and pi_z, this will be ignored due to presence of vfoptions.z_grid_J and vfoptions.pi_z_J
[V_init, Policy_init]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
toc

%% Let's take a quick look at what we have calculated, namely V and Policy

% The value function V depends on the state, so now it depends on both asset holdings and age.
if Params.scenario<3
    % We can plot V as a 3d plot (surf is matlab command for 3d plot)
    zind=ceil(n_z.household/2);
    if ~small_z_no_e
        eind=ceil(vfoptions.n_e.household/2);
    end
    figure(1)
    subplot(2,1,1);
    if small_z_no_e
        surf(a_grid.household*ones(1,Params.J),ones(n_a.household,1)*(1:1:Params.J),reshape(V_init.household(:,zind,:),[n_a.household,Params.J]))
    else
        surf(a_grid.household*ones(1,Params.J),ones(n_a.household,1)*(1:1:Params.J),reshape(V_init.household(:,zind,eind,:),[n_a.household,Params.J]))
    end
    title('Value function: median value of z')
    xlabel('Assets (a)')
    ylabel('Age j')
    subplot(2,1,2);
    if small_z_no_e
        surf(a_grid.household*ones(1,Params.J),ones(n_a.household,1)*(agejshifter+(1:1:Params.J)),reshape(V_init.household(:,zind,:),[n_a.household,Params.J]))
    else
        surf(a_grid.household*ones(1,Params.J),ones(n_a.household,1)*(agejshifter+(1:1:Params.J)),reshape(V_init.household(:,zind,eind,:),[n_a.household,Params.J]))
    end
    title('Value function: median value of z')
    xlabel('Assets (a)')
    ylabel('Age in Years')
end


%% Test
disp('Test StationaryDist')
StationaryDist_init=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy_init,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);

%% Test
% Note: Because we used simoptions we must include this as an input
disp('Test AggVars')
if Params.scenario<3
    FnsToEvaluate2=FnsToEvaluate2_12;
elseif Params.scenario<4
    FnsToEvaluate2=FnsToEvaluate2_3;
else
    FnsToEvaluate2=FnsToEvaluate2_4;
end
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_init, Policy_init, FnsToEvaluate2, Params, n_d, n_a, n_z,N_j,Names_i, d_grid, a_grid, z_grid,simoptions);

% Next few lines were used to try a few parameter values so as to get a
% decent initial guess before actually solving the general equilbrium
fprintf('Check: L_h, L_f, K \n')
[AggVars.L_h.Mean,AggVars.L_f.Mean,AggVars.K.Mean]
fprintf('Check: K/L_f (should be about 2.03) \n')
AggVars.K.Mean/AggVars.L_f.Mean
if Params.scenario<3
    fprintf('Check: S \n')
    [AggVars.S.Mean]
elseif Params.scenario<4
    fprintf('Check: S, A, H, PV_h\n')
    [AggVars.S.Mean,AggVars.A.Mean,AggVars.H.Mean,AggVars.PV_h.Mean]
else
    fprintf('Check: S, A, H, PV_h, PV_f \n')
    [AggVars.S.Mean,AggVars.A.Mean,AggVars.H.Mean,AggVars.PV_h.Mean,AggVars.PV_f.Mean]
end
fprintf('Check: ShareIssuance GE condition \n')
Params.P0-((((1-Params.tau_cg)*Params.P0 + (1-Params.tau_d)*Params.D_pp)/(1+Params.r_pp-Params.tau_cg))-AggVars.S.Mean)
end

solve_GE_temp=solve_GE; clear solve_GE
solve_TPath_temp=solve_TPath; clear solve_TPath
save tpathElectrify0.mat
solve_GE=solve_GE_temp; solve_TPath=solve_TPath_temp;
else
    load tpathElectrify0.mat
end % solve_setup
clear solve_GE_temp solve_TPath_temp

%% Solve for the General Equilibrium
if mod(solve_GE,2)==1
    % heteroagentoptions.fminalgo=4 % CMA-ES algorithm 
    
    heteroagentoptions.verbose=1;
    if Params.scenario<3
        FnsToEvaluate=FnsToEvaluate_12;
        heteroagentoptions.toleranceGEprices=10^(-4);
        heteroagentoptions.toleranceGEcondns=10^(-4); % This is the hard one
        if solve_TPath
            % heteroagentoptions.maxiter=200;
        end
    else
        heteroagentoptions.toleranceGEprices=10^(-2);
        heteroagentoptions.toleranceGEcondns=10^(-1); % This is the hard one
        heteroagentoptions.maxiter=35*(1+logical(small_z_no_e)+logical(small_model));                % About 3 hours for 35 iterations

        if Params.scenario<4
            FnsToEvaluate=FnsToEvaluate_3;
            heteroagentoptions.CustomModelStats=@(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions) ...
                Electrify_CustomModelStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions);
        else
            FnsToEvaluate=FnsToEvaluate_4;
            heteroagentoptions.CustomModelStats=@(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions) ...
                Electrify_4CustomModelStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions);
        end
    end

    [p_eqm_init,GEcondns_init]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j,Names_i,[],pi_z,d_grid,a_grid,z_grid,jequaloneDist,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,AgeWeightsParamNames,PTypeDistParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
    % p_eqm contains the general equilibrium parameter values
    % Put this into Params so we can calculate things about the initial equilibrium
    % GEcondns tells us the values of the GeneralEqmEqns, should be near zero
    Params.pension=p_eqm_init.pension;
    Params.AccidentBeqS_pp=p_eqm_init.AccidentBeqS_pp;
    % To use a '_tminus1' variable we must include its initial value
    transpathoptions.initialvalues.BeqleftS_pp=p_eqm_init.AccidentBeqS_pp;
    if Params.scenario>2
        Params.AccidentBeqAH_pp=p_eqm_init.AccidentBeqAH_pp;
        transpathoptions.initialvalues.BeqleftAH_pp=p_eqm_init.AccidentBeqAH_pp;
    end
    Params.G_pp=p_eqm_init.G_pp;
    Params.w=p_eqm_init.w;
    % Params.firmbeta=p_eqm_init.firmbeta;
    Params.P0=p_eqm_init.P0;

    % Re-Calculate a few things related to the general equilibrium.
    [V_init, Policy_init]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    StationaryDist_init=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy_init,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);

    % Calculate various stats
    AllStats_init=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_init,Policy_init,FnsToEvaluate2,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
    % Calculate the life-cycle profiles
    AgeConditionalStats_init=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_init,Policy_init,FnsToEvaluate2,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
    
    % Note: Only part of this initial stationary general eqm we actually 'need'
    % is the agent distribution. Rest is just out of interest.
    
    AgentDist_init=StationaryDist_init; % Just to emphasize that there is no need for the
       % initial agent distribution to be a stationary dist (it is in this
       % example, but does not need to be for transition paths)
    
    % Just to see it...
    GEcondns_init
    
    %%
    ParamPath_temp=ParamPath; clear ParamPath
    solve_GE_temp=solve_GE; clear solve_GE
    solve_TPath_temp=solve_TPath; clear solve_TPath
    if solve_GE_temp<2
        % Keep these parameters in case we want to graph results later
        save tpathElectrifyA.mat
    else
        % We will let the solution to GE_final have final say on Params
        Params_temp=Params; clear Params
        save tpathElectrifyA.mat
        Params=Params_temp; clear Params_temp
    end
    ParamPath=ParamPath_temp; solve_GE=solve_GE_temp; solve_TPath=solve_TPath_temp;
    % load tpathElectrifyA.mat
else
    load tpathElectrifyA.mat
end
clear ParamPath_temp solve_GE_temp solve_TPath_temp

%% Solve for final stationary general eqm with Params at time T
% Must ensure that our T does not conflict with any other dimensions
if small_T==1
    T_end=length(Names_i)+1;
else
    T_end=T;
end
last_n_a_dim=@(n_a_field) n_a_field(end);
last_n_a_dims=structfun(last_n_a_dim, n_a);
while any(ismember(last_n_a_dims,T_end))
    T_end=T_end+1;
end
if T_end>T
    error("impossible dimensions for T and ")
end

% ParamPath on Ek (Energy Use) and ek (Energy Efficiency)
% More transitions down in the demographics section
ParamPath.Ek=linspace(1,1.01,T); Params.Ek=ParamPath.Ek(1);
ParamPath.ek=linspace(1,1.01,T); Params.ek=ParamPath.ek(1);
ParamPath.carbon_tax=linspace(35,2450,T); Params.carbon_tax=ParamPath.carbon_tax(1);
ParamPath.energy_pct_brown=linspace(0.80,0.05,T); Params.energy_pct_brown=ParamPath.carbon_tax(1);

% Model inflation as a series of 10-year supply-side shocks across 100 year transition period
% These are shocks above "normal" cpi inflation
shock_period=10;
shock_years=(1+shock_period:shock_period:max_age+1);
% Exponentially increasing every shock_period years from from ~1% to ~3% after initial shock-free period
shock_pct=[zeros(1,shock_period), repelem(cumsum(exp((shock_years-(shock_period+1))/100)/100),1,shock_period)];

% Translate shock years into periods and periods into transition periods
ParamPath.cpi=shock_pct(1:Params.ypp:Params.J*Params.ypp); % CPI per period j
ParamPath.cpi(end+1:T*jpT)=ParamPath.cpi(end); % CPI extended to the jth period implied by final T
ParamPath.cpi=ParamPath.cpi(1:T); % CPI on a per transition period basis
Params.cpi=ParamPath.cpi(1);

% Steady increase of fossil costs above "normal" cpi inflation
ParamPath.cpi_energy=1.001.^((0:Params.J-1)*Params.ypp)-1; % Params.J periods of energy cost increases
% Translate energy periods (j) into transition periods
ParamPath.cpi_energy(end+1:T*jpT)=ParamPath.cpi_energy(end); % Energy cost increases extended to the jth period implied by final T
ParamPath.cpi_energy=ParamPath.cpi_energy(1:T); % Energy cost increases on a per transition period basis
Params.cpi_energy=ParamPath.cpi_energy(1);

if solve_GE>=2
    % 40 years of changing demographics
    % 60 years in final demographic state (to allow time to converge to final stationary general eqm)
    % Conditional survival probabilities
    ParamPath.sj=[Params.sj_init+(Params.sj_final-Params.sj_init).*linspace(0,1,ceil(40/(Params.ypp*jpT)))'; Params.sj_final.*ones(T-ceil(40/(Params.ypp*jpT)),1)];
    % T-by-N_j (whether this or N_j-by_T, toolkit understands both)
    % Calculate the implied mewj from the sj
    ParamPath.mewj=cumprod([ones(T,1), ParamPath.sj(:,1:end-1)], 2); % mass of age jj is the mass of jj-1 that survive
    % Factor in population growth; In N_j dimension, older people are from earlier (smaller) populations
    % ...in the T dimension, we see overall population growth as T increases
    ParamPath.mewj=ParamPath.mewj./((1+Params.n_pp).^((1:Params.J)-1)); % Population shrinks in the N_j dimension
    ParamPath.mewj=ParamPath.mewj.*((1+Params.n_pp).^(jpT*((1:T)-1)))'; % Population grows in the T dimension
    ParamPath.mewj=ParamPath.mewj./sum(ParamPath.mewj,2); % normalize age-masses to sum to one
    % Looking at ParamPath.mewj you can see that as tt increases, the mass at older ages increases

    Params.Ek=ParamPath.Ek(T_end);
    Params.ek=ParamPath.ek(T_end);
    Params.carbon_tax=ParamPath.carbon_tax(T_end);
    Params.energy_pct_brown=ParamPath.energy_pct_brown(T_end);
    Params.sj=ParamPath.sj(T_end,:); % conditional survival probabilities
    Params.mewj=ParamPath.mewj(T_end,:);
    Params.cpi=ParamPath.cpi(T_end);
    Params.cpi_energy=ParamPath.cpi_energy(T_end);
    Params.Lhscale=ParamPath.Lhscale(T_end);

    %% Let's take a quick look at what we have calculated, namely V and Policy

    % Evaluate the final stationary general eqm
    disp('Test ValueFnIter')
    [V_final, Policy_final]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
    disp('Test StationaryDist')
    StationaryDist_final=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy_final,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
    disp('Test AggVars')
    AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_final, Policy_final, FnsToEvaluate2, Params, n_d, n_a, n_z,N_j,Names_i, d_grid, a_grid, z_grid,simoptions);

    % Next few lines were used to try a few parameter values so as to get a
    % decent initial guess before actually solving the general equilbrium
    fprintf('Check: L_h, L_f, K \n')
    [AggVars.L_h.Mean,AggVars.L_f.Mean,AggVars.K.Mean]
    fprintf('Check: K/L_f (should be about 2.03) \n')
    AggVars.K.Mean/AggVars.L_f.Mean
    if Params.scenario<3
        fprintf('Check: S \n')
        [AggVars.S.Mean]
    elseif Params.scenario<4
        fprintf('Check: S, A, H, PV_h\n')
        [AggVars.S.Mean,AggVars.A.Mean,AggVars.H.Mean,AggVars.PV_h.Mean]
    else
        fprintf('Check: S, A, H, PV_h, PV_f \n')
        [AggVars.S.Mean,AggVars.A.Mean,AggVars.H.Mean,AggVars.PV_h.Mean,AggVars.PV_f.Mean]
    end
    fprintf('Check: ShareIssuance GE condition \n')
    Params.P0-((((1-Params.tau_cg)*Params.P0 + (1-Params.tau_d)*Params.D_pp)/(1+Params.r_pp-Params.tau_cg))-AggVars.S.Mean)

    % And now, the GE for the final conditions!
    [p_eqm_final,GEcondns_final]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,[],pi_z,d_grid,a_grid,z_grid,jequaloneDist,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,AgeWeightsParamNames,PTypeDistParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
    % Done, the general eqm prices are in p_eqm
    % GEcondns tells us the values of the GeneralEqmEqns, should be near zero
    Params.pension=p_eqm_final.pension;
    Params.AccidentBeqS_pp=p_eqm_final.AccidentBeqS_pp;
    if Params.scenario>2
        Params.AccidentBeqAH_pp=p_eqm_final.AccidentBeqAH_pp;
    end
    Params.G_pp=p_eqm_final.G_pp;
    Params.w=p_eqm_final.w;
    % Params.firmbeta=p_eqm_final.firmbeta;
    Params.P0=p_eqm_final.P0;

    % Calculate various stats
    AllStats_final=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_final, Policy_final, FnsToEvaluate2, Params, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid, z_grid,simoptions);
    % Calculate the life-cycle profiles
    AgeConditionalStats_final=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_final,Policy_final, FnsToEvaluate2,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

    % Note: Only part of this final stationary general eqm we actually 'need'
    % is the value fn (although we likely want p_eqm_final for initial guess of PricePath0). 
    % Rest is just out of interest.
    
    % Double-check that the general eqm is accurate before we start the
    % transition path, because if it is not then it won't solve
    GEcondns_final

    %%
    solve_TPath_temp=solve_TPath; clear solve_TPath
    save tpathElectrifyB.mat
    solve_TPath=solve_TPath_temp;
else
    load tpathElectrifyB.mat
end % solve_GE_final
clear solve_TPath_temp

if solve_TPath
    %% Setup for the transition path
    % T=100; % number of periods for transition path
    
    % Already created ParamPath.sj and ParamPath.mewj
    
    % Initial guess for general eqm parameters; if small_T, V_final is V_init at T=T_end
    if small_T==1
        T_eq=ceil(T_end/2);
        paramnames=fieldnames(ParamPath);
        for nn=1:length(paramnames)
            if size(ParamPath.(paramnames{nn}),1)==1
                ParamPath0.(paramnames{nn})=ParamPath.(paramnames{nn})(1,1:T_end);
            else
                ParamPath0.(paramnames{nn})=ParamPath.(paramnames{nn})(1:T_end,:);
            end
        end
    else
        T_eq=ceil(0.618*T); % Demographic change has stopped and T_eq begins period of transition equilibrium-finding
        ParamPath0=ParamPath;
    end
    PricePath0.w=[linspace(p_eqm_init.w, p_eqm_final.w,T_eq), p_eqm_final.w*ones(1,T_end-T_eq)];
    % PricePath0.firmbeta=[linspace(p_eqm_init.firmbeta, p_eqm_final.firmbeta,T_eq), p_eqm_final.firmbeta*ones(1,T_end-T_eq)];
    PricePath0.P0=[linspace(p_eqm_init.P0, p_eqm_final.P0,T_eq), p_eqm_final.P0*ones(1,T_end-T_eq)];
    PricePath0.pension=[linspace(p_eqm_init.pension, p_eqm_final.pension,T_eq), p_eqm_final.pension*ones(1,T_end-T_eq)];
    PricePath0.AccidentBeqS_pp=[linspace(p_eqm_init.AccidentBeqS_pp,p_eqm_final.AccidentBeqS_pp,T_eq), p_eqm_final.AccidentBeqS_pp*ones(1,T_end-T_eq)];
    if Params.scenario>2
        PricePath0.AccidentBeqAH_pp=[linspace(p_eqm_init.AccidentBeqAH_pp,p_eqm_final.AccidentBeqAH_pp,T_eq), p_eqm_final.AccidentBeqAH_pp*ones(1,T_end-T_eq)];
    end
    PricePath0.G_pp=[linspace(p_eqm_init.G_pp, p_eqm_final.G_pp,T_eq), p_eqm_final.G_pp*ones(1,T_end-T_eq)];
    % PricePath0.TargetKdivL=2.03*ones(1,T_end);

    % General eqm eqns, same idea as with the stationary general eqm
    % GeneralEqmEqns_Transition.capitalmarket=@(r_pp,alpha_k,alpha_l,delta,K,L,ypp) r_pp-(alpha_k*(K^(alpha_k-1))*(L^(alpha_l))-((delta+1)^ypp-1)); % r=marginal product of capital
    GeneralEqmEqns_Transition.labormarket=@(w,alpha_k,alpha_l,K,L_f) w-(alpha_l)*(K^alpha_k)*(L_f^(alpha_l-1)); % w=marginal product of labor
    % GeneralEqmEqns_Transition.firmdiscounting=GeneralEqmEqns.firmdiscounting;
    % GeneralEqmEqns_Transition.dividends=GeneralEqmEqns.dividends;
    GeneralEqmEqns_Transition.ShareIssuance=GeneralEqmEqns.ShareIssuance;
    GeneralEqmEqns_Transition.pensions=GeneralEqmEqns.pensions;
    GeneralEqmEqns_Transition.govbudgetbalance=GeneralEqmEqns.govbudget;
    % Note: bequests are left in t-1 and received in t
    GeneralEqmEqns_Transition.bequestsS_pp=@(BeqleftS_pp_tminus1,AccidentBeqS_pp,n_pp) BeqleftS_pp_tminus1/(1+n_pp)-AccidentBeqS_pp; % Accidental share bequests received equal accidental share bequests left
    if Params.scenario>2
        GeneralEqmEqns_Transition.bequestsAH_pp=@(BeqleftAH_pp_tminus1,AccidentBeqAH_pp,n_pp) BeqleftAH_pp_tminus1/(1+n_pp)-AccidentBeqAH_pp; % Accidental asset+house bequests received equal accidental asset+house bequests left
    end
    
    % Note: in this example these are actually identical to the general eqm
    % eqns for the stationary general eqm, but that is not often the case.
    
    % Set up the shooting algorithm
    transpathoptions.GEnewprice=3;
    % Need to explain to transpathoptions how to use the GeneralEqmEqns to update the general eqm transition prices (in PricePath).
    transpathoptions.GEnewprice3.howtoupdate=... % a row is: GEcondn, price, add, factor
        {'labormarket','w',0,0.03;... % labormarket GE condition will be positive if w is too big, so subtract
        ... % 'firmdiscounting','firmbeta',0,0.03;... % firmdiscounting GE condition will be positive if firmbeta is too big, so subtract
        'ShareIssuance','P0',0,0.03;... % ShareIssuance GE condition will be positive if P0 is too big, so subtract
        'pensions','pension',0,0.03;... % pensions GE condition will be positive if pension is too big, so subtract
        'govbudgetbalance','G_pp',0,0.03;... % govbudget GE condition will be positive if G_pp is too big, so subtract
        'bequestsS_pp','AccidentBeqS_pp',1,0.03;... % bequests GE condition will be negative if BeqS_pp is too big, so add
        'bequestsAH_pp','AccidentBeqAH_pp',1,0.03;... % bequests GE condition will be negative if BeqAH_pp is too big, so add
        };
    if Params.scenario<3
        mask=strcmp(transpathoptions.GEnewprice3.howtoupdate(:,1),'bequestsAH_pp');
        transpathoptions.GEnewprice3.howtoupdate(mask,:)=[];
    elseif Params.scenario==4
        for pp=1:size(transpathoptions.GEnewprice3.howtoupdate,1)
            transpathoptions.GEnewprice3.howtoupdate{pp,4}=0.008;
        end
    end

    % Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
    % Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
    % A small 'factor' will make the convergence to solution take longer, but too large a value will make it 
    % unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.

    % TESTING -- pensions doesn't depend on PType
    % transpathoptions.GEptype={'pensions'};
    
    %% Solve the transition path
    % Setup the options relating to the transition path
    transpathoptions.verbose=1;
    transpathoptions.maxiter=100; % default is 1000
    transpathoptions.fastOLG=0; % PTypes will force this on `simoptions`; must we match that energy?
    transpathoptions.graphpricepath=1; % plots of the ParamPath that get updated every interation
    transpathoptions.graphaggvarspath=1; % plots of the AggVarsPath that get updated every iteration
    
    % Running, it was about stuck iterating around 2 or 3*10^(-4) but had clearly solved. So
    transpathoptions.tolerance=4*10^(-3); % default is 10^(-4), which is a very demanding accuracy

    %%

    save tpathElectrifyC.mat
    % load tpathElectrifyC.mat

    % And go! (with FnsToEvaluate2)
    [PricePath,GECondnsPath]=TransitionPath_Case1_FHorz_PType(PricePath0, ParamPath0, T_end, V_final, AgentDist_init, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate2, GeneralEqmEqns_Transition, Params, DiscountFactorParamNames, AgeWeightsParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);

    %%
    solve_TPath_temp=solve_TPath; clear solve_TPath
    save tpathElectrifyD.mat
    solve_TPath=solve_TPath_temp;
else
    load tpathElectrifyD.mat
end % solve_TPath
clear solve_TPath_temp

    %% Now calculate some things about the transition path (path for Value fn, Policy fn, Agent Distribution)
    % You can calculate the value and policy functions for the transition path
    [VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath0, T_end, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
    
    % You can then use these to calculate the agent distribution for the transition path
    AgentDistPath=AgentDistOnTransPath_Case1_FHorz_PType(StationaryDist_init, jequaloneDist, PricePath, ParamPath0, PolicyPath, AgeWeightsParamNames,n_d,n_a,n_z,N_j,Names_i,pi_z,T_end, Params, transpathoptions, simoptions);
    
    %% Analyse the transition path
    % And then we can calculate AggVars for the path
    AggVarsPath=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath,PolicyPath, PricePath, ParamPath0, Params, T_end, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptions, simoptions);
    
    %% Plot some paths
    figure(1)
    % Plot of K and w
    % Note: include periods -3 to 0 (the initial stationary eqm) so can see any jump in period 1
    subplot(2,1,1); plot(1:1:T_end,AggVarsPath.K.Mean)
    hold on
    plot(-3:1:0,AllStats_init.K.Mean*ones(1,4),'w')
    hold off
    xlim([-3,T_end])
    title('Path of aggregate capital (K)')
    subplot(2,1,2); plot(1:1:T_end,PricePath.w)
    hold on
    plot(-3:1:0,p_eqm_init.w*ones(1,4),'w')
    hold off
    xlim([-3,T_end])
    title('Path of wage rate (w)')

% Can just use the same FnsToEvaluate as before
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_init,Policy_init,FnsToEvaluate2,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

if max(AgeConditionalStats.S.Maximum)==share_grid(end)
    warning("share_grid maximum reached")
end
if Params.scenario>2
    if max(AgeConditionalStats.A.Maximum)==asset_grid(end)
        warning("asset_grid maximum reached")
    end
    if max(AgeConditionalStats.H.Maximum)==house_grid(end)
        warning("house_grid maximum reached")
    end
    if max(AgeConditionalStats.PV_h.Maximum)==pv_grid_hh(end)
        warning("pv_grid_hh maximum reached")
    end
end

%% Plot the life cycle profiles of capital and labour for the inital and final eqm.

figure_c=figure(10);
if Params.scenario<3
    subplot(3,1,1); plot(1:1:Params.J,AgeConditionalStats.L_h.Mean)
    title('Life Cycle Profile: Effective Labour Supply')
    subplot(3,1,2); plot(1:1:Params.J,AgeConditionalStats.S.Mean)
    title('Life Cycle Profile: Share holdings')
    subplot(3,1,3); plot(1:1:Params.J,Params.kappa_j)
    title('Life Cycle Profile: kappa_j')
else
    subplot(3,2,1); plot(1:1:Params.J,AgeConditionalStats.L_h.Mean)
    title('Life Cycle Profile: Effective Labour Supply')
    subplot(3,2,3); plot(1:1:Params.J,AgeConditionalStats.S.Mean)
    title('Life Cycle Profile: Share holdings')
    subplot(3,2,5); plot(1:1:Params.J,Params.kappa_j)
    title('Life Cycle Profile: kappa_j')
end
if Params.scenario>2
    subplot(3,2,2); plot(1:1:Params.J,AgeConditionalStats.A.Mean)
    title('Life Cycle Profile: Asset holdings')
    subplot(3,2,4); plot(1:1:Params.J,AgeConditionalStats.H.Mean)
    title('Life Cycle Profile: House holdings')
    subplot(3,2,6); plot(1:1:Params.J,AgeConditionalStats.PV.Mean)
    title('Life Cycle Profile: Solar PV installed')
end
saveas(figure_c,'./SavedOutput/Graphs/Electrify_LifeCycleProfiles','pdf')

%% Calculate some aggregates and print findings about them

% Add consumption to FnsToEvaluate2
if Params.scenario<3
    FnsToEvaluate2.Consumption.household=@( ...
            labor,sprime,s,z,e, ...
            pension,AccidentBeqS_pp,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,cpi,rentprice,energy_pct_cost ...
        ) Electrify_HouseholdConsumptionFn( ...
            labor,0,sprime,0,0,s,0,0,0,z,e, ...
            pension,AccidentBeqS_pp,0,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,0,cpi,rentprice,0,0,energy_pct_cost);
    FnsToEvaluate2.Income.household=@( ...
            labor,sprime,s,z,e, ...
            pension,AccidentBeqS_pp,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,cpi,energy_pct_cost ...
        ) Electrify_HouseholdIncomeFn( ...
            labor,0,sprime,0,0,s,0,0,0,z,e, ...
            pension,AccidentBeqS_pp,0,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,cpi,energy_pct_cost);
elseif Params.scenario<4
    FnsToEvaluate2.Consumption.household=@( ...
            labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
            pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,r_wedge_pp,f_htc,rentprice,cpi_energy,pv_pct_cost,energy_pct_cost ...
        ) Electrify_HouseholdConsumptionFn( ...
            labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
            pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,r_wedge_pp,f_htc,rentprice,cpi_energy,pv_pct_cost,energy_pct_cost);
    FnsToEvaluate2.Income.household=@( ...
            labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
            pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,cpi_energy,energy_pct_cost ...
        ) Electrify_HouseholdIncomeFn( ...
            labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
            pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,cpi_energy,energy_pct_cost);
else
    FnsToEvaluate2.Consumption.household=@( ...
            labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e, ...
            pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,r_wedge_pp,f_htc,rentprice,cpi_energy,pv_pct_cost,energy_pct_cost,energy_pct_brown,carbon_tax ...
        ) Electrify_4HouseholdConsumptionFn( ...
            labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e, ...
            pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,r_wedge_pp,f_htc,rentprice,cpi_energy,pv_pct_cost,energy_pct_cost,energy_pct_brown,carbon_tax);
    FnsToEvaluate2.Income.household=@( ...
            labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e, ...
            pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,cpi_energy,energy_pct_cost,energy_pct_brown,carbon_tax ...
        ) Electrify_4HouseholdIncomeFn( ...
            labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e, ...
            pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
            kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
            r_pp,cpi_energy,energy_pct_cost,energy_pct_brown,carbon_tax);
end

AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_init, Policy_init, FnsToEvaluate2, Params, n_d, n_a, n_z,N_j, Names_i, d_grid, a_grid, z_grid,simoptions);

Y=AggVars.Output.Mean;

P=((1-Params.tau_cg)*Params.P0 + (1-Params.tau_d)*Params.D_pp)/(1+Params.r_pp-Params.tau_cg);

% Calculate the aggregate TFP as output/((capital^alpha_k)*(labor^alpha_l))
AggregateTFP=Y/((AggVars.K.Mean^Params.alpha_k)*(AggVars.L_f.Mean^Params.alpha_l));

% Total value of firms
temp=V_init.firm.*StationaryDist_init.firm;
temp(StationaryDist_init.firm==0)=0; % Get rid of points that have V=-inf but zero mass which would give nan
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
