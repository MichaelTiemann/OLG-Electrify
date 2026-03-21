%% OLG Electrification (based on OLGModels14: Heterogenous households and heterogeneous firms
%% and also Life-Cycle Model 35: Portfolio-Choice with Housing)
% See https://www.vfitoolkit.com/updates-blog/2021/an-introduction-to-life-cycle-models/
% OLGModel14.m in the repo https://github.com/vfitoolkit/IntroToOLGModels
% and LifeCycleModel35.m in the repo https://github.com/vfitoolkit/IntroToLifeCycleModels

% A line some need for running on the Server
addpath(genpath('./MatlabToolkits/'))

%% Basic statistical abstract (NZD)
% NZ GDP: $440B ($80K per capita, $152K per employed worker)
% NZ Wages: $60K living, $70K median, $80K average * 2.9M workers = $232B wages
% NZ Energy:    - 525 PJ/year
%   Oil         - 270 PJ
%   Electricity - 144 PJ
%   Gas         -  58 PJ
%   Biomass     -  40 PJ
%   Coal        -  18 PJ
% NZ Electricity retail: $350/MWh
% NZ HH Energy: 20 kWh/day electricity =>  7 MWh/year =>  $2500/year => 3.5% wages
% NZ HH Energy: 73 kWh/day overall     => 27 MWh/year => $10000/year => 14.0% wages
% NZ Firm Energy retail: $150/MWh
%   Transport    - 200 PJ
%   Industrial   - 160 PJ
%   Commercial   -  55 PJ
%   Ag,Forest,Fish- 30 PJ
%   Total: 445 PJ => 125 TWh => $20B energy costs => 4.5% of 440B GDP
% Energy is 8.6% of Labor costs
% Net capital stocks of NZ $1,329B less $690B real estate = $630B
% K/L = $630B/232B = 2.72

solve_setup=true;
solve_GE_init=true;
solve_GE_final=true;

solve_TPath=true;
% If true, shrink n_z down to 3 (the min for discretization)
% and make e parameter always zero (no e_grid).firm
small_z_no_e=false;
solve_demographic_change=true;

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
max_age=[100,100,100,100];
Params.agejshifter=19; % Age 20 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=ceil((max_age(Params.scenario)-Params.agejshifter)/Params.ypp); % =60/ypp, Number of period in life-cycle
N_j.household=Params.J; % Number of periods in finite horizon

jpT=1; % Default: one transition period=1 time period; Could have multiple j's per T
T=ceil(Params.J*1.4/jpT);
if T==Params.J
    % The toolkit thinks that T and J must be different (T larger to reach equilibrium post J)
    T=T+1;
end

% ParamPath on Ek (Energy Use) and ek (Energy Efficiency)
% More transitions down in the demographics section
ParamPath.Ek=linspace(1,1.2,T); Params.Ek=ParamPath.Ek(1);
ParamPath.ek=linspace(1,1.5,T); Params.ek=ParamPath.ek(1);
ParamPath.carbon_tax=linspace(35,2450,T); Params.carbon_tax=ParamPath.carbon_tax(1);
ParamPath.energy_pct_brown=linspace(80,20,T); Params.energy_pct_brown=ParamPath.carbon_tax(1);

% Model inflation as a series of 10-year supply-side shocks across 100 year transition period
% These are shocks above "normal" cpi inflation
shock_period=10;
shock_years=(1+shock_period:shock_period:max_age(Params.scenario)+1);
% Exponentially increasing every shock_period years from from ~2% to ~7% after initial shock-free period
shock_pct=[zeros(1,shock_period), repelem(cumsum(exp((shock_years-(shock_period+1))/64)/50),1,shock_period)];

% Translate shock years into periods and periods into transition periods
ParamPath.cpi=shock_pct(1:Params.ypp:Params.J*Params.ypp); % CPI per period j
ParamPath.cpi(end+1:T*jpT)=ParamPath.cpi(end); % CPI extended to the jth period implied by final T
ParamPath.cpi=ParamPath.cpi(1:jpT:T); % CPI on a per transition period basis
Params.cpi=ParamPath.cpi(1);

% Steady increase of fossil costs above "normal" cpi inflation
ParamPath.cpi_energy=1.01.^((0:Params.J-1)*Params.ypp)-1; % Params.J periods of energy cost increases
% Translate energy periods (j) into transition periods
ParamPath.cpi_energy(end+1:T*jpT)=ParamPath.cpi_energy(end); % Energy cost increases extended to the jth period implied by final T
ParamPath.cpi_energy=ParamPath.cpi_energy(1:jpT:T); % Energy cost increases on a per transition period basis
Params.cpi_energy=ParamPath.cpi_energy(1);

%% Grid sizes to use for household
if Params.scenario<3
    n_d.household=101;
    n_a.household=201;
    n_z.household=3+2*floor(1.7*log(min(Params.J,60))); % AR(1) with age-dependent params = 15 with 60 periods
    vfoptions.lowmemory.household=0;
else
    n_d.household=[21,5]; % Decisions: labor, buyhouse (5)
    if Params.scenario<4
        % 21,33,4,5 => labor strike
        % 15,23,4,5 => ok
        n_a.household=[5,31,4,5]; % Endogenous shares, assets (>=6), housing (>=2), and solarpv (>=2) assets (0-60 kW generation)
    else
        n_a.household=[5,31,3,4,3]; % Endogenous shares, assets (>=6), car (3), housing (>=2), and solarpv (>=2) assets (0-60 kW generation)
    end
    n_z.household=1+2*floor(1.2*log(min(Params.J,60))); % AR(1) with age-dependent params = 7 with 60 periods
    if Params.scenario<4
        if small_z_no_e
            vfoptions.lowmemory.household=0;
        else
            vfoptions.lowmemory.household=2;
        end
    else
        vfoptions.lowmemory.household=3;
    end
end
if small_z_no_e
    n_z.household=1;
    Params.e=0;
else
    % Exogenous labor productivity units shocks (next two lines)
    vfoptions.n_e.household=3; % iid
end

%% Grids to use for firm
if Params.scenario<4
    n_d.firm=101; % Dividend payment
    n_a.firm=201; % Capital holdings
else
    n_d.firm=0; % Not Yet Used: Electrification investment
    n_a.firm=[51,42]; % Capital holdings and PV assets
end
if small_z_no_e
    n_z.firm=1;
else
    n_z.firm=3+2*floor(log(min(Params.J,60))); % Productivity shock; scaled to model, not firm horizon
end
N_j.firm=Inf; % Infinite horizon
vfoptions.lowmemory.firm=logical(Params.scenario==4 && ~small_z_no_e);

%% Grids to use for energy
if Params.scenario<4
    n_d.energy=0; % What decisions?
    n_a.energy=1; % What assets?
else
    n_d.energy=101; % Invest in PV
    n_a.energy=202; % PV assets
end
if small_z_no_e
    n_z.energy=1;
else
    n_z.energy=3+2*floor(log(min(Params.J,60))); % Productivity shock; scaled to model, not firm horizon
end
N_j.energy=Inf; % Infinite horizon
vfoptions.lowmemory.energy=0;

%% Global parameters (applies to household and firm)
% Note: with w=1, labor tax=20%, kappa_j(1)=0.5, agents have 0.4 budget to start
% If rent=0.3 energy=0.1, they live, but cannot save; they strike

% Annual risk-free rate of return
r=0.05; % We will discover risk-free rate of return per period in GE
r_wedge=0.05; Params.r_wedge_pp=(1+r_wedge)^Params.ypp-1;

Lhscale=[0.25,0.25,0.21,0.21]; % Scaling the household labor supply; we scale model and GE finds its own equilibrium
if Params.scenario>2 && Params.ypp>1
    if Params.ypp<8
        Lhscale(3:4)=0.38*ones(1,2);
    elseif Params.ypp<11
        Lhscale(3:4)=0.28*ones(1,2);
    else
        Lhscale(3:4)=0.23*ones(1,2);
    end
    if Params.scenario==4 && Params.ypp==5
        if small_z_no_e
            Lhscale(4)=1.3;
        else
            Lhscale(4)=0.5;
        end
    end
end

ParamPath.Lhscale=linspace(Lhscale(Params.scenario),2,T);
Params.Lhscale=ParamPath.Lhscale(1);

%% Parameters for households
% Discount rate; Changed to get S to increase nearer to 1 given r=0.05
% (ran it with beta=0.99, got S=0.3, so increased this; note that it interacts with sj to give the actual discount factor)
beta=[0.95,0.95,0.99,0.99];
Params.beta_pp = beta(Params.scenario)^Params.ypp;

% Housing
% Params.minhouse % set below, is the minimum value of house that can be purchased
rentprice=[0,0.3,0.3,0.3]; % I figured setting rent a decent fraction of income is sensible
Params.rentprice=rentprice(Params.scenario);
houseservices=[0,0.5,0.5,0.5]; % housing services as a fraction of house value
Params.houseservices=houseservices(Params.scenario);
energy_pct_cost=[0,0.07,0.07,0.05]; % Electricity: 3%; Gas: 1-2%; Petrol: 1-2%; Scenario 4 disaggregates petrol from this cost
Params.energy_pct_cost=energy_pct_cost(Params.scenario);
if Params.scenario>2
    Params.f_htc=0.05; % transaction cost of buying/selling house (is a percent of h+hprime)
    f_coll=[0,0,0.5,0.5]; % collateral contraint (fraction of house value that can be borrowed)
    Params.f_coll=f_coll(Params.scenario);
    pv_pct_cost=[0,0,0.05,0.05]; % modeling a $30K install for a $600K house
    Params.pv_pct_cost=pv_pct_cost(Params.scenario);
end

% Preferences
Params.sigma = 2; % Coeff of relative risk aversion (curvature of consumption)
sigma_h=[0,0,0.5,0.2];
Params.sigma_h=sigma_h(Params.scenario); % Relative importance of housing services (vs consumption) in utility
sigma_c=[0,0,0.5,0.7];
Params.sigma_c=sigma_c(Params.scenario); % Relative importance of housing services (vs consumption) in utility
Params.eta = 1.5; % Curvature of leisure (This will end up being 1/Frisch elasty)
psi = [2, 1, 1, 1]; % Weight on leisure
Params.psi=psi(Params.scenario);

% Labor productivity at start, peak, and end of working life
k_j1 = [0.5, 0.5, 0.5, 0.5];
k_j2 = [2, 2, 2, 2];
k_j2_length = [0,0,5,5];
k_j3 = [1, 1, 1, 1];

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

Params.carservices_j=0.1*ones(Params.J,1);
% Cars start to be useful as people ramp up family life
age1=ceil((28-Params.agejshifter)/Params.ypp);
age2=ceil((32-Params.agejshifter)/Params.ypp);
Params.carservices_j(age1:age2)=linspace(0.5,2,age2-age1+1);
age1=ceil((32-Params.agejshifter)/Params.ypp);
age2=ceil((44-Params.agejshifter)/Params.ypp);
Params.carservices_j(age1:age2)=2*ones(age2-age1+1,1);
age1=ceil((44-Params.agejshifter)/Params.ypp);
age2=ceil((65-Params.agejshifter)/Params.ypp);
Params.carservices_j(age1:age2)=linspace(2,1,age2-age1+1);
age1=ceil((65-Params.agejshifter)/Params.ypp);
age2=ceil((80-Params.agejshifter)/Params.ypp);
Params.carservices_j(age1:age2)=linspace(1,0,age2-age1+1);
Params.carservices_j(age2+1:end)=0;
Params.carservices_j=Params.carservices_j*2;

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
% Note: when Params.ypp==1, the product over the reshaped array is over a single year period (i.e. trivial)
sj_init=prod(1-reshape(dj(1:Params.ypp*Params.J),[Params.ypp,Params.J]),1); % p5-year survival rates
sj_init(end)=0; % In the present model the last period (j=J) value of sj is actually irrelevant

% Add 5 years of life expectancy...age sj(65) in the future will be sj(60) by today's statistics
% Part of this is achieved by improving early childhood survival as well...feeding two birds with one worm
sj_final=prod(1-reshape([dj(1:2:10), repelem(dj(11:15), 3), dj(16:Params.ypp*Params.J-5)],[Params.ypp,Params.J]),1);
sj_final(end)=0; % In the present model the last period (j=J) value of sj is actually irrelevant

if ~solve_demographic_change
    sj_final=sj_init;
end

%% Setup for sj and mewj transitions (T-by-N_j)
% We defer doing transition maths until we calculate GE final
Params.sj=sj_init;
Params.mewj=cumprod([1,Params.sj(1:end-1)],2); % mass of age jj is the mass of jj-1 that survive
Params.mewj=Params.mewj./((1+Params.n_pp).^(Params.ypp*((1:Params.J)-1))); % Population shrinks in the N_j dimension
Params.mewj=Params.mewj./sum(Params.mewj); % normalize age-masses to sum to one

% Note: This is rather incomplete, as really you should also have
% population growth rate n. But this does not change any thing in terms of
% the 'objects to compute'. Instead you need to renormalize the model for
% the population growth, and this just means you get a 'n' appearing in
% some equations below. But other than 'n' in some equations, the way you
% do this with the toolkit does not change.

% Warm glow of bequest
Params.warmglow1=0.3; % (relative) importance of bequests
Params.warmglow2=3; % bliss point of bequests (essentially, the target amount)
Params.warmglow3=Params.sigma; % By using the same curvature as the utility of consumption it makes it much easier to guess appropraite parameter values for the warm glow

% The warmglow parameters will help us find the GE solution to actual bequest rates/values
AccidentBeqS=[0.02,0.02,0.02,0.02]; % Accidental bequests (this is the lump sum transfer of shares)
Params.AccidentBeqS_pp=AccidentBeqS(Params.scenario)*Params.ypp;
if Params.scenario>2
    AccidentBeqAH=[0,0,0.02,0.02]; % Accidental bequests (this is the lump sum transfer of assets+house value)
    Params.AccidentBeqAH_pp=AccidentBeqAH(Params.scenario)*Params.ypp;
end

% Taxes
Params.tau_l = 0.2; % Tax rate on labour income

%% Parameters for firm
% Production
Params.alpha_k=0.311; % diminishing returns to capital and energy inputs
Params.alpha_l=0.650; % diminishing returns to labor input
Params.delta=0.054; % Annual depreciation of physical capital
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

%% Parameters for energy
% Idiosyncatic productivity shocks
Params.rho_z_energy=0.767;
Params.sigma_z_e_energy=0.211;

% Set the firm discount factor below (as it is determined in general eqm)
% Params.firmbeta=1/(1+Params.r_pp/(1-Params.tau_cg)); % 1/(1+r_pp) but returns net of capital gains tax

%% Remaining Parameters will be set in GE below

%% Grids for household

% Grid for labour choice
labor_grid=linspace(0,1,n_d.household(1))'; % Notice that it is imposing the 0<=h<=1 condition implicitly

% Grid for share holdings, always > 0
% For later scenarios, shrink the grid for more accuracy
s_grid_cubed=linspace(0,1,ceil(n_a.household(1)/2)).^3; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.
s_grid_linear=linspace(1,10,floor(n_a.household(1)/2)+1);
% share_grid=[s_grid_cubed, s_grid_linear(2:end)]';
share_grid=16*linspace(0,1,n_a.household(1))';

% Set up d for VFI Toolkit (is the two decision variables)
if Params.scenario<3
    d_grid.household=labor_grid;
    a_grid.household=share_grid;
    Params.minhouse=1;
else
    % Grid for bank account; a negative balance implies a mortgage
    a_grid_cubed=linspace(-1,0,ceil(n_a.household(2)/2)-1).^3;
    a_grid_linear=linspace(0,16,floor(n_a.household(2)/2)+2);
    asset_grid=[a_grid_cubed, a_grid_linear(2:end)]';
    
    % Make it so that there is a zero assets
    % Find closest to zero assets
    [~,zeroassetindex]=min(abs(asset_grid));
    asset_grid(zeroassetindex)=0;
    
    % age20avgincome=Params.w*Params.kappa_j(1);
    % house_grid=[0; logspace(2*age20avgincome, 12*age20avgincome, 5)'];
    if Params.scenario<4
        house_grid=(0:1:n_a.household(3)-1)';
    else
        house_grid=(0:1:n_a.household(4)-1)';
    end

    % Note, we can see from w*kappa_j*z and the values of these, that average
    % income is going to be around one, so will just use this simpler house grid
    % [We can think about the values of the house_grid as being relative the average income (or specifically average at a given age)]
    Params.minhouse=house_grid(2); % first is zero (no house)
    
    if Params.scenario<4
        car_grid=zeros(0);
    else
        car_grid=(0:1:n_a.household(3)-1)'; % car assets: no car; petrol car; EV car
    end
    
    % buyhouse decisions
    %  0=no house
    %  1=buy house w/o pv this period
    %  2=buy house w/ pv this period
    %  3=keep house; no pv upgrade
    %  4=keep house; pv upgrade (if possible)
    %  5=testing (not used)
    buyhouse_grid=(0:1:n_d.household(2)-1)';
    
    % kW of solar generation installed, 10 kW per grid element
    if Params.scenario<4
        pv_grid_hh=(0:1:n_a.household(4)-1)';
    else
        pv_grid_hh=(0:1:n_a.household(5)-1)';
    end
    
    d_grid.household=[labor_grid; buyhouse_grid];
    a_grid.household=[share_grid; asset_grid; car_grid; house_grid; pv_grid_hh];
    
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
    if Params.scenario<4
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
end

% z_grid and pi_z for household
z_grid.household=z_grid_J;
pi_z.household=pi_z_J;

if ~small_z_no_e
    % Any (iid) e variable always has to go into vfoptions and simoptions
    vfoptions.e_grid.household=e_grid_J;
    vfoptions.pi_e.household=pi_e_J;
    simoptions.n_e.household=vfoptions.n_e.household;
    simoptions.e_grid.household=e_grid_J;
    simoptions.pi_e.household=pi_e_J;
end


%% Grids for firm
if Params.scenario<4
    d_grid.firm=linspace(0,1+floor(log(Params.ypp)),n_d.firm)'; % Notice that it is imposing the d>=0 condition implicitly
    % k_max=10 replicates OLGModel14; K>4=infeasible when ypp=1, but need more as ypp increases
    k_max=[10,6+ceil(log(Params.ypp)),10+ceil(log(Params.ypp)),10+ceil(log(Params.ypp))];
    k_grid_cubed=linspace(0,1,ceil(n_a.firm/2)).^3; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.
    k_grid_linear=linspace(1,k_max(Params.scenario),floor(n_a.firm/2)+1);
    k_grid=[k_grid_cubed, k_grid_linear(2:end)];
    a_grid.firm=k_grid';
else
    d_grid.firm=linspace(0,1,n_d.firm(1))'; % Electrification investment
    % k_max=10 replicates OLGModel14; K>4=infeasible when ypp=1, but need more as ypp increases
    k_max=6+ceil(log(Params.ypp));
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
if Params.scenario < 4
    d_grid.energy=0; % Notice that it is imposing the d>=0 condition implicitly
    a_grid.energy=linspace(0,1,n_a.energy)'; % Nothing in particular
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


%% Now, create the return function

% For households
DiscountFactorParamNames.household={'beta_pp','sj'};

if Params.scenario<3
    % Hardwire buyhouse, hprime, aprime, h, a, and solarpv to zero
    ReturnFn.household=@( ...
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
elseif Params.scenario<4
    % Notice we use 'Electrify_HouseholdReturnFn'
    ReturnFn.household=@( ...
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
else
    ReturnFn.household=@( ...
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
end

% For firms
DiscountFactorParamNames.firm={'firmbeta'};
if Params.scenario<4
    % Notice we use 'Electrify_FirmReturnFn'
    ReturnFn.firm=@( ...
            d,kprime,k,z, ...
            w, ...
            ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg ...
        ) Electrify_FirmReturnFn( ...
            d,kprime,k,z, ...
            w, ...
            ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg ...
        );
else
    % Notice we use 'Electrify_4FirmReturnFn'
    ReturnFn.firm=@( ...
            kprime,pvprime,k,pv,z, ...
            w, ...
            ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg,Ek,ek,pv_max_firm,carbon_tax ...
        ) Electrify_4FirmReturnFn( ...
            0,kprime,pvprime,k,pv,z, ...
            w, ...
            ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg,Ek,ek,pv_max_firm,carbon_tax ...
        );
end

% For energy
DiscountFactorParamNames.energy={};
if Params.scenario<4
    % Notice we use 'Electrify_EnergyReturnFn'
    ReturnFn.energy=@( ...
            aprime,a,z ...
        ) Electrify_EnergyReturnFn( ...
            aprime,a,z ...
        );
else
    % Notice we use 'Electrify_EnergyReturnFn'
    ReturnFn.energy=@( ...
            d,aprime,a,z ...
        ) Electrify_4EnergyReturnFn( ...
            d,aprime,a,z ...
        );
end

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
Params.r_pp=(1+r)^Params.ypp-1;
Params.pv_max_firm=pv_grid_firm(end);
Params.firmbeta=1/(1+Params.r_pp/(1-Params.tau_cg)); % 1/(1+r_pp) but returns net of capital gains tax

% Solved by GE

% Some initial values/guesses for variables that will be determined in general eqm
Params.P0=1;
Params.w=1;
Params.pension=0.4; % Initial guess (this will be determined in general eqm)
Params.G_pp=0.1*Params.ypp; % Government expenditure

% And some initial values/guesses for AggVar values that will be calculated while calculating the general eqm
Params.D_pp=(1+0.10)^Params.ypp-1; % The dividends paid by the firm per period
Params.EnergyCosts_h=0.3; % Energy used by households
Params.EnergyCosts_f=0.7; % Energy used by firms
Params.CarbonCosts_h=0.05; % Carbon tax paid by households
Params.CarbonCosts_f=0.3; % Cabron tax by firms

%% General eqm variables
if Params.scenario<3
    GEPriceParamNames={'w','P0','pension','G_pp','AccidentBeqS_pp'};
elseif Params.scenario<4
    GEPriceParamNames={'w','P0','pension','G_pp','AccidentBeqS_pp','AccidentBeqAH_pp'};
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

% Stationary Distribution Aggregates from households (important that ordering of Names and Functions is the same)
if Params.scenario<3
    FnsToEvaluate.L_h.household=@(labor,sprime,s,z,e,kappa_j,Lhscale) ...
        labor*kappa_j*exp(z+e)*Lhscale;  % Aggregate labour supply in efficiency units, not scaled by ypp
    FnsToEvaluate.S.household=@(labor,sprime,s,z,e) s; % Aggregate share holdings
    FnsToEvaluate.PensionSpending.household=@(labor,sprime,s,z,e,pension,ypp,agej,Jr) ...
        (agej>=Jr)*pension*ypp; % Total spending on pensions
    FnsToEvaluate.PayrollTaxRevenue.household=@(labor,sprime,s,z,e,ypp,agej,Jr,tau_l,w,kappa_j,Lhscale) ...
        (agej<Jr)*tau_l*labor*w*kappa_j*exp(z+e)*ypp*Lhscale; % Total spending on payroll taxes
    FnsToEvaluate.CapitalGainsTaxRevenue.household=@(labor,sprime,s,z,e,tau_cg,P0,D_pp,tau_d,r_pp) ...
        tau_cg*(P0-(((1-tau_cg)*P0 + (1-tau_d)*D_pp)/(1+r_pp-tau_cg)))*s; % tau_cg*(P0-Plag)*s, but substitute P=Plag, and then substitute for P
    FnsToEvaluate.BeqleftS_pp.household=@(labor,sprime,s,z,e,sj) ...
        sprime*(1-sj); % Accidental share bequests left by people who die
elseif Params.scenario<4
    FnsToEvaluate.L_h.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,kappa_j,Lhscale) ...
        labor*kappa_j*exp(z+e)*Lhscale;  % Aggregate labour supply in efficiency units, not scaled by ypp
    FnsToEvaluate.S.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) s; % Aggregate share holdings
    FnsToEvaluate.PensionSpending.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,pension,ypp,agej,Jr) ...
        (agej>=Jr)*pension*ypp; % Total spending on pensions
    FnsToEvaluate.PayrollTaxRevenue.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,ypp,agej,Jr,tau_l,w,kappa_j,Lhscale) ...
        (agej<Jr)*tau_l*labor*w*kappa_j*exp(z+e)*Lhscale*ypp; % Total spending on payroll taxes
    FnsToEvaluate.CapitalGainsTaxRevenue.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,tau_cg,P0,D_pp,tau_d,r_pp) ...
        tau_cg*(P0-(((1-tau_cg)*P0 + (1-tau_d)*D_pp)/(1+r_pp-tau_cg)))*s+(1-tau_d)*r_pp*max(a,0); % tau_cg*(P0-Plag)*s + deposit interest, but substitute P=Plag, and then substitute for P
    FnsToEvaluate.BeqleftS_pp.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,sj) ...
        sprime*(1-sj); % Accidental share bequests left by people who die
    % AccidentalBeqAHLeft is zero (if in debt) or accidental asset+house bequests left by people who die
    FnsToEvaluate.BeqleftAH_pp.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,scenario,sj,cpi) ...
        max(0,(aprime+(1+cpi)*hprime)*(1-sj));
    % BadDebt is the debt somebody accidentally leaves behind, or zero if net worth is positive
else
    FnsToEvaluate.L_h.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,kappa_j,Lhscale) ...
        labor*kappa_j*exp(z+e)*Lhscale;  % Aggregate labour supply in efficiency units, not scaled by ypp
    FnsToEvaluate.S.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e) s; % Aggregate share holdings
    FnsToEvaluate.PensionSpending.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,pension,ypp,agej,Jr) ...
        (agej>=Jr)*pension*ypp; % Total spending on pensions
    FnsToEvaluate.PayrollTaxRevenue.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,ypp,agej,Jr,tau_l,w,kappa_j,Lhscale) ...
        (agej<Jr)*tau_l*labor*w*kappa_j*exp(z+e)*Lhscale*ypp; % Total spending on payroll taxes
    FnsToEvaluate.CapitalGainsTaxRevenue.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,tau_cg,P0,D_pp,tau_d,r_pp) ...
        tau_cg*(P0-(((1-tau_cg)*P0 + (1-tau_d)*D_pp)/(1+r_pp-tau_cg)))*s+(1-tau_d)*r_pp*max(a,0); % tau_cg*(P0-Plag)*s + deposit interest, but substitute P=Plag, and then substitute for P
    FnsToEvaluate.BeqleftS_pp.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,sj) ...
        sprime*(1-sj); % Accidental share bequests left by people who die
    % AccidentalBeqAHLeft is zero (if in debt) or accidental asset+house bequests left by people who die
    FnsToEvaluate.BeqleftAH_pp.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,scenario,sj,cpi) ...
        max(0,(aprime+(1+cpi)*hprime)*(1-sj));
    % BadDebt is the debt somebody accidentally leaves behind, or zero if net worth is positive
    FnsToEvaluate.EnergyCosts_h.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,w,ypp,cpi_energy,energy_pct_cost,energy_pct_brown,carbon_tax) ...
        Electrify_4HouseholdEnergyCosts(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,w,ypp,cpi_energy,energy_pct_cost,energy_pct_brown,carbon_tax);
    FnsToEvaluate.CarbonCosts_h.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,w,ypp,cpi_energy,energy_pct_cost,energy_pct_brown,carbon_tax) ...
        Electrify_4HouseholdCarbonCosts(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,w,ypp,cpi_energy,energy_pct_cost,energy_pct_brown,carbon_tax);
end

% From firms
if Params.scenario<4
    FnsToEvaluate.L_f.firm=@(d,kprime,k,z,w,alpha_k,alpha_l) ...
        (w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % (effective units of) labor demanded by firm, not scaled by ypp
    FnsToEvaluate.K.firm=@(d,kprime,k,z,w,alpha_k,alpha_l) k; % physical capital
    FnsToEvaluate.D_pp.firm=@(d,kprime,k,z,ypp) (1+d)^ypp-1; % dividend paid by firm
    FnsToEvaluate.Sissued.firm=@(d,kprime,k,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi) ...
        Electrify_FirmShareIssuance(d,kprime,k,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi); % Share issuance
    FnsToEvaluate.CorpTaxRevenue.firm=@(d,kprime,k,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi) ...
        Electrify_FirmCorporateTaxRevenue(d,kprime,k,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi); % revenue from the corporate profits tax
else
    FnsToEvaluate.L_f.firm=@(kprime,pvprime,k,pv,z,w,alpha_k,alpha_l) ...
        (w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % (effective units of) labor demanded by firm, not scaled by ypp
    FnsToEvaluate.K.firm=@(kprime,pvprime,k,pv,z,w,alpha_k,alpha_l) k; % physical capital
    FnsToEvaluate.PV_f.firm=@(kprime,pvprime,k,pv,z,w,alpha_k,alpha_l) pv; % firm's solarPV generation capacity
    FnsToEvaluate.D_pp.firm=@(kprime,pvprime,k,pv,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,Ek,ek,pv_max_firm,carbon_tax) ...
        Electrify_4FirmDividend(0,kprime,pvprime,k,pv,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,Ek,ek,pv_max_firm,carbon_tax); % dividend paid by firm
    FnsToEvaluate.Sissued.firm=@(kprime,pvprime,k,pv,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,Ek,ek,pv_max_firm,carbon_tax) ...
        Electrify_4FirmShareIssuance(0,kprime,pvprime,k,pv,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,Ek,ek,pv_max_firm,carbon_tax); % Share issuance
    FnsToEvaluate.CorpTaxRevenue.firm=@(kprime,pvprime,k,pv,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,Ek,ek,pv_max_firm,carbon_tax) ...
        Electrify_4FirmCorporateTaxRevenue(0,kprime,pvprime,k,pv,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,Ek,ek,pv_max_firm,carbon_tax); % revenue from the corporate profits tax
    FnsToEvaluate.EnergyCosts_f.firm=@(kprime,pvprime,k,pv,z,w,ypp,alpha_k,alpha_l,Ek,ek,pv_max_firm,carbon_tax) ...
        Electrify_4FirmEnergyCosts(0,kprime,pvprime,k,pv,z,w,ypp,alpha_k,alpha_l,Ek,ek,pv_max_firm,carbon_tax);
    FnsToEvaluate.CarbonCosts_f.firm=@(kprime,pvprime,k,pv,z,w,ypp,alpha_k,alpha_l,Ek,ek,pv_max_firm,carbon_tax) ...
        Electrify_4FirmCarbonCosts(0,kprime,pvprime,k,pv,z,w,ypp,alpha_k,alpha_l,Ek,ek,pv_max_firm,carbon_tax);
end

% From energy -- there must be at least one
if Params.scenario<4
    FnsToEvaluate.EnergyRevenue.energy=@(aprime,a,z) 0;
else
    FnsToEvaluate.EnergyRevenue.energy=@(invest,aprime,a,z,EnergyCosts_h,EnergyCosts_f) EnergyCosts_h+EnergyCosts_f;
    FnsToEvaluate.TransitionInvestment.energy=@(invest,aprime,a,z,CarbonCosts_h,CarbonCosts_f) CarbonCosts_h+CarbonCosts_f;
end

% General Equilibrium conditions (these should evaluate to zero in general equilbrium)
GeneralEqmEqns.sharemarket=@(S) S-1; % mass of all shares equals one
GeneralEqmEqns.labormarket=@(L_h,L_f) (Params.scenario+1)*(L_h-L_f)*Params.ypp; % labor supply of households equals labor demand of firms (scaled by ypp)
GeneralEqmEqns.pensions=@(PensionSpending,PayrollTaxRevenue) PensionSpending-PayrollTaxRevenue; % Retirement benefits equal Payroll tax revenue: pension*fractionretired-tau*w*H
GeneralEqmEqns.govbudget=@(G_pp,tau_d,D_pp,CapitalGainsTaxRevenue,CorpTaxRevenue) G_pp-tau_d*D_pp-CapitalGainsTaxRevenue-CorpTaxRevenue; % G is equal to the target, GdivYtarget*Y
% GeneralEqmEqns.firmdiscounting=@(firmbeta,r_pp,tau_cg) firmbeta-1/(1+r_pp/(1-tau_cg)); % Firms discount rate is related to market return rate
% GeneralEqmEqns.dividends=@(D_pp,D_pp) (Params.scenario+1)*(D_pp-D_pp); % That the dividend households receive equals that which firms give
GeneralEqmEqns.ShareIssuance=@(Sissued,P0,D_pp,tau_cg,tau_d,r_pp) ...
    P0-((((1-tau_cg)*P0 + (1-tau_d)*D_pp)/(1+r_pp-tau_cg))-Sissued); % P0=P-S, but substitute for P (see derivation inside the return fn)
GeneralEqmEqns.CapitalOutputRatio=@(K,L_f,TargetKdivL) (K/L_f-TargetKdivL)/100; % Ratio not based on ypp
GeneralEqmEqns.bequestsS_pp=@(BeqleftS_pp,AccidentBeqS_pp,n_pp) BeqleftS_pp/(1+n_pp)-AccidentBeqS_pp; % Accidental share bequests received equal accidental share bequests left
if Params.scenario>2
    GeneralEqmEqns.bequestsAH_pp=@(BeqleftAH_pp,AccidentBeqAH_pp,n_pp) BeqleftAH_pp/(1+n_pp)-AccidentBeqAH_pp; % Accidental asset+house bequests received equal accidental asset+house bequests left
end

% For analysing the model
FnsToEvaluate2=FnsToEvaluate;
if Params.scenario<3
    FnsToEvaluate2.earnings.household=@(labor,aprime,a,z,e,w,kappa_j,Lhscale) w*kappa_j*labor*exp(z+e)*Lhscale; % w*kappa_j is the labor earnings
    FnsToEvaluate2.A.household=@(labor,aprime,a,z,e,w,kappa_j,Lhscale) a; % w*kappa_j is the labor earnings
    FnsToEvaluate2.BeqleftS_pp.household=@(labor,aprime,a,z,e,sj) aprime*(1-sj); % Accidental asset bequests left by people who die
elseif Params.scenario<4
    FnsToEvaluate2.earnings.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,w,kappa_j,Lhscale) w*kappa_j*labor*exp(z+e)*Lhscale; % w*kappa_j is the labor earnings
    FnsToEvaluate2.A.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) a; % Aggregate asset/mortgage holdings
    FnsToEvaluate2.S.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) s; % Aggregate share holdings
    FnsToEvaluate2.H.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) h; % Aggregate house holdings
    FnsToEvaluate2.PV.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) solarpv; % Aggregate solarpv holdings
    FnsToEvaluate2.BeqleftS_pp.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,sj) sprime*(1-sj); % Accidental share bequests left by people who die
    FnsToEvaluate2.BeqleftAH_pp.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,scenario,sj,cpi) max(0,(aprime+(1+cpi)*hprime)*(1-sj));
    FnsToEvaluate2.BadDebt_pp.household=@(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e,scenario,sj,cpi) ...
        min(0,(aprime+(1+cpi)*hprime)*(1-sj));
else
    FnsToEvaluate2.earnings.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,w,kappa_j,Lhscale) w*kappa_j*labor*exp(z+e)*Lhscale; % w*kappa_j is the labor earnings
    FnsToEvaluate2.A.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e) a; % Aggregate asset/mortgage holdings
    FnsToEvaluate2.S.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e) s; % Aggregate share holdings
    FnsToEvaluate2.Car.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e) car; % Aggregate house holdings
    FnsToEvaluate2.H.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e) h; % Aggregate house holdings
    FnsToEvaluate2.PV_h.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e) solarpv; % Aggregate solarpv holdings
    FnsToEvaluate2.BeqleftS_pp.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,sj) sprime*(1-sj); % Accidental share bequests left by people who die
    FnsToEvaluate2.BeqleftAH_pp.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,scenario,sj,cpi) max(0,(aprime+(1+cpi)*hprime)*(1-sj));
    FnsToEvaluate2.BadDebt_pp.household=@(labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e,scenario,sj,cpi) ...
        min(0,(aprime+(1+cpi)*hprime)*(1-sj));
end
if Params.scenario<4
    FnsToEvaluate2.Output.firm=@(d,kprime,k,z,w,ypp,alpha_k,alpha_l) ...
        z*(k^alpha_k)*((w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)))^alpha_l*ypp; % Production function z*(k^alpha_k)*(l^alpha_l) (substituting for l)
else
    FnsToEvaluate2.Output.firm=@(kprime,pvprime,k,pv,z,w,ypp,alpha_k,alpha_l,Ek,ek) ...
        (Ek*ek)*z*(k^alpha_k)*((w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)))^alpha_l*ypp; % Production function z*(k^alpha_k)*(l^alpha_l) (substituting for l)
    FnsToEvaluate2.CarbonCosts_f.firm=@(kprime,pvprime,k,pv,z,w,ypp,alpha_k,alpha_l,Ek,ek,pv_max_firm,carbon_tax) ...
        Electrify_4FirmCarbonCosts(0,kprime,pvprime,k,pv,z,w,ypp,alpha_k,alpha_l,Ek,ek,pv_max_firm,carbon_tax);
end

% Note: I keep the FnsToEvaluate use in general eqm to a minimum (to reduce
% runtimes) and then use FnsToEvaluate2 to analyse model with more stats.
% Note: FnsToEvaluate may need 'e' grids, but AggVars and other stats use a
% a joint ze grid (which reads as z in their parameter lists).

if true
    [V_f, Policy_f]=ValueFnIter_Case1_PType(n_d,n_a,n_z, {'firm'}, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
    % We can plot V as a 3d plot (surf is matlab command for 3d plot)
    figure(1)
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
        surf(a_grid.household*ones(1,Params.J),ones(n_a.household,1)*(Params.agejshifter+(1:1:Params.J)),reshape(V_init.household(:,zind,:),[n_a.household,Params.J]))
    else
        surf(a_grid.household*ones(1,Params.J),ones(n_a.household,1)*(Params.agejshifter+(1:1:Params.J)),reshape(V_init.household(:,zind,eind,:),[n_a.household,Params.J]))
    end
    title('Value function: median value of z')
    xlabel('Assets (a)')
    ylabel('Age in Years')
end

if Params.scenario<3
    checkfeasible_household_case12(V_init,Policy_init,a_grid.household,1,n_d,n_a,n_z,N_j,d_grid,a_grid,Params,vfoptions);
    checkfeasible_firm_case12(V_init,Policy_init,n_d,n_a,n_z,d_grid,a_grid,Params,vfoptions);
elseif Params.scenario<4
    checkfeasible_household_case3(V_init,Policy_init,buyhouse_grid,asset_grid,zeroassetindex,house_grid,n_d,n_a,n_z,N_j,d_grid,a_grid,Params,vfoptions);
end

%% Initial distribution of agents at birth (j=1)
% Before we plot the life-cycle profiles we have to define how agents are
% at age j=1. We will give them all zero shares (and possibly zero assets, no house, no solarpv).
if small_z_no_e
    jequaloneDist.household=zeros([n_a.household,n_z.household],'gpuArray'); % Put no households anywhere on grid
    if Params.scenario<3
        % All agents start with zero shares, and the median shocks
        jequaloneDist.household(1,floor((n_z.household+1)/2))=1;
    elseif Params.scenario<4
        % All agents start with zero shares, assets, houses, solarpv, and median shocks
        jequaloneDist.household(1,zeroassetindex,1,1,floor((n_z.household+1)/2))=1;
    else
        % All agents start with zero shares, assets, cars, houses, solarpv, and median shocks
        jequaloneDist.household(1,zeroassetindex,1,1,1,floor((n_z.household+1)/2))=1;
    end
else
    jequaloneDist.household=zeros([n_a.household,n_z.household,vfoptions.n_e.household],'gpuArray'); % Put no households anywhere on grid
    if Params.scenario<3
        % All agents start with zero shares, and the median shocks
        jequaloneDist.household(1,floor((n_z.household+1)/2),floor((simoptions.n_e.household+1)/2))=1;
    elseif Params.scenario<4
        % All agents start with zero shares, assets, houses, solarpv, and median shocks
        jequaloneDist.household(1,zeroassetindex,1,1,floor((n_z.household+1)/2),floor((simoptions.n_e.household+1)/2))=1;
    else
        % All agents start with zero shares, assets, cars, houses, solarpv, and median shocks
        jequaloneDist.household(1,zeroassetindex,1,1,1,floor((n_z.household+1)/2),floor((simoptions.n_e.household+1)/2))=1;
    end
end

% Note that because the firms are infinite horizon they do not have an age=1 distribution

%% Agents age distribution
AgeWeightsParamNames=struct('household',{{'mewj'}}); % So VFI Toolkit knows which parameter is the mass of agents of each age

%% Test
disp('Test StationaryDist')
StationaryDist_init=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy_init,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);

%% Test
% Note: Because we used simoptions we must include this as an input
disp('Test AggVars')
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_init, Policy_init, FnsToEvaluate2, Params, n_d, n_a, n_z,N_j,Names_i, d_grid, a_grid, z_grid,simoptions);

% Next few lines were used to try a few parameter values so as to get a
% decent initial guess before actually solving the general equilbrium
fprintf('Check: L_h, L_f, K \n')
[AggVars.L_h.Mean,AggVars.L_f.Mean,AggVars.K.Mean]
fprintf('Check: K/L_f (should be about 2.03) \n')
AggVars.K.Mean/AggVars.L_f.Mean
if Params.scenario<3
    fprintf('Check: S, D_pp \n')
    [AggVars.S.Mean,AggVars.D_pp.Mean]
elseif Params.scenario<4
    fprintf('Check: S, A, H, PV_h\n')
    [AggVars.S.Mean,AggVars.A.Mean,AggVars.H.Mean,AggVars.PV_h.Mean]
else
    fprintf('Check: S, A, H, PV_h, PV_f \n')
    [AggVars.S.Mean,AggVars.A.Mean,AggVars.H.Mean,AggVars.PV_h.Mean,AggVars.PV_f.Mean]
end
fprintf('Check: ShareIssuance GE condition \n')
Params.P0-((((1-Params.tau_cg)*Params.P0 + (1-Params.tau_d)*Params.D_pp)/(1+Params.r_pp-Params.tau_cg))-AggVars.S.Mean)

end % solve_setup

%% Solve for the General Equilibrium
if solve_GE_init
    % heteroagentoptions.fminalgo=4 % CMA-ES algorithm 
    
    heteroagentoptions.verbose=1;
    if Params.scenario<3
        heteroagentoptions.toleranceGEprices=10^(-4);
        heteroagentoptions.toleranceGEcondns=10^(-4); % This is the hard one
        if solve_TPath
            % heteroagentoptions.maxiter=200;
        end
    else
        heteroagentoptions.toleranceGEprices=10^(-2);
        heteroagentoptions.toleranceGEcondns=10^(-1); % This is the hard one
        heteroagentoptions.maxiter=50;                % About 3 hours for 35 iterations
    end
    if Params.scenario>3
        % heteroagentoptions.useCustomModelStats=1;
        heteroagentoptions.household.CustomModelStats=@( ...
            V,Policy,StationaryDist,Parameters,FnsToEvaluate, ...
            n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,heteroagentoptions,vfoptions,simoptions ...
            ) Electrify_4HouseholdCustomModelStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate, ...
            n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,heteroagentoptions,vfoptions,simoptions);
    elseif Params.scenario>2
        % heteroagentoptions.useCustomModelStats=1;
        heteroagentoptions.household.CustomModelStats=@( ...
            V,Policy,StationaryDist,Parameters,FnsToEvaluate, ...
            n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,heteroagentoptions,vfoptions,simoptions ...
            ) Electrify_HouseholdCustomModelStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate, ...
            n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,heteroagentoptions,vfoptions,simoptions);
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
    if solve_TPath
        clear solve_TPath
        save tpathElectrifyA.mat
        solve_TPath=true;
    else
        clear solve_TPath
        save tpathElectrifyA.mat
        solve_TPath=false;
    end
    % load tpathElectrifyA.mat
else
    if solve_GE_final
        load tpathElectrifyA.mat
        solve_GE_final=true;
    else
        load tpathElectrifyA.mat
        solve_GE_final=false;
    end
end

if solve_GE_final
    % 40 years of changing demographics
    % 60 years in final demographic state (to allow time to converge to final stationary general eqm)
    % Conditional survival probabilities
    ParamPath.sj=[sj_init+(sj_final-sj_init).*linspace(0,1,ceil(40/(Params.ypp*jpT)))'; sj_final.*ones(T-ceil(40/(Params.ypp*jpT)),1)];
    % T-by-N_j (whether this or N_j-by_T, toolkit understands both)
    % Calculate the implied mewj from the sj
    ParamPath.mewj=cumprod([ones(T,1), ParamPath.sj(:,1:end-1)], 2); % mass of age jj is the mass of jj-1 that survive
    % Factor in population growth; In N_j dimension, older people are from earlier (smaller) populations
    % ...in the T dimension, we see overall population growth as T increases
    ParamPath.mewj=ParamPath.mewj./((1+Params.n_pp).^(Params.ypp*((1:Params.J)-1))); % Population shrinks in the N_j dimension
    ParamPath.mewj=ParamPath.mewj.*((1+Params.n_pp).^(Params.ypp*jpT*((1:T)-1)))'; % Population grows in the T dimension
    ParamPath.mewj=ParamPath.mewj./sum(ParamPath.mewj,2); % normalize age-masses to sum to one
    % Looking at ParamPath.mewj you can see that as tt increases, the mass at older ages increases

    %% Solve for final stationary general eqm with Params at time T
    Params.Ek=ParamPath.Ek(T);
    Params.ek=ParamPath.ek(T);
    Params.carbon_tax=ParamPath.carbon_tax(T);
    Params.energy_pct_brown=ParamPath.energy_pct_brown(T);
    Params.sj=ParamPath.sj(T,:); % conditional survival probabilities
    Params.mewj=ParamPath.mewj(T,:);
    Params.cpi=ParamPath.cpi(T);
    Params.cpi_energy=ParamPath.cpi_energy(T);
    Params.Lhscale=ParamPth.Lhscale(T);

if true
    ParamPath.Ek=linspace(1,1.2,T); Params.Ek=ParamPath.Ek(1);
    ParamPath.ek=linspace(1,1.2,T); Params.ek=ParamPath.ek(1);
    Tx=T;
    Params.ek=ParamPath.ek(Tx);
    Params.carbon_tax=ParamPath.carbon_tax(Tx);
    Params.energy_pct_brown=ParamPath.energy_pct_brown(Tx);
    Params.sj=ParamPath.sj(Tx,:); % conditional survival probabilities
    Params.mewj=ParamPath.mewj(Tx,:);
    Params.cpi=ParamPath.cpi(Tx);
    Params.cpi_energy=ParamPath.cpi_energy(Tx);
    [V_f, Policy_f]=ValueFnIter_Case1_PType(n_d,n_a,n_z, {'firm'}, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
    % We can plot V as a 3d plot (surf is matlab command for 3d plot)
    figure(1)
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
        fprintf('Check: S, D_pp \n')
        [AggVars.S.Mean,AggVars.D_pp.Mean]
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

    save tpathElectrifyB.mat
else
    if solve_TPath
        load tpathElectrifyB.mat
        solve_TPath=true;
    else
        load tpathElectrifyB.mat
        solve_TPath=false;
    end
end % solve_GE_final

    if ~solve_demographic_change
        % TESTING!  This resets our GEqm to initial rather than final state
        p_eqm_final=p_eqm_init;
        Params.pension=p_eqm_final.pension;
        Params.AccidentBeqS_pp=p_eqm_final.AccidentBeqS_pp;
        if Params.scenario>2
            Params.AccidentBeqAH_pp=p_eqm_final.AccidentBeqAH_pp;
        end
        Params.G_pp=p_eqm_final.G_pp;
        Params.w=p_eqm_final.w;
        % Params.firmbeta=p_eqm_final.firmbeta;
        Params.P0=p_eqm_final.P0;
    
        % Evaluate the final stationary general eqm
        [V_final, Policy_final]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
        StationaryDist_final=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy_final,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
        % Calculate various stats
        AllStats_final=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_final, Policy_final, FnsToEvaluate2, Params, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid, z_grid,simoptions);
        % Calculate the life-cycle profiles
        AgeConditionalStats_final=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_final,Policy_final, FnsToEvaluate2,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
    end

if solve_TPath
    %% Setup for the transition path
    % T=100; % number of periods for transition path
    
    % Already created ParamPath.sj and ParamPath.mewj
    
    % Initial guess for general eqm parameters
    T_eq=ceil(T/2); % Demographic change has stopped and T_eq begins period of transition equilibrium-finding
    PricePath0.w=[linspace(p_eqm_init.w, p_eqm_final.w,T_eq), p_eqm_final.w*ones(1,T-T_eq)];
    % PricePath0.firmbeta=[linspace(p_eqm_init.firmbeta, p_eqm_final.firmbeta,T_eq), p_eqm_final.firmbeta*ones(1,T-T_eq)];
    PricePath0.P0=[linspace(p_eqm_init.P0, p_eqm_final.P0,T_eq), p_eqm_final.P0*ones(1,T-T_eq)];
    PricePath0.pension=[linspace(p_eqm_init.pension, p_eqm_final.pension,T_eq), p_eqm_final.pension*ones(1,T-T_eq)];
    PricePath0.AccidentBeqS_pp=[linspace(p_eqm_init.AccidentBeqS_pp,p_eqm_final.AccidentBeqS_pp,T_eq), p_eqm_final.AccidentBeqS_pp*ones(1,T-T_eq)];
    if Params.scenario>2
        PricePath0.AccidentBeqAH_pp=[linspace(p_eqm_init.AccidentBeqAH_pp,p_eqm_final.AccidentBeqAH_pp,T_eq), p_eqm_final.AccidentBeqAH_pp*ones(1,T-T_eq)];
    end
    PricePath0.G_pp=[linspace(p_eqm_init.G_pp, p_eqm_final.G_pp,T_eq), p_eqm_final.G_pp*ones(1,T-T_eq)];
    % PricePath0.TargetKdivL=2.03*ones(1,T);
    % Just some reasonable guesses I made up.
    
    % General eqm eqns, same idea as with the stationary general eqm
    % GeneralEqmEqns_Transition.capitalmarket=@(r_pp,alpha_k,alpha_l,delta,K,L,ypp) r_pp-(alpha_k*(K^(alpha_k-1))*(L^(alpha_l))-((delta+1)^ypp-1)); % r=marginal product of capital
    GeneralEqmEqns_Transition.labormarket=@(w,alpha_k,alpha_l,K,L_f) w-(alpha_l)*(K^alpha_k)*(L_f^(alpha_l-1)); % w=marginal product of labor
    GeneralEqmEqns_Transition.firmdiscounting=GeneralEqmEqns.firmdiscounting;
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
            transpathoptions.GEnewprice3.howtoupdate{pp,4}=0.01;
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
    [PricePath,GECondnsPath]=TransitionPath_Case1_FHorz_PType(PricePath0, ParamPath, T, V_final, AgentDist_init, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate2, GeneralEqmEqns_Transition, Params, DiscountFactorParamNames, AgeWeightsParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);
    
    %%
    save tpathElectrifyD.mat
    % load tpathElectrifyD.mat

    %% Now calculate some things about the transition path (path for Value fn, Policy fn, Agent Distribution)
    % You can calculate the value and policy functions for the transition path
    [VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
    
    % You can then use these to calculate the agent distribution for the transition path
    AgentDistPath=AgentDistOnTransPath_Case1_FHorz_PType(StationaryDist_init, jequaloneDist, PricePath, ParamPath, PolicyPath, AgeWeightsParamNames,n_d,n_a,n_z,N_j,Names_i,pi_z,T, Params, transpathoptions, simoptions);
    
    %% Analyse the transition path
    % And then we can calculate AggVars for the path
    AggVarsPath=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath,PolicyPath, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptions, simoptions);
    
    %% Plot some paths
    figure(1)
    % Plot of K and w
    % Note: include periods -3 to 0 (the initial stationary eqm) so can see any jump in period 1
    subplot(2,1,1); plot(1:1:T,AggVarsPath.K.Mean)
    hold on
    plot(-3:1:0,AllStats_init.K.Mean*ones(1,4),'w')
    hold off
    xlim([-3,T])
    title('Path of aggregate capital (K)')
    subplot(2,1,2); plot(1:1:T,PricePath.w)
    hold on
    plot(-3:1:0,p_eqm_init.w*ones(1,4),'w')
    hold off
    xlim([-3,T])
    title('Path of wage rate (w)')


end % solve_TPath

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

function feasible=checkfeasible_household_case12(V,Policy,asset_grid,zeroassetindex,n_d,n_a,n_z,N_j,d_grid,a_grid,Params,vfoptions)

return
feasible=true;
fhh_options.tolerance=vfoptions.tolerance;
fhh_options.lowmemory=vfoptions.lowmemory.household;
z_idx=ceil(n_z.household/2);
if isfield(vfoptions,'n_e')
    fhh_options.n_e=vfoptions.n_e.household;
    fhh_options.e_grid=vfoptions.e_grid.household;
    fhh_options.pi_e=vfoptions.pi_e.household;
    e_idx=ceil(fhh_options.n_e/2);
    fhh=PolicyInd2Val_FHorz(Policy.household,n_d.household,n_a.household,n_z.household,N_j.household,d_grid.household,a_grid.household,fhh_options);
    fhh=squeeze(fhh(:,:,:,:,z_idx,e_idx,:));
else
    fhh=PolicyInd2Val_FHorz(Policy.household,n_d.household,n_a.household,n_z.household,N_j.household,d_grid.household,a_grid.household,fhh_options);
    fhh=squeeze(fhh(:,:,:,:,z_idx,:));
end

z_idx=ceil(n_z.household/2);
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
        Params.r_pp,0,0,Params.energy_pct_cost);
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

end

function feasible=checkfeasible_firm_case12(V,Policy,n_d,n_a,n_z,d_grid,a_grid,Params,vfoptions)
return
feasible=true;
fxx_options.tolerance=vfoptions.tolerance;
fxx_options.lowmemory=vfoptions.lowmemory.firm;
z_idx=ceil(n_z.firm/2);
fxx=PolicyInd2Val_Case1(Policy.firm,n_d.firm,n_a.firm,n_z.firm,d_grid.firm,a_grid.firm,fxx_options);
fxx=squeeze(fxx(:,:,z_idx));

z_idx=ceil(n_z.firm/2);
kprime_index_last=1;
[~,d_index]=min(abs(d_grid.firm-0.2));
[k_value,k_index]=min(abs(a_grid.firm-0.35));
for fxx_agej=1:5
    fxx_next=squeeze(fxx(:,d_index));
    kprime_index_next=find(k_value==fxx_next(2,kprime_index_last));
    if fxx_next(1,kprime_index_next)==0
        feasible=false;
        error("capital strike")
    end

    income = Electrify_FirmCorporateTaxRevenue( ...
        1,0,asset_grid(aprime_index_next),0,0,asset_grid(aprime_index_last),0,0,0,0,0, ...
        Params.pension,Params.AccidentBeqS_pp,0,Params.w,Params.P0,Params.D_pp, ...
        Params.kappa_j(fxx_agej),Params.tau_l,Params.tau_d,Params.tau_cg,Params.ypp,fxx_agej,Params.Jr, ...
        Params.r_pp,0,0,Params.energy_pct_cost);
    if income<0
        feasible=false;
        error("income negative")
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

function feasible=checkfeasible_household_case3(V,Policy,buyhouse_grid,asset_grid,zeroassetindex,house_grid,n_d,n_a,n_z,N_j,d_grid,a_grid,Params,vfoptions)

feasible=true;
fhh_options.tolerance=vfoptions.tolerance;
fhh_options.lowmemory=vfoptions.lowmemory.household;
fhh_options.experienceasset=vfoptions.experienceasset.household;
fhh_options.aprimeFn=vfoptions.aprimeFn.household;
fhh_options.refine_d=vfoptions.refine_d.household;
z_idx=ceil(n_z.household/2);
if isfield(vfoptions,'n_e')
    fhh_options.n_e=vfoptions.n_e.household;
    fhh_options.e_grid=vfoptions.e_grid.household;
    fhh_options.pi_e=vfoptions.pi_e.household;
    e_idx=ceil(fhh_options.n_e/2);
    fhh=PolicyInd2Val_FHorz(Policy.household,n_d.household,n_a.household,n_z.household,N_j.household,d_grid.household,a_grid.household,fhh_options);
    fhh=squeeze(fhh(:,:,:,:,:,z_idx,e_idx,:));
else
    fhh=PolicyInd2Val_FHorz(Policy.household,n_d.household,n_a.household,n_z.household,N_j.household,d_grid.household,a_grid.household,fhh_options);
    fhh=squeeze(fhh(:,:,:,:,:,z_idx,:));
end

buyhouse_index_last=1;
aprime_index_last=zeroassetindex;
hprime_index_last=1;
for fhh_agej=1:5
    % We have already removed z_idx and possibly e_idx above
    fhha_next=squeeze(fhh(:,1,:,1,1,fhh_agej));
    aprime_index_next=find(asset_grid==fhha_next(4,aprime_index_last));
    % hprime_index_next=find(house_grid==fhh_next(2,hprime_index_last));
    if fhha_next(1,aprime_index_next)==0
        feasible=false;
        warning("labor strike")
    end
    income = Electrify_HouseholdIncomeFn( ...
        1,0,0,asset_grid(aprime_index_next),0,0,asset_grid(aprime_index_last),0,0,0,0, ...
        Params.pension,Params.AccidentBeqS_pp,Params.AccidentBeqAH_pp,Params.w,Params.P0,Params.D_pp, ...
        Params.kappa_j(fhh_agej),Params.tau_l,Params.tau_d,Params.tau_cg,Params.ypp,fhh_agej,Params.Jr, ...
        Params.r_pp,Params.cpi_energy,Params.energy_pct_cost);
    if income<0
        feasible=false;
        error("income negative")
    end
    consumption=Electrify_HouseholdConsumptionFn( ...
        1,0,0,asset_grid(aprime_index_next),0,0,asset_grid(aprime_index_last),0,0,0,0, ...
        Params.pension,Params.AccidentBeqS_pp,Params.AccidentBeqAH_pp,Params.w,Params.P0,Params.D_pp, ...
        Params.kappa_j(fhh_agej),Params.tau_l,Params.tau_d,Params.tau_cg,Params.ypp,fhh_agej,Params.Jr, ...
        Params.r_pp,Params.r_wedge_pp,Params.f_htc,Params.rentprice,Params.cpi_energy,Params.pv_pct_cost,Params.energy_pct_cost);
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
            Params.scenario,Params.r_pp,Params.r_wedge_pp,Params.f_htc,Params.minhouse,Params.rentprice,Params.f_coll,Params.houseservices,Params.cpi_energy,Params.pv_pct_cost,Params.energy_pct_cost);
        if isfinite(F)
            fprintf("F = %.2f \n", F);
        else
            feasible=false;
            warning("infeasible last->next")
        end
    end
    aprime_index_last=aprime_index_next;
end

% all(squeeze(Policy.household(:,1,1:zeroassetindex-2,1,1,4:6,:,2:20))==1,[3 4 5])
% if all(squeeze(Policy.household(:,1,1:zeroassetindex-2,1,1,3:7,:,2:20))==1,'all')
%     warning("infeasible Stationary Distribution")
% end

end
