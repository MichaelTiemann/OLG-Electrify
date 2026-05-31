%% OLG Electrification (based on OLGModels14: Heterogeneous households and heterogeneous firms
%% and also Life-Cycle Model 35: Portfolio-Choice with Housing)
% See https://www.vfitoolkit.com/updates-blog/2021/an-introduction-to-life-cycle-models/
% OLGModel14.m in the repo https://github.com/vfitoolkit/IntroToOLGModels
% and LifeCycleModel35.m in the repo https://github.com/vfitoolkit/IntroToLifeCycleModels
% See https://github.com/MichaelTiemann/OLG-Electrify/blob/main/README.md for more info

% A line some need for running on the Server
addpath(genpath('./MatlabToolkits/'))

solve_setup=true;
solve_GE=3; % 0: skip GE; 1: solve initial, 2: solve final, 3: solve both
solve_TPath=true;
small_z_no_e=false; % n_z=1; n_e=0
small_model=true; % Minimal vs. maximal grid sizes
small_T=2; % small_T==1 means just do T=1, T=2 (or smallest not-to-be-confused-with-dimension); small_T==2 means use jpT

if solve_setup

Names_i={'firm','household','energy'};
PTypeDistParamNames={'ptypemass'};
Params.ptypemass=[0.5,0.4,0.1]; % Mass of households, firms, and energy sum to one

%% Parameters for household (4 scenarios)
% Scenario 1: no housing, no assets, no inflation
% Scenario 2: add rental+energy costs, but no housing/assets/inflation
% Scenario 3: add housing/assets/pv/inflation
% Scenario 4: add cars/detailed energy
Params.scenario=4;

% To be able to solve such a big problem, I switched to 5 year model period.
% Note that ypp (years-per-period) must be at most 15 (for kappa_j labor productivity evolution).
% Discounting parameters (beta and sj) defined in terms of ypp
Params.ypp=2; % model period, in years (just used this to modify some parameters from annual to model period)

% Lets model agents from age 20 to age 100, so 81 periods (or 61 for scenario 3)
max_age=80;
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

energy_pct_cost=[0,0.07,0.07,0.05]; % Electricity: 3%; Gas: 1-2%; Petrol: 1-2%; Scenario 4 dis-aggregates petrol from this cost

% Demographics
% Population growth rate
n=[0.02,0.02,0.01,0.01]; % percentage rate (expressed as fraction) of population growth per period

% Age-dependent labor productivity units
% Stage 1: starting out (typ. first 25-30 years)
% Stage 2: peak earnings (typ. years 25-30 (meaning ages 45-50))
% Stage 3: winding down (typ. last 14 years before retirement (ages 50-64))
% Stage r: retirement
% Labor productivity at start, peak, and end of working life
k_j1 = [0.5, 0.5, 0.5, 0.5];
k_j2 = [2, 2, 2, 2];
k_j2_length = [0,0,5,5]; % years...that will be scaled by YPP if/when needed
k_j3 = [1, 1, 1, 1];

% Note: These iid shocks will interact with the endogenous labor so the final labor
% earnings process will not equal that of Karahan & Ozkan (2013)
% Note: Karahan & Ozkan (2013) also have a fixed effect (which they call alpha) and which I ignore here.

% Warm glow of bequest
Params.warmglow1=0.3; % (relative) importance of bequests
Params.warmglow2=3; % bliss point of bequests (essentially, the target amount)
Params.warmglow3=Params.sigma; % By using the same curvature as the utility of consumption it makes it much easier to guess appropriate parameter values for the warm glow

AccidentBeqS=[0.02,0.02,0.02,0.02]; % Accidental bequests (this is the lump sum transfer of shares)
AccidentBeqAH=[0,0,0.02,0.02]; % Accidental bequests (this is the lump sum transfer of assets+house value)

% Preferences
% Relative importance of housing services (vs consumption) in utility
sigma_h=[0,0,0.5,0.5];
% Relative importance of car services (vs consumption) in utility
sigma_c=[0,0,0.5,0.3];
Params.eta=1.5; % Curvature of leisure (This will end up being 1/Frisch elasticity)
psi = [2, 1, 1, 1]; % Weight on leisure

%% Energy parameters
transition_target_firm=0.90; % What fraction of electrification transition do we target (zero to one)?
transition_target_energy=0.90; % What fraction of electrification transition do we target (zero to one)?
if Params.scenario==4
    Params.pvinstalled_firm=0; Params.pvinstalled_energy=0;
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

vfoptions=struct(); simoptions=struct();
Params=Electrify_Scenario_YPP_Setup(Params,Params.scenario,Params.ypp,small_z_no_e,max_age,agejshifter,r,r_wedge,beta,n,k_j1,k_j2,k_j2_length,k_j3,sigma_h,sigma_c,psi,Params.tau_cg,energy_pct_cost,G,D,AccidentBeqS,AccidentBeqAH);
[ReturnFn,FnsToEvaluate,FnsToEvaluate2,FnsToEvaluate3,vfoptions,simoptions]=Electrify_Scenario_Fn_Setup(Params,vfoptions,simoptions);

% Housing (ignored/overwritten if no housing in scenario)
% Params.minhouse % set below, is the minimum value of house that can be purchased
Params.rentprice=0.3; % To make real fraction of income, must be multiplied by kappa_j in scenarios 3 & 4
% housing services as a fraction of house value (ignored in scenarios w/o housing)
Params.houseservices=0.5;
Params.f_htc=0.05; % transaction cost of buying/selling house (is a percent of h+hprime)
Params.f_coll=0.5; % collateral contraint (fraction of house value that can be borrowed)
Params.pv_pct_cost=0.033; % modeling a $15K install for a 5kW unit install
%% Parameters for firm
% Production
Params.alpha_k=0.311; % diminishing returns to capital and energy inputs
Params.alpha_l=0.650; % diminishing returns to labor input
% Capital adjustment costs
Params.capadjconstant=1.21; % term in the capital adjustment cost

% Idiosyncratic productivity shocks
Params.rho_z_firm=0.767;
Params.sigma_z_e_firm=0.211;

%% Parameters for energy
% Idiosyncratic productivity shocks
Params.rho_z_energy=0.767;
Params.sigma_z_e_energy=0.211;

% Set the firm discount factor below (as it is determined in general eqm)
% Params.firmbeta=1/(Params.r_plus1_ypp/(1-Params.tau_cg)); % 1/(1+r)^ypp but returns net of capital gains tax

%% Create our Grids from Scenario and Parameters
[n_d,n_a,n_z,N_j,vfoptions]=Electrify_GridSizeSetup(Params.scenario, Params.J, small_z_no_e, small_model, vfoptions);
[d_grid,a_grid,z_grid,pi_z,jequaloneDist,share_asset_grid,house_grid,pv_grid_hh,k_grid,pv_grid_firm,pv_grid_energy,Params,vfoptions,simoptions]=Electrify_GridSetup(Params.scenario, n_d, n_a, n_z, small_z_no_e, Params, vfoptions, simoptions);

% Set up Transition Path control parameters
if small_T==2
    jpT=3;
else
    jpT=1; % Default: one transition period=1 time period; Could have multiple j's per T when ypp>1
end

last_n_a_dim=@(n_a_field) n_a_field(end);
last_n_a_dims=structfun(last_n_a_dim, n_a);
T=ceil(Params.J*1.4/jpT)+1;
if T==length(Names_i)
    T=T+1;
end
if T==Params.J
    % The toolkit thinks that T and J must be different (T larger to reach equilibrium post J)
    T=T+1;
end
while any(ismember(last_n_a_dims,T))
    T=T+1;
end

%% Remaining Parameters will be set in GE below

%% Now, create the return function

% For households
DiscountFactorParamNames.household={'beta','sj'};

% For firms
DiscountFactorParamNames.firm={'firmbeta'};

% For energy
DiscountFactorParamNames.energy={'energybeta'};

%% Begin setting up to use VFI Toolkit to solve
% vfoptions.howardsgreedy=0;
% vfoptions.howards=80;
% vfoptions.maxhowards=200;
if Params.scenario<3 && small_model==false
    vfoptions.tolerance=10^(-9);
else
    vfoptions.tolerance=10^(-6);
end
% Note that simoptions.tolerance is used very differently than vfoptions.tolerance

% The user can experiment with gridinterplayer=0 (pure discretization) or gridinterplayer=1 (linear interpolation b/w grid points).
% If gridinterplayer=1, then you must set vfoptions.divideandconquer=1 (required for transition).
vfoptions.gridinterplayer.household  = 0;
vfoptions.level1n.household          = 5;
vfoptions.divideandconquer.household = 0; % Divide and conquer presently works with at most 2 endogenous states; we have 4+Experience Asset
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

if small_model
    P0=[2,2,4.8,1.9];
else
    P0=[2,2,4.8,2.6];
end
Params.P0=P0(Params.scenario); % This price is not 1 because we need price for older and younger agents to balance
% We build a simple model of acquiring and disposing of stock over a lifetime
Params.S_agej_first=ceil(20/Params.ypp); % the age at which we start acquiring more stock than noise
Params.S_agej_peak_first=Params.Jr-1; % the age of first peak acquisition
Params.S_agej_peak_last=Params.Jr+ceil(5/Params.ypp); % the age of last peak acquisition
Params.S_agej_last=Params.J-ceil(5/Params.ypp); % the age of final disposal

%% Agents age distribution
AgeWeightsParamNames=struct('household',{{'mewj'}}); % So VFI Toolkit knows which parameter is the mass of agents of each age

% Solved by GE

% Some initial values/guesses for variables that will be determined in general eqm
Params.w=1;
Params.pension=0.4; % Initial guess (this will be determined in general eqm)
Params.max_benefit=0.4; % Initial guess (this will be determined in general eqm)
% Params.G=0.1; % Government expenditure

% And some initial values/guesses for AggVar values that will be calculated while calculating the general eqm
Params.EnergyCosts_h=0.3; % Energy used by households
Params.EnergyCosts_f=0.7; % Energy used by firms
Params.CarbonCosts_h=0.05; % Carbon tax paid by households
Params.CarbonCosts_f=0.3; % Carbon tax by firms

%% General eqm variables
if Params.scenario<3
    GEPriceParamNames={'w','D','P0','pension'};
elseif Params.scenario<4
    GEPriceParamNames={'w','D','P0','pension'};
else
    GEPriceParamNames={'w','P0','pension', 'max_benefit'};
end
heteroagentoptions.constrainpositive=GEPriceParamNames;

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
% Note also we must differentiate based on Scenarios...

% General Equilibrium conditions (these should evaluate to zero in general equilibrium)
GeneralEqmEqns.sharemarket=@(S) S-1; % mass of all shares equals one
GeneralEqmEqns.labormarket=@(L_h,L_f) (L_h-L_f)*max(2,Params.ypp); % labor supply of households equals labor demand of firms (scaled by ypp)
GeneralEqmEqns.pensions=@(PensionSpending,PayrollTaxRevenue,BenefitSpending) PensionSpending-(PayrollTaxRevenue-BenefitSpending); % Retirement benefits equal Payroll tax revenue (pension*fractionretired-tau*w*H) less benefit
GeneralEqmEqns.benefits=@(PensionSpending,PayrollTaxRevenue,BenefitSpending) BenefitSpending-(PayrollTaxRevenue-PensionSpending); % Welfare benefits equal Payroll tax revenue (benefit-tau*w*H) less pendsions
% GeneralEqmEqns.firmdiscounting=@(firmbeta,r,tau_cg) firmbeta-1/(1+r/(1-tau_cg)); % Firms discount rate is related to market return rate
if Params.scenario<4
    GeneralEqmEqns.dividends=@(dividend,D) dividend-D; % That the dividend households receive equals that which firms give
    GeneralEqmEqns.ShareIssuance=@(Sissued,P0,D,tau_cg,tau_d,r) ...
        P0-((((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg))-Sissued); % P0=P-S, but substitute for P (see derivation inside the return fn)
end
GeneralEqmEqns.CapitalOutputRatio=@(K,L_f,TargetKdivL) (K/L_f-TargetKdivL)/100; % Ratio not based on ypp

Params_Lhscale=Params;
Params_Lhscale.Lhscale=1;
[V_Lhscale, Policy_Lhscale]=ValueFnIter_MixHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid, a_grid, z_grid, pi_z,ReturnFn, Params_Lhscale, DiscountFactorParamNames, vfoptions);
StationaryDist_Lhscale=StationaryDist_MixHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames, Policy_Lhscale,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_Lhscale,simoptions);

%% Test
% Note: Because we used simoptions we must include this as an input
FnsToEvaluate_final.L_h=FnsToEvaluate.L_h;
FnsToEvaluate_final.L_f=FnsToEvaluate.L_f;
AggVars_final=EvalFnOnAgentDist_AggVars_MixHorz_Case1_PType(StationaryDist_Lhscale,Policy_Lhscale, FnsToEvaluate_final, Params_Lhscale, n_d, n_a, n_z,N_j,Names_i,d_grid, a_grid, z_grid,simoptions);
Params.Lhscale=Params_Lhscale.Lhscale*AggVars_final.L_f.Mean/AggVars_final.L_h.Mean;
fprintf("Setting Lhscale to %.2f (with Lhscale==%.2f, L_h was %.2f, L_f was %.2f) \n", Params.Lhscale, Params_Lhscale.Lhscale, AggVars_final.L_h.Mean, AggVars_final.L_f.Mean);
clear FnsToEvaluate_final
FnsToEvaluate_final.S=FnsToEvaluate.S;
FnsToEvaluate_final.D=FnsToEvaluate.D;
[Params.P0,Params.D,V_init,Policy_init,StationaryDist_init]=Calibrate_P0(Params,AgeWeightsParamNames,PTypeDistParamNames,DiscountFactorParamNames,FnsToEvaluate_final,ReturnFn,Policy_Lhscale,jequaloneDist,StationaryDist_Lhscale,n_d, n_a, n_z,N_j,{'household'},Names_i,d_grid,a_grid,z_grid,pi_z,vfoptions,simoptions);

clear Params_Lhscale V_Lhscale Policy_Lhscale StationaryDist_Lhscale FnsToEvaluate_final AggVars_final

%% Now solve the whole value function iteration problem with Lhscale set, just to check that things are working before we go to General Equilbrium
% disp('Test ValueFnIter')
% tic;
% Note: z_grid and pi_z, this will be ignored due to presence of vfoptions.z_grid_J and vfoptions.pi_z_J
% [V_init, Policy_init]=ValueFnIter_MixHorz_PType(n_d,n_a,n_z,N_j,Names_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
% toc

% Plot some things from the firm perspective
if Params.scenario==3
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
    % We can plot V as a 3d plot (surf is matlab command for 3d plot)
    figure(5)
    subplot(2,1,1);
    % Plot F as K vs PV
    surf(pv_grid_firm,k_grid,max(sum(V_init.firm.*reshape(pi_z.firm(:,ceil(n_z.firm/2)),1,1,[]),3),0))
    title('Value function: F as K vs PV')
    xlabel('PV')
    ylabel('K')
end

if false && Params.scenario>=3 % shares vs. assets not possible in Scenarios 1 and 2
    % household = [S+A, Car, House, PV, z, agej]
    [~,unitassetindex]=min(abs(share_asset_grid-1));
    asset_grid=share_asset_grid(share_asset_grid<=1);
    share_grid=share_asset_grid(share_asset_grid>=1);
    % We can plot V as a 3d plot (surf is matlab command for 3d plot)
    figure(6)
    for row=1:3
        for col=1:3
            subplot(3,3,(row-1)*3+col);
            if Params.ypp<6
                agej=(row-1)*6+col*2-1;
            else
                agej=(row-1)*3+col;
            end
            if agej>Params.J
                % This happens when max_age < 100
                break
            end
            % Plot F as S vs H
            if small_z_no_e
                zval=1;
            else
                zval=4;
            end
            if Params.scenario<4
                z_assets_top=V_init.household(1:length(asset_grid),2,1,zval,agej); z_assets_bot=V_init.household(1:length(asset_grid),1,1,zval,agej);
                z_shares_top=V_init.household(length(asset_grid):end-1,2,1,zval,agej); z_shares_bot=V_init.household(length(asset_grid):end-1,1,1,zval,agej);
            else
                z_assets_top=V_init.household(1:length(asset_grid),2-small_model,2,2,zval,agej); z_assets_bot=V_init.household(1:length(asset_grid),1,1,1,zval,agej);
                z_shares_top=V_init.household(length(asset_grid):end-1,2-small_model,2,2,zval,agej); z_shares_bot=V_init.household(length(asset_grid):end-1,1,1,1,zval,agej);
            end
            % z_mat=[(z_assets_top(ceil(linspace(1,length(asset_grid)/2,8)))-z_assets_bot(ceil(linspace(1,length(asset_grid)/2,8))))'; ...
            %     (z_assets_top(ceil(linspace(length(asset_grid)/2,length(asset_grid),8)))-z_assets_bot(ceil(linspace(length(asset_grid)/2,length(asset_grid),8))))'; ...
            %     reshape(z_shares_top,[14,8])-reshape(z_shares_bot,[14,8]) ];
            if small_model
                stride=5;
            else
                stride=8;
            end
            z_mat=[(z_assets_top(ceil(linspace(1,length(asset_grid)/2,stride))))'; ...
                (z_assets_top(ceil(linspace(length(asset_grid)/2,length(asset_grid),stride))))'; ...
                reshape(z_shares_top,[14,stride]) ];
            surf(linspace(asset_grid(ceil(end/2)),asset_grid(end),stride),-1:14,z_mat)
            title(sprintf('Value function: F as S vs A at age %d', (agej-1)*Params.ypp+agejshifter+1))
            xlabel('A')
            ylabel('S')
        end
    end
end


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


%% Test -- Already set when calibrating P0
% disp('Test StationaryDist')
% StationaryDist_init=StationaryDist_MixHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy_init,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);

% Calculate the life-cycle profiles
AgeConditionalStats_init=LifeCycleProfiles_MixHorz_PType(StationaryDist_init,Policy_init,FnsToEvaluate2.S,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
Params.S_agej_first=max(find(AgeConditionalStats_init.household.Mean>0.1,1,'first')-1,1);
Params.S_agej_last=min(find(AgeConditionalStats_init.household.Mean>0.5,1,'last')+1,length(AgeConditionalStats_init.household.Mean)); % warmglow creates extra long tail we want to ignore
[~,S_agej_peak]=max(AgeConditionalStats_init.household.Mean);
S_peak_inflection_value=0.95*AgeConditionalStats_init.household.Mean(S_agej_peak);
for S_agej_peak_first=S_agej_peak:-1:Params.S_agej_first
    if AgeConditionalStats_init.household.Mean(S_agej_peak_first)<S_peak_inflection_value
        break
    end
end
Params.S_agej_peak_first=S_agej_peak_first;
for S_agej_peak_last=S_agej_peak:Params.S_agej_last
    if AgeConditionalStats_init.household.Mean(S_agej_peak_last)<S_peak_inflection_value
        break
    end
end
Params.S_agej_peak_last=S_agej_peak_last;

%% Test
% Note: Because we used simoptions we must include this as an input
disp('Test AggVars')
AggVars=EvalFnOnAgentDist_AggVars_MixHorz_Case1_PType(StationaryDist_init, Policy_init, FnsToEvaluate2, Params, n_d, n_a, n_z,N_j,Names_i, d_grid, a_grid, z_grid,simoptions);

% Next few lines were used to try a few parameter values so as to get a
% decent initial guess before actually solving the general equilibrium
fprintf('Check: L_h, L_f, K \n')
[AggVars.L_h.Mean,AggVars.L_f.Mean,AggVars.K.Mean]
fprintf('Check: K/L_f (should be about 2.03) \n')
AggVars.K.Mean/AggVars.L_f.Mean
if Params.scenario<3
    fprintf('Check: S \n')
    [AggVars.S.Mean]
    % fprintf('Check: ShareIssuance GE condition \n')
    Params.P0-((((1-Params.tau_cg)*Params.P0 + (1-Params.tau_d)*Params.D)/(1+Params.r-Params.tau_cg))-AggVars.S.Mean)
elseif Params.scenario<4
    fprintf('Check: S, A, H, PV_h\n')
    [AggVars.S.Mean,AggVars.A.Mean,AggVars.H.Mean,AggVars.PV_h.Mean]
else
    fprintf('Check: S, A, H, PV_h, Benefit, pvnew_f, pv_f, Params.pvinstalled_firm\n')
    [AggVars.S.Mean,AggVars.A.Mean,AggVars.H.Mean,AggVars.PV_h.Mean,AggVars.Benefit.Mean,AggVars.pvnew_f.Mean,AggVars.pv_f.Mean,Params.pvinstalled_firm]
end


if Params.scenario==3
    Electrify_CustomModelStats(V_init,Policy_init,StationaryDist_init,Params,FnsToEvaluate2,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,[],vfoptions,simoptions)
elseif Params.scenario==4
    Electrify_4CustomModelStats(V_init,Policy_init,StationaryDist_init,Params,FnsToEvaluate2,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,[],vfoptions,simoptions)
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
        heteroagentoptions.toleranceGEprices=10^(-4);
        heteroagentoptions.toleranceGEcondns=10^(-4); % This is the hard one
        if solve_TPath
            % heteroagentoptions.maxiter=200;
        end
    else
        heteroagentoptions.toleranceGEprices=10^(-3);
        heteroagentoptions.toleranceGEcondns=10^(-2); % This is the hard one
        heteroagentoptions.maxiter=15*(1+logical(small_z_no_e)+logical(small_model));                % About 3 hours for 35 iterations

        if Params.scenario<4
            heteroagentoptions.CustomModelStats=@(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions) ...
                Electrify_CustomModelStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions);
        else
            heteroagentoptions.CustomModelStats=@(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions) ...
                Electrify_4CustomModelStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions);
        end
    end

    [p_eqm_init,GEcondns_init]=HeteroAgentStationaryEqm_MixHorz_PType(n_d, n_a, n_z, N_j,Names_i,[],pi_z,d_grid,a_grid,z_grid,jequaloneDist,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,AgeWeightsParamNames,PTypeDistParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
    % p_eqm contains the general equilibrium parameter values
    % Put this into Params so we can calculate things about the initial equilibrium
    % GEcondns tells us the values of the GeneralEqmEqns, should be near zero
    Params.pension=p_eqm_init.pension;
    Params.max_benefit=p_eqm_init.max_benefit;
    % Plot twist: we are going to use the good AggVar value of bequests as elements of the GEqm so that we can transition from init to final across demographic changes
    p_eqm_init.AccidentBeqS=Params.AccidentBeqS;
    p_eqm_init.G=Params.tau_d*Params.D*Params.ypp+AggVars.CapitalGainsTaxRevenue.household.Mean+AggVars.CorpTaxRevenue.firm.Mean;
    % To use a '_tminus1' variable we must include its initial value
    % transpathoptions.initialvalues.BeqleftS= Params.AccidentBeqS;
    if Params.scenario>2
        % Plot twist as above...
        p_eqm_init.AccidentBeqAH=Params.AccidentBeqAH;
        % transpathoptions.initialvalues.BeqleftAH= Params.AccidentBeqAH;
    end
    transpathoptions.initialvalues.pvinstalled_firm=0;
    transpathoptions.initialvalues.pvinstalled_energy=0;

    Params.P0=p_eqm_init.P0;
    % Params.G=p_eqm_init.G;
    Params.w=p_eqm_init.w;
    % Params.firmbeta=p_eqm_init.firmbeta;

    % Re-Calculate a few things related to the general equilibrium.
    [V_init, Policy_init]=ValueFnIter_MixHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    StationaryDist_init=StationaryDist_MixHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy_init,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);

    % Calculate various stats
    AllStats_init=EvalFnOnAgentDist_AllStats_MixHorz_PType(StationaryDist_init,Policy_init,FnsToEvaluate2,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
    % Calculate the life-cycle profiles
    AgeConditionalStats_init=LifeCycleProfiles_MixHorz_PType(StationaryDist_init,Policy_init,FnsToEvaluate2,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
    
    if abs(1-AllStats_init.L_h.household.Mean/AllStats_init.L_f.firm.Mean)>0.05
        warning("L_h and L_f have diverged; check Params.Lhscale")
    end

    % Note: Only part of this initial stationary general eqm we actually 'need'
    % is the agent distribution. Rest is just out of interest.
    
    AgentDist_init=StationaryDist_init; % Just to emphasize that there is no need for the
       % initial agent distribution to be a stationary dist (it is in this
       % example, but does not need to be for transition paths)
    
    % Just to see it...
    GEcondns_init


    %% Plot the life cycle profiles of capital and labour for the initial eqm.
    % Can just use the same FnsToEvaluate as before
    AgeConditionalStats_init=LifeCycleProfiles_MixHorz_PType(StationaryDist_init,Policy_init,FnsToEvaluate2,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

    figure_c=figure(10);
    if Params.scenario==4
        rows=4;
    else
        rows=3;
    end
    if Params.scenario<3
        subplot(rows,1,1); plot(1:1:Params.J,AgeConditionalStats_init.L_h.Mean)
        title('Life Cycle Profile: Effective Labour Supply')
        subplot(rows,1,2); plot(1:1:Params.J,AgeConditionalStats_init.S.Mean)
        title('Life Cycle Profile: Share holdings')
        subplot(rows,1,3); plot(1:1:Params.J,Params.kappa_j)
        title('Life Cycle Profile: kappa_j')
    else
        subplot(rows,2,1); plot(1:1:Params.J,AgeConditionalStats_init.L_h.Mean)
        title('Life Cycle Profile: Effective Labour Supply')
        subplot(rows,2,3); plot(1:1:Params.J,AgeConditionalStats_init.S.Mean)
        title('Life Cycle Profile: Share holdings')
        subplot(rows,2,5); plot(1:1:Params.J,Params.kappa_j)
        title('Life Cycle Profile: kappa_j')
    end
    if Params.scenario>2
        subplot(rows,2,2); plot(1:1:Params.J,AgeConditionalStats_init.A.Mean)
        title('Life Cycle Profile: Asset holdings')
        subplot(rows,2,4); plot(1:1:Params.J,AgeConditionalStats_init.H.Mean)
        title('Life Cycle Profile: House holdings')
        subplot(rows,2,6); plot(1:1:Params.J,AgeConditionalStats_init.PV_h.Mean)
        title('Life Cycle Profile: Solar PV installed')
    end
    if Params.scenario>3
        subplot(4,2,7); plot(1:1:Params.J,Params.carservices_j)
        title('Life Cycle Profile: carservices_j')
        subplot(4,2,8); plot(1:1:Params.J,AgeConditionalStats_init.Car_none.Mean)
        title('Life Cycle Profile: Car (blue=none,red=petrol,yellow=ev)')
        hold on
        plot(1:1:Params.J,AgeConditionalStats_init.Car_petrol.Mean)
        plot(1:1:Params.J,AgeConditionalStats_init.Car_ev.Mean)
        hold off
    end

    saveas(figure_c,'./SavedOutput/Graphs/Electrify_LifeCycleProfiles_init','pdf')

    if max(AgeConditionalStats_init.S.Maximum)==share_asset_grid(end)
        warning("share_grid maximum reached")
    end
    if Params.scenario>2
        if max(AgeConditionalStats_init.H.Maximum)==house_grid(end)
            warning("house_grid maximum reached")
        end
        if max(AgeConditionalStats_init.PV_h.Maximum)==pv_grid_hh(end)
            warning("pv_grid_hh maximum reached")
        end
    end

    %% Calculate some aggregates and print findings about them

    AggVars=EvalFnOnAgentDist_AggVars_MixHorz_Case1_PType(StationaryDist_init, Policy_init, FnsToEvaluate3, Params, n_d, n_a, n_z,N_j, Names_i, d_grid, a_grid, z_grid,simoptions);

    Y=AggVars.Output_f.Mean;

    P=((1-Params.tau_cg)*Params.P0 + (1-Params.tau_d)*Params.D)/(1+Params.r-Params.tau_cg);

    G=Params.tau_d*Params.D+AggVars.CapitalGainsTaxRevenue.household.Mean+AggVars.CorpTaxRevenue.firm.Mean; % If G influenced any GEqm equations, we'd need it to equilibrate with them

    % Calculate the aggregate TFP as output/((capital^alpha_k)*(labor^alpha_l))
    AggregateTFP=Y/((AggVars.K.Mean^Params.alpha_k)*(AggVars.L_f.Mean^Params.alpha_l));

    % Total value of firms
    temp=V_init.firm.*StationaryDist_init.firm;
    temp(StationaryDist_init.firm==0)=0; % Get rid of points that have V=-inf but zero mass which would give nan
    TotalValueOfFirms=sum(temp(isfinite(temp)));

    fileID = fopen('SavedOutput\aggs_init.txt','w');
    fprintf(fileID,'Following are some aggregates of the model economy (Scenario %d): \n', Params.scenario);
    fprintf(fileID,'Output: Y=%8.2f \n',AggVars.Output_f.Mean);
    fprintf(fileID,'Aggregate TFP: Y=%8.2f \n',AggregateTFP);
    fprintf(fileID,'Capital-Output ratio (firm side): K/Y=%8.2f \n',AggVars.K.Mean/Y);
    if Params.scenario<3
        fprintf(fileID,'Total share value (HH side): P*S (%.2f) = %8.2f\n',P*AggVars.S.Mean,P*AggVars.S.Mean);
    else
        fprintf(fileID,'Total share+asset value (HH side): P*S (%.2f) + A (%.2f) = %8.2f\n',P*AggVars.S.Mean,AggVars.A.Mean,P*AggVars.S.Mean+AggVars.A.Mean);
        fprintf(fileID,'Total house value (HH side): H=%8.2f \n',AggVars.H.Mean);
        fprintf(fileID,'Total bad debt (HH side): P*S=%8.2f \n',AggVars.BadDebt.Mean);
    end
    fprintf(fileID,'Total firm value (firm side): Value of firm=%8.2f \n',TotalValueOfFirms);
    fprintf(fileID,'Consumption-Output ratio: C/Y=%8.2f \n',AggVars.Consumption.Mean/(Y*Params.ypp));
    fprintf(fileID,'Government-to-Output ratio: G/Y=%8.2f \n', G/Y);
    fprintf(fileID,'Wage: w=%8.2f \n',Params.w);
    fclose(fileID);

    type 'SavedOutput\aggs_init.txt'

    %%
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
    solve_GE=solve_GE_temp; solve_TPath=solve_TPath_temp;
    % load tpathElectrifyA.mat
else
    load tpathElectrifyA.mat
end
clear solve_GE_temp solve_TPath_temp

%% Solve for final stationary general eqm with Params at time T
% Must ensure that our T does not conflict with any other dimensions
if small_T==1
    T_end=length(Names_i)+1;
else
    T_end=T;
end
while any(ismember(last_n_a_dims,T_end))
    T_end=T_end+1;
end
if T_end>T
    error("impossible dimensions for T and T_end")
end

% ParamPath on Ek (Energy Use) and ek (Energy Efficiency)
% More transitions down in the demographics section
ParamPath.pvinstalled_firm=floor(transition_target_firm*linspace(0,Params.pvmax_firm,T)); Params.pvinstalled_firm=ParamPath.pvinstalled_firm(1);
ParamPath.pvinstalled_energy=floor(transition_target_energy*linspace(0,Params.pvmax_energy,T));Params.pvinstalled_energy=ParamPath.pvinstalled_energy(1);
ParamPath.Ek=linspace(1,1.01,T); Params.Ek=ParamPath.Ek(1);
ParamPath.ek=linspace(1,1.01,T); Params.ek=ParamPath.ek(1);
ParamPath.carbon_tax=linspace(42,2450,T); Params.carbon_tax=ParamPath.carbon_tax(1);
ParamPath.energy_pct_brown=linspace(0.55,0.05,T); Params.energy_pct_brown=ParamPath.carbon_tax(1);

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
ParamPath.cpi_energy=1.01.^((0:Params.J-1)*Params.ypp)-1; % Params.J periods of energy cost increases
% Translate energy periods (j) into transition periods
ParamPath.cpi_energy(end+1:T*jpT)=ParamPath.cpi_energy(end); % Energy cost increases extended to the jth period implied by final T
ParamPath.cpi_energy=ParamPath.cpi_energy(1:T); % Energy cost increases on a per transition period basis
Params.cpi_energy=ParamPath.cpi_energy(1);

% Params.P0=2.05;

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
    ParamPath.mewj=ParamPath.mewj./((1+Params.n).^(Params.ypp*(1:Params.J)-1)); % Population shrinks in the N_j dimension
    ParamPath.mewj=ParamPath.mewj.*((1+Params.n).^(Params.ypp*jpT*((1:T)-1)))'; % Population grows in the T dimension
    ParamPath.mewj=ParamPath.mewj./sum(ParamPath.mewj,2); % normalize age-masses to sum to one
    % Looking at ParamPath.mewj you can see that as tt increases, the mass at older ages increases

    Params.pvinstalled_firm=ParamPath.pvinstalled_firm(T_end);
    Params.pvinstalled_energy=ParamPath.pvinstalled_energy(T_end);
    Params.Ek=ParamPath.Ek(T_end);
    Params.ek=ParamPath.ek(T_end);
    Params.carbon_tax=ParamPath.carbon_tax(T_end);
    Params.energy_pct_brown=ParamPath.energy_pct_brown(T_end);
    Params.sj=ParamPath.sj(T_end,:); % conditional survival probabilities
    Params.mewj=ParamPath.mewj(T_end,:);
    Params.cpi=ParamPath.cpi(T_end);
    Params.cpi_energy=ParamPath.cpi_energy(T_end);

    %% TEST
    labor=0.85;
    buyhouse=0;
    saprime=0.25;
    cprime=0;
    hprime=0;
    sa=1.0;
    car=0;
    h=0;
    solarpv=0;
    z=0;
    e=0;
    agej=6;
    rentprice=Params.rentprice;
    houseservices=Params.houseservices;
    carservices_j=Params.carservices_j;
    cpi_energy=Params.cpi_energy;
    energy_pct_cost=Params.energy_pct_cost;
    energy_pct_brown=Params.energy_pct_brown;
    carbon_tax=Params.carbon_tax;
    v1=Electrify_4HouseholdReturnFn( ...
        labor,buyhouse,saprime,cprime+2,hprime,sa,car,h,solarpv,z,e, ...
        Params.pension,Params.max_benefit,Params.AccidentBeqS,Params.AccidentBeqAH,Params.w,Params.P0,Params.D,Params.sigma,Params.psi,Params.eta,Params.sigma_h,Params.sigma_c,Params.kappa_j(agej),Params.warmglow1,Params.warmglow2,Params.tau_l,Params.tau_d,Params.tau_cg,Params.S_agej_first,Params.S_agej_peak_first,Params.S_agej_peak_last,Params.S_agej_last, ...
        Params.ypp,agej,Params.Jr,Params.J,Params.r,Params.r_wedge,Params.f_htc,Params.minhouse,rentprice,Params.f_coll,houseservices,carservices_j(agej),cpi_energy,Params.pv_pct_cost,energy_pct_cost,energy_pct_brown,carbon_tax); % Level=0, Refine=0
    v2=Electrify_4HouseholdReturnFn( ...
        labor,buyhouse,saprime+0.45,cprime+1,hprime,sa,car,h,solarpv,z,e, ...
        Params.pension,Params.max_benefit,Params.AccidentBeqS,Params.AccidentBeqAH,Params.w,Params.P0,Params.D,Params.sigma,Params.psi,Params.eta,Params.sigma_h,Params.sigma_c,Params.kappa_j(agej),Params.warmglow1,Params.warmglow2,Params.tau_l,Params.tau_d,Params.tau_cg,Params.S_agej_first,Params.S_agej_peak_first,Params.S_agej_peak_last,Params.S_agej_last, ...
        Params.ypp,agej,Params.Jr,Params.J,Params.r,Params.r_wedge,Params.f_htc,Params.minhouse,rentprice,Params.f_coll,houseservices,carservices_j(agej),cpi_energy,Params.pv_pct_cost,energy_pct_cost,energy_pct_brown,carbon_tax); % Level=0, Refine=0
    fprintf("v1-v2: %.2f - %.2f = %.2f \n", v1, v2, v1-v2);

    %% Let's take a quick look at what we have calculated, namely V and Policy

    % Evaluate the final stationary general eqm
    disp('Test ValueFnIter')
    [V_final, Policy_final]=ValueFnIter_MixHorz_PType(n_d,n_a,n_z,N_j,Names_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
    disp('Test StationaryDist')
    StationaryDist_final=StationaryDist_MixHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy_final,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);

    FnsToEvaluate_final.L_h=FnsToEvaluate.L_h;
    FnsToEvaluate_final.S=FnsToEvaluate.S;
    FnsToEvaluate_final.L_f=FnsToEvaluate.L_f;
    FnsToEvaluate_final.D=FnsToEvaluate.D;
    AggVars_final=EvalFnOnAgentDist_AggVars_MixHorz_Case1_PType(StationaryDist_final,Policy_final, FnsToEvaluate_final, Params, n_d, n_a, n_z,N_j,Names_i,d_grid, a_grid, z_grid,simoptions);
    Params.D=AggVars_final.D.firm.Mean;
    Lhscale_final_ratio=AggVars_final.L_f.Mean/AggVars_final.L_h.Mean;
    while abs(Lhscale_final_ratio-1)>0.1
        % Make these asymmetric in size so they don't oscillate too much
        if Lhscale_final_ratio<1
            Params.w=Params.w*0.98;
        else
            Params.w=Params.w*1.06;
        end
        if AggVars_final.S.Mean<0.1
            Params.P0=Params.P0/2;
        elseif AggVars_final.S.Mean<0.5
            Params.P0-Params.P0*.9;
        elseif AggVars_final.S.Mean>4
            Params.P0=Params.P0*2;
        elseif AggVars_final.S.Mean>2
            Params.P0=Params.P0*1.4;
        elseif AggVars_final.S.Mean>1.4
            Params.P0=Params.P0*1.2;
        end
        [V_final, Policy_final]=ValueFnIter_MixHorz_PType(n_d,n_a,n_z,N_j,Names_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
        StationaryDist_final=StationaryDist_MixHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy_final,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
        AggVars_final=EvalFnOnAgentDist_AggVars_MixHorz_Case1_PType(StationaryDist_final,Policy_final, FnsToEvaluate_final, Params, n_d, n_a, n_z,N_j,Names_i,d_grid, a_grid, z_grid,simoptions);
        Params.D=AggVars_final.D.firm.Mean;
        Lhscale_final_ratio=AggVars_final.L_f.Mean/AggVars_final.L_h.Mean;
    end
    fprintf("Setting w to %.2f, P0 to %.2f, D to %.2f \n", Params.w, Params.P0, Params.D);
    % Re-calculate P0 as we have changed some important parameters...
    FnsToEvaluate_final=rmfield(FnsToEvaluate_final,{'L_h','L_f'});
    Plag=Params.P0;
    [Params.P0,Params.D,V_final,Policy_final,StationaryDist_final]=Calibrate_P0(Params,AgeWeightsParamNames,PTypeDistParamNames,DiscountFactorParamNames,FnsToEvaluate_final,ReturnFn,Policy_final,jequaloneDist,StationaryDist_final,n_d, n_a, n_z,N_j,{'household'},Names_i,d_grid,a_grid,z_grid,pi_z,vfoptions,simoptions);

    % Calculate the life-cycle profiles
    AgeConditionalStats_final=LifeCycleProfiles_MixHorz_PType(StationaryDist_final,Policy_final,FnsToEvaluate2.S,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
    Params.S_agej_first=max(find(AgeConditionalStats_final.household.Mean>0.1,1,'first')-1,1);
    Params.S_agej_last=min(find(AgeConditionalStats_final.household.Mean>0.5,1,'last')+1,length(AgeConditionalStats_final.household.Mean)); % warmglow creates extra long tail we want to ignore
    [~,S_agej_peak]=max(AgeConditionalStats_final.household.Mean);
    S_peak_inflection_value=0.95*AgeConditionalStats_final.household.Mean(S_agej_peak);
    for S_agej_peak_first=S_agej_peak:-1:Params.S_agej_first
        if AgeConditionalStats_final.household.Mean(S_agej_peak_first)<S_peak_inflection_value
            break
        end
    end
    Params.S_agej_peak_first=S_agej_peak_first;
    for S_agej_peak_last=S_agej_peak:Params.S_agej_last
        if AgeConditionalStats_final.household.Mean(S_agej_peak_last)<S_peak_inflection_value
            break
        end
    end
    % Could shift S_agej_peak_last to S_agej_peak if S_agej_peak < S_agej_last and S_agej_peak_last==S_agej_last 
    Params.S_agej_peak_last=S_agej_peak_last;

    disp('Test AggVars')
    % [V_final, Policy_final]=ValueFnIter_MixHorz_PType(n_d,n_a,n_z,N_j,Names_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
    % StationaryDist_final=StationaryDist_MixHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy_final,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
    AggVars=EvalFnOnAgentDist_AggVars_MixHorz_Case1_PType(StationaryDist_final, Policy_final, FnsToEvaluate2, Params, n_d, n_a, n_z,N_j,Names_i, d_grid, a_grid, z_grid,simoptions);
    Params.D=AggVars.D.firm.Mean;

    % Next few lines were used to try a few parameter values so as to get a
    % decent initial guess before actually solving the general equilibrium
    fprintf('Check: L_h, L_f, K, D \n')
    [AggVars.L_h.Mean,AggVars.L_f.Mean,AggVars.K.Mean,AggVars.D.Mean]
    fprintf('Check: K/L_f (should be about 2.03) \n')
    AggVars.K.Mean/AggVars.L_f.Mean
    if Params.scenario<3
        fprintf('Check: S \n')
        [AggVars.S.Mean]
        fprintf('Check: ShareIssuance GE condition \n')
        Params.P0-((((1-Params.tau_cg)*Params.P0 + (1-Params.tau_d)*Params.D)/(1+Params.r-Params.tau_cg))-AggVars.S.Mean)
    elseif Params.scenario<4
        fprintf('Check: S, A, H, PV_h\n')
        [AggVars.S.Mean,AggVars.A.Mean,AggVars.H.Mean,AggVars.PV_h.Mean]
    else
        fprintf('Check: S, A, H, PV_h, Benefit, pvnew_f, pv_f, Params.pvinstalled_firm \n')
        [AggVars.S.Mean,AggVars.A.Mean,AggVars.H.Mean,AggVars.PV_h.Mean,AggVars.Benefit.Mean,AggVars.pvnew_f.Mean,AggVars.pv_f.Mean,Params.pvinstalled_firm]
    end

    % And now, the GE for the final conditions!
    [p_eqm_final,GEcondns_final]=HeteroAgentStationaryEqm_MixHorz_PType(n_d,n_a,n_z,N_j,Names_i,[],pi_z,d_grid,a_grid,z_grid,jequaloneDist,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,AgeWeightsParamNames,PTypeDistParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
    % Done, the general eqm prices are in p_eqm
    % GEcondns tells us the values of the GeneralEqmEqns, should be near zero
    Params.pension=p_eqm_final.pension;
    Params.max_benefit=p_eqm_final.max_benefit;
    p_eqm_final.AccidentBeqS=Params.AccidentBeqS;
    if Params.scenario>2
        p_eqm_final.AccidentBeqAH=Params.AccidentBeqAH;
    end
    % Params.firmbeta=p_eqm_final.firmbeta;
    Params.P0=p_eqm_final.P0;
    p_eqm_final.G=Params.tau_d*Params.D*Params.ypp+AggVars.CapitalGainsTaxRevenue.household.Mean+AggVars.CorpTaxRevenue.firm.Mean;
    Params.w=p_eqm_final.w;


    % Calculate various stats
    AllStats_final=EvalFnOnAgentDist_AllStats_MixHorz_PType(StationaryDist_final, Policy_final, FnsToEvaluate2, Params, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid, z_grid,simoptions);
    % Calculate the life-cycle profiles
    AgeConditionalStats_final=LifeCycleProfiles_MixHorz_PType(StationaryDist_final,Policy_final, FnsToEvaluate2,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

    % Note: Only part of this final stationary general eqm we actually 'need'
    % is the value fn (although we likely want p_eqm_final for initial guess of PricePath0). 
    % Rest is just out of interest.
    
    % Double-check that the general eqm is accurate before we start the
    % transition path, because if it is not then it won't solve
    GEcondns_final

    %% Plot the life cycle profiles of capital and labour for the final eqm.
    AgeConditionalStats_final=LifeCycleProfiles_MixHorz_PType(StationaryDist_final,Policy_final,FnsToEvaluate2,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

    figure_d=figure(11);
    if Params.scenario==4
        rows=4;
    else
        rows=3;
    end
    if Params.scenario<3
        subplot(rows,1,1); plot(1:1:Params.J,AgeConditionalStats_final.L_h.Mean)
        title('Life Cycle Profile: Effective Labour Supply')
        subplot(rows,1,2); plot(1:1:Params.J,AgeConditionalStats_final.S.Mean)
        title('Life Cycle Profile: Share holdings')
        subplot(rows,1,3); plot(1:1:Params.J,Params.kappa_j)
        title('Life Cycle Profile: kappa_j')
    else
        subplot(rows,2,1); plot(1:1:Params.J,AgeConditionalStats_final.L_h.Mean)
        title('Life Cycle Profile: Effective Labour Supply')
        subplot(rows,2,3); plot(1:1:Params.J,AgeConditionalStats_final.S.Mean)
        title('Life Cycle Profile: Share holdings')
        subplot(rows,2,5); plot(1:1:Params.J,Params.kappa_j)
        title('Life Cycle Profile: kappa_j')
    end
    if Params.scenario>2
        subplot(rows,2,2); plot(1:1:Params.J,AgeConditionalStats_final.A.Mean)
        title('Life Cycle Profile: Asset holdings')
        subplot(rows,2,4); plot(1:1:Params.J,AgeConditionalStats_final.H.Mean)
        title('Life Cycle Profile: House holdings')
        subplot(rows,2,6); plot(1:1:Params.J,AgeConditionalStats_final.PV_h.Mean)
        title('Life Cycle Profile: Solar PV installed')
    end
    if Params.scenario>3
        subplot(4,2,7); plot(1:1:Params.J,Params.carservices_j)
        title('Life Cycle Profile: carservices_j')
        subplot(4,2,8); plot(1:1:Params.J,AgeConditionalStats_final.Car_none.Mean)
        title('Life Cycle Profile: Car (blue=none,red=petrol,yellow=ev)')
        hold on
        plot(1:1:Params.J,AgeConditionalStats_final.Car_petrol.Mean)
        plot(1:1:Params.J,AgeConditionalStats_final.Car_ev.Mean)
        hold off
    end
    saveas(figure_d,'./SavedOutput/Graphs/Electrify_LifeCycleProfiles_final','pdf')

    if max(AgeConditionalStats_final.S.Maximum)==share_asset_grid(end)
        warning("share_grid maximum reached")
    end
    if Params.scenario>2
        if max(AgeConditionalStats_final.H.Maximum)==house_grid(end)
            warning("house_grid maximum reached")
        end
        if max(AgeConditionalStats_final.PV_h.Maximum)==pv_grid_hh(end)
            warning("pv_grid_hh maximum reached")
        end
    end

    %% Calculate some aggregates and print findings about them

    AggVars=EvalFnOnAgentDist_AggVars_MixHorz_Case1_PType(StationaryDist_init, Policy_init, FnsToEvaluate3, Params, n_d, n_a, n_z,N_j, Names_i, d_grid, a_grid, z_grid,simoptions);

    Y=AggVars.Output_f.Mean;

    P=((1-Params.tau_cg)*Params.P0 + (1-Params.tau_d)*Params.D)/(1+Params.r-Params.tau_cg);

    G=Params.tau_d*Params.D+AggVars.CapitalGainsTaxRevenue.household.Mean+AggVars.CorpTaxRevenue.firm.Mean; % If G influenced any GEqm equations, we'd need it to equilibrate with them

    % Calculate the aggregate TFP as output/((capital^alpha_k)*(labor^alpha_l))
    AggregateTFP=Y/((AggVars.K.Mean^Params.alpha_k)*(AggVars.L_f.Mean^Params.alpha_l));

    % Total value of firms
    temp=V_init.firm.*StationaryDist_init.firm;
    temp(StationaryDist_init.firm==0)=0; % Get rid of points that have V=-inf but zero mass which would give nan
    TotalValueOfFirms=sum(temp(isfinite(temp)));

    fileID = fopen('SavedOutput\aggs_final.txt','w');
    fprintf(fileID,'Following are some aggregates of the model economy (Scenario %d): \n', Params.scenario);
    fprintf(fileID,'Output: Y=%8.2f \n',AggVars.Output_f.Mean);
    fprintf(fileID,'Aggregate TFP: Y=%8.2f \n',AggregateTFP);
    fprintf(fileID,'Capital-Output ratio (firm side): K/Y=%8.2f \n',AggVars.K.Mean/Y);
    if Params.scenario<3
        fprintf(fileID,'Total share value (HH side): P*S (%.2f) = %8.2f\n',P*AggVars.S.Mean,P*AggVars.S.Mean);
    else
        fprintf(fileID,'Total share+asset value (HH side): P*S (%.2f) + A (%.2f) = %8.2f\n',P*AggVars.S.Mean,AggVars.A.Mean,P*AggVars.S.Mean+AggVars.A.Mean);
        fprintf(fileID,'Total house value (HH side): H=%8.2f \n',AggVars.H.Mean);
        fprintf(fileID,'Total bad debt (HH side): P*S=%8.2f \n',AggVars.BadDebt.Mean);
    end
    fprintf(fileID,'Total firm value (firm side): Value of firm=%8.2f \n',TotalValueOfFirms);
    fprintf(fileID,'Consumption-Output ratio: C/Y=%8.2f \n',AggVars.Consumption.Mean/(Y*Params.ypp));
    fprintf(fileID,'Government-to-Output ratio: G/Y=%8.2f \n', G/Y);
    fprintf(fileID,'Wage: w=%8.2f \n',Params.w);
    fclose(fileID);

    type 'SavedOutput\aggs_final.txt'
    %%
    solve_TPath_temp=solve_TPath; clear solve_TPath
    save tpathElectrifyB.mat
    solve_TPath=solve_TPath_temp;
else
    load tpathElectrifyB.mat
end % solve_GE_final
clear solve_TPath_temp

ParamPath.Lhscale=linspace(Params.Lhscale,Params.Lhscale*AggVars.L_f.Mean/AggVars.L_h.Mean,T);
Params.Lhscale=ParamPath.Lhscale(T_end);

if solve_TPath
    if Params.scenario==4
        vfoptions.lowmemory.household=3;
    end
    %% Setup for the transition path
    % T=100; % number of periods for transition path
    
    % Already created ParamPath.sj and ParamPath.mewj

    if small_T==1
        T_eq=ceil(T_end/2);
    else
        T_eq=ceil(0.618*T); % Demographic change has stopped and T_eq begins period of transition equilibrium-finding
    end

    ParamPath.AccidentBeqS=[linspace(p_eqm_init.AccidentBeqS,p_eqm_final.AccidentBeqS,T_eq), p_eqm_final.AccidentBeqS*ones(1,T_end-T_eq)];
    if Params.scenario>2
        ParamPath.AccidentBeqAH=[linspace(p_eqm_init.AccidentBeqAH,p_eqm_final.AccidentBeqAH,T_eq), p_eqm_final.AccidentBeqAH*ones(1,T_end-T_eq)];
    end
    ParamPath.G=[linspace(p_eqm_init.G, p_eqm_final.G,T_eq), p_eqm_final.G*ones(1,T_end-T_eq)];

    % Initial guess for general eqm parameters; if small_T, V_final is V_init at T=T_end
    if small_T==1
        paramnames=fieldnames(ParamPath);
        for nn=1:length(paramnames)
            if size(ParamPath.(paramnames{nn}),1)==1
                ParamPath0.(paramnames{nn})=ParamPath.(paramnames{nn})(1,1:T_end);
            else
                ParamPath0.(paramnames{nn})=ParamPath.(paramnames{nn})(1:T_end,:);
            end
        end
    else
        ParamPath0=ParamPath;
    end

    PricePath0.w=[linspace(p_eqm_init.w, p_eqm_final.w,T_eq), p_eqm_final.w*ones(1,T_end-T_eq)];
    % PricePath0.firmbeta=[linspace(p_eqm_init.firmbeta, p_eqm_final.firmbeta,T_eq), p_eqm_final.firmbeta*ones(1,T_end-T_eq)];
    if isfield(p_eqm_init, 'P0')
        PricePath0.P0=[linspace(p_eqm_init.P0, p_eqm_final.P0,T_eq), p_eqm_final.P0*ones(1,T_end-T_eq)];
    end
    PricePath0.pension=[linspace(p_eqm_init.pension, p_eqm_final.pension,T_eq), p_eqm_final.pension*ones(1,T_end-T_eq)];
    PricePath0.max_benefit=[linspace(p_eqm_init.max_benefit, p_eqm_final.max_benefit,T_eq), p_eqm_final.max_benefit*ones(1,T_end-T_eq)];
    % PricePath0.TargetKdivL=2.03*ones(1,T_end);

    % General eqm eqns, same idea as with the stationary general eqm
    % GeneralEqmEqns_Transition.capitalmarket=@(r,alpha_k,alpha_l,delta,K,L) r-(alpha_k*(K^(alpha_k-1))*(L^(alpha_l))-delta); % r=marginal product of capital
    GeneralEqmEqns_Transition.labormarket=@(w,alpha_k,alpha_l,K,L_f) w-(alpha_l)*(K^alpha_k)*(L_f^(alpha_l-1)); % w=marginal product of labor
    % GeneralEqmEqns_Transition.firmdiscounting=GeneralEqmEqns.firmdiscounting;
    % GeneralEqmEqns_Transition.dividends=GeneralEqmEqns.dividends;
    if Params.scenario<4
        GeneralEqmEqns_Transition.ShareIssuance=GeneralEqmEqns.ShareIssuance;
    else
        GeneralEqmEqns_Transition.ShareIssuance=@(Sissued,P0,D,tau_cg,tau_d,r) ...
            P0-((((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg))-Sissued); % P0=P-S, but substitute for P (see derivation inside the return fn)
    end
    GeneralEqmEqns_Transition.pensions=GeneralEqmEqns.pensions;
    GeneralEqmEqns_Transition.benefits=GeneralEqmEqns.benefits;
    % GeneralEqmEqns_Transition.govbudgetbalance=GeneralEqmEqns.govbudget;
    % Note: bequests are left in t-1 and received in t
    % GeneralEqmEqns_Transition.bequestsS=@(BeqleftS_tminus1,AccidentBeqS,n,ypp) BeqleftS_tminus1/(1+n)^ypp-AccidentBeqS; % Accidental share bequests received equal accidental share bequests left
    % if Params.scenario>2
    %     GeneralEqmEqns_Transition.bequestsAH=@(BeqleftAH_tminus1,AccidentBeqAH,n,ypp) BeqleftAH_tminus1/(1+n)^ypp-AccidentBeqAH; % Accidental asset+house bequests received equal accidental asset+house bequests left
    % end
    
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
        'benefits','max_benefit',0,0.03;... % benefits GE condition will be positive if pension is too big, so subtract
        ... % 'govbudgetbalance','G',0,0.03;... % govbudget GE condition will be positive if G is too big, so subtract
        ... % 'bequestsS','AccidentBeqS',1,0.03;... % bequests GE condition will be negative if BeqS is too big, so add
        ... % 'bequestsAH','AccidentBeqAH',1,0.03;... % bequests GE condition will be negative if BeqAH is too big, so add
        };
    % if Params.scenario<3
    %     mask=strcmp(transpathoptions.GEnewprice3.howtoupdate(:,1),'bequestsAH');
    %     transpathoptions.GEnewprice3.howtoupdate(mask,:)=[];
    % elseif Params.scenario==4
    %     for pp=1:size(transpathoptions.GEnewprice3.howtoupdate,1)
    %         transpathoptions.GEnewprice3.howtoupdate{pp,4}=0.008;
    %     end
    % end

    % Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
    % Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
    % A small 'factor' will make the convergence to solution take longer, but too large a value will make it 
    % unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.

    % TESTING -- pensions doesn't depend on PType
    % transpathoptions.GEptype={'pensions'};
    
    %% Solve the transition path
    vfoptions.fastOLG.household=1; simoptions.fastOLG.household=1; % Needs to be set up for transition paths
    % Setup the options relating to the transition path
    transpathoptions.verbose=1;
    transpathoptions.maxiter=2; % default is 1000
    transpathoptions.fastOLG=0; % PTypes will force this on `simoptions`; must we match that energy?
    transpathoptions.graphpricepath=1; % plots of the ParamPath that get updated every interation
    transpathoptions.graphaggvarspath=1; % plots of the AggVarsPath that get updated every iteration
    
    % Running, it was about stuck iterating around 2 or 3*10^(-4) but had clearly solved. So
    transpathoptions.tolerance=4*10^(-2); % default is 10^(-4), which is a very demanding accuracy
    transpathoptions.updateaccuracycutoff=transpathoptions.tolerance/(2*length(GEPriceParamNames));

    transpathoptions.verbose=1;
    transpathoptions.graphpricepath=1; % 1: creates a graph of the 'current' price path which updates each iteration.
    transpathoptions.graphaggvarspath=1; % 1: creates a graph of the 'current' aggregate variables which updates each iteration.
    transpathoptions.graphGEcondns=1;  % 1: creates a graph of the 'current' general eqm conditions which updates each iteration.
    transpathoptions.historyofpricepath=1;
    %%

    save tpathElectrifyC.mat
    % load tpathElectrifyC.mat

    % And go! (with FnsToEvaluate2)
    vfoptions.refine_d.firm=[1,0,1];
    vfoptions.refine_d.energy=[1,0,1];
    [PricePath,GECondnsPath]=TransitionPath_MixHorz_PType(PricePath0, ParamPath0, T_end, V_final, AgentDist_init, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate2, GeneralEqmEqns_Transition, Params, DiscountFactorParamNames, AgeWeightsParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);

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
[VPath,PolicyPath]=ValueFnOnTransPath_MixHorz_PType(PricePath, ParamPath0, T_end, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);

% You can then use these to calculate the agent distribution for the transition path
AgentDistPath=AgentDistOnTransPath_MixHorz_PType(StationaryDist_init, jequaloneDist, PricePath, ParamPath0, PolicyPath, AgeWeightsParamNames,n_d,n_a,n_z,N_j,Names_i,pi_z,T_end, Params, transpathoptions, simoptions);

%% Analyse the transition path
% And then we can calculate AggVars for the path
AggVarsPath=EvalFnOnTransPath_AggVars_MixHorz_PType(FnsToEvaluate, AgentDistPath,PolicyPath, PricePath, ParamPath0, Params, T_end, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptions, simoptions);

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


function [P0,D,V,Policy,StationaryDist]=Calibrate_P0(Params_S,AgeWeightsParamNames,PTypeDistParamNames,DiscountFactorParamNames,FnsToEvaluate,ReturnFn,Policy,jequaloneDist,StationaryDist,n_d,n_a,n_z,N_j,FHorz_names,Names_i,d_grid,a_grid,z_grid,pi_z,vfoptions,simoptions)
AggVars=EvalFnOnAgentDist_AggVars_MixHorz_Case1_PType(StationaryDist,Policy, FnsToEvaluate, Params_S, n_d, n_a, n_z,N_j,Names_i,d_grid, a_grid, z_grid,simoptions);
S=sum_S_FHorz(AggVars.S, FHorz_names);
D=AggVars.D.firm.Mean;
while S>1.4
    Params_S.P0=Params_S.P0*1.2;
    [V, Policy]=ValueFnIter_MixHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid, a_grid, z_grid, pi_z,ReturnFn, Params_S, DiscountFactorParamNames,vfoptions);
    StationaryDist=StationaryDist_MixHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames, Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_S,simoptions);
    AggVars=EvalFnOnAgentDist_AggVars_MixHorz_Case1_PType(StationaryDist,Policy, FnsToEvaluate, Params_S, n_d, n_a, n_z,N_j,Names_i,d_grid, a_grid, z_grid,simoptions);
    S=sum_S_FHorz(AggVars.S, FHorz_names);
    D=AggVars.D.firm.Mean;
end
while S<0.5
    Params_S.P0=Params_S.P0*0.90;
    [V, Policy]=ValueFnIter_MixHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid, a_grid, z_grid, pi_z,ReturnFn, Params_S, DiscountFactorParamNames,vfoptions);
    StationaryDist=StationaryDist_MixHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames, Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_S,simoptions);
    AggVars=EvalFnOnAgentDist_AggVars_MixHorz_Case1_PType(StationaryDist,Policy, FnsToEvaluate, Params_S, n_d, n_a, n_z,N_j,Names_i,d_grid, a_grid, z_grid,simoptions);
    S=sum_S_FHorz(AggVars.S, FHorz_names);
    D=AggVars.D.firm.Mean;
end
if S>1.1
    P0=Params_S.P0*1.05;
elseif S<0.9
    P0=Params_S.P0*0.98;
else
    P0=Params_S.P0;
end
Params_S.P0=P0;
Params_S.D=D;
fprintf("Setting P0 to %.2f, D to %.2f \n", P0, D);

[V, Policy]=ValueFnIter_MixHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid, a_grid, z_grid, pi_z,ReturnFn, Params_S, DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_MixHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames, Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_S,simoptions);
% AggVars=EvalFnOnAgentDist_AggVars_MixHorz_Case1_PType(StationaryDist,Policy, FnsToEvaluate, Params_S, n_d, n_a, n_z,N_j,Names_i,d_grid, a_grid, z_grid,simoptions);
% S=sum_S_FHorz(AggVars.S, FHorz_names);
% D=AggVars.D.firm.Mean;

end

function S=sum_S_FHorz(AggVars_S, FHorz_names)
S=0;
for ii=1:length(FHorz_names)
    S=S+AggVars_S.(FHorz_names{ii}).Mean;
end


end