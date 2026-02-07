%% Rooftop Solar Life-Cycle Model (based on Life-Cycle Model 35: Portfolio-Choice with Housing)
% Modify Life-Cycle Model 35, with solarpv as an experienceasset:

% In terms of code, using 'riskyasset' alongside a standard asset means:
% The return function and functions to evaluate have first inputs (d,hprime,h,a,z,...) [no aprime like in Case1]
% Notice that we have hprime and h (for the standard asset), but only a for
% the risky asset (no aprime).
% We need to define aprime(d,u)

% There is not much agreement on how to handle mortality risk with Epstein-Zin preferences
% We can treat them as a risk
% vfoptions.survivalprobability='sj';
% DiscountFactorParamNames={'beta'};
% Or we could just treat them as another discount factor
DiscountFactorParamNames={'beta','sj'};


%% How does VFI Toolkit think about this?
%
% Two decision variables: labor and solarpv
% Three endogenous state variables: a, h, and solarpv (assets, housing, and solarpv)
% One stochastic exogenous state variable: z, an AR(1) process (in logs), idiosyncratic shock to labor productivity units
% One between-period i.i.d. variable: u, the return to the risky asset
% Age: j

%% Begin setting up to use VFI Toolkit to solve
% Lets model agents from age 20 to age 100, so 81 periods

Params.agejshifter=19; % Age 20 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=100-Params.agejshifter; % =81, Number of period in life-cycle

% Grid sizes to use
n_d=[51,2]; % Decisions: labor, install solarpv
n_a=[31,5,6]; % Endogenous asset, housing, and solarpv assets (0-40kW generation)
n_z=7; % Exogenous labor productivity units shock
% n_u=5; % Between period i.i.d. shock
N_j=Params.J; % Number of periods in finite horizon

% When there is more than one endogenous state, the experienceasset is the last one


% Specify Epstein-Zin preferences
% vfoptions.exoticpreferences='EpsteinZin';
% vfoptions.EZpositiveutility=0; % Epstein-Zin preferences in utility-units have to be handled differently depending on whether the utility funciton is positive or negative valued (this is all done internally, you just need to use vfoptions to specify which)
% vfoptions.EZriskaversion='phi'; % additional risk-aversion
% Params.phi is set below

%% To speed up the use of experienceasset we use 'refine_d', which requires us to set the decision variables in a specific order
vfoptions.refine_d=[1,0,1]; % tell the code how many d1, d2, and d3 there are
% Idea is to distinguish three categories of decision variable:
%  d1: decision is in the ReturnFn but not in aprimeFn
%  d2: decision is in the aprimeFn but not in ReturnFn
%  d3: decision is in both ReturnFn and in aprimeFn
% Note: ReturnFn must use inputs (d1,d3,..) 
%       aprimeFn must use inputs (d2,d3,..)
% n_d must be set up as n_d=[n_d1, n_d2, n_d3]
% d_grid must be set up as d_grid=[d1_grid; d2_grid; d3_grid];
% It is possible to solve models without any d2, as is the case here.
simoptions.refine_d=vfoptions.refine_d;

vfoptions.gridinterplayer=0;
vfoptions.ngridinterp=5;
simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;

%% Parameters

% Housing
Params.f_htc=0.05; % transaction cost of buying/selling house (is a percent of h+hprime)
% Params.minhouse % set below, is the minimum value of house that can be purchased
Params.rentprice=0.3; % I figured setting rent a decent fraction of income is sensible
Params.f_coll=0.5; % collateral contraint (fraction of house value that can be borrowed)
Params.houseservices=0.5; % housing services as a fraction of house value
Params.pv_pct_cost=0.05; % modeling a $30K install for a $600K house
Params.energy_pct_cost=0.07; % Electricity: 3%; Gas: 2%; Petrol: 2%

% Discount rate
Params.beta = 0.96;
% Preferences
Params.sigma=2; % Coeff of relative risk aversion (curvature of consumption)
% Params.phi=10; % Additional risk aversion (from Epstein-Zin preferences)
Params.sigma_h=0.25; % Relative importance of housing services (vs consumption) in utility
Params.eta = 1.5; % Curvature of leisure (This will end up being 1/Frisch elasty)
Params.psi = 0.5; % Weight on leisure

% Prices
Params.w=1; % Wage

% Asset returns
Params.r=0.05; % Rate of return on risk free asset
Params.r_wedge=0.02; % Rate difference between deposits and loans
% u is the stochastic component of the excess returns to the risky asset
% Params.rp=0.03; % Mean excess returns to the risky asset (so the mean return of the risky asset will be r+rp)
% Params.sigma_u=0.025; % Standard deviation of innovations to the risky asset
% Params.rho_u=0; % Asset return risk component is modeled as iid (if you regresse, e.g., the percent change in S&P500 on it's one year lag you get a coefficient of essentially zero)
% [u_grid, pi_u]=discretizeAR1_FarmerToda(Params.rp,Params.rho_u,Params.sigma_u,n_u);
% pi_u=pi_u(1,:)'; % This is iid

% Demographics
Params.agej=1:1:Params.J; % Is a vector of all the agej: 1,2,3,...,J
Params.Jr=65-Params.agejshifter;

% Pensions
Params.pension=0.4; % Increased to be greater than rental costs

% Age-dependent labor productivity units
Params.kappa_j=[linspace(0.5,2,Params.Jr-15),linspace(2,1,14),zeros(1,Params.J-Params.Jr+1)];
% Exogenous shock process: AR1 on labor productivity units
Params.rho_z=0.9;
Params.sigma_epsilon_z=0.03;

% Conditional survival probabilities: sj is the probability of surviving to be age j+1, given alive at age j
% Most countries have calculations of these (as they are used by the government departments that oversee pensions)
% In fact I will here get data on the conditional death probabilities, and then survival is just 1-death.
% Here I just use them for the US, taken from "National Vital Statistics Report, volume 58, number 10, March 2010."
% I took them from first column (qx) of Table 1 (Total Population)
% Conditional death probabilities
Params.dj=[0.006879, 0.000463, 0.000307, 0.000220, 0.000184, 0.000172, 0.000160, 0.000149, 0.000133, 0.000114, 0.000100, 0.000105, 0.000143, 0.000221, 0.000329, 0.000449, 0.000563, 0.000667, 0.000753, 0.000823,...
    0.000894, 0.000962, 0.001005, 0.001016, 0.001003, 0.000983, 0.000967, 0.000960, 0.000970, 0.000994, 0.001027, 0.001065, 0.001115, 0.001154, 0.001209, 0.001271, 0.001351, 0.001460, 0.001603, 0.001769, 0.001943, 0.002120, 0.002311, 0.002520, 0.002747, 0.002989, 0.003242, 0.003512, 0.003803, 0.004118, 0.004464, 0.004837, 0.005217, 0.005591, 0.005963, 0.006346, 0.006768, 0.007261, 0.007866, 0.008596, 0.009473, 0.010450, 0.011456, 0.012407, 0.013320, 0.014299, 0.015323,...
    0.016558, 0.018029, 0.019723, 0.021607, 0.023723, 0.026143, 0.028892, 0.031988, 0.035476, 0.039238, 0.043382, 0.047941, 0.052953, 0.058457, 0.064494,...
    0.071107, 0.078342, 0.086244, 0.094861, 0.104242, 0.114432, 0.125479, 0.137427, 0.150317, 0.164187, 0.179066, 0.194979, 0.211941, 0.229957, 0.249020, 0.269112, 0.290198, 0.312231, 1.000000]; 
% dj covers Ages 0 to 100
Params.sj=1-Params.dj(21:101); % Conditional survival probabilities
Params.sj(end)=0; % In the present model the last period (j=J) value of sj is actually irrelevant

%% Grids

% Grid for labour choice
labor_grid=linspace(0,1,n_d(1))'; % Notice that it is imposing the 0<=h<=1 condition implicitly

% The ^3 means that there are more points near 0 and near 10. We know from
% theory that the value function will be more 'curved' near zero assets,
% and putting more points near curvature (where the derivative changes the most) increases accuracy of results.
asset_grid=-3+13*(linspace(0,1,n_a(1)))'; % Note, I use equal spacing (normally would put most points near zero)
% note: will go from -3 to 13-3
% Make it so that there is a zero assets
% Find closest to zero assets
[~,zeroassetindex]=min(abs(asset_grid));
asset_grid(zeroassetindex)=0;

% age20avgincome=Params.w*Params.kappa_j(1);
% house_grid=[0; logspace(2*age20avgincome, 12*age20avgincome, 5)'];
house_grid=(0:1:n_a(2)-1)';
% Note, we can see from w*kappa_j*z and the values of these, that average
% income is going to be around one, so will just use this simpler house grid
% [We can think about the values of the house_grid as being relative the average income (or specifically average at a given age)]
Params.minhouse=house_grid(2); % first is zero (no house)

% kW of solar generation installed, 10kW per grid element
solarpv_grid=10*(0:1:n_a(3)-1)';

% First, the AR(1) process z
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,Params.rho_z,Params.sigma_epsilon_z,n_z);
z_grid=exp(z_grid); % Take exponential of the grid
[mean_z,~,~,~]=MarkovChainMoments(z_grid,pi_z); % Calculate the mean of the grid so as can normalise it
z_grid=z_grid./mean_z; % Normalise the grid on z (so that the mean of z is exactly 1)

% Share of assets invested in the risky asset
% riskyshare_grid=linspace(0,1,n_d(1))'; % Share of assets, from 0 to 1

% installpv is a binary choice
installpv_grid=[0; 1];

% Set up d for VFI Toolkit (is the two decision variables)
d_grid=[labor_grid; installpv_grid];

a_grid=[asset_grid; house_grid; solarpv_grid];

%% Solar PV is an experience asset
vfoptions.experienceasset=1;
vfoptions.lowmemory=0;
simoptions.experienceasset=1;

%% Define aprime function used for the experience asset

% experienceasset: aprime_val=aprimeFn(d,a)
% vfoptions.refine_d: the decision variables input to aprimeFn are d1,d3
aprimeFn=@(installpv, solarpv) ElectrifyHousing_aprimeFn(installpv, solarpv); % Will return the value of aprime (solarpv)

%% Put the experience asset into vfoptions and simoptions
vfoptions.aprimeFn=aprimeFn;
% vfoptions.n_u=n_u;
% vfoptions.u_grid=u_grid;
% vfoptions.pi_u=pi_u;
simoptions.aprimeFn=aprimeFn;
% simoptions.n_u=n_u;
% simoptions.u_grid=u_grid;
% simoptions.pi_u=pi_u;
% Because a_grid and d_grid are involved in experience assets, but are not
% normally needed for agent distriubiton simulation, we have to also
% include these in simoptions
simoptions.a_grid=a_grid;
simoptions.d_grid=d_grid;

%% Now, create the return function 
% DiscountFactorParamNames={'beta','sj'};

% Use 'ElectrifyHousing_ReturnFn'
ReturnFn=@(labor,installpv,aprime,hprime,a,h,solarpv,z, ...
        r,r_wedge,w,sigma,agej,Jr,pension,kappa_j,psi,eta,sigma_h,f_htc,minhouse,rentprice,f_coll,houseservices,pv_pct_cost,energy_pct_cost) ...
    ElectrifyHousing_ReturnFn(labor,installpv,aprime,hprime,a,h,solarpv,z, ...
        r,r_wedge,w,sigma,agej,Jr,pension,kappa_j,psi,eta,sigma_h,f_htc,minhouse,rentprice,f_coll,houseservices,pv_pct_cost,energy_pct_cost)
% vfoptions.refine_d: only (d1,d3,..) are input to ReturnFn

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium
disp('Test ValueFnIter')
vfoptions.verbose=1;
tic;
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
toc

% V is now (a,z,j). This was already true, just that previously z was trivial (a single point) 
% Compare
size(V)
% with
% [n_a,n_z,N_j]
% there are the same.
% Policy is
size(Policy)
% which is the same as
[length(n_d)+1+1,n_a,n_z,N_j]
% The n_a,n_z,N_j represent the state on which the decisions/policys
% depend, and there is one decision for each decision variable 'd' plus one
% more for the experience asset and one more for the standard asset


%% Let's take a quick look at what we have calculated, namely V and Policy

% Default asset index value is at 33rd %-ile
a_ii=floor(n_a(1)/3);
% Default house index value is h=1
h_ii=2;
% Default SolarPV index value is at 10kW
pv_ii=2;
% Ages to plot
agej_to_plot = [1 20 45 46 81];

% The value function V depends on the state, so now it depends on both asset holdings and age.

% We can plot V as a 3d plot (surf is matlab command for 3d plot)
% Which z value should we plot? I will plot the median
zind=floor(n_z+1)/2; % This will be the median
figure(1)
subplot(2,1,1); surf(asset_grid*ones(1,Params.J),ones(n_a(1),1)*(1:1:Params.J),exp(reshape(V(:,h_ii,pv_ii,zind,:),[n_a(1),Params.J])))
title('Value function: median value of exp(z)')
xlabel(sprintf('Assets (a), h=%d, solarpv=%d', house_grid(h_ii), solarpv_grid(pv_ii)))
ylabel('Age j')
subplot(2,1,2); surf(asset_grid*ones(1,Params.J),ones(n_a(1),1)*(Params.agejshifter+(1:1:Params.J)),exp(reshape(V(:,h_ii,pv_ii,zind,:),[n_a(1),Params.J])))
title('Value function: median value of exp(z)')
xlabel(sprintf('Assets (a), h=%d, solarpv=%d', house_grid(h_ii), solarpv_grid(pv_ii)))
ylabel('Age in Years')

% Do another plot of V, this time as a function (of assets) for a given age (I do a few for different ages)
figure(2)
for jj=1:5
    subplot(5,3,jj*3-2);
    plot(asset_grid,V(:,h_ii,pv_ii,1,agej_to_plot(jj)), ...
        asset_grid,V(:,h_ii,pv_ii,zind,agej_to_plot(jj)), ...
        asset_grid,V(:,h_ii,pv_ii,end,agej_to_plot(jj)))
    if agej_to_plot(jj)==Params.Jr
        title('Value fn at age j=46 (first year of retirement)')
    else
        title(sprintf('Value fn at age j=%d', agej_to_plot(jj)))
        if jj==1
            legend('min z','median z','max z') % Just include the legend once in the top subplot
        end
    end
end
xlabel(sprintf('Assets (a), h=%d, solarpv=%d kW', house_grid(h_ii), solarpv_grid(pv_ii)))

for jj=1:5
    subplot(5,3,jj*3-1);
    plot(house_grid,V(a_ii,:,pv_ii,1,agej_to_plot(jj)), ...
        house_grid,V(a_ii,:,pv_ii,zind,agej_to_plot(jj)), ...
        house_grid,V(a_ii,:,pv_ii,end,agej_to_plot(jj)))
    if agej_to_plot(jj)==Params.Jr
        title('Value fn at age j=46 (first year of retirement)')
    else
        title(sprintf('Value fn at age j=%d', agej_to_plot(jj)))
        if jj==1
            legend('min z','median z','max z') % Just include the legend once in the top subplot
        end
    end
end
xlabel(sprintf('Housing (h), a=%.2f, solarpv=%d kW', asset_grid(a_ii), solarpv_grid(pv_ii)))

for jj=1:5
    subplot(5,3,jj*3);
    plot(solarpv_grid,squeeze(V(a_ii,pv_ii,:,1,agej_to_plot(jj))), ...
        solarpv_grid,squeeze(V(a_ii,pv_ii,:,zind,agej_to_plot(jj))), ...
        solarpv_grid,squeeze(V(a_ii,pv_ii,:,end,agej_to_plot(jj))))
    if agej_to_plot(jj)==Params.Jr
        title('Value fn at age j=46 (first year of retirement)')
    else
        title(sprintf('Value fn at age j=%d', agej_to_plot(jj)))
        if jj==1
            legend('min z','median z','max z') % Just include the legend once in the top subplot
        end
    end
end
xlabel(sprintf('Solar PV (solarpv), a=%.2f, h=%d', asset_grid(a_ii), house_grid(h_ii)))

% Convert the policy function to values (rather than indexes).
% Policy(1,:,:,:) is labor, Policy(2,:,:,:) is installpv [as function of (a,z,j)]
% Policy(3,:,:,:) is aprime (assets), Policy(4,:,:,:) is hprime (housing)
% Plot both as a 3d plot, again I arbitrarily choose the median value of z
figure(3)
PolicyVals=PolicyInd2Val_Case1_FHorz(Policy,n_d,n_a,n_z,N_j,d_grid,a_grid,simoptions);
subplot(3,1,1); surf(asset_grid*ones(1,Params.J),ones(n_a(1),1)*(1:1:Params.J),reshape(PolicyVals(1,:,h_ii,pv_ii,zind,:),[n_a(1),Params.J]))
title(sprintf('Policy function: labor, median z, h=%d, solarpv=%d kW', house_grid(h_ii), solarpv_grid(pv_ii)))
xlabel('Assets (a)')
ylabel('Age j')
zlabel('Labor')
subplot(3,1,2); surf(house_grid*ones(1,Params.J),ones(n_a(2),1)*(1:1:Params.J),reshape(PolicyVals(4,a_ii,:,pv_ii,zind,:),[n_a(2),Params.J]))
title(sprintf('Policy function: hprime, median z, a=%.2f, solarpv=10kW', asset_grid(a_ii), solarpv_grid(pv_ii)))
xlabel('Housing (h)')
ylabel('Age j')
zlabel('Housing')
subplot(3,1,3); surf(solarpv_grid*ones(1,Params.J),ones(n_a(3),1)*(1:1:Params.J),reshape(PolicyVals(2,a_ii,h_ii,:,zind,:),[n_a(3),Params.J]))
title(sprintf('Policy function: solarpv, median z, a=%.2f, h=%d', asset_grid(a_ii), house_grid(h_ii)))
xlabel('Solar PV (solarpv)')
ylabel('Age j')
zlabel('Solar PV Installed')

for ii=1:5
    % % Again, plot both policies (h and aprime), this time as a function (of assets) for a given age  (I do a few for different ages)
    figure(3+ii)
    a_ii=ii*floor(n_a(1)/5);
    h_ii=ii;

    agej_to_plot = [1 20 45 46 81];
    for jj=1:5
        subplot(5,2,jj*2-1);
        plot(asset_grid,squeeze(PolicyVals(1,:,h_ii,pv_ii,1,agej_to_plot(jj))), ...
            asset_grid,squeeze(PolicyVals(1,:,h_ii,pv_ii,zind,agej_to_plot(jj))), ...
            asset_grid,squeeze(PolicyVals(1,:,h_ii,pv_ii,end,agej_to_plot(jj))))
        if agej_to_plot(jj)==Params.Jr
            title('Policy for labor at age j=46 (first year of retirement)')
        else
            title(sprintf('Policy for labor at age j=%d', agej_to_plot(jj)))
        end
    end
    xlabel(sprintf('Assets (a), h=%d, solarpv=%d kW', house_grid(h_ii), solarpv_grid(pv_ii)))

    for jj=1:5
        subplot(5,2,jj*2);
        plot(solarpv_grid,squeeze(PolicyVals(2,a_ii,h_ii,:,1,agej_to_plot(jj))), ...
            solarpv_grid,squeeze(PolicyVals(2,a_ii,h_ii,:,zind,agej_to_plot(jj))), ...
            solarpv_grid,squeeze(PolicyVals(2,a_ii,h_ii,:,end,agej_to_plot(jj))))
        if agej_to_plot(jj)==Params.Jr
            title('Policy for labor at age j=46 (first year of retirement)')
        else
            if jj==1
                legend('min z','median z','max z') % Just include the legend once in the top-right subplot
            end
            title(sprintf('Policy for labor at age j=%d', agej_to_plot(jj)))
        end
    end
    xlabel(sprintf('Solar PV (solarpv), a=%.2f, h=%d', asset_grid(a_ii), house_grid(h_ii)))
end

%% Now, we want to graph Life-Cycle Profiles

%% Initial distribution of agents at birth (j=1)
% Before we plot the life-cycle profiles we have to define how agents are
% at age j=1. We will give them all zero assets.
jequaloneDist=zeros([n_a,n_z],'gpuArray'); % Put no households anywhere on grid
jequaloneDist(zeroassetindex,1,1,floor((n_z+1)/2))=1; % All agents start with zero assets, no house, no solarpv, and the median shock

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
% Like with return function, we have to include (h,aprime,a,z) as first
% inputs, then just any relevant parameters.
% FnsToEvaluate.riskyshare=@(savings,riskyshare,hprime,h,a,z) riskyshare; % riskyshare, is the fraction of savings invested in the risky asset
FnsToEvaluate.earnings=@(labor,installpv,aprime,hprime,a,h,solarpv,z,w,kappa_j) labor*w*kappa_j*z; % labor earnings
FnsToEvaluate.assets=@(labor,installpv,aprime,hprime,a,h,solarpv,z) a; % a is the current asset holdings
FnsToEvaluate.housing=@(labor,installpv,aprime,hprime,a,h,solarpv,z) h; % h is the current house holdings
FnsToEvaluate.solarpv=@(labor,installpv,aprime,hprime,a,h,solarpv,z) solarpv; % solarpv is the current solarpv holdings

% notice that we have called these riskyshare, earnings and assets

%% Calculate the life-cycle profiles
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);

% For example
% AgeConditionalStats.earnings.Mean
% There are things other than Mean, but in our current deterministic model
% in which all agents are born identical the rest are meaningless.

%% Plot the life cycle profiles of fraction-of-time-worked, earnings, and assets

figure(10)
% subplot(4,1,1); plot(1:1:Params.J,AgeConditionalStats.riskyshare.Mean)
% title('Life Cycle Profile: Share of savings invested in risky asset (riskyshare)')
subplot(4,1,1); plot(1:1:Params.J,AgeConditionalStats.earnings.Mean)
title('Life Cycle Profile: Labor Earnings (w kappa_j z)')
subplot(4,1,2); plot(1:1:Params.J,AgeConditionalStats.assets.Mean)
title('Life Cycle Profile: Assets (a)')
subplot(4,1,3); plot(1:1:Params.J,AgeConditionalStats.housing.Mean)
title('Life Cycle Profile: Housing (h)')
subplot(4,1,4); plot(1:1:Params.J,AgeConditionalStats.solarpv.Mean)
title('Life Cycle Profile: SolarPV (solarpv)')