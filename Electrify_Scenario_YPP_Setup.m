function [Params]=Electrify_Scenario_YPP_Setup(Params,scenario,ypp,small_z_no_e,max_age,agejshifter,r,r_wedge,beta,n,k_j1,k_j2,k_j2_length,k_j3,sigma_h,sigma_c,psi,tau_cg,energy_pct_cost,G,D,AccidentBeqS,AccidentBeqAH)

if small_z_no_e
    Params.e=0;
end

Params.scenario=scenario;
Params.ypp=ypp;

Params.J=ceil((max_age-agejshifter)/ypp); % =60/ypp, Number of period in life-cycle
Params.Jr=round((65-agejshifter)/ypp); % Age 65 (period 10 is ages 65-69 in the 5 year case)
Params.agej=1:1:Params.J; % Is a vector of all the periods: 1,2,3,...,J

if Params.Jr>5
    kappa_j12=linspace(k_j1(scenario),k_j2(scenario),Params.Jr-round((15+k_j2_length(scenario))/ypp));
    kappa_j2s=k_j2(scenario)*ones(1,ceil(k_j2_length(scenario)/ypp));
    kappa_j23=linspace(k_j2(scenario),k_j3(scenario),ceil(14/ypp));
else
    kappa_j12=linspace(k_j1(scenario),k_j2(scenario),Params.Jr-1-min(k_j2_length(scenario),1));
    kappa_j2s=k_j2(scenario)*ones(1,min(k_j2_length(scenario),1)); % At most one period of max wage
    kappa_j23=k_j3(scenario)*ones(1,1); % One period of "pre-retirement" work
end
kappa_jr=zeros(1,Params.J-Params.Jr+1);
kappa_j=[kappa_j12, kappa_j2s, kappa_j23, kappa_jr];

% If Params.J is rounded up, don't add extra zeros
Params.kappa_j=kappa_j(1:Params.J);

Params.r_wedge=r_wedge;
Params.beta = beta(scenario);
Params.n=n(scenario); % percentage rate (expressed as fraction) of population growth per period

Params.carservices_j=0.1*ones(Params.J,1);
% Cars start to be useful as people ramp up family life
age1=ceil((28-agejshifter)/ypp);
age2=ceil((32-agejshifter)/ypp);
Params.carservices_j(age1:age2)=linspace(0.5,2,age2-age1+1);
age1=age2; age2=ceil((44-agejshifter)/ypp);
Params.carservices_j(age1:age2)=2*ones(age2-age1+1,1);
age1=age2; age2=ceil((65-agejshifter)/ypp);
Params.carservices_j(age1:age2)=linspace(2,1,age2-age1+1);
age1=age2; age2=ceil((80-agejshifter)/ypp);
Params.carservices_j(age1:age2)=linspace(1,0,age2-age1+1);
Params.carservices_j(age2+1:end)=0;

% Life-cycle AR(1) process z, on (log) labor productivity units
% Chosen following Karahan & Ozkan (2013) [as used by Fella, Gallipoli & Pan (2019)]
% Note that 37 covers 24 to 60 inclusive (as in the original)
% Now repeat the first and last values to fill in working age, and put zeros for retirement (where it is anyway irrelevant)
ones_pp4y=ones(1,ceil(4/ypp));
rho_z=0.7596+0.2039*((1:ypp:37)/10)-0.0535*((1:ypp:37)/10).^2+0.0028*((1:ypp:37)/10).^3; % Chosen following Karahan & Ozkan (2013) [as used by Fella, Gallipoli & Pan (2019)]
sigma_epsilon_z=0.0518-0.0405*((1:ypp:37)/10)+0.0105*((1:ypp:37)/10).^2-0.0002*((1:ypp:37)/10).^3; % Chosen following Karahan & Ozkan (2013) [as used by Fella, Gallipoli & Pan (2019)]

% Here we allow one period each at the start and end of working age, followed by retirement
Params.rho_z=[rho_z(1)*ones_pp4y, ...
    rho_z, ...
    rho_z(end)*ones_pp4y, ...
    zeros(1,Params.J-Params.Jr+1)];
Params.rho_z=Params.rho_z(1:Params.J);
Params.sigma_epsilon_z=[sigma_epsilon_z(1)*ones_pp4y, ...
    sigma_epsilon_z, ...
    sigma_epsilon_z(end)*ones_pp4y, ...
    sigma_epsilon_z(end)*ones(1,Params.J-Params.Jr+1)];
Params.sigma_epsilon_z=Params.sigma_epsilon_z(1:Params.J);

% Transitory iid shock
sigma_e=0.0410+0.0221*((24:ypp:60)/10)-0.0069*((24:ypp:60)/10).^2+0.0008*((24:ypp:60)/10).^3;
Params.sigma_e=[sigma_e(1)*ones_pp4y, ...
    sigma_e, ...
    sigma_e(end)*ones_pp4y, ...
    sigma_e(end)*ones(1,Params.J-Params.Jr+1)];
Params.sigma_e=Params.sigma_e(1:Params.J);

% Note: These iid shocks will interact with the endogenous labor so the final labor
% earnings process will not equal that of Karahan & Ozkan (2013)
% Note: Karahan & Ozkan (2013) also have a fixed effect (which they call alpha) and which I ignore here.

% Conditional survival probabilities: sj is the probability of surviving to be age j+1, given alive at age j
% Most countries have calculations of these (as they are used by the government departments that oversee pensions)
% In fact I will here get data on the conditional death probabilities, and then survival is just 1-death.
% Here I just use them for the US, taken from "National Vital Statistics Report, volume 58, number 10, March 2010."
% I took them from first column (qx) of Table 1 (Total Population)
% Conditional death probabilities
dj=[0.006879, 0.000463, 0.000307, 0.000220, 0.000184, 0.000172, 0.000160, 0.000149, 0.000133, 0.000114, 0.000100, 0.000105, 0.000143, 0.000221, 0.000329, 0.000449, 0.000563, 0.000667, 0.000753, 0.000823,... % Ages 1-20
    0.000894, 0.000962, 0.001005, 0.001016, 0.001003, 0.000983, 0.000967, 0.000960, 0.000970, 0.000994, 0.001027, 0.001065, 0.001115, 0.001154, 0.001209, 0.001271, 0.001351, 0.001460, 0.001603, 0.001769,... % Ages 21-40
    0.001943, 0.002120, 0.002311, 0.002520, 0.002747, 0.002989, 0.003242, 0.003512, 0.003803, 0.004118, 0.004464, 0.004837, 0.005217, 0.005591, 0.005963, 0.006346, 0.006768, 0.007261, 0.007866, 0.008596,... % Ages 41-60
    0.009473, 0.010450, 0.011456, 0.012407, 0.013320, 0.014299, 0.015323,...                                                                                                                                   % Ages 61-67
    0.016558, 0.018029, 0.019723, 0.021607, 0.023723, 0.026143, 0.028892, 0.031988, 0.035476, 0.039238, 0.043382, 0.047941, 0.052953, 0.058457, 0.064494,...                                                   % Ages 68-82
    0.071107, 0.078342, 0.086244, 0.094861, 0.104242, 0.114432, 0.125479, 0.137427, 0.150317, 0.164187, 0.179066, 0.194979, 0.211941, 0.229957, 0.249020, 0.269112, 0.290198, 0.312231, 1.000000];             % Ages 83-101
dj=resize(dj,101+ypp,FillValue=1);
% dj covers Ages 0-100, plus extras at end to make it period-friendly
% Note: when ypp==1, the product over the reshaped array is over a single year period (i.e. trivial)
sj_init=prod(1-reshape(dj(1:ypp*Params.J),[ypp,Params.J]),1); % p5-year survival rates
sj_init(end)=0; % In the present model the last period (j=J) value of sj is actually irrelevant
Params.sj_init=sj_init;

% Add 5 years of life expectancy...age sj(65) in the future will be sj(60) by today's statistics
% Part of this is achieved by improving early childhood survival as well...feeding two birds with one worm
sj_final=prod(1-reshape([dj(1:2:10), repelem(dj(11:15), 3), dj(16:ypp*Params.J-5)],[ypp,Params.J]),1);
sj_final(end)=0; % In the present model the last period (j=J) value of sj is actually irrelevant
Params.sj_final=sj_final;

%% Setup for sj and mewj transitions (T-by-N_j)
% We defer doing transition maths until we calculate GE final
Params.sj=sj_init;
Params.mewj=cumprod([1,Params.sj(1:end-1)],2); % mass of age jj is the mass of jj-1 that survive
Params.mewj=Params.mewj./((1+Params.n).^(ypp*(1:Params.J)-1)); % Population shrinks in the N_j dimension
Params.mewj=Params.mewj./sum(Params.mewj); % normalize age-masses to sum to one

% Note: This is rather incomplete, as really you should also have
% population growth rate n. But this does not change any thing in terms of
% the 'objects to compute'. Instead you need to renormalize the model for
% the population growth, and this just means you get a 'n' appearing in
% some equations below. But other than 'n' in some equations, the way you
% do this with the toolkit does not change.

% The warmglow parameters will help us find the GE solution to actual bequest rates/values
Params.AccidentBeqS=AccidentBeqS(scenario);
if scenario>2
    Params.AccidentBeqAH=AccidentBeqAH(scenario);
end
Params.r=r;
Params.firmbeta=1/(1+r/(1-tau_cg)); % 1/(1+r) but returns net of capital gains tax
Params.energybeta=1/(1+r/(1-tau_cg)); % 1/(1+r) but returns net of capital gains tax

Params.sigma_h=sigma_h(scenario);
Params.sigma_c=sigma_c(scenario);
Params.psi=psi(scenario);
Params.energy_pct_cost=energy_pct_cost(scenario);

if scenario==4
    Params.Ek=1; Params.ek=1;
    Params.carbon_tax=35;
    Params.energy_pct_brown=0.8;
end

Params.delta=0.054; % Depreciation of physical capital per period
Params.G=G; % Government expenditure
Params.D=D; % The dividends paid by the firm per period
if scenario>=3
    Params.pv_delta=0.02; % Depreciation of PVs per period
end

end