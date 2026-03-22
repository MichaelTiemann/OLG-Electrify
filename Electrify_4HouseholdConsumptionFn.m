function c=Electrify_4HouseholdConsumptionFn( ...
    labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e, ...
    pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
    kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
    r_pp,r_wedge_pp,f_htc,rentprice,energy_cpi,pv_pct_cost,energy_pct_cost,energy_pct_brown,carbon_tax)

% Implement depreciation model:
%   Car services A(t) = (1-delta_a)*A(t-1) + I(a,t)
%   Housing services H(t) = (1-delta_h)*H(t-1) + I(h,t)
%   delta_a is high depreciation; delta_h is low depreciation

% Note: experienceasset, so first inputs are (d,a,z,e,...)
% vfoptions.refine_d: only decisions d1,d3 are input to ReturnFn

carcost=0;
rentalcosts=rentprice*w*ypp;
htc=0; % house transaction cost
pvinstallcost=0;
% A Tally of energy costs, which will be deducted at the end
energy_cost_pp=0;

% Housing services (based on housing stock)
if h==0
    hs=0.5*houseservices*minhouse;
else
    hs=houseservices*h;
    rentalcosts=0;
end
% Houses start at 4x annual wage
hcost=4*h*w;
hprimecost=4*hprime*w;
% Make buying/selling a house costly/illiquid
if hprime~=h
    htc=f_htc*(hcost+hprimecost);
end

% buyhouse 2 and 4 are install/upgrade PV options
if buyhouse==2 || buyhouse==4
    if (h+hprime)==0
        % No house -> no solar
        pvinstallcost=Inf;
    elseif h==hprime
        % Pay the retrofit penalty
        pvinstallcost=1.1*pv_pct_cost*hcost;
    else % Changing house
        % PV costs approximately 5% of new house ($30K system for $600K house)
        pvinstallcost=pv_pct_cost*hprimecost;
    end
end

%% Car matters
% Car costs 50% annual wage, or can trade at 25% annual wage
if cprime==0
    if car>0
        carcost=-0.25*w; % Selling a car: get back 1/2 of what was paid for it
    end
else
    if car==0
        if cprime==1
            carcost=0.25*w; % Buying from scratch; cheap petrol car
        else
            carcost=0.5*w; % Buying from scratch; pay full price (50% of w)
        end
    elseif cprime~=car
        carcost=0.25*w; % Trading cars; pay half price (25%) with trade-in
    end
    % annual insurance, maintenance, WOF, etc.
    carcost=carcost+0.02*w*ypp;
end

% We can get P (share price) from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have
P=((1-tau_cg)*P0 + (1-tau_d)*D_pp)/(1+r_pp-tau_cg);

Plag=P; % As stationary general eqm

if agej<Jr % If working age
    %consumption = labor income + "other income" below
    c=(1-tau_l)*labor*w*kappa_j*exp(z+e)*ypp; 
else % Retirement
    c=pension*ypp;
end
% Other income: accidental share bequest + share holdings (including dividend) - dividend tax + accidental asset+house bequest + net housing assets
c=c+((1-tau_d)*D_pp+P0)*(s+AccidentBeqS_pp)+AccidentBeqS_pp+AccidentBeqAH_pp+(hcost-hprimecost);
if a<0 % In both cases, resulting `a` is added to consumption, then `aprime` subtracted
    % Subtract loan interest by adding diminishing assets
    c=c+(1+r_pp+r_wedge_pp)*a;
else
    % Deposit interest included in augmented assets
    c=c+(1+r_pp)*a;
end
% ...subtract capital gains tax and next period share, asset holdings
c=c-tau_cg*(P0-Plag)*(s+AccidentBeqS_pp)-P*sprime-aprime;
% ...subtract housing-related costs: transaction costs, rental or home maintenance costs, pv installation
c=c-htc-rentalcosts-hcost*0.02*ypp-pvinstallcost;

% ...subtract car costs (purchase, sale, and/or maintenance)
if carcost~=0
    c=c-carcost;
    % Energy costs...
    if car==1
        energy_cost_pp=energy_cost_pp+0.02*w*ypp;
    else
        if solarpv>0.5
            solarpv=solarpv-0.5;
        else
            energy_pct_cost=energy_pct_cost+0.02;
        end
    end
else
    % Public transportation cost...
    c=c-0.2*w;
end

% Add cost of housing
energy_cost_pp=energy_cost_pp+(1+energy_cpi)*energy_pct_cost*max(h^1.5,1)*ypp;

if car~=2
    % car batteries make solarpv more effective...
    solarpv=solarpv/2;
end

% PV generation: 30kW (2 solar units) meets h==1 energy needs
energy_cost_pp=energy_cost_pp+(1+energy_cpi)*energy_pct_cost*(max(h^1.5,1)-solarpv/2)*ypp;
carbon_tax_pp=energy_cost_pp*energy_pct_brown*carbon_tax/2000;

c=c-energy_cost_pp-carbon_tax_pp;


end
