function c_pp=Electrify_4HouseholdConsumptionFn( ...
    labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e, ...
    pension,AccidentBeqS,AccidentBeqAH,w,P0,D, ...
    kappa_j,tau_l,tau_d,tau_cg,S_agej_first,S_agej_peak_first,S_agej_peak_last,S_agej_last, ...
    ypp,agej,Jr,r,r_wedge,f_htc,rentprice,cpi_energy,pv_pct_cost,energy_pct_cost,energy_pct_brown,carbon_tax)

% Implement depreciation model:
%   Car services A(t) = (1-delta_a)*A(t-1) + I(a,t)
%   Housing services H(t) = (1-delta_h)*H(t-1) + I(h,t)
%   delta_a is high depreciation; delta_h is low depreciation

% Note: experienceasset, so first inputs are (d,a,z,e,...)
% vfoptions.refine_d: only decisions d1,d3 are input to ReturnFn

[sprime,aprime,s,a]=decode_sa(saprime,sa);

carcost=0;
if h==0
    rentalcosts_pp=rentprice*sqrt(kappa_j)*w*ypp;
else
    rentalcosts_pp=0;
end
htc=0; % house transaction cost
pvinstallcost=0;

% Houses start at 4x annual wage
hcost=4*h*w;
hprimecost=4*hprime*w;
% Make buying/selling a house costly/illiquid
if hprime~=h
    htc=f_htc*(hcost+hprimecost);
end

% buyhouse 3 and 4 are install/upgrade PV options
if buyhouse==3 || buyhouse==4
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
    carcost_pp=0;
    if car~=0
        if car==1 % Selling a car: get back <= 1/2 of what was paid for it
            carcost=-0.2*w;
        elseif car==2
            carcost=-0.4*w;
        end
        carcost_pp=carcost+0.02*w*ypp;
    end
else
    if car==0
        if cprime==1
            carcost=0.4*w; % Buying from scratch; cheap petrol car
        else
            carcost=0.8*w; % Buying from scratch; pay full price (50% of w)
        end
    elseif car<cprime
        carcost=0.5*w; % Minuscule trade-in value of petrol car
    else
        carcost=0.0*w; % Minimal extra cost to move backwards
    end
    % annual insurance, maintenance, WOF, etc.
    carcost_pp=carcost+0.02*w*ypp;
end

if sprime>=s
    cg=0; % We are holding or buying, so no capital gains
else
    if agej<=S_agej_peak_first
        Plag=P0*(1-2*r)^ypp; % Dispose of shares presumably acquired recently
    elseif agej<S_agej_peak_last
        % We have been holding since peak acquisition
        agej_bought=S_agej_peak_first;
        Plag=P0*(1-2*r)^(ypp*(agej-agej_bought));
    elseif S_agej_peak_last==S_agej_last % Bulk liquidation
        % Sell all remaining shares from first acquisition to buy-point (using geometric mean to average acquisition cost)
        agej_bought=S_agej_peak_first-sqrt(S_agej_peak_first-S_agej_first);
        Plag=P0*(1-2*r)^(ypp*(agej-agej_bought));
    else
        % Estimate where we are past peak accumulation and mirror around to
        % proportional acquisition point
        agej_selling_pct=(agej-S_agej_peak_last)/(S_agej_last-S_agej_peak_last);
        agej_bought=S_agej_peak_first-agej_selling_pct*(S_agej_peak_first-S_agej_first);
        Plag=P0*(1-2*r)^(ypp*(agej-agej_bought));
    end
    cg=tau_cg*(P0-Plag)*(s+AccidentBeqS-sprime);
end

if agej<Jr % If working age
    %consumption = labor income + "other income" below
    c_pp=(1-tau_l)*labor*w*kappa_j*exp(z+e)*ypp; 
else % Retirement
    c_pp=pension*ypp;
end

% Other income: accidental share bequest + share holdings (including dividend) - dividend tax + accidental asset+house bequest + net housing assets
c_pp=c_pp+((1-tau_d)*D*ypp+P0)*(s+AccidentBeqS)+AccidentBeqAH+(hcost-hprimecost);
if a<0 % In both cases, resulting `a` is added to consumption, then `aprime` subtracted
    % Subtract loan interest by adding diminishing assets
    c_pp=c_pp+(1+r+r_wedge)^ypp*a;
else
    % Deposit interest included in augmented assets
    c_pp=c_pp+(1+r)^ypp*a;
end
% ...subtract capital gains, next period share, asset holdings
c_pp=c_pp-cg-P0*sprime-aprime;
% ...subtract housing-related costs: transaction costs, rental or home maintenance costs, pv installation
c_pp=c_pp-htc-rentalcosts_pp-hcost*0.01*ypp-pvinstallcost;

% ...subtract car costs (purchase, sale, and/or maintenance)
if carcost_pp~=0
    c_pp=c_pp-carcost_pp;
else
    % Public transportation cost...
    c_pp=c_pp-0.1*w*ypp;
end

% Add energy cost of housing, less PV generation: 30kW (2 solar units) meets h==1 energy needs
energy_cost_pp=Electrify_4HouseholdEnergyCosts(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e,w,ypp,cpi_energy,energy_pct_cost);
carbon_tax_pp=energy_cost_pp*energy_pct_brown*carbon_tax*ypp/200;

c_pp=c_pp-energy_cost_pp-carbon_tax_pp;


end

function [sprime,aprime,s,a]=decode_sa(saprime,sa)

if saprime<1
    sprime=0;
    aprime=saprime;
else
    sprime=floor(saprime);
    aprime=rem(saprime,1);
end

if sa<1
    s=0;
    a=sa;
else
    s=floor(sa);
    a=rem(sa,1);
end


end
