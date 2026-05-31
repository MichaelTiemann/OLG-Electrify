function F=Electrify_4HouseholdReturnFn( ...
    labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e, ...
    pension,max_benefit,AccidentBeqS,AccidentBeqAH,w,P0,D,sigma,psi,eta,sigma_h,sigma_c,kappa_j,warmglow1,warmglow2,tau_l,tau_d,tau_cg,S_agej_first,S_agej_peak_first,S_agej_peak_last,S_agej_last, ...
    ypp,agej,Jr,J,r,r_wedge,f_htc,minhouse,rentprice,f_coll,houseservices,carservices_j,cpi_energy,pv_pct_cost,energy_pct_cost,energy_pct_brown,carbon_tax ...
    )
% Implement depreciation model:
%   Car services A(t) = (1-delta_a)*A(t-1) + I(a,t)
%   Housing services H(t) = (1-delta_h)*H(t-1) + I(h,t)
%   delta_a is high depreciation; delta_h is low depreciation

% Note: experienceasset, so first inputs are (d,a,z,e,...)
% vfoptions.refine_d: only decisions d1,d3 are input to ReturnFn

F=-Inf;

%% Housing matters
% buyhouse decisions--needed for solarpv experience asset
%  0=no house/sell house
%  1=buy house w/o pv this period
%  2=keep house; no pv upgrade
%  3=buy house w/ pv this period
%  4=keep house; pv upgrade (if possible)
if buyhouse==0
    if hprime~=0
        % Forbid owning house when buyhouse=0
        return
    end
elseif mod(buyhouse,2)==0
    if hprime==0 || hprime~=h
        % Forbid selling/changing house we say we are keeping
        return
    end
end

% Houses start at 4x household wage ($114K across 2M NZ households)
[sprime,aprime,s,a]=decode_sa(saprime,sa);
hcost=4*h*w;
hprimecost=4*hprime*w;

%% Allow/Disallow some trivial agent decisions
if (sprime-s>0 && aprime+hcost<0 ...                  % Cannot buy shares with negative net worth
    || agej*ypp>=11 && aprime<-f_coll*hprimecost ...  % Collateral constraint on borrowing (for older buyers that earn real money)
    || hprime<h && aprime<0 ...                       % Cannot sell down a house that is collateralized
    || agej>=Jr && hprime==0 && aprime<0)             % Ban pensioners from negative assets (if they don't own houses)
    return 
end

% Housing services (based on housing stock)
if h==0
    hs=0.5*houseservices*minhouse;
else
    hs=houseservices*h;
end

% Calculate car services value as we see it
if car~=0
    carservices=carservices_j;
else
    carservices=0.5;
end

c_pp=Electrify_4HouseholdConsumptionFn(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e, ...
    pension,AccidentBeqS,AccidentBeqAH,w,P0,D, ...
    kappa_j,tau_l,tau_d,tau_cg,S_agej_first,S_agej_peak_first,S_agej_peak_last,S_agej_last, ...
    ypp,agej,Jr,r,r_wedge,f_htc,rentprice,cpi_energy,pv_pct_cost,energy_pct_cost,energy_pct_brown,carbon_tax);

% If we are aiming for a starter loan, what loan can we afford?  Car not included
net_worth_prime=P0*sprime+aprime+hprimecost;
if aprime<0 && agej*ypp<11
    maxloan=-0.5*((10+ypp)-agej*ypp)/10;
    if net_worth_prime<maxloan
        if net_worth_prime+c_pp>maxloan
            % We could have put this into aprime ...
            % ... but asset_grid might be too small
            % This keeps state feasible, but disfavored
            c_pp=exp(-500);
        else
            % Limit starter loan needed to get people going
            return
        end
    end
end

benefit_used=0;
if c_pp<=0 && saprime<1 && cprime==0 && hprime==0
    % Agent can't make ends meet and sold down what they can sell: take the benefit
    if c_pp+max_benefit>0.2
        benefit_used=0.2-c_pp;
        c_pp=0.2;
    end
end

if c_pp>0
    % Adding one to all valid solutions doesn't alter results of searching
    % for optima, but does make output more legible when debugging.
    F=1+(((c_pp^(1-sigma_h-sigma_c))*(hs^sigma_h)*(carservices^sigma_c))^(1-sigma))/(1-sigma) -psi*(labor^(1+eta))/(1+eta); % The utility function
    % Disfavor using benefit...forces max labor participation
    F=F-100*benefit_used;
end

% Warm-glow bequest; must handle aprime<0
if agej==J % Final period
    if net_worth_prime<0
        % Died too far in debt...shouldn't happen
        F=-Inf;
    else
        % Our warmglow includes selling our next period house assets
        warmglow=warmglow1*(net_worth_prime^(1-warmglow2))/(1-warmglow2);
        F=F+warmglow;
    end
end


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
