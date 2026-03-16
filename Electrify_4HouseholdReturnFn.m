function F=Electrify_4HouseholdReturnFn( ...
    labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e, ...
    pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
    sigma,psi,eta,sigma_h,kappa_j,tau_l,tau_d,tau_cg,warmglow1,warmglow2,ypp,agej,Jr,J,...
    scenario,r_pp,r_wedge_pp,f_htc,minhouse,rentprice,f_coll,houseservices,carservices_j,cpi_cost,pv_pct_cost,energy_pct_cost ...
    )
% Get rid of progressive taxes
% Implement depreciation model:
%   Car services A(t) = (1-delta_a)*A(t-1) + I(a,t)
%   Housing services H(t) = (1-delta_h)*H(t-1) + I(h,t)
%   delta_a is high depreciation; delta_h is low depreciation

% Note: experienceasset, so first inputs are (d,a,z,e,...)
% vfoptions.refine_d: only decisions d1,d3 are input to ReturnFn

F=-Inf;

%% Housing matters
if buyhouse==0
    if hprime~=0
        % Forbid owning house when buyhouse=0
        return
    end
elseif buyhouse>=3
    if hprime==0 || hprime~=h
        % Forbid selling/changing house we say we are keeping
        return
    end
end

rentalcosts=rentprice*ypp;
hs=1; % Housing services (based on housing stock)
htc=0; % house transaction cost
hcost=0;
hprimecost=0;
pvinstallcost=0;
if scenario==3
    if h==0
        hs=0.5*houseservices*minhouse;
    else
        hs=houseservices*h;
        rentalcosts=0;
    end
    % Houses start at 2x annual wage
    hcost=2*h*(1+cpi_cost);
    hprimecost=2*hprime*(1+cpi_cost);
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
end

%% Allow/Disallow some trivial agent decisions
if (sprime-s>0 && aprime+hprimecost<0 ...             % Cannot buy shares with negative net worth
    || agej*ypp>=11 && aprime<-f_coll*hprimecost ...  % Collateral constraint on borrowing (for older buyers that earn real money)
    || hprime<h && aprime<0 ...                       % Cannot sell down a house that is collateralized
    || agej>=Jr && aprime<0)                          % Ban pensioners from negative assets (even if they own houses)
    return 
end

%% Car matters
sigma_c=sigma_h/3;
% Car costs 50% annual wage, or can trade at 25% annual wage
if cprime==0
    if car>0
        carcost=-0.25*(1+cpi_cost);
    else
        carcost=0;
    end
elseif cprime==3
    carcost=3;
else
    if car==0
        carcost=0.5*(1+cpi_cost);
    else
        carcost=0.25*(1+cpi_cost);
    end
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
% Other income: accidental share bequest + share holdings (including dividend) - dividend tax + accidental asset+house bequest + (inflation-shock adjusted) net housing assets
c=c+((1-tau_d)*D_pp+P0)*(s+AccidentBeqS_pp)+AccidentBeqS_pp+AccidentBeqAH_pp+(hcost-hprimecost);
% PV generation: 30kW (2 solar units) meets h==1 energy needs
c=c+(1+cpi_cost)*energy_pct_cost*(solarpv/2)*ypp;
if a<0 % In both cases, resulting `a` is added to consumption, then `aprime` subtracted
    % Subtract loan interest by adding diminishing assets
    c=c+(1+r_pp+r_wedge_pp)*a;
else
    % Deposit interest included in augmented assets
    c=c+(1+r_pp)*a;
end
% ...subtract capital gains tax and next period share, asset holdings
c=c-tau_cg*(P0-Plag)*(s+AccidentBeqS_pp)-P*sprime-aprime;
% ...subtract housing-related costs: transaction costs, rental or home maintenance costs, pv installation, and scaled energy costs
c=c-htc-rentalcosts-hcost*0.02*ypp-pvinstallcost-(1+cpi_cost)*energy_pct_cost*max(h^1.5,1)*ypp;
% ...subtract car costs (purchase, sale, and/or maintenance)
if carcost~=0
    c=c-carcost;
if cprime>0
    % annual insurance, maintenance, WOF, fuel etc.
    % TODO: sort petrol vs. electric
    c=c-0.10*ypp;
end

% If we are aiming for a starter loan, what loan can we afford?  Car not included
net_worth_prime=P*sprime+aprime+hprimecost;
if aprime<0 && agej*ypp<11
    maxloan=-0.5*((10+ypp)-agej*ypp)/10;
    if net_worth_prime<maxloan
        if net_worth_prime+c>maxloan
            % We could have put this into aprime ...
            % ... but asset_grid might be too small
            % This keeps state feasible, but disfavored
            c=exp(-500);
        else
            % Limit starter loan needed to get people going
            return
        end
    end
end

if c>0
    F=(((c^(1-sigma_h-sigma_c))*(hs^sigma_h)*(carservices_j^sigma_c))^(1-sigma))/(1-sigma) -psi*(labor^(1+eta))/(1+eta); % The utility function
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
