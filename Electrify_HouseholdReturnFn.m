function F=Electrify_HouseholdReturnFn( ...
    labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
    pension,AccidentBeqS,AccidentBeqAH,w,P0,D, ...
    sigma,psi,eta,sigma_h,kappa_j,tau_l,tau_d,tau_cg,warmglow1,warmglow2,agej,Jr,J,...
    scenario,r,r_wedge,f_htc,minhouse,rentprice,f_coll,houseservices,cpi,pv_pct_cost,energy_pct_cost ...
    )

% Note: experienceasset, so first inputs are (d,a,z,e,...)
% vfoptions.refine_d: only decisions d1,d3 are input to ReturnFn

F=-Inf;

% buyhouse decisions
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

% Housing matters
rentalcosts=0; % Overwrite if housing in scenario
hs=1; % Housing services (based on housing stock)
htc=0; % house transaction cost
hcost=0;
hprimecost=0;
pvinstallcost=0;
if scenario==3
    rentalcosts=rentprice*sqrt(kappa_j);
    if h==0
        hs=0.5*houseservices*minhouse;
    else
        hs=houseservices*h;
        rentalcosts=0;
    end
    % Houses start at 4x annual wage
    hcost=4*h*(1+cpi);
    hprimecost=4*hprime*(1+cpi);
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
end

%% Allow/Disallow some trivial agent decisions
if (sprime-s>0 && aprime+hprimecost<0 ...             % Cannot buy shares with negative net worth
    || agej>=11 && aprime<-f_coll*hprimecost ...  % Collateral constraint on borrowing (for older buyers that earn real money)
    || hprime<h && aprime<0 ...                       % Cannot sell down a house that is collateralized
    || agej>=Jr && aprime<0)                          % Ban pensioners from negative assets (even if they own houses)
    return 
end

% We can get P (share price) from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have
P=((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg);

Plag=P; % As stationary general eqm

if agej<Jr % If working age
    %consumption = labor income + "other income" below
    c=(1-tau_l)*labor*w*kappa_j*exp(z+e); 
else % Retirement
    c=pension;
end
% Other income: accidental share bequest + share holdings (including dividend) - dividend tax + accidental asset+house bequest + (inflation-shock adjusted) net housing assets
c=c+((1-tau_d)*D+P0)*(s+AccidentBeqS)+AccidentBeqAH+(hcost-hprimecost);
% PV generation: 30kW (2 solar units) meets h==1 energy needs
%%% WTF c=c+(1+cpi)*energy_pct_cost*(solarpv/2);
if a<0 % In both cases, resulting `a` is added to consumption, then `aprime` subtracted
    % Subtract loan interest by adding diminishing assets
    c=c+(1+r+r_wedge)*a;
else
    % Deposit interest included in augmented assets
    c=c+(1+r)*a;
end
% ...subtract capital gains tax and next period share, asset holdings
c=c-tau_cg*(P0-Plag)*(s+AccidentBeqS)-P*sprime-aprime;
% ...subtract housing-related costs: transaction costs, rental or home maintenance costs, pv installation, and scaled energy costs
c=c-htc-rentalcosts-hcost*0.02-pvinstallcost-(1+cpi)*energy_pct_cost*max(h^1.5,1);

% If we are aiming for a starter loan, what loan can we afford?
net_worth_prime=P*sprime+aprime+hprimecost;
if aprime<0 && agej<11
    maxloan=-0.5*(11-agej)/10;
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
    F=(((c^(1-sigma_h))*(hs^sigma_h))^(1-sigma))/(1-sigma) -psi*(labor^(1+eta))/(1+eta); % The utility function
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
