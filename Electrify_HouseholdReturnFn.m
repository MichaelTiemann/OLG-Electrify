function F=Electrify_HouseholdReturnFn( ...
    labor,buyhouse,saprime,hprime,sa,h,solarpv,z,e, ...
    pension,AccidentBeqS,AccidentBeqAH,w,P0,D, ...
    sigma,psi,eta,sigma_h,kappa_j,tau_l,tau_d,tau_cg,S_agej_first,S_agej_peak_first,S_agej_peak_last,S_agej_last,warmglow1,warmglow2,ypp,agej,Jr,J,...
    scenario,r,r_wedge,f_htc,minhouse,rentprice,f_coll,houseservices,cpi,pv_pct_cost,energy_pct_cost ...
    )

% Note: experienceasset, so first inputs are (d,a,z,e,...)
% vfoptions.refine_d: only decisions d1,d3 are input to ReturnFn

F=-Inf;

[sprime,aprime,s,a]=decode_sa(saprime,sa);

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
rentalcosts_pp=0; % Overwrite if housing in scenario
hs=1; % Housing services (based on housing stock)
htc=0; % house transaction cost
hcost=0;
hprimecost=0;
pvinstallcost=0;
if scenario==3
    rentalcosts_pp=rentprice*sqrt(kappa_j)*ypp;
    if h==0
        hs=0.5*houseservices*minhouse;
    else
        hs=houseservices*h;
        rentalcosts_pp=0;
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
    || agej*ypp>=11 && aprime<-f_coll*hprimecost ...  % Collateral constraint on borrowing (for older buyers that earn real money)
    || hprime<h && aprime<0 ...                       % Cannot sell down a house that is collateralized
    || agej>=Jr && aprime<0)                          % Ban pensioners from negative assets (even if they own houses)
    return 
end

% We can get P (share price) from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have P=((1-tau_cg)*P0 + (1-tau_d)*D_pp)/(1+r_pp-tau_cg);

% But in fact the price does not meaningfully represent the acquisition
% cost by a younger generation that is now older in this stationary
% distribution.  However, we can use the history of acquisition and
% disposals to impute when agents are buying and selling, and thus what
% capital gains they should pay.  We imagine that stocks earn 2x the
% risk-free rate of return (i.e., 2*r_pp) and that if we are selling before
% they peak, we are selling recently acquired stocks, whereas if we are
% selling at or after the peak of acquisition, we are selling long-term
% gains in a LIFO fashion.

% We take P0 as the price of the current stationary distribution, and we
% back-calculate what the price Plag may have been in the past.

P=P0;
if sprime>=s
    cg=0; % We are holding or buying, so no capital gains
else
    if agej<=S_agej_peak_first
        Plag=P0*(1-r)^ypp; % Dispose of shares presumably acquired recently
    elseif S_agej_peak_last==S_agej_last % Bulk liquidation
        % Sell all remaining shares from first acquisition to buy-point (using geometric mean to average acquisition cost)
        agej_bought=S_agej_peak_first-sqrt(S_agej_peak_first-S_agej_first);
        Plag=P0*(1-r)^(ypp*(agej-agej_bought));
    else
        % Estimate where we are past peak accumulation and mirror around to
        % proportional acquisition point
        agej_selling_pct=(agej-S_agej_peak_last)/(S_agej_last-S_agej_peak_last);
        agej_bought=S_agej_peak_first-agej_selling_pct*(S_agej_peak_first-S_agej_first);
        Plag=P0*(1-r)^(ypp*(agej-agej_bought));
    end
    cg=tau_cg*(P0-Plag)*(s+AccidentBeqS-sprime);
end

if agej<Jr % If working age
    %consumption = labor income + "other income" below
    c=(1-tau_l)*labor*w*kappa_j*exp(z+e)*ypp; 
else % Retirement
    c=pension*ypp;
end
% Other income: accidental share bequest + share holdings (including dividend) - dividend tax + accidental asset+house bequest + (inflation-shock adjusted) net housing assets
c=c+((1-tau_d)*D*ypp+P0)*(s+AccidentBeqS)+AccidentBeqAH+(hcost-hprimecost);
% PV generation: 30kW (2 solar units) meets h==1 energy needs
%%% WTF c=c+(1+cpi)*energy_pct_cost*(solarpv/2)*ypp;
if a<0 % In both cases, resulting `a` is added to consumption, then `aprime` subtracted
    % Subtract loan interest by adding diminishing assets
    c=c+(1+r+r_wedge)^ypp*a;
else
    % Deposit interest included in augmented assets
    c=c+(1+r)^ypp*a;
end
% ...subtract capital gains tax and next period share, asset holdings
c=c-cg-P*sprime-aprime;
% ...subtract housing-related costs: transaction costs, rental or home maintenance costs, pv installation, and scaled energy costs
c=c-htc-rentalcosts_pp-hcost*0.02*ypp-pvinstallcost-(1+cpi)*energy_pct_cost*max(h^1.5,1)*ypp;

% If we are aiming for a starter loan, what loan can we afford?
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
