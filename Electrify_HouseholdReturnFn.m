function F=Electrify_HouseholdReturnFn( ...
    labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
    r,pension,AccidentBeqS,AccidentBeqAH,w,P0,D,Lhscale, ...
    sigma,psi,eta,sigma_h,kappa_j,tau_l,tau_d,tau_cg,warmglow1,warmglow2,agej,Jr,J,...
    r_wedge,f_htc,minhouse,rentprice,f_coll,houseservices,agej_pct_cost,pv_pct_cost,energy_pct_cost ...
    )
% Get rid of progressive taxes
% Add Lhnormalize

% Note: experienceasset, so first inputs are (d,a,z,e,...)
% vfoptions.refine_d: only decisions d1,d3 are input to ReturnFn

F=-Inf;

if hprime>0 || aprime~=0
    % Not buying houses or assets right now
    return
    % Later...if buyhouse=0, forbid hprime>h
end

if buyhouse==0
    if hprime~=0
        % Forbid owning house when buyhouse=0
        return
    end
elseif hprime==0 || hprime~=h
    if buyhouse>=3
        % Forbid selling/changing house we say we are keeping
        return
    end
end

% We can get P (share price) from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have
P=((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg);

%% Allow/Disallow some trivial agent decisions
if sprime>0 && aprime+(1+agej_pct_cost)*hprime<0
    % Cannot buy shares with negative net worth
    return
end
net_worth_prime=P*sprime+aprime+(1+agej_pct_cost)*hprime;
if agej<6
    if net_worth_prime<-0.55*(6-agej)/5
        % Starter loan needed to get people going
        return
    end
elseif (net_worth_prime < 0 ...                          % Don't allow net debt
        || aprime<-f_coll*(1+agej_pct_cost)*hprime ...   % Collateral constraint on borrowing
        || agej>=Jr && aprime<0)                         % Ban pensioners from negative assets (even if they own houses)
    return 
end

Plag=P; % As stationary general eqm

% Make buying/selling a house costly/illiquid
htc=0; % house transaction cost
if hprime~=h
    htc=f_htc*(1+agej_pct_cost)*(h+hprime);
end

pvinstallcost=0;
% buyhouse 2 and 4 are install/upgrade PV options
if buyhouse==2 || buyhouse==4
    if (h+hprime)==0
        % No house -> no solar
        pvinstallcost=Inf;
    elseif h==hprime
        % Pay the retrofit penalty
        pvinstallcost=1.1*pv_pct_cost*(1+agej_pct_cost)*h;
    else % Changing house
        % PV costs approximately 5% of new house ($30K system for $600K house)
        pvinstallcost=pv_pct_cost*(1+agej_pct_cost)*hprime;
    end
end

% Housing services (based on housing stock)
if h==0
    hs=0.5*houseservices*minhouse;
    rentalcosts=rentprice;
else
    hs=houseservices*sqrt(h);
    rentalcosts=0;
end

if agej<Jr % If working age
    %consumption = labor income plus other "other income" below
    c=(1-tau_l)*labor*w*kappa_j*exp(z+e)*Lhscale; 
else % Retirement
    c=pension;
end
% Other income: accidental share bequest + share holdings (including dividend) - dividend tax + accidental asset+house bequest + (inflation-shock adjusted) net housing assets
c=c+((1-tau_d)*D+P0)*(s+AccidentBeqS)+AccidentBeqAH+(1+agej_pct_cost)*(h-hprime);
if a<0
    % Subtract loan interest by adding a negative number
    c=c+(1+r+r_wedge)*a;
else
    % Add deposit interest
    c=c+(1+r)*a;
end
% ...subtract the rest of the things:
% house transaction costs - rental - pvinstall - energy costs (offset by pv generation) - capital gains tax - next period share, asset holdings
c=c-htc-rentalcosts-pvinstallcost-(1+agej_pct_cost)*(energy_pct_cost*max(h,1)^2*(1-solarpv/2))-tau_cg*(P0-Plag)*(s+AccidentBeqS)-P*sprime-aprime;

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
