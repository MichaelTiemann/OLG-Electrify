function F=Electrify_HouseholdReturnFn( ...
    labor,buyhouse,aprime,hprime,a,h,solarpv,z,e, ...
    r,pension,AccidentBeq,w,P0,D,Lhscale, ...
    sigma,psi,eta,sigma_h,kappa_j,tau_l,tau_d,tau_cg,warmglow1,warmglow2,agej,Jr,J,...
    r_wedge,f_htc,minhouse,rentprice,f_coll,houseservices,agej_pct_cost,pv_pct_cost,energy_pct_cost ...
    )
% Replace assets with 'share holdings'
% Get rid of progressive taxes
% Add Lhnormalize

% Note: experienceasset, so first inputs are (d,a,z,e,...)
% vfoptions.refine_d: only decisions d1,d3 are input to ReturnFn

% Make buying/selling a house costly/illiquid
htc=0; % house transaction cost
if hprime~=h
    htc=f_htc*(h+hprime);
end

pvinstallcost=0;
% buyhouse 2 and 4 are install/upgrade PV options
if buyhouse==2 || buyhouse==4
    if (h+hprime)==0
        % No house -> no solar
        pvinstallcost=Inf;
    elseif h==hprime
        % Pay the retrofit penalty
        pvinstallcost=1.1*pv_pct_cost*h;
    else % Changing house
        % PV costs approximately 5% of new house ($30K system for $600K house)
        pvinstallcost=pv_pct_cost*hprime;
    end
end

% Housing services (based on housing stock)
if h==0
    s=0.5*houseservices*minhouse;
    rentalcosts=rentprice;
else
    s=houseservices*sqrt(h);
    rentalcosts=0;
end

% Interest rate (loans vs. deposits)
r_rate=r;
if a<0
    r_rate=r+r_wedge;
end

% We can get P from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have
P=((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg);

Plag=P; % As stationary general eqm

F=-Inf;
if agej<Jr % If working age
    %consumption = labor income + accidental bequest + share holdings
    %(including dividend) + net house holdings...less a few things below
    c=(1-tau_l)*labor*w*kappa_j*exp(z+e)*Lhscale+((1-tau_d)*D+P0)*(a+AccidentBeq)+(1+agej_pct_cost)*(h-hprime); 
else % Retirement
    c=pension+((1-tau_d)*D+P0)*(a+AccidentBeq)+(1+agej_pct_cost)*(h-hprime);
end
% ...subtract the rest of the things:
% - house transaction costs - rental - pvinstall - energy costs (offset by pv generation) - capital gains tax - next period share holdings
c=c-htc-rentalcosts-(1+agej_pct_cost)*(pvinstallcost+energy_pct_cost*min(h,1)^2*(1-solarpv/3))-tau_cg*(P0-Plag)*(a+AccidentBeq)-P*aprime;

if c>0
    F=(((c^(1-sigma_h))*(s^sigma_h))^(1-sigma))/(1-sigma) -psi*(labor^(1+eta))/(1+eta); % The utility function
end

% Warm-glow bequest; must handle aprime<0
if agej==J % Final period
    if aprime+hprime<0
        % Died too far in debt...shouldn't happen
        F=-Inf;
    else
        % Our warmglow includes selling our next period house assets
        warmglow=warmglow1*((aprime+hprime)^(1-warmglow2))/(1-warmglow2);
        F=F+warmglow;
    end
end

if aprime<-f_coll*hprime
    F=-Inf; % Collateral constraint on borrowing
end

%% Ban pensioners from negative assets
if agej>=Jr && aprime<0
    F=-Inf;
end


end
