function c=Electrify_HouseholdConsumptionFn( ...
    labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
    pension,AccidentBeqS,AccidentBeqAH,w,P0,D, ...
    kappa_j,tau_l,tau_d,tau_cg,agej,Jr,Lhscale,...
    r,r_wedge,f_htc,rentprice,agej_pct_cost,pv_pct_cost,energy_pct_cost)
% Get rid of progressive taxes
% Add Lhnormalize

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
        pvinstallcost=1.1*pv_pct_cost*(1+agej_pct_cost)*h;
    else % Changing house
        % PV costs approximately 5% of new house ($30K system for $600K house)
        pvinstallcost=pv_pct_cost*(1+agej_pct_cost)*hprime;
    end
end

% Housing services (based on housing stock)
if h==0
    rentalcosts=rentprice;
else
    rentalcosts=0;
end

% We can get P from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have
P=((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg);

Plag=P; % As stationary general eqm

if agej<Jr % If working age
    %consumption = labor income + "other income" below
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
% house transaction costs - rental - pvinstall - energy costs (offset by pv generation) - capital gains tax - next period share, asset, and house holdings
c=c-htc-rentalcosts-pvinstallcost-(1+agej_pct_cost)*(energy_pct_cost*max(h,1)^2*(1-solarpv/2))-tau_cg*(P0-Plag)*(s+AccidentBeqS)-P*sprime-aprime;

end
