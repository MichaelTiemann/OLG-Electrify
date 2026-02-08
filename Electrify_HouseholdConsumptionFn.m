function c=Electrify_HouseholdConsumptionFn( ...
    labor,buyhouse,aprime,hprime,a,h,solarpv,z,e, ...
    r,pension,AccidentBeq,w,P0,D,Lhscale, ...
    kappa_j,tau_l,tau_d,tau_cg,agej,Jr,...
    r_wedge,f_htc,rentprice,agej_pct_cost,pv_pct_cost,energy_pct_cost)
% Replace assets with 'share holdings'
% Get rid of progressive taxes
% Add Lhnormalize
% Sort out asset return vs credit wedge

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
    rentalcosts=rentprice;
else
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
    % don't forget to add +e to z in exp(z)
    c=(1-tau_l)*labor*w*kappa_j*exp(z+e)*Lhscale+((1-tau_d)*D+P0)*(a+AccidentBeq)+(1+agej_pct_cost)*(h-hprime); 
else % Retirement
    c=pension+((1-tau_d)*D+P0)*(a+AccidentBeq)+(1+agej_pct_cost)*(h-hprime);
end
% ...subtract the rest of the things:
% - house transaction costs - rental - pvinstall - energy costs (offset by pv generation) - capital gains tax - next period share holdings
c=c-htc-rentalcosts-(1+agej_pct_cost)*(pvinstallcost+energy_pct_cost*min(h,1)^2*(1-solarpv/3))-tau_cg*(P0-Plag)*(a+AccidentBeq)-P*aprime;


end
