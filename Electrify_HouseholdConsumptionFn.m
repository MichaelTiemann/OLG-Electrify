function c=Electrify_HouseholdConsumptionFn( ...
    labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
    pension,AccidentBeqS,AccidentBeqAH,w,P0,D, ...
    kappa_j,tau_l,tau_d,tau_cg,agej,Jr, ...
    r,r_wedge,f_htc,rentprice,cpi,pv_pct_cost,energy_pct_cost)

% Housing matters
rentalcosts=0;
htc=0; % house transaction cost
hcost=0;
hprimecost=0;
pvinstallcost=0;
if h+hprime>0
    % Houses start at 4x annual wage
    hcost=4*h*(1+cpi);
    hprimecost=4*hprime*(1+cpi);
elseif h==0
    rentalcosts=rentprice;
end

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
        pvinstallcost=1.1*pv_pct_cost*(1+cpi)*h;
    else % Changing house
        % PV costs approximately 5% of new house ($30K system for $600K house)
        pvinstallcost=pv_pct_cost*(1+cpi)*hprime;
    end
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
    c=(1-tau_l)*labor*w*kappa_j*exp(z+e); 
else % Retirement
    c=pension;
end
% Other income: accidental share bequest + share holdings (including dividend) - dividend tax + accidental asset+house bequest + (inflation-shock adjusted) net housing assets
c=c+((1-tau_d)*D+P0)*(s+AccidentBeqS)+AccidentBeqS+AccidentBeqAH+(hcost-hprimecost);
% PV generation: 30kW (2 solar units) meets h==1 energy needs
c=c+(1+cpi)*energy_pct_cost*(solarpv/2);
if a<0
    % Subtract loan interest by adding a negative number
    c=c+(1+r+r_wedge)*a;
else
    % Add deposit interest
    c=c+(1+r)*a;
end
% ...subtract capital gains tax and next period share, asset holdings
c=c-tau_cg*(P0-Plag)*(s+AccidentBeqS)-P*sprime-aprime;
% ...subtract housing-related costs:  pv installation/upgrade, house transaction costs, rental or home maintenance costs, and scaled energy costs
c=c-htc-rentalcosts-hcost*0.02-pvinstallcost-(1+cpi)*energy_pct_cost*max(h^1.5,1);

end
