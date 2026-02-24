function c=Electrify_HouseholdConsumptionFn( ...
    labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
    pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
    kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
    r_pp,r_r_wedge_pp,f_htc,rentprice,agej_pct_cost,pv_pct_cost,energy_pct_cost)

% Housing matters
hcost=0;
hprimecost=0;
rentalcosts=0;
if h==0
    rentalcosts=rentprice*ypp;
end
if h+hprime>0
    % Houses start at 2x annual wage
    hcost=2*h*(1+agej_pct_cost);
    hprimecost=2*hprime*(1+agej_pct_cost);
end

% We can get P from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have
P=((1-tau_cg)*P0 + (1-tau_d)*D_pp)/(1+r_pp-tau_cg);

Plag=P; % As stationary general eqm
htc=0; % house transaction cost
pvinstallcost=0; % solarpv installation cost

% Make buying/selling a house costly/illiquid
if hprime~=h
    htc=f_htc*(h+hprime);
end

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

if agej<Jr % If working age
    %consumption = labor income + "other income" below
    c=(1-tau_l)*labor*w*kappa_j*exp(z+e)*ypp; 
else % Retirement
    c=pension*ypp;
end
% Other income: accidental share bequest + share holdings (including dividend) - dividend tax + accidental asset+house bequest + (inflation-shock adjusted) net housing assets
c=c+((1-tau_d)*D_pp+P0)*(s+AccidentBeqS_pp)+AccidentBeqAH_pp+(hcost-hprimecost);
% PV generation: 30kW (2 solar units) meets h==1 energy needs
c=c+(1+agej_pct_cost)*energy_pct_cost*(solarpv/2)*ypp;
if a<0
    % Subtract loan interest by adding a negative number
    c=c+(1+r_r_wedge_pp)*a;
else
    % Add deposit interest
    c=c+(1+r_pp)*a;
end
% ...subtract capital gains tax and next period share, asset holdings
c=c-tau_cg*(P0-Plag)*(s+AccidentBeqS_pp)-P*sprime-aprime;
% ...subtract housing-related costs:  pv installation/upgrade, house transaction costs, rental or home maintenance costs, and scaled energy costs
c=c-pvinstallcost-htc-rentalcosts-hcost*1.02*ypp-(1+agej_pct_cost)*energy_pct_cost*max(h^2,1)*ypp;

end
