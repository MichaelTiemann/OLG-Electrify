function c_pp=Electrify_HouseholdConsumptionFn( ...
    labor,buyhouse,saprime,hprime,sa,h,solarpv,z,e, ...
    pension,AccidentBeqS,AccidentBeqAH,w,P0,D, ...
    kappa_j,tau_l,tau_d,tau_cg,S_agej_first,S_agej_peak_first,S_agej_peak_last,S_agej_last,ypp,agej,Jr, ...
    r,r_wedge,f_htc,rentprice,cpi,pv_pct_cost,energy_pct_cost)

[sprime,aprime,s,a]=decode_sa(saprime,sa);

% Housing matters
rentalcosts_pp=0;
htc=0; % house transaction cost
hcost=0;
hprimecost=0;
pvinstallcost=0;
if h+hprime>0
    % Houses start at 4x annual wage
    hcost=4*h*(1+cpi);
    hprimecost=4*hprime*(1+cpi);
elseif h==0
    rentalcosts_pp=rentprice*sqrt(kappa_j)*ypp;
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
% And thus we have P=((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r_pp-tau_cg);

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
    c_pp=(1-tau_l)*labor*w*kappa_j*exp(z+e)*ypp; 
else % Retirement
    c_pp=pension*ypp;
end
% Other income: accidental share bequest + share holdings (including dividend) - dividend tax + accidental asset+house bequest + (inflation-shock adjusted) net housing assets
c_pp=c_pp+((1-tau_d)*D*ypp+P0)*(s+AccidentBeqS)+AccidentBeqS+AccidentBeqAH+(hcost-hprimecost);
% PV generation: 30kW (2 solar units) meets h==1 energy needs
c_pp=c_pp+(1+cpi)*energy_pct_cost*(solarpv/2)*ypp;
if a<0
    % Subtract loan interest by adding a negative number
    c_pp=c_pp+(1+r+r_wedge)^ypp*a;
else
    % Add deposit interest
    c_pp=c_pp+(1+r)^ypp*a;
end
% ...subtract capital gains tax and next period share, asset holdings
c_pp=c_pp-cg-P*sprime-aprime;
% ...subtract housing-related costs:  pv installation/upgrade, house transaction costs, rental or home maintenance costs, and scaled energy costs
c_pp=c_pp-htc-rentalcosts_pp-hcost*0.02*ypp-pvinstallcost-(1+cpi)*energy_pct_cost*max(h^1.5,1)*ypp;


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
