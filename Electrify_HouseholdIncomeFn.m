function income_pp=Electrify_HouseholdIncomeFn( ...
    labor,buyhouse,saprime,hprime,sa,h,solarpv,z,e, ...
    pension,AccidentBeqS,AccidentBeqAH,w,P0,D, ...
    kappa_j,tau_l,tau_d,tau_cg,S_agej_first,S_agej_peak_first,S_agej_peak_last,S_agej_last, ...
    ypp,agej,Jr,r,cpi,energy_pct_cost)

if saprime<1
    sprime=0;
else
    sprime=floor(saprime);
end

if sa<1
    s=0;
    a=sa;
else
    s=floor(sa);
    a=rem(sa,1);
end

hcost=0;
hprimecost=0;
if h+hprime>0
    % Houses start at 4x annual wage
    hcost=4*h*(1+cpi);
    hprimecost=4*hprime*(1+cpi);
end

% We can get P from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have P=((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg);

P=P0;
if saprime>=s
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
    % income just is consumption but without subtracting the term for next period share holdings (-P*sprime) or asset holdings (aprime)
    income_pp=(1-tau_l)*labor*w*kappa_j*exp(z+e)*ypp;
else % Retirement
    income_pp=pension*ypp;
end
% Other income: accidental share bequest + share holdings (including dividend) - capital gains + accidental asset+house bequest + net housing assets
income_pp=income_pp+((1-tau_d)*D*ypp+P0)*(s+AccidentBeqS)-cg+AccidentBeqAH+(hcost-hprimecost);
% PV generation: 30kW (2 solar units) meets h==1 energy needs
income_pp=income_pp+(1+cpi)*energy_pct_cost*(solarpv/2)*ypp;
if a>0
    % Add deposit interest
    income_pp=income_pp+((1+r)^ypp-1)*a;
end

end
