function income=Electrify_HouseholdIncomeFn( ...
    labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
    pension,AccidentBeqS,AccidentBeqAH,w,P0,D, ...
    kappa_j,tau_l,tau_d,tau_cg,agej,Jr, ...
    r,cpi,energy_pct_cost)

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
% And thus we have
P=((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg);

Plag=P; % As stationary general eqm

if agej<Jr % If working age
    %consumption = labor income + "other income" below
    % income just is consumption but without subtracting the term for next period share holdings (-P*sprime) or asset holdings (aprime)
    income=(1-tau_l)*labor*w*kappa_j*exp(z+e);
else % Retirement
    income=pension;
end
% Other income: accidental share bequest + share holdings (including dividend) - capital gains + accidental asset+house bequest + net housing assets
income=income+((1-tau_d)*D+P0)*(s+AccidentBeqS)-tau_cg*(P0-Plag)*(s+AccidentBeqS)+AccidentBeqAH+(hcost-hprimecost);
% PV generation: 30kW (2 solar units) meets h==1 energy needs
income=income+(1+cpi)*energy_pct_cost*(solarpv/2);
if a>0
    % Add deposit interest
    income=income+r*a;
end

end
