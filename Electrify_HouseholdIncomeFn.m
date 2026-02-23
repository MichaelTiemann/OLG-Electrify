function income=Electrify_HouseholdIncomeFn( ...
    labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
    pension,AccidentBeqS_pp,AccidentBeqAH_pp,w,P0,D_pp, ...
    kappa_j,tau_l,tau_d,tau_cg,ypp,agej,Jr, ...
    r_pp,r_r_wedge_pp,agej_pct_cost)
% Replace assets with 'share holdings'
% Get rid of progressive taxes
% Add Lhnormalize
% Figure out real estate transactions

% We can get P from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have
P=((1-tau_cg)*P0 + (1-tau_d)*D_pp)/(1+r_pp-tau_cg);

Plag=P; % As stationary general eqm

if agej<Jr % If working age
    %consumption = labor income + "other income" below
    % income just is consumption but without subtracting the term for next period share holdings (-P*sprime) or asset holdings (aprime)
    income=(1-tau_l)*labor*w*kappa_j*exp(z+e)*ypp;
else % Retirement
    income=pension*ypp;
end
% Other income: accidental share bequest + share holdings (including dividend) - capital gains + accidental asset+house bequest + (inflation-shock adjusted) net housing assets
income=income+((1-tau_d)*D_pp+P0)*(s+AccidentBeqS_pp)-tau_cg*(P0-Plag)*(s+AccidentBeqS_pp)+AccidentBeqAH_pp+(1+agej_pct_cost)*(h-hprime);
if a<0
    % Subtract loan interest by adding a negative number
    income=income+r_r_wedge_pp*a;
else
    % Add deposit interest
    income=income+r_pp*a;
end

end
