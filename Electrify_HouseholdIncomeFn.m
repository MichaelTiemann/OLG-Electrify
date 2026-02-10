function income=Electrify_HouseholdIncomeFn( ...
    labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e, ...
    r,pension,AccidentBeqS,AccidentBeqAH,w,P0,D,Lhscale, ...
    kappa_j,tau_l,tau_d,tau_cg,agej,Jr, ...
    r_wedge,agej_pct_cost)
% Replace assets with 'share holdings'
% Get rid of progressive taxes
% Add Lhnormalize
% Figure out real estate transactions

% We can get P from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have
P=((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg);

Plag=P; % As stationary general eqm

if agej<Jr % If working age
    %consumption = labor income + "other income" below
    % income just is consumption but without subtracting the term for next period share holdings (-P*sprime)
    income=(1-tau_l)*labor*w*kappa_j*exp(z+e)*Lhscale;
else % Retirement
    income=pension;
end
% Other income: accidental share bequest + share holdings (including dividend) + accidental asset+house bequest + (inflation-shock adjusted) net housing assets
c=c+((1-tau_d)*D+P0)*(s+AccidentBeqS)+AccidentBeqAH+(1+agej_pct_cost)*(h-hprime);
if a<0
    % Subtract loan interest by adding a negative number
    c=c+(1+r+r_wedge)*a;
else
    % Add deposit interest
    c=c+(1+r)*a;
end

end
