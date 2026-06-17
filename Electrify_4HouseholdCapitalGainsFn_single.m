function cg=Electrify_4HouseholdCapitalGainsFn_single( ...
    labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e, ...
    ypp,agej,at_death,P0,AccidentBeqS,r,tau_cg,S_agej_first,S_agej_peak_first,S_agej_peak_last,S_agej_last)
% Replace assets with 'share holdings'
% Get rid of progressive taxes
% Add Lhnormalize

% We can get P from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have P=((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg);

single_0=single(0); single_1=single(1);

sprime=single_0;
s=single_0;
if saprime>=single_1
    sprime=single(floor(saprime));
end
if sa>=single_1
    s=single(floor(sa));
end

maybe_liquidate=single_1;
if sprime>=s
    agej_bought=agej;
    cg=single_0; % We are holding or buying, so no capital gains
else
    if agej<=S_agej_peak_first
        agej_bought=agej-single_1;
        Plag=P0*(single_1-2*r)^ypp; % Dispose of shares presumably acquired recently
    elseif agej<S_agej_peak_last
        % We have been holding since peak acquisition
        agej_bought=S_agej_peak_first;
        Plag=P0*(single_1-2*r)^(ypp*(agej-agej_bought));
    elseif S_agej_peak_last==S_agej_last % Bulk liquidation
        % Sell all remaining shares from first acquisition to buy-point (using geometric mean to average acquisition cost)
        agej_bought=S_agej_peak_first-sqrt(S_agej_peak_first-S_agej_first);
        Plag=P0*(single_1-2*r)^(ypp*(agej-agej_bought));
        maybe_liquidate=single_0;
    else
        % Estimate where we are past peak accumulation and mirror around to
        % proportional acquisition point
        agej_selling_pct=(agej-S_agej_peak_last)/(S_agej_last-S_agej_peak_last);
        agej_bought=S_agej_peak_first-agej_selling_pct*(S_agej_peak_first-S_agej_first);
        Plag=P0*(single_1-2*r)^(ypp*(agej-agej_bought));
    end
    cg=tau_cg*(P0-Plag)*(s+AccidentBeqS-sprime);
end

if at_death && agej_bought>=S_agej_first && maybe_liquidate>0
    % Sell all remaining shares from first acquisition to buy-point (using geometric mean to average acquisition cost)
    Plag=P0*(single_1-2*r)^(ypp*(agej_bought-sqrt(agej_bought-S_agej_first)));
    cg=cg+tau_cg*(P0-Plag)*sprime;
end

end
