function energy_cost=Electrify_4HouseholdEnergyCosts( ...
    labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e, ...
    w, ...
    energy_cpi,energy_pct_cost ...
    )
% Implement depreciation model:
%   Car services A(t) = (1-delta_a)*A(t-1) + I(a,t)
%   Housing services H(t) = (1-delta_h)*H(t-1) + I(h,t)
%   delta_a is high depreciation; delta_h is low depreciation

% Note: experienceasset, so first inputs are (d,a,z,e,...)
% vfoptions.refine_d: only decisions d1,d3 are input to ReturnFn

energy_cost=0;

% Car energy costs...
if car==1
    energy_cost=energy_cost+0.041*w;
elseif car==2
    if solarpv>0.5
        solarpv=solarpv-0.5;
    else
        energy_cost=energy_cost+0.02*w;
    end
end

if car~=2
    % car batteries make solarpv more effective...
    solarpv=solarpv/2;
end

% Add cost of housing energy; PV generation: 30kW (2 solar units) meets h==1 energy needs
energy_cost=energy_cost+(1+energy_cpi)*energy_pct_cost*(max(h^1.5,1)-solarpv/2);

end
