function energy_cost_pp=Electrify_4HouseholdEnergyCosts( ...
    labor,buyhouse,sprime,aprime,cprime,hprime,s,a,car,h,solarpv,z,e, ...
    w, ...
    ypp,energy_cpi,energy_pct_cost,energy_pct_brown,carbon_tax ...
    )
% Implement depreciation model:
%   Car services A(t) = (1-delta_a)*A(t-1) + I(a,t)
%   Housing services H(t) = (1-delta_h)*H(t-1) + I(h,t)
%   delta_a is high depreciation; delta_h is low depreciation

% Note: experienceasset, so first inputs are (d,a,z,e,...)
% vfoptions.refine_d: only decisions d1,d3 are input to ReturnFn

% Add cost of housing
energy_cost_pp=energy_pct_cost*max(h^1.5,1)*ypp;

% Add car energy costs
if car>0
    % Energy costs...
    if car==1
        energy_cost_pp=energy_cost_pp+0.02*w*ypp;
    elseif solarpv>0.5 % car==2
        solarpv=solarpv-0.5; % Deduct car charging from energy demand
    end
end

if car~=2
    % car batteries make solarpv more effective...
    solarpv=solarpv/2;
end

% PV generation: 10kW/30kWh (2 solar units) meets h==1 energy needs
energy_cost_pp=energy_cost_pp-energy_pct_cost*(solarpv/2)*ypp;

energy_cost_pp=(1+energy_cpi)*energy_cost_pp;

end
