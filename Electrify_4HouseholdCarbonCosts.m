function carbon_cost_pp=Electrify_4HouseholdCarbonCosts( ...
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

energy_cost_pp=0;

% ...subtract car costs (purchase, sale, and/or maintenance)
if car>0
    % Energy costs...
    if car==1
        energy_cost_pp=energy_cost_pp+0.02*w*ypp;
    else
        if solarpv>0.5
            solarpv=solarpv-0.5;
        else
            energy_pct_cost=energy_pct_cost+0.02;
        end
    end
end

% Add cost of housing
energy_cost_pp=energy_cost_pp+(1+energy_cpi)*energy_pct_cost*max(h^1.5,1)*ypp;

if car~=2
    % car batteries make solarpv more effective...
    solarpv=solarpv/2;
end

% PV generation: 30kW (2 solar units) meets h==1 energy needs
carbon_cost_pp=energy_pct_cost*(max(h^1.5,1)-solarpv/2)*ypp*energy_pct_brown*carbon_tax/2000;


end
