function carbon_cost_pp=Electrify_4HouseholdCarbonCosts( ...
    labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e, ...
    w, ...
    ypp,cpi_energy,energy_pct_cost,energy_pct_brown,carbon_tax ...
    )
% Implement depreciation model:
%   Car services A(t) = (1-delta_a)*A(t-1) + I(a,t)
%   Housing services H(t) = (1-delta_h)*H(t-1) + I(h,t)
%   delta_a is high depreciation; delta_h is low depreciation

% Note: experienceasset, so first inputs are (d,a,z,e,...)
% vfoptions.refine_d: only decisions d1,d3 are input to ReturnFn

energy_cost_pp=Electrify_4HouseholdEnergyCosts(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e,w,ypp,cpi_energy,energy_pct_cost);
carbon_cost_pp=energy_cost_pp*energy_pct_brown*carbon_tax*ypp/200; % Magic divisor to hit 0.7% hh income at $42/t CO2e


end
