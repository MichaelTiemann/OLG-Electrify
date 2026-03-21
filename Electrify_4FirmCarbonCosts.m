function y_carbon_cost_pp=Electrify_4FirmCarbonCosts( ...
    electrification,kprime,pvprime,k,pv,z, ...
    w, ...
    ypp,alpha_k,alpha_l,Ek,ek,pv_max,carbon_tax)
% Ek is the energy required by capital allocation k
% ek is the energy efficiency (more is better, like TPF)

% We can solve a static problem to get the firm labor input
l=(w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % This is just w=Marg. Prod. Labor, but rearranged

% Output.  See https://profstevekeen.substack.com/p/the-role-of-energy-in-economics
% We could use (Ek*ek)^alpha_k or (Ek*ek) as part of the TFP multiplier
y_pp=(Ek*ek)*z*(k^alpha_k)*(l^alpha_l)*ypp;

% If Y is full GDP ($440B), then Ek=125 TWh and ek=$440B/125TWh=$3.52 GDP/kWh
% 69 TWh to be electrified (56 TWh already renewable); need 46,000 MW generation
y_carbon_tax=76.4e6*carbon_tax/440e9; % Energy sector emitted 76.4 Mt CO2e; cost of carbon = NZD $35-$2450 / tCO2e
y_carbon_cost_pp=y_carbon_tax*(1-pv/pv_max)*y_pp;

end
