function y_carbon_cost=Electrify_4FirmCarbonCosts( ...
    pvnew,kprime,k,pv,z, ...
    w, ...
    pvinstalled_firm,pvmax_firm,delta,pv_delta,alpha_k,alpha_l,Ek,ek,energy_pct_brown,carbon_tax)
% Ek is the energy required by capital allocation k
% ek is the energy efficiency (more is better, like TPF)

% We can solve a static problem to get the firm labor input
l=(w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % This is just w=Marg. Prod. Labor, but rearranged

% Output.  See https://profstevekeen.substack.com/p/the-role-of-energy-in-economics
% We could use (Ek*ek)^alpha_k or (Ek*ek) as part of the TFP multiplier
y=(Ek*ek)*z*(k^alpha_k)*(l^alpha_l);

y_energy_cost=0.131*y/ek; % Assume energy cost is 13.1% of firm production

% $30M/year revenues / $10B = 0.003 units of K per year (not period)
pv_cost_offset=min(pv*0.003,y_energy_cost);

% Energy sector emitted 76.4 Mt CO2e; 49/51 HH/firm split; cost of carbon = NZD $42-$2450 / tCO2e
% We use a magic number to get cost of carbon tax to be $1.6B @ $42/ton,
% which is 16% of a "unit of K", thus 0.16
y_carbon_cost=76.4e6*0.51*carbon_tax*energy_pct_brown*(y_energy_cost-pv_cost_offset)/7e9;


end
