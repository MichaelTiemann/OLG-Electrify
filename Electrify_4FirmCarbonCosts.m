function y_carbon_cost=Electrify_4FirmCarbonCosts( ...
    pvnew,kprime,k,pv,z, ...
    w, ...
    ypp,pvinstalled_firm,pvmax_firm,delta,pv_delta,alpha_k,alpha_l,Ek,ek,energy_pct_brown,carbon_tax)

y_energy_cost=Electrify_4FirmEnergyCosts(pvnew,kprime,k,pv,z,w,ypp,pvinstalled_firm,pvmax_firm,delta,pv_delta,alpha_k,alpha_l,Ek,ek,carbon_tax);

% $30M/year revenues / $10B = 0.003 units of K per year (not period)
pv_cost_offset=min(pv*0.003,y_energy_cost);

% Energy sector emitted 76.4 Mt CO2e; 49/51 HH/firm split; cost of carbon = NZD $42-$2450 / tCO2e
% We use a magic number to get cost of carbon tax to be $1.6B @ $42/ton,
% which is 16% of a "unit of K", thus 0.16
y_carbon_cost=76.4e6*0.51*carbon_tax*energy_pct_brown*(y_energy_cost-pv_cost_offset)/7e9;


end
