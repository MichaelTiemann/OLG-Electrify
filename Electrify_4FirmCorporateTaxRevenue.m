function revenue=Electrify_4FirmCorporateTaxRevenue(electrification,kprime,pvprime,k,pv,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,Ek,ek)
% Whether we set it up so that dividends or equity issuance is the decision
% variable is unimportant, here I use dividends as the decision variable.

% Note: r is not needed anywhere here, it is relevant to the firm via the discount factor.

% We can solve a static problem to get the firm labor input
l=(w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % This is just w=Marg. Prod. Labor, but rearranged

% Output.  See https://profstevekeen.substack.com/p/the-role-of-energy-in-economics
% We could use (Ek*ek)^alpha_k or (Ek*ek) as part of the TFP multiplier
y=(Ek*ek)*z*(k^alpha_k)*(l^alpha_l)*ypp;

% If Y is full GDP ($440B), then Ek=125 TWh and ek=$440B/125TWh=$3.52/kWh
% 69 TWh to be electrified (56 TWh already renewable); need 46,000 MW generation
y_energy_cost=0.045*y; % Assume energy cost is 4.5% of firm production
% 200GWh PV/year * 1000 MWh/GWh * $150/MWh = $30M PV/year (vs $630B)
pv_cost_offset=pv*k/21000;

% 200GWh/year = 133MW*1500h/yr = $220M cost @ $1.65M/MW; $220M/$630B = 0.00035
new_pv_cost=(pvprime-pv)*k/3000;

% Profit
profit=y-(w*l+y_energy_cost-pv_cost_offset)*ypp;

% Investment
delta_pp=(1+delta)^ypp-1;
invest=kprime-(1-delta)^ypp*k+new_pv_cost;

% Capital-adjustment costs
capitaladjcost=(capadjconstant/2)*((invest/k-delta_pp)^2) *k; 

% Taxable corporate income
T=profit-delta_pp*k-phi*capitaladjcost;
% -delta*k: investment expensing
% phi is the fraction of capitaladjcost that can be deducted from corporate taxes

revenue=tau_corp*T;

end
