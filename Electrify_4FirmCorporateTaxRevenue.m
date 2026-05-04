function revenue_pp=Electrify_4FirmCorporateTaxRevenue(pvnew,kprime,k,pv,z,w,ypp,pvinstalled_firm,pvmax_firm,delta_pp,pv_delta_pp,alpha_k,alpha_l,capadjconstant,tau_corp,phi,Ek,ek,energy_pct_brown,carbon_tax)
% Whether we set it up so that dividends or equity issuance is the decision
% variable is unimportant, here I use dividends as the decision variable.

% Note: r is not needed anywhere here, it is relevant to the firm via the discount factor.

% We can solve a static problem to get the firm labor input
l=(w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % This is just w=Marg. Prod. Labor, but rearranged

% Output.  See https://profstevekeen.substack.com/p/the-role-of-energy-in-economics
% We could use (Ek*ek)^alpha_k or (Ek*ek) as part of the TFP multiplier
y_pp=(Ek*ek)*z*(k^alpha_k)*(l^alpha_l)*ypp;

y_energy_cost_pp=0.131*y_pp/ek; % Assume energy cost is 13.1% of firm production

% For sake of argument, say K = $10B (so we need 11 K to get full Y)
% 200GWh/year = 133MW*1500h/yr = $220M cost @ $1.65M/MW; $220M/$10B = 0.022 units of K per PV
new_pv_cost=pvnew*0.022;
% $30M/year revenues / $10B = 0.003 units of K per year (not period)
pv_cost_offset_pp=min((pvinstalled_firm+pv)*ypp*0.003,y_energy_cost_pp);

% If Y is $110B, then Ek=96 TWh and ek=$110B/96TWh=$1146 Y/MWh
% 53 TWh to be electrified (43 TWh already renewable); need 35,333 MW generation
% 200GWh PV/year * 1000 MWh/GWh * $150/MWh = $30M/PV/year (vs $110B)
% 53 TWh to electrify = $58.3B total costs

% Energy sector emitted 76.4 Mt CO2e; 49/51 HH/firm split; cost of carbon = NZD $42-$2450 / tCO2e
% We use a magic number to get cost of carbon tax to be $1.6B @ $42/ton,
% which is 16% of a "unit of K", thus 0.16
y_carbon_tax_pp=76.4e6*0.51*carbon_tax*energy_pct_brown*(1-pv_cost_offset_pp/y_energy_cost_pp)/7e9;

% Profit
profit_pp=y_pp-w*l*ypp-y_energy_cost_pp+pv_cost_offset_pp-y_carbon_tax_pp;

% Investment
invest_pp=kprime-(1-delta_pp)*k;

% Capital-adjustment costs (k>0 always)
if invest_pp>=0
    k_pp=k*ypp;
    capitaladjcost_pp=(capadjconstant/2)*((invest_pp/k_pp-delta_pp)^2)*k_pp;
else
    capitaladjcost_pp=0;
end

% Taxable corporate income
T_pp=max(profit_pp-delta_pp*k-pv_delta_pp*pv-phi*capitaladjcost_pp,0);
% -delta_pp*k: investment expensing; -pv_delta_pp*pv: pv expensing
% phi is the fraction of capitaladjcost that can be deducted from corporate taxes

revenue_pp=tau_corp*T_pp;

end
