function s=Electrify_4FirmShareIssuance(electrification,kprime,pvprime,k,pv,z,w,ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,Ek,ek,carbon_tax)
% Whether we set it up so that dividends or equity issuance is the decision
% variable is unimportant, here I use dividends as the decision variable.

% Note: r is not needed anywhere here, it is relevant to the firm via the discount factor.

% We can solve a static problem to get the firm labor input
l=(w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % This is just w=Marg. Prod. Labor, but rearranged

% Output.  See https://profstevekeen.substack.com/p/the-role-of-energy-in-economics
% We could use (Ek*ek)^alpha_k or (Ek*ek) as part of the TFP multiplier
y_pp=(Ek*ek)*z*(k^alpha_k)*(l^alpha_l)*ypp;

% If Y is full GDP ($440B), then Ek=125 TWh and ek=$440B/125TWh=$3.52 GDP/kWh
% 69 TWh to be electrified (56 TWh already renewable); need 46,000 MW generation
y_carbon_tax=76.4e6*carbon_tax/440e9; % Energy sector emitted 76.4 Mt CO2e; social cost of carbon = NZD $2450 / tCO2e
y_energy_cost_pp=(0.045+y_carbon_tax)*y_pp; % Assume energy cost is 4.5% of firm production
% 200GWh PV/year * 1000 MWh/GWh * $150/MWh = $30M PV/year (vs $440B)
% 69 TWh to electrify = $10350M total costs
pv_cost_offset_pp=min(pv*ypp*30/10350,y_energy_cost_pp);

% 200GWh/year = 133MW*1500h/yr = $220M cost @ $1.65M/MW; $220M/$440B = 0.0005 max GDP
new_pv_cost=(pvprime-pv)*1/2000;

% Profit
profit_pp=y_pp-w*l*ypp-y_energy_cost_pp+pv_cost_offset_pp-new_pv_cost;

% Investment
delta_pp=(1+delta)^ypp-1;
invest_pp=kprime-(1-delta)^ypp*k+new_pv_cost;

% Capital-adjustment costs
capitaladjcost_pp=(capadjconstant/2)*((invest_pp/k-delta_pp)^2) *k; 

% Taxable corporate income
T=profit_pp-delta_pp*k-phi*capitaladjcost_pp;
% -delta*k: investment expensing
% phi is the fraction of capitaladjcost that can be deducted from corporate taxes

% Firms financing constraint gives the new equity issuance
s=0;
mid_dividend_pp=1.2^ypp-1;
dividend_pp=s+(profit_pp-tau_corp*T)-invest_pp-capitaladjcost_pp;
if dividend_pp<0
    % We will issue new shares and provide a discounted dividend
    low_dividend_pp=1.1^ypp-1;
    s=low_dividend_pp-dividend_pp;
    dividend_pp=low_dividend_pp;
elseif dividend_pp<=0.2
    % We will issue new shares and provide a full dividend
    s=mid_dividend_pp-dividend_pp;
    dividend_pp=mid_dividend_pp;
else
    % We don't need to issue shares and can pay rich dividend
    s=0.12345;
end

end
