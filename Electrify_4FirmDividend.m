function dividend=Electrify_4FirmDividend( ...
    pvnew,kprime,k,pv,z, ...
    w, ...
    pvinstalled_firm,pvmax_firm,delta,pv_delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,Ek,ek,energy_pct_brown,carbon_tax)
% Whether we set it up so that dividends or equity issuance is the decision
% variable is unimportant, here I use dividends as the decision variable.

% Note: r is not needed anywhere here, it is relevant to the firm via the discount factor.

% We can solve a static problem to get the firm labor input
l=(w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % This is just w=Marg. Prod. Labor, but rearranged

% Output.  See https://profstevekeen.substack.com/p/the-role-of-energy-in-economics
% We could use (Ek*ek)^alpha_k or (Ek*ek) as part of the TFP multiplier
y=(Ek*ek)*z*(k^alpha_k)*(l^alpha_l);

y_energy_cost=0.131*y/ek; % Assume energy cost is 13.1% of firm production

% For sake of argument, say K = $10B (so we need 11 K to get full Y)
% 200GWh/year = 133MW*1500h/yr = $220M cost @ $1.65M/MW; $220M/$10B = 0.022 units of K per PV
new_pv_cost=pvnew*0.022;
% $30M/year revenues / $10B = 0.003 units of K per year (not period)
pv_cost_offset=min(pv*0.003,y_energy_cost);

% If Y is $110B, then Ek=96 TWh and ek=$110B/96TWh=$1146 Y/MWh
% 53 TWh to be electrified (43 TWh already renewable); need 35,333 MW generation
% 200GWh PV/year * 1000 MWh/GWh * $150/MWh = $30M/PV/year (vs $110B)
% 53 TWh to electrify = $58.3B total costs

% Energy sector emitted 76.4 Mt CO2e; 49/51 HH/firm split; cost of carbon = NZD $42-$2450 / tCO2e
% We use a magic number to get cost of carbon tax to be $1.6B @ $42/ton,
% which is 16% of a "unit of K", thus 0.16
y_carbon_tax=76.4e6*0.51*carbon_tax*energy_pct_brown*(1-pv_cost_offset/y_energy_cost)/7e9;

% Profit
profit=y-w*l-y_energy_cost+pv_cost_offset-y_carbon_tax;

% Investment
pv_delta=0.02;
invest=kprime+new_pv_cost-(1-delta)*k-pv_delta*pv;

% Capital-adjustment costs (k>0 always)
if invest>=0
    capitaladjcost=(capadjconstant/2)*((invest/k-delta-pv_delta)^2)*k;
else
    capitaladjcost=0;
end

% Taxable corporate income
T=max(profit-delta*k-pv_delta*pv-phi*capitaladjcost,0);
% -delta*k: investment expensing; -pv_delta*pv: pv expensing
% phi is the fraction of capitaladjcost that can be deducted from corporate taxes


% This is the marginal dividend payable without allocating new shares
s=0;
mid_dividend=0.2;
dividend=s+(profit-tau_corp*T)-invest-capitaladjcost;
if dividend<0
    % We will issue new shares and provide a discounted dividend
    low_dividend=mid_dividend/2;
    s=low_dividend-dividend;
    dividend=low_dividend;
elseif dividend<=0.2
    % We will issue new shares and provide a full dividend
    s=mid_dividend-dividend;
    dividend=mid_dividend;
else % We don't need to issue shares and can pay rich dividend
end


end
