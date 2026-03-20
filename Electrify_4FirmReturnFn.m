function F=Electrify_4FirmReturnFn( ...
    electrification,kprime,pvprime,k,pv,z, ...
    w, ...
    ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg,Ek,ek)
% Whether we set it up so that dividends or equity issuance is the decision
% variable is unimportant, here I use dividends as the decision variable.

% Note: r is not needed anywhere here, it is relevant to the firm via the discount factor.

F=-Inf;

% Cannot uninstall PVs
if pvprime < pv || pvprime-pv > 2
    return
end

% We can solve a static problem to get the firm labor input
l=(w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % This is just w=Marg. Prod. Labor, but rearranged

% Output.  See https://profstevekeen.substack.com/p/the-role-of-energy-in-economics
% We could use (Ek*ek)^alpha_k or (Ek*ek) as part of the TFP multiplier
y=(Ek*ek)*z*(k^alpha_k)*(l^alpha_l)*ypp;

% If Y is full GDP ($440B), then Ek=125 TWh and ek=$440B/125TWh=$3.52/kWh
% 69 TWh to be electrified (56 TWh already renewable); need 46,000 MW generation
y_energy_cost=0.045*y; % Assume energy cost is 4.5% of firm production
% 200GWh PV/year * 1000 MWh/GWh * $150/MWh = $30M PV/year (vs $630B)
pv_cost_offset=pv*1/21000;

% 200GWh/year = 133MW*1500h/yr = $220M cost @ $1.65M/MW; $220M/$630B = 0.00035
new_pv_cost=(pvprime-pv)*1000/3000;

% Profit
profit_pp=y-w*l*ypp-y_energy_cost+pv_cost_offset*ypp;

% Investment
delta_pp=(1+delta)^ypp-1;
invest_pp=kprime-(1-delta)^ypp*k+new_pv_cost;

% Capital-adjustment costs
capitaladjcost_pp=(capadjconstant/2)*((invest_pp/(k*ypp)-delta_pp)^2) *(k*ypp); 

% Taxable corporate income
T=profit_pp-delta_pp*k-phi*capitaladjcost_pp;
% -delta_pp*k: investment expensing
% phi is the fraction of capitaladjcost that can be deducted from corporate taxes

% Firms financing constraint gives the new equity issuance
% dividend_pp=(1+dividend)^ypp-1;
% s=dividend_pp+invest+capitaladjcost-(profit-tau_corp*T);

% This is the marginal dividend payable without allocating new shares
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
end

% Firms per-period objective
if s>=0 % enforce that 'no share repurchases allowed'
    % When tau_d==tau_cg, F=profit-invest-capitaladjcost-tau_corp*(profit-delta_pp*k-phi*capitaladjcost)
    % Add term to prefer greater Y and dividends closer to 20%
    F=(((1-tau_d)/(1-tau_cg))*dividend_pp-s)+y*(1-(dividend_pp-mid_dividend_pp)^2)/10;
end

% Note: dividend payments cannot be negative is enforced by the grid on
% dividends which has a minimum value of zero

end
