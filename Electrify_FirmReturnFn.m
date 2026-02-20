function F=Electrify_FirmReturnFn( ...
    dividend,kprime,k,z, ...
    w,D_pp, ...
    ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg)
% Whether we set it up so that dividends or equity issuance is the decision
% variable is unimportant, here I use dividends as the decision variable.

% Note: r is not needed anywhere here, it is relevant to the firm via the discount factor.

F=-Inf;

% We can solve a static problem to get the firm labor input
l=(w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % This is just w=Marg. Prod. Labor, but rearranged

% Output
y=z*(k^alpha_k)*(l^alpha_l)*ypp;

% Profit
profit=y-w*l*ypp;

% Investment
delta_pp=(1+delta)^ypp-1;
invest=kprime-(1-delta)^ypp*k;

% Capital-adjustment costs
capitaladjcost=(capadjconstant/2)*((invest/k-delta_pp)^2) *k; 

% Taxable corporate income
T=profit-delta_pp*k-phi*capitaladjcost;
% -delta_pp*k: investment expensing
% phi is the fraction of capitaladjcost that can be deducted from corporate taxes

% Firms financing constraint gives the new equity issuance
dividend_pp=(1+dividend)^ypp-1;
s=dividend_pp+invest+capitaladjcost-(profit-tau_corp*T);

% Firms per-period objective
if s>=0 % enforce that 'no share repurchases allowed'
    F=((1-tau_d)/(1-tau_cg))*dividend_pp-s;
    % Disfavor discrepencies between dividends paid and expected (D)
    F=F-(D_pp-dividend_pp)^2;
end

% Note: dividend payments cannot be negative is enforced by the grid on
% dividends which has a minimum value of zero

end
