function F=Electrify_FirmReturnFn( ...
    dividend,kprime,k,z, ...
    w, D_target, ...
    delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg)
% Whether we set it up so that dividends or equity issuance is the decision
% variable is unimportant, here I use dividends as the decision variable.

% Note: r is not needed anywhere here, it is relevant to the firm via the discount factor.

F=-Inf;

% We can solve a static problem to get the firm labor input
% When alpha_k=0.311, alpha_l=0.65, then
%  @ k=0.2: l=0.07 (0.35 multiplier)
%  @ k=1.2: l=0.34 (0.29 multiplier)
%  @ k=2.2: l=0.59 (0.27 multiplier)
%  @ k=4.2: l=1.05 (0.25 multiplier)
%  @ k=8.2: l=1.89 (0.23 multiplier)
l=(w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % This is just w=Marg. Prod. Labor, but rearranged

% Output; y ~ 1.54 * labor (close to 55% GM)
% [0.2, 1.2, 2.2, 4.2, 8.2].^0.311 .* [0.07, 0.34, 0.59, 1.05, 1.89].^0.65 => 0.1076    0.5249    0.9069    1.6129    2.9100
y=z*(k^alpha_k)*(l^alpha_l);

% Profit; profit ~ 50% of labor
% [0.11, 0.52, 0.91, 1.61,2.91] - [0.07, 0.34, 0.59, 1.05, 1.89] => 0.0376    0.1849    0.3169    0.5629    1.0200
profit=y-w*l;

% Investment
invest=kprime-(1-delta)*k; % largely tracks delta*k = 0.054*k

% Capital-adjustment costs (k>0 always)
capitaladjcost=(capadjconstant/2)*((invest/k-delta)^2)*k; 

% Taxable corporate income
T=profit-delta*k-phi*capitaladjcost;
% -delta*k: investment expensing
% phi is the fraction of capitaladjcost that can be deducted from corporate taxes

% Firms financing constraint gives the new equity issuance
% Larger dividends means more shares must be sold
s=dividend+invest+capitaladjcost-(profit-tau_corp*T);

% Firms per-period objective
if s>=0 % enforce that 'no share repurchases allowed'
    % When tau_d==tau_cg, F=profit-invest-capitaladjcost-tau_corp*(profit-delta*k-phi*capitaladjcost)
    % Thus, while larger dividends gets us here, they net out and don't influence F
    % Add term to prefer greater Y and dividends closer to 20%
    F=(((1-tau_d)/(1-tau_cg))*dividend-s)+(1-(dividend-D_target)^2);
end

% Note: dividend payments cannot be negative is enforced by the grid on
% dividends which has a minimum value of zero

end
