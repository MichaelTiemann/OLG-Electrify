function F=Electrify_4FirmReturnFn( ...
    electrification,kprime,pvprime,k,pv,z, ...
    w, ...
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
invest=kprime-(1-delta)^ypp*k+electrification;

% Capital-adjustment costs
capitaladjcost=(capadjconstant/2)*((invest/k-delta_pp)^2) *k; 

% Taxable corporate income
T=profit-delta_pp*k-phi*capitaladjcost;
% -delta_pp*k: investment expensing
% phi is the fraction of capitaladjcost that can be deducted from corporate taxes

% Firms financing constraint gives the new equity issuance
% dividend_pp=(1+dividend)^ypp-1;
% s=dividend_pp+invest+capitaladjcost-(profit-tau_corp*T);

% This is the marginal dividend payable without allocating new shares
s=0;
dividend_pp=s+(profit-tau_corp*T)-invest-capitaladjcost;
if dividend_pp<0
    % We will issue new shares and provide a discounted dividend
    s=0.1-dividend_pp;
    dividend_pp=0.1;
elseif dividend_pp<=0.2
    % We will issue new shares and provide a full dividend
    s=0.2-dividend_pp;
    dividend_pp=0.2;
else
    % We are earning too much and cannot buy back shares.
    s=-1;
end

% Firms per-period objective
if s>=0 % enforce that 'no share repurchases allowed'
    % When tau_d==tau_cg, F=profit-invest-capitaladjcost-tau_corp*(profit-delta_pp*k-phi*capitaladjcost)
    % Add term to prefer greater Y and dividends closer to 20%
    F=(((1-tau_d)/(1-tau_cg))*dividend_pp-s)+y*(1-(dividend_pp-0.2)^2)/10;
end

% Note: dividend payments cannot be negative is enforced by the grid on
% dividends which has a minimum value of zero

end
