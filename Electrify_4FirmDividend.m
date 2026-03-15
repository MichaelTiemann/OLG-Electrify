function dividend_pp=Electrify_4FirmDividend( ...
    electrification,kprime,pvprime,k,pv,z, ...
    w, ...
    ypp,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi)
% Whether we set it up so that dividends or equity issuance is the decision
% variable is unimportant, here I use dividends as the decision variable.

% Note: r is not needed anywhere here, it is relevant to the firm via the discount factor.

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

% This is the marginal dividend payable without allocating new shares
s=0;
dividend_pp=s+(profit-tau_corp*T)-invest-capitaladjcost;
if dividend_pp<0
    % We will issue new shares and provide a discounted dividend
    dividend_pp=0.1;
elseif dividend_pp<=0.2
    % We will issue new shares and provide a full dividend
    dividend_pp=0.2;
end

end
