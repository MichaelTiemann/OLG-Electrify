function revenue_pp=Electrify_FirmCorporateTaxRevenue(dividend,kprime,k,z,w,ypp,delta_pp,alpha_k,alpha_l,capadjconstant,tau_corp,phi)
% Whether we set it up so that dividends or equity issuance is the decision
% variable is unimportant, here I use dividends as the decision variable.

% Note: r is not needed anywhere here, it is relevant to the firm via the discount factor.

% We can solve a static problem to get the firm labor input
l=(w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % This is just w=Marg. Prod. Labor, but rearranged

% Output
y=z*(k^alpha_k)*(l^alpha_l);

% Profit
profit_pp=(y-w*l)*ypp;

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
T_pp=profit_pp-delta_pp*k-phi*capitaladjcost_pp;
% -delta*k: investment expensing
% phi is the fraction of capitaladjcost that can be deducted from corporate taxes

revenue_pp=tau_corp*T_pp;

end
