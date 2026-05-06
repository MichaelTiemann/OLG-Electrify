function y_energy_cost_pp=Electrify_4FirmEnergyCosts( ...
    pvnew,kprime,k,pv,z, ...
    w, ...
    ypp,pvinstalled_firm,pvmax_firm,delta,pv_delta,alpha_k,alpha_l,Ek,ek,carbon_tax)
% Ek is the energy required by capital allocation k
% ek is the energy efficiency (more is better, like TPF)

% We can solve a static problem to get the firm labor input
l=(w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % This is just w=Marg. Prod. Labor, but rearranged

% Output.  See https://profstevekeen.substack.com/p/the-role-of-energy-in-economics
% We could use (Ek*ek)^alpha_k or (Ek*ek) as part of the TFP multiplier
y_pp=(Ek*ek)*z*(k^alpha_k)*(l^alpha_l)*ypp;

y_energy_cost_pp=0.131*y_pp/ek; % Assume energy cost is 13.1% of firm production


end
