function energy_cost_pp=Electrify_4FirmEnergyCosts( ...
    electrification,kprime,pvprime,k,pv,z, ...
    w, ...
    ypp,alpha_k,alpha_l,Ek,ek)
% Ek is the energy required by capital allocation k
% ek is the energy efficiency (more is better, like TPF)

% We can solve a static problem to get the firm labor input
l=(w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % This is just w=Marg. Prod. Labor, but rearranged

% Output.  See https://profstevekeen.substack.com/p/the-role-of-energy-in-economics
% We could use (Ek*ek)^alpha_k or (Ek*ek) as part of the TFP multiplier
y_pp=(Ek*ek)*z*(k^alpha_k)*(l^alpha_l)*ypp;

% If Y is full GDP ($440B), then Ek=125 TWh and ek=$440B/125TWh=$3.52/kWh
% 69 TWh to be electrified (56 TWh already renewable); need 46,000 MW generation
energy_cost_pp=0.045*y_pp; % Assume energy cost is 4.5% of firm production

end
