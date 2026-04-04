function aprime=ElectrifyEnergy_aprimeFn(installpv,pv,ypp)

% Decision function limits the state space of `installpv`
aprime=min(pv+installpv,200);


end
