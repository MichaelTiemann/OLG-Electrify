function aprime=ElectrifyEnergy_aprimeFn(pvnew,pv,ypp,pvinstalled_energy,pvmax_energy,pv_delta_pp)

% Decision function limits the state space of `pvnew`
aprime=min(pv+pvnew,pvmax_energy);

aprime=aprime*pv_delta_pp;

end
