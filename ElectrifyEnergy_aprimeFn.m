function aprime=ElectrifyEnergy_aprimeFn(pvnew,pv,pvinstalled_energy,pvmax_energy,pv_delta)

% Decision function limits the state space of `pvnew`
aprime=min(pv+pvnew,pvmax_energy);

aprime=aprime*(1-pv_delta);

end
