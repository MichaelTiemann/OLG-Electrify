function aprime=ElectrifyFirm_aprimeFn(pvnew,pv,ypp,pvinstalled_firm,pvmax_firm,pv_delta_pp)

% Decision function limits the state space of `pvnew`
aprime=min(pv+pvnew,pvmax_firm);

aprime=aprime*pv_delta_pp;

end
