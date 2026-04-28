function aprime=ElectrifyFirm_aprimeFn(pvnew,pv,pvinstalled_firm,pvmax_firm,pv_delta)

% Decision function limits the state space of `pvnew`
aprime=min(pv+pvnew,pvmax_firm);

aprime=aprime*(1-pv_delta);

end
