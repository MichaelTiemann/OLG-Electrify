function aprime=ElectrifyFirm_aprimeFn(installpv,pv,ypp)

if installpv<3
    aprime=min(pv+installpv,100);
else
    aprime=0.99*pv;
end

end
