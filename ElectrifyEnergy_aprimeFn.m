function aprime=ElectrifyEnergy_aprimeFn(installpv,pv,ypp)

if installpv<11
    aprime=min(pv+installpv,200);
else
    aprime=0.99*pv;
end

end
