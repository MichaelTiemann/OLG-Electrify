function energy_cost_pp=Electrify_4HouseholdEnergyCosts( ...
    labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e, ...
    w, ...
    ypp,energy_cpi,energy_pct_cost ...
    )

energy_cost_pp=0;

% Car energy costs...
if car==1
    energy_cost_pp=0.04*w*(1+energy_cpi)*ypp;
elseif car==2
    if solarpv>0.5
        solarpv=solarpv-0.5;
    else
        energy_cost_pp=0.02*w*(1+energy_cpi)^0.5*ypp;
    end
end

if car~=2
    % car batteries make solarpv more effective...
    solarpv=solarpv/2;
end

% Add cost of housing energy; PV generation: 30kW (2 solar units) meets h==1 energy needs
energy_cost_pp=energy_cost_pp+(1+energy_cpi)*energy_pct_cost*(max(h^1.5,1)-solarpv/2)*ypp;

end

