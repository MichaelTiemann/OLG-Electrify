function energy_cost_pp=Electrify_4HouseholdEnergyCosts_single( ...
    labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e, ...
    w, ...
    ypp,cpi_energy,energy_pct_cost ...
    )

single_1=single(1);

energy_cost_pp=single(0);

% Car energy costs...
if car==1
    energy_cost_pp=0.04*w*(single_1+cpi_energy)*ypp;
elseif car==2
    if solarpv>0.5
        solarpv=solarpv-single(0.5);
    else
        energy_cost_pp=0.02*w*(single_1+cpi_energy)^0.5*ypp;
    end
end

if car~=2
    % car batteries make solarpv more effective...
    solarpv=solarpv/2;
end

% Add cost of housing energy; PV generation: 30kW (2 solar units) meets h==1 energy needs
energy_cost_pp=energy_cost_pp+(single_1+cpi_energy)*energy_pct_cost*(max(h^2,single_1)-solarpv/2)*ypp;

end

