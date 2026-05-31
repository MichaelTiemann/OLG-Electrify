%% Calculate benefit needed to meet consumption needs
function benefit=Electrify_4HouseholdBenefitNeededFn(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e, ...
    pension,AccidentBeqS,AccidentBeqAH,w,P0,D, ...
    kappa_j,tau_l,tau_d,tau_cg,S_agej_first,S_agej_peak_first,S_agej_peak_last,S_agej_last, ...
    ypp,agej,Jr,r,r_wedge,f_htc,rentprice,cpi_energy,pv_pct_cost,energy_pct_cost,energy_pct_brown,carbon_tax)

benefit=0;
if saprime>=1 || cprime>0 || hprime>0
    % No benefit if agent has meaningful assets to sell down
    return
end

c_pp=Electrify_4HouseholdConsumptionFn(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e, ...
    pension,AccidentBeqS,AccidentBeqAH,w,P0,D, ...
    kappa_j,tau_l,tau_d,tau_cg,S_agej_first,S_agej_peak_first,S_agej_peak_last,S_agej_last, ...
    ypp,agej,Jr,r,r_wedge,f_htc,rentprice,cpi_energy,pv_pct_cost,energy_pct_cost,energy_pct_brown,carbon_tax)

if c_pp<=0
    if agej<Jr
        benefit=0.2-c_pp;
    else
        benefit=0.4-c_pp; % Big sorry to pensioners who cannot afford to live
    end
end

end
