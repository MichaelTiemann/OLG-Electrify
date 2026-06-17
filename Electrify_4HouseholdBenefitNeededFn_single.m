%% Calculate benefit needed to meet consumption needs
function benefit=Electrify_4HouseholdBenefitNeededFn_single(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e, ...
    pension,AccidentBeqS,AccidentBeqAH,w,P0,D, ...
    kappa_j,tau_l,tau_d,tau_cg,S_agej_first,S_agej_peak_first,S_agej_peak_last,S_agej_last, ...
    ypp,agej,Jr,r,r_wedge,f_htc,rentprice,cpi_energy,pv_pct_cost,energy_pct_cost,energy_pct_brown,carbon_tax)

benefit=single(0);
if saprime>=1 || cprime>0 || hprime>0 || labor==0
    % No benefit if agent has meaningful assets to sell down and no labor to offer
    return
end

c_pp=Electrify_4HouseholdConsumptionFn_single(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e, ...
    pension,AccidentBeqS,AccidentBeqAH,w,P0,D, ...
    kappa_j,tau_l,tau_d,tau_cg,S_agej_first,S_agej_peak_first,S_agej_peak_last,S_agej_last, ...
    ypp,agej,Jr,r,r_wedge,f_htc,rentprice,cpi_energy,pv_pct_cost,energy_pct_cost,energy_pct_brown,carbon_tax);

if c_pp<=0
    if agej<Jr
        % Scale benefit by labor to favor those who work most
        benefit=0.1*(single(1)+labor)-c_pp;
    else
        benefit=single(0.2)-c_pp; % Big sorry to pensioners who cannot afford to live
    end
end

end
