function F=Electrify_4EnergyReturnFn(pvnew,kprime,k,pv,z,pvinstalled_energy,pvmax_energy,delta,pv_delta,EnergyCosts_h,EnergyCosts_f)
% Choosing unit industrial PV to be 200GWh/year generation...
% 200GWh/year = 133MW*1500h/yr = $220M cost @ $1.65M/MW; $220M/$440B = 0.0005 max GDP
% 200GWh PV/year * 1000 MWh/GWh * $150/MWh = Firm PV Energy Cost Offset $30M/PV/year (vs $440B)
% 69 TWh/year to electrify = 345*200GWh/year * $220M/200GWh/year = $75900M cost of 100% conversion

F=-Inf;

% Convert household demand to energy supply
% Starting point: 0.46 over 5 years = 0.092/a = $7360 per person * 2.9M = $21B
% $21B / $350/MWh = 60 TWh; 60 TWh => 40000 MW gen
MW_hh=(EnergyCosts_h*80000)*2.9e6/350/1500;
revenue_hh=EnergyCosts_h*80000*2.9e6; % for the period

% Starting point: 0.64 over 5 years = 0.128/a /0.045 = 2.844 * 155B = $440B/a
% 0.128 * 155B = $19.8B
% $19.8B / $150/MWh = 132 TWh (call it 125 TWh)
% 125 TWh / 1500 h / year = 83,333 MW gen
Y_firm=(EnergyCosts_f/0.045)*155e9;
MW_firm=Y_firm*125e12/440e9/1500e6;
revenue_firm=EnergyCosts_f*155e9; % for the period

% Each PV is 133MW generation
pv_req = (MW_hh+MW_firm)/133;

pv_maint=(pv*220e6)*delta;

% 200GWh/year = 133MW*1500h/yr = $220M cost @ $1.65M/MW; $220M/$440B = 0.0005 max GDP per PV
pvnew_cost=pvnew*220e6;

% Can't buy what we don't have money to buy
if pvnew_cost>kprime
    return
end

if pv_req<pv
    F=(revenue_hh+revenue_firm)-pv_maint;
else
    if MW_hh<pv*133
        MW_remainder=pv*133-MW_hh;
        MW_unmet=MW_firm-MW_remainder;
        F=revenue_hh+revenue_firm*MW_remainder/MW_firm-pv_maint-0.1*revenue_firm*(MW_unmet/MW_firm);
    else
        MW_unmet=MW_hh-pv*133+MW_firm;
        F=revenue_hh*pv*133/MW_hh-pv_maint-0.1*(revenue_hh*(1-pv*133/MW_hh)+revenue_firm);
    end
    F=F/440e9;
end


end
