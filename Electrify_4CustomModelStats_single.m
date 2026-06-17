function CustomStats=Electrify_4CustomModelStats_single(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions)
CustomStats=struct();
% Just use the median values; 1st is mean; 5th is min/max
simoptions.whichstats=zeros(1,7,'single');
simoptions.whichstats([1,2,5])=1;

CFnsToEvaluate.H1buy = @(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e) h~=hprime && hprime==1;
CFnsToEvaluate.H2buy = @(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e) h~=hprime && hprime==2;
CFnsToEvaluate.H3buy = @(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e) h~=hprime && hprime==3;
CFnsToEvaluate.H4buy = @(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e) h~=hprime && hprime==4;
CFnsToEvaluate.H1sell = @(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e) h~=hprime && h==1;
CFnsToEvaluate.H2sell = @(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e) h~=hprime && h==2;
CFnsToEvaluate.H3sell = @(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e) h~=hprime && h==3;
CFnsToEvaluate.H4sell = @(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e) h~=hprime && h==4;
CFnsToEvaluate.H1 = @(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e) h==1; % Total H1 house holdings
CFnsToEvaluate.H2 = @(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e) h==2; % Total H2 house holdings
CFnsToEvaluate.H3 = @(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e) h==3; % Total H3 house holdings
CFnsToEvaluate.H_u = @(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e) (h>0); % Unit house holdings
CFnsToEvaluate.PV_h = @(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e) solarpv; % total solarpv holdings
CFnsToEvaluate.PV_u = @(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e) (solarpv>0); % Unit solarpv holdings
CFnsToEvaluate.petrol_car = @(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e) (car==1);
CFnsToEvaluate.ev_car = @(labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e) (car==2);
CFnsToEvaluate.BenefitNeeded = @( ...
    labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e, ...
    pension,AccidentBeqS,AccidentBeqAH,w,P0,D, ...
    kappa_j,tau_l,tau_d,tau_cg,S_agej_first,S_agej_peak_first,S_agej_peak_last,S_agej_last, ...
    ypp,agej,Jr,r,r_wedge,f_htc,rentprice,cpi_energy,pv_pct_cost,energy_pct_cost,energy_pct_brown,carbon_tax...
    ) Electrify_4HouseholdBenefitNeededFn_single( ...
    labor,buyhouse,saprime,cprime,hprime,sa,car,h,solarpv,z,e, ...
    pension,AccidentBeqS,AccidentBeqAH,w,P0,D, ...
    kappa_j,tau_l,tau_d,tau_cg,S_agej_first,S_agej_peak_first,S_agej_peak_last,S_agej_last, ...
    ypp,agej,Jr,r,r_wedge,f_htc,rentprice,cpi_energy,pv_pct_cost,energy_pct_cost,energy_pct_brown,carbon_tax);

[~,ii] = find(strcmp(Names_i, 'household'));
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist.household,Policy.household,CFnsToEvaluate,Parameters,[],n_d.household,n_a.household,n_z.household,N_j.household,d_grid.household,a_grid.household,z_grid.household,PType_Options(simoptions,Names_i{ii}));

mean_buyers=[AgeConditionalStats.H1buy.Mean; AgeConditionalStats.H2buy.Mean; AgeConditionalStats.H3buy.Mean; AgeConditionalStats.H4buy.Mean];
mean_sellers=[AgeConditionalStats.H1sell.Mean; AgeConditionalStats.H2sell.Mean; AgeConditionalStats.H3sell.Mean; AgeConditionalStats.H4sell.Mean];
[max_benefit,max_benefit_J]=max(AgeConditionalStats.BenefitNeeded.Maximum);

hbuyers=1000*mean_buyers.*Parameters.mewj;
hsellers=1000*mean_sellers.*Parameters.mewj;
hunits=(single(1):4)';

hdemand=sum(round(hbuyers-hsellers),2).*hunits;

hdemand_total=sum(hdemand);

CustomStats.hbuyers_total=sum(sum(round(hbuyers),2).*hunits);
CustomStats.hsellers_total=sum(sum(round(hsellers),2).*hunits);
CustomStats.hdemand_total=hdemand_total;
CustomStats.H1=sum(1000*AgeConditionalStats.H1.Mean.*Parameters.mewj);
CustomStats.H2=sum(1000*AgeConditionalStats.H2.Mean.*Parameters.mewj);
CustomStats.H3=sum(1000*AgeConditionalStats.H3.Mean.*Parameters.mewj);
CustomStats.PV_h=sum(1000*AgeConditionalStats.PV_h.Mean.*Parameters.mewj);
CustomStats.petrol_car=sum(1000*AgeConditionalStats.petrol_car.Mean.*Parameters.mewj);
CustomStats.ev_car=sum(1000*AgeConditionalStats.ev_car.Mean.*Parameters.mewj);
if max_benefit>Parameters.max_benefit
    CustomStats.max_benefit=max_benefit;
    CustomStats.max_benefit_J=max_benefit_J;
end

FnsToEvaluate2.pvnew_f=@(pvnew,kprime,k,pv,z,pvinstalled_firm,pvmax_firm) pvnew;
FnsToEvaluate2.pv_f=@(pvnew,kprime,k,pv,z,pvinstalled_firm,pvmax_firm) pv;
simoptions_temp=PType_Options(simoptions,Names_i{strcmp(Names_i, 'firm')});
simoptions_temp.parallel=2;
AggVars=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist.firm, Policy.firm, FnsToEvaluate2, Parameters, [], n_d.firm, n_a.firm, n_z.firm, d_grid.firm, a_grid.firm, z_grid.firm, simoptions_temp);
CustomStats.pvnew_f=AggVars.pvnew_f.Mean;
CustomStats.pv_f=AggVars.pv_f.Mean;

return

end
