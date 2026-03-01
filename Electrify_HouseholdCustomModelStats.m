function CustomStats=Electrify_HouseholdCustomModelStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate, ...
        n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,heteroagentoptions,vfoptions,simoptions)
CustomStats=struct();
% Just use the median values; 1st is mean; 5th is min/max
simoptions.whichstats=zeros(1,7);
simoptions.whichstats([1,2,5])=1;

CFnsToEvaluate.H1buy = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) h~=hprime && hprime==1;
CFnsToEvaluate.H2buy = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) h~=hprime && hprime==2;
CFnsToEvaluate.H3buy = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) h~=hprime && hprime==3;
CFnsToEvaluate.H4buy = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) h~=hprime && hprime==4;
CFnsToEvaluate.H1sell = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) h~=hprime && h==1;
CFnsToEvaluate.H2sell = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) h~=hprime && h==2;
CFnsToEvaluate.H3sell = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) h~=hprime && h==3;
CFnsToEvaluate.H4sell = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) h~=hprime && h==4;
CFnsToEvaluate.H_u = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) (h>0); % Unit house holdings
CFnsToEvaluate.PV_u = @(labor,buyhouse,sprime,aprime,hprime,s,a,h,solarpv,z,e) (solarpv>0); % Unit solarpv holdings

AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,CFnsToEvaluate,Parameters,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions);

mean_buyers=[AgeConditionalStats.H1buy.Mean; AgeConditionalStats.H2buy.Mean; AgeConditionalStats.H3buy.Mean; AgeConditionalStats.H4buy.Mean];
mean_sellers=[AgeConditionalStats.H1sell.Mean; AgeConditionalStats.H2sell.Mean; AgeConditionalStats.H3sell.Mean; AgeConditionalStats.H4sell.Mean];

hbuyers=1000*mean_buyers.*Parameters.mewj;
hsellers=1000*mean_sellers.*Parameters.mewj;
hunits=(1:4)';

hdemand=sum(round(hbuyers-hsellers),2).*hunits;

hdemand_total=sum(hdemand);

CustomStats.hbuyers_total=sum(sum(round(hbuyers),2).*hunits);
CustomStats.hsellers_total=sum(sum(round(hsellers),2).*hunits);
CustomStats.hdemand_total=hdemand_total;

if false
    pvals=PolicyInd2Val_FHorz(Policy,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions);
    
    % buyers have 1-3 shares; any assets; no house; no PV; z; younger side
    uu_s123_b=squeeze(Policy(5,1:3,:,2:5,:,1,20:40));
    logical_index=uu_s12_b>1;
    buyers=uu_s12_b(logical_index); % hprime > 1 is a buy decision
    
    % sellers have 1-3 shares; any assets; a house; any PV; z; older side
    uu_s123_s=squeeze(Policy(5,1:3,:,2:5,:,1,20:40));
    sellers=uu_s12_s(uu_s12_s==1); % hprime == 1 is a sell decision
    
    % keepers can have any shares, any assets; house doesn't change
    % we change the house index to be first
    uu_s12_k=squeeze(Policy(5,:,:,2:5,:,1,20:60));
    keepers = uu_s12_k(uu_s12_k(:,3) == uu_s12_k);
end

return
