function F=ElectrifyHousing_ReturnFn(savings,installpv,aprime,hprime,a,h,solarpv,z,r,w,sigma,agej,Jr,pension,kappa_j,sigma_h,f_htc,minhouse,rentprice,f_coll,houseservices,pv_pct_cost,energy_pct_cost)
% Note: experienceasset, so first inputs are (d,a,z,...)
% vfoptions.refine_d: only decisions d1,d3 are input to ReturnFn (and this model has no d1)

% Ban people from overinstalling solar
if installpv && solarpv>0
    F=-Inf;
    return;
end

% Make buying/selling a house costly/illiquid
htc=0; % house transaction cost
if hprime~=h
    htc=f_htc*(h+hprime);
end

pvinstallcost=0;
if installpv
    if (h+hprime)==0
        pvinstallcost=Inf;
    elseif h==hprime
        % Pay the retrofit penalty
        pvinstallcost=1.1*pv_pct_cost*h;
    else % Changing house
        % PV costs approximately 5% of house ($30K system for $600K house)
        pvinstallcost=pv_pct_cost*hprime;
    end
end

% Housing services
if h==0
    s=0.5*houseservices*minhouse;
    rentalcosts=rentprice;
else
    s=houseservices*h;
    rentalcosts=0;
end

F=-Inf;
if agej<Jr % If working age
    c=w*kappa_j*z+(1+r)*a-savings+(h-hprime)-htc-rentalcosts-pvinstallcost-(energy_pct_cost*(1-solarpv/30)); % Note: +h and -hprime
else % Retirement
    % give a rent subsidy to elderly for no good reason-rentalcosts;
    c=pension+(1+r)*a-savings+(h-hprime)-htc-rentalcosts-pvinstallcost-(energy_pct_cost*(1-solarpv/30));
end

if c>0 && s>0
    F=(((c^(1-sigma_h))*(s^sigma_h))^(1-sigma))/(1-sigma); % The utility function
end

if savings<-f_coll*hprime
    F=-Inf; % Collateral constraint on borrowing
end

% Negative savings is only allowed in the form of a safe mortgage. This is dealt with via the savingsFn.

%% Ban pensioners from negative assets
if agej>=Jr && savings<0
    F=-Inf;
end


end
