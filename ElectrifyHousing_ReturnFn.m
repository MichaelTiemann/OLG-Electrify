function F=ElectrifyHousing_ReturnFn(labor,installpv,aprime,hprime,a,h,solarpv,z, ...
    r,r_wedge,w,sigma,agej,Jr,pension,kappa_j,psi,eta,sigma_h,f_htc,minhouse,rentprice,f_coll,houseservices,pv_pct_cost,energy_pct_cost)
% Note: experienceasset, so first inputs are (d,a,z,...)
% vfoptions.refine_d: only decisions d1,d3 are input to ReturnFn (and this model has no d1)

% Ban people from overinstalling solar
% if installpv && solarpv>0
%     F=-Inf;
%     return;
% end

% Make buying/selling a house costly/illiquid
htc=0; % house transaction cost
if hprime~=h
    htc=f_htc*(h+hprime);
end

pvinstallcost=0;
if installpv
    if (h+hprime)==0
        % No house -> no solar
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
    s=houseservices*sqrt(h);
    rentalcosts=0;
end

% Interest rate (loans vs. deposits)
r_rate=r;
if a<0
    r_rate=r+r_wedge;
end

F=-Inf;


if agej<Jr % If working age
    % Note: +h and -hprime; if a is negative, interest is negative
    c=(1+r_rate)*a+labor*w*kappa_j*z+h-hprime;
else % Retirement
    c=(1+r_rate)*a+pension+h-hprime;
end
c=c-htc-rentalcosts-pvinstallcost-(energy_pct_cost*min(h,1)^3*(1-solarpv/30))-aprime;

if c>0
    F=(((c^(1-sigma_h))*(s^sigma_h))^(1-sigma))/(1-sigma) -psi*(labor^(1+eta))/(1+eta); % The utility function
end

if aprime<-f_coll*hprime
    F=-Inf; % Collateral constraint on borrowing
end

%% Ban pensioners from negative assets
if agej>=Jr && a<0
    F=-Inf;
end


end
