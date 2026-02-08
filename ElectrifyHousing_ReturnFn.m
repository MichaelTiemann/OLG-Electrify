function F=ElectrifyHousing_ReturnFn(labor,buyhouse,aprime,hprime,a,h,solarpv,z, ...
    r,r_wedge,w,sigma,agej,Jr,pension,kappa_j,psi,eta,sigma_h,f_htc,minhouse,rentprice,f_coll,houseservices,agej_pct_cost,pv_pct_cost,energy_pct_cost)
% Note: experienceasset, so first inputs are (d,a,z,...)
% vfoptions.refine_d: only decisions d1,d3 are input to ReturnFn

% Make buying/selling a house costly/illiquid
htc=0; % house transaction cost
if hprime~=h
    htc=f_htc*(h+hprime);
end

pvinstallcost=0;
% buyhouse 2 and 4 are install/upgrade PV options
if buyhouse==2 || buyhouse==4
    if (h+hprime)==0
        % No house -> no solar
        pvinstallcost=Inf;
    elseif h==hprime
        % Pay the retrofit penalty
        pvinstallcost=1.1*pv_pct_cost*h;
    else % Changing house
        % PV costs approximately 5% of new house ($30K system for $600K house)
        pvinstallcost=pv_pct_cost*hprime;
    end
end

% Housing services (based on housing stock)
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
    c=(1+r_rate)*a+labor*w*kappa_j*z+(1+agej_pct_cost)*(h-hprime);
else % Retirement
    c=(1+r_rate)*a+pension+(1+agej_pct_cost)*(h-hprime);
end
% 3x10kW is full solar for single house; we allow up to 5x10kW
c=c-htc-rentalcosts-(1+agej_pct_cost)*(pvinstallcost+energy_pct_cost*min(h,1)^2*(1-solarpv/3))-aprime;

if c>0
    F=(((c^(1-sigma_h))*(s^sigma_h))^(1-sigma))/(1-sigma) -psi*(labor^(1+eta))/(1+eta); % The utility function
end

if aprime<-f_coll*hprime
    F=-Inf; % Collateral constraint on borrowing
end

%% Ban pensioners from negative assets
if agej>=Jr && aprime<0
    F=-Inf;
end


end
