function F=ElectrifyHousing_ReturnFn(labor,buyhouse,aprime,hprime,a,h,solarpv,z, ...
    r,r_wedge,w,sigma,agej,Jr,pension,kappa_j,psi,eta,sigma_h,f_htc,minhouse,rentprice,f_coll,houseservices,agej_pct_cost,pv_pct_cost,energy_pct_cost)
% Note: experienceasset, so first inputs are (d,a,z,...)
% vfoptions.refine_d: only decisions d1,d3 are input to ReturnFn

F=-Inf;

%% Allow/Disallow some trivial agent decisions
net_worth_prime=aprime+(1+agej_pct_cost)*hprime;
if agej<6
    if net_worth_prime<-0.55*(6-agej)/5
        % Starter loan needed to get people going
        return
    end
elseif (net_worth_prime < 0 ...                          % Don't allow net debt
        || aprime<-f_coll*(1+agej_pct_cost)*hprime ...   % Collateral constraint on borrowing
        || agej>=Jr && aprime<0)                         % Ban pensioners from negative assets (even if they own houses)
    return 
end

% Make buying/selling a house costly/illiquid
htc=0; % house transaction cost
if hprime~=h
    htc=f_htc*(1+agej_pct_cost)*(h+hprime);
end

pvinstallcost=0;
% buyhouse 2 and 4 are install/upgrade PV options
if buyhouse==2 || buyhouse==4
    if (h+hprime)==0
        % No house -> no solar
        pvinstallcost=Inf;
    elseif h==hprime
        % Pay the retrofit penalty
        pvinstallcost=1.1*pv_pct_cost*(1+agej_pct_cost)*h;
    else % Changing house
        % PV costs approximately 5% of new house ($30K system for $600K house)
        pvinstallcost=pv_pct_cost*(1+agej_pct_cost)*hprime;
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

if agej<Jr % If working age
    % Note: +h and -hprime; if a is negative, interest is negative
    c=labor*w*kappa_j*z+(1+agej_pct_cost)*(h-hprime);
else % Retirement
    c=pension+(1+agej_pct_cost)*(h-hprime);
end
c=c+(1+r_rate)*a-aprime;
% 3x10kW is full solar for single house; we allow up to 5x10kW
c=c-htc-rentalcosts-(1+agej_pct_cost)*(pvinstallcost+energy_pct_cost*max(h,1)^2*(1-solarpv/3));

if c>0
    F=(((c^(1-sigma_h))*(s^sigma_h))^(1-sigma))/(1-sigma) -psi*(labor^(1+eta))/(1+eta); % The utility function
end


end
