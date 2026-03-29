function [n_d,n_a,n_z,N_j,e,vfoptions]=Electrify_GridSizeSetup(scenario, J, small_z_no_e, small_model, vfoptions)

%% Grid sizes to use for household
if scenario<3
    n_d.household=101;
    n_a.household=201;
    n_z.household=3+2*floor(1.7*log(min(J,60))); % AR(1) with age-dependent params = 15 with 60 periods
    vfoptions.lowmemory.household=0;
else
    if scenario<4
        % Endogenous shares, assets (>=6), housing (>=2), and solarpv (>=2) x10kW PV
        if small_model
            n_d.household=[21,3]; % Decisions: labor, buyhouse (3 w/o PV; 5 w/PV)
            n_a.household=[5,31,2,2]; 
        else
            n_d.household=[21,5]; % Decisions: labor, buyhouse (3 w/o PV; 5 w/PV)
            n_a.household=[5,31,5,7];
        end
        if small_z_no_e
            vfoptions.lowmemory.household=0;
        else
            vfoptions.lowmemory.household=2;
        end
    else
        % Endogenous shares, assets (>=6), car (3), housing (>=2), and solarpv (>=2) x15kW PV
        if small_model
            n_d.household=[21,3]; % Decisions: labor, buyhouse (3 w/o PV; 5 w/PV)
            n_a.household=[9,31,1,2,2];
        else
            n_d.household=[21,5]; % Decisions: labor, buyhouse (3 w/o PV; 5 w/PV)
            n_a.household=[5,31,3,4,5];
        end
        vfoptions.lowmemory.household=3;
    end
    n_z.household=1+2*floor(1.2*log(min(J,60))); % AR(1) with age-dependent params = 7 with 60 periods
end
if small_z_no_e
    n_z.household=1;
    e=0;
else
    % Exogenous labor productivity units shocks (next two lines)
    vfoptions.n_e.household=3; % iid
end
N_j.household=J; % Number of periods in finite horizon

%% Grids to use for firm
if scenario<4
    n_d.firm=101; % Dividend payment
    n_a.firm=201; % Capital holdings
else
    n_d.firm=0; % Not Yet Used: Electrification investment
    n_a.firm=[51,42]; % Capital holdings and PV assets
end
if small_z_no_e
    n_z.firm=1;
else
    n_z.firm=3+2*floor(log(min(J,60))); % Productivity shock; scaled to model, not firm horizon
end
N_j.firm=Inf; % Infinite horizon
vfoptions.lowmemory.firm=logical(scenario==4 && ~small_z_no_e);

%% Grids to use for energy
if scenario<4
    n_d.energy=0; % What decisions?
    n_a.energy=1; % What assets?
else
    n_d.energy=101; % Invest in PV
    n_a.energy=202; % PV assets
end
if small_z_no_e
    n_z.energy=1;
else
    n_z.energy=3+2*floor(log(min(J,60))); % Productivity shock; scaled to model, not firm horizon
end
N_j.energy=Inf; % Infinite horizon
vfoptions.lowmemory.energy=0;


end