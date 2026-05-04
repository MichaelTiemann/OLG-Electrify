function [n_d,n_a,n_z,N_j,vfoptions]=Electrify_GridSizeSetup(scenario, J, small_z_no_e, small_model, vfoptions)

%% Grid sizes to use for household
if scenario<3
    n_d.household=101;
    n_a.household=201;
    n_z.household=3+2*floor(1.7*log(min(J,60))); % AR(1) with age-dependent params = 15 with 60 periods
    vfoptions.lowmemory.household=0;
else
    if scenario<4
        % Endogenous shares+assets (>=6), housing (>=2), and solarpv (>=2) x5kW PV
        if small_model
            n_d.household=[21,3]; % Decisions: labor, buyhouse (3 w/o PV; 5 w/PV)
            n_a.household=[151,2,2]; 
        else
            n_d.household=[21,5]; % Decisions: labor, buyhouse (3 w/o PV; 5 w/PV)
            n_a.household=[151,5,7];
        end
        if small_z_no_e
            vfoptions.lowmemory.household=0;
        elseif small_model
            vfoptions.lowmemory.household=1;
        else
            vfoptions.lowmemory.household=2;
        end
    else
        % Endogenous shares, assets (>=6), car (3), housing (>=2), and solarpv (>=2) x5kW PV
        if small_model
            n_d.household=[21,3]; % Decisions: labor, buyhouse (3 w/o PV; 5 w/PV)
            n_a.household=[95,1,2,2]; % note: car fixed at zero
        else
            n_d.household=[21,5]; % Decisions: labor, buyhouse (3 w/o PV; 5 w/PV)
            n_a.household=[151,3,4,5];
        end
        vfoptions.lowmemory.household=3;
    end
    n_z.household=1+2*floor(1.2*log(min(J,60))); % AR(1) with age-dependent params = 7 with 60 periods
    vfoptions.experienceasset.household=1;

end
if small_z_no_e
    n_z.household=1;
else
    % Exogenous labor productivity units shocks (next two lines)
    vfoptions.n_e.household=3; % iid
end
N_j.household=J; % Number of periods in finite horizon

%% Grids to use for firm
if scenario<4
    n_d.firm=101; % Dividend payment
    n_a.firm=201; % Capital holdings
    vfoptions.experienceasset.firm=0;
else
    n_d.firm=5; % Per-period PV investment (0, 1, 2, 3, or 4 PV arrays per period)
    n_a.firm=[51,n_d.firm(end)]; % Capital holdings and incremental PV (experience) assets
    vfoptions.experienceasset.firm=1;
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
    vfoptions.experienceasset.energy=0;
else
    n_d.energy=11; % Per-period PV investment (0-10 PV arrays per period)
    n_a.energy=[n_a.firm(1),n_d.energy(end)]; % Capital holdings and incremental PV (experience) assets
    vfoptions.experienceasset.energy=1;
end
if small_z_no_e
    n_z.energy=1;
else
    n_z.energy=3+2*floor(log(min(J,60))); % Productivity shock; scaled to model, not firm horizon
end
N_j.energy=Inf; % Infinite horizon
vfoptions.lowmemory.energy=0;

end