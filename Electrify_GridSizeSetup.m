function [n_d,n_a,n_z,N_j,vfoptions]=Electrify_GridSizeSetup(scenario, J, small_z_no_e, small_model, vfoptions)

%% Grid sizes to use for household
if small_z_no_e
    n_z.household=1;
    n_e.household=1; % not put into vfoptions, so just scales lowmemory calculations
else
    if scenario<3
        n_z.household=3+2*floor(1.7*log(min(J,60))); % AR(1) with age-dependent params = 15 with 60 periods
    else
        n_z.household=1+2*floor(1.2*log(min(J,60))); % AR(1) with age-dependent params = 7 with 60 periods
    end
    % Exogenous labor productivity units shocks (next two lines)
    n_e.household=3;
    vfoptions.n_e.household=n_e.household; % iid
end
if scenario<3
    n_d.household=101;
    n_a.household=201;
    vfoptions.lowmemory.household=0;
else
    if scenario<4
        % Endogenous shares+assets (>=6), housing (>=2), and solarpv (>=2) x5kW PV
        if small_model
            n_d.household=[21,3]; % Decisions: labor, buyhouse (3 w/o PV; 5 w/PV)
            n_a.household=[95,2,2]; 
        else
            n_d.household=[21,5]; % Decisions: labor, buyhouse (3 w/o PV; 5 w/PV)
            n_a.household=[151,5,7];
        end
    else
        % Endogenous shares+assets (>=6), car (3), housing (>=2), and solarpv (>=2) x5kW PV
        if small_model
            n_d.household=[41,5]; % Decisions: labor, buyhouse (3 w/o PV; 5 w/PV)
            n_a.household=[41,3,3,3];
        else
            n_d.household=[61,5]; % Decisions: labor, buyhouse (3 w/o PV; 5 w/PV)
            n_a.household=[60,3,5,6]; % 31, 46, 60, 75, 90, 104, 119, 134, or 148
        end
    end
    vfoptions.experienceasset.household=1;
    numel_lowmem3=prod(n_d.household)*prod(n_a.household(1:end-1))^2;
    if numel_lowmem3 < 2^31
        numel_lowmem2=numel_lowmem3*n_a.household(end);
        if numel_lowmem2 < 2^31
            numel_lowmem1=numel_lowmem2*prod(n_z.household);
            if numel_lowmem1 < 2^31
                numel_lowmem0=numel_lowmem1*prod(n_e.household);
                if numel_lowmem0 < 2^31
                    vfoptions.lowmemory.household=0;
                else
                    vfoptions.lowmemory.household=1;
                end
            else
                vfoptions.lowmemory.household=2;
            end
        else
            vfoptions.lowmemory.household=3;
        end
        fprintf("vfoptions.lowmemory.household = %d \n", vfoptions.lowmemory.household)
    else
        error("Model size exceeds GPU maximum variable size");
    end
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
    n_z.firm=3;
else
    n_z.firm=11; % Productivity shock
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
    n_z.energy=11; % Productivity shock
end
N_j.energy=Inf; % Infinite horizon
vfoptions.lowmemory.energy=0;

end