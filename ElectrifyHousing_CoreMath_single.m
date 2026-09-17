function F = ElectrifyHousing_CoreMath_single(installpv, buyhouse, aprime, hprime, a, h, solarpv, ...
    pbefore, pafter, yearsowned, olddownpayment, z, ...
    w, r, sigma, agej, Jr, pension, kappa_j, sigma_h, f_htc, minhouse, rentprice, houseservices, mortgageduration, pv_pct_cost, energy_pct_cost)

% --- 0. Strictly Typed Constants ---
single_1 = single(1);
single_minf = single(-Inf);
single_0 = single(0);
single_02 = single(0.2);

% --- 1. Invalid State Bans (Short-Circuit to save GPU compute!) ---
if (installpv == 1) && (solarpv > 0)
    F = single_minf; return;
end
if (agej >= Jr) && (a < 0)
    F = single_minf; return;
end
if (hprime == 0) && (buyhouse ~= 0)
    F = single_minf; return;
end
if (hprime ~= 0) && (hprime == h) && (buyhouse ~= 4)
    F = single_minf; return;
end
if (hprime ~= 0) && (hprime ~= h) && ~(buyhouse > 0 && buyhouse < 4)
    F = single_minf; return;
end
if (installpv == 1) && (buyhouse == 0)
    F = single_minf; return;
end

% --- 2. Old House & Equity Cashed Out ---
if h == single_0
    old_house_purchase_value = single_0;
    old_house_current_value  = single_0;
    outstandingdebt          = single_0;
    equity_from_old          = single_0;
else
    old_house_purchase_value = h * pbefore;
    old_house_current_value  = h * pbefore * pafter;

    old_originalmortgage = (single_1 - olddownpayment) * old_house_purchase_value;
    rate_factor = (single_1 + r)^mortgageduration;
    pmt_factor  = (r * rate_factor) / (rate_factor - single_1);

    if yearsowned < mortgageduration
        outstandingdebt = old_originalmortgage * (rate_factor - (single_1 + r)^(yearsowned + single_1)) / (rate_factor - single_1);
    else
        outstandingdebt = single_0;
    end

    % If moving, liquidate the old house
    if hprime ~= h
        equity_from_old = old_house_current_value - outstandingdebt;
    else
        equity_from_old = single_0;
    end
end

% --- 3. New House & Current Mortgage ---
if buyhouse == 4 % Holding current house
    relevantdownpayment = olddownpayment;
    current_house_purchase_value = old_house_purchase_value;
    current_yearsowned = yearsowned;
elseif buyhouse > 0 % Buying a new house (buyhouse = 1, 2, or 3)
    relevantdownpayment = single_02 * buyhouse;
    current_house_purchase_value = hprime * pbefore * pafter;
    current_yearsowned = single_0; % Reset years owned
else % Renting
    relevantdownpayment = single_0;
    current_house_purchase_value = single_0;
    current_yearsowned = single_0;
end

if buyhouse > 0
    current_originalmortgage = (single_1 - relevantdownpayment) * current_house_purchase_value;
    if current_yearsowned < mortgageduration
        % Recalculate pmt_factor in case we bypassed it above
        rate_factor = (single_1 + r)^mortgageduration;
        pmt_factor  = (r * rate_factor) / (rate_factor - single_1);
        mortgagepayment = current_originalmortgage * pmt_factor;
    else
        mortgagepayment = single_0;
    end
else
    mortgagepayment = single_0;
end

% --- 4. Transactions, PV, and Housing Services ---
if hprime ~= h && buyhouse > 0 && buyhouse < 4
    cash_for_new_downpayment = relevantdownpayment * hprime * pbefore * pafter;
else
    cash_for_new_downpayment = single_0;
end

costofnewhouse = cash_for_new_downpayment - equity_from_old;

if hprime ~= h
    htc = f_htc * pafter * hprime;
else
    htc = single_0;
end

if installpv == 1
    if buyhouse > 0 && buyhouse < 4
        % PV cost scales with the NEW house
        pvinstallcost = pv_pct_cost * hprime * pbefore * pafter;
    elseif buyhouse == 4
        pvinstallcost = single(1.1) * pv_pct_cost * h * max(pbefore, pafter);
    else
        pvinstallcost = single(Inf); 
    end
else
    pvinstallcost = single_0;
end

% Housing services and rent based strictly on the house you live in THIS period
if hprime == single_0
    s = single(0.5) * houseservices * minhouse;
    rentalcosts = rentprice;
else
    s = houseservices * hprime;
    rentalcosts = single_0;
end

% --- 5. Budget Constraint (Consumption) ---
energy_cost = energy_pct_cost * (single_1 - solarpv / single(30));

if agej < Jr
    c = w * kappa_j * z + (single_1 + r)*a - aprime - costofnewhouse - htc - rentalcosts - mortgagepayment - pvinstallcost - energy_cost;
else
    c = pension + (single_1 + r)*a - aprime - costofnewhouse - htc - rentalcosts - mortgagepayment - pvinstallcost - energy_cost;
end

% --- 6. Utility Assembly ---
if c <= 0
    F = single_minf;
else
    F = (((c^(single_1 - sigma_h)) * (s^sigma_h))^(single_1 - sigma)) / (single_1 - sigma);
end


end