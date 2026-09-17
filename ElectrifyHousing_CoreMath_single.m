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

% --- 2. House and Mortgage Setup ---
if buyhouse == 4
    relevantdownpayment = olddownpayment;
    housevalueatpurchase = h * pbefore;
elseif buyhouse > 0 && buyhouse < 4
    relevantdownpayment = single_02 * buyhouse;
    housevalueatpurchase = h * pbefore * pafter;
else
    relevantdownpayment = single_0;
    housevalueatpurchase = single_0;
end

if buyhouse > 0
    originalmortgage = (single_1 - relevantdownpayment) * housevalueatpurchase;
else
    originalmortgage = single_0;
end

if (buyhouse > 0) && (yearsowned < 20)
    rate_factor = (single_1 + r)^mortgageduration;
    pmt_factor  = (r * rate_factor) / (rate_factor - single_1);
    mortgagepayment = originalmortgage * pmt_factor;

    debt_factor = (rate_factor - (single_1 + r)^(yearsowned + single_1)) / (rate_factor - single_1);
    outstandingdebt = originalmortgage * debt_factor;
else
    mortgagepayment = single_0;
    outstandingdebt = single_0;
end

% --- 3. Transactions and Costs ---
if hprime ~= h
    costofnewhouse = relevantdownpayment * pbefore * pafter * hprime - outstandingdebt;
    htc = f_htc * pafter * hprime;
else
    costofnewhouse = single_0;
    htc = single_0;
end

if installpv == 1
    if buyhouse > 0 && buyhouse < 4
        pvinstallcost = pv_pct_cost * h * pbefore;
    elseif buyhouse == 4
        pvinstallcost = single(1.1) * pv_pct_cost * h * max(pbefore, pafter);
    else
        pvinstallcost = single(Inf); % Catch-all (though mostly blocked by rules above)
    end
else
    pvinstallcost = single_0;
end

if h == 0
    s = single(0.5) * houseservices * minhouse;
    rentalcosts = rentprice;
else
    s = houseservices * h;
    rentalcosts = single_0;
end

% --- 4. Budget Constraint (Consumption) ---
energy_cost = energy_pct_cost * (single_1 - solarpv / single(30));

if agej < Jr
    c = w * kappa_j * z + (single_1 + r)*a - aprime - costofnewhouse - htc - rentalcosts - mortgagepayment - pvinstallcost - energy_cost;
else
    c = pension + (single_1 + r)*a - aprime - costofnewhouse - htc - rentalcosts - mortgagepayment - pvinstallcost - energy_cost;
end

% --- 5. Utility Assembly ---
if c <= 0
    F = single_minf;
else
    F = (((c^(single_1 - sigma_h)) * (s^sigma_h))^(single_1 - sigma)) / (single_1 - sigma);
end


end