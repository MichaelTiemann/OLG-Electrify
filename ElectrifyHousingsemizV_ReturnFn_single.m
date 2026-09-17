function F = ElectrifyHousingsemizV_ReturnFn_single(installpv, buyhouse, aprime, hprime, a, h, solarpv, ...
    pbefore, pafter, yearsowned, olddownpayment, z, ...
    w, r, sigma, agej, Jr, pension, kappa_j, sigma_h, f_htc, minhouse, rentprice, houseservices, mortgageduration, pv_pct_cost, energy_pct_cost)

% --- 0. Master Tensor Setup ---
master_tensor = a + h + hprime + aprime + solarpv + z + buyhouse + installpv + pbefore + pafter + yearsowned + olddownpayment;
single_1=single(1);
single_minf=single(-Inf);

% --- 1. House and Mortgage Setup ---
relevantdownpayment = 0.2 * buyhouse;
relevantdownpayment(buyhouse == 4) = olddownpayment(buyhouse == 4);

housevalueatpurchase = zeros(size(master_tensor), 'like', master_tensor);

mask_buy4 = (buyhouse == 4);
housevalueatpurchase(mask_buy4) = h(mask_buy4) .* pbefore(mask_buy4);

mask_buy13 = (buyhouse > 0) & (buyhouse < 4);
housevalueatpurchase(mask_buy13) = h(mask_buy13) .* pbefore(mask_buy13) .* pafter(mask_buy13);

outstandingdebt = zeros(size(master_tensor), 'like', master_tensor);
mortgagepayment = zeros(size(master_tensor), 'like', master_tensor);

own_mask = (buyhouse > 0);
originalmortgage = zeros(size(master_tensor), 'like', master_tensor);
originalmortgage(own_mask) = (1 - relevantdownpayment(own_mask)) .* housevalueatpurchase(own_mask);

% Apply fixed exponential compounding (.^ instead of .*)
paying_mask = own_mask & (yearsowned < 20);
rate_factor = (single_1+r).^mortgageduration;
pmt_factor = (r .* rate_factor) ./ (rate_factor - single_1);

mortgagepayment(paying_mask) = originalmortgage(paying_mask) .* pmt_factor;
outstandingdebt(paying_mask) = originalmortgage(paying_mask) .* (rate_factor - (single_1+r).^(yearsowned(paying_mask) + single_1)) ./ (rate_factor - single_1);

% --- 2. Transactions and Costs ---
costofnewhouse = zeros(size(master_tensor), 'like', master_tensor);
move_mask = (hprime ~= h);
costofnewhouse(move_mask) = relevantdownpayment(move_mask) .* pbefore(move_mask) .* pafter(move_mask) .* hprime(move_mask) - outstandingdebt(move_mask);

htc = zeros(size(master_tensor), 'like', master_tensor);
htc(move_mask) = f_htc * pafter(move_mask) .* hprime(move_mask);

pvinstallcost = zeros(size(master_tensor), 'like', master_tensor);
inst_mask = (installpv == 1);

pv_buy0 = inst_mask & (buyhouse == 0);
pvinstallcost(pv_buy0) = single(Inf);

pv_buy13 = inst_mask & (buyhouse > 0) & (buyhouse < 4);
pvinstallcost(pv_buy13) = pv_pct_cost * h(pv_buy13) .* pbefore(pv_buy13);

pv_buy4 = inst_mask & (buyhouse == 4);
pvinstallcost(pv_buy4) = 1.1 * pv_pct_cost * h(pv_buy4) .* max(pbefore(pv_buy4), pafter(pv_buy4));

s = houseservices * h;
rentalcosts = zeros(size(master_tensor), 'like', master_tensor);

h_zero = (h == 0);
s(h_zero) = 0.5 * houseservices * minhouse;
rentalcosts(h_zero) = rentprice;

% --- 3. Budget Constraint (Consumption) ---
if agej < Jr
    c = w * kappa_j * z + (single_1+r)*a - aprime - costofnewhouse - htc - rentalcosts - mortgagepayment - pvinstallcost - (energy_pct_cost * (single_1 - solarpv/30));
else
    c = pension + (single_1+r)*a - aprime - costofnewhouse - htc - rentalcosts - mortgagepayment - pvinstallcost - (energy_pct_cost * (single_1 - solarpv/30));
end

% --- 4. Utility Assembly ---
F = -inf(size(master_tensor), 'like', master_tensor);

valid_c = (c > 0);
F(valid_c) = (((c(valid_c).^(single_1 - sigma_h)) .* (s(valid_c).^sigma_h)).^(single_1 - sigma)) ./ (single_1 - sigma);

% --- 5. Invalid State Bans ---
F((installpv == 1) & (solarpv > 0)) = single_minf;

if agej >= Jr
    F(a < 0) = single_minf;
end

F((hprime == 0) & (buyhouse ~= 0)) = single_minf;
F((hprime ~= 0) & (hprime == h) & (buyhouse ~= 4)) = single_minf;
F((hprime ~= 0) & (hprime ~= h) & ~(buyhouse > 0 & buyhouse < 4)) = single_minf;


end