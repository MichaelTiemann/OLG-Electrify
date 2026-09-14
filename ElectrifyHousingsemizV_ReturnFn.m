function F = ElectrifyHousingsemizV_ReturnFn(installpv, buyhouse, aprime, hprime, a, h, solarpv, ...
    pbefore, pafter, yearsowned, olddownpayment, z, ...
    w, r, sigma, agej, Jr, pension, kappa_j, sigma_h, f_htc, minhouse, rentprice, houseservices, mortgageduration, pv_pct_cost, energy_pct_cost)

% --- 0. Master Tensor Setup ---
% 1. Add all inputs together to force MATLAB to find the max broadcast size
master_tensor = a + h + hprime + aprime + solarpv + z + buyhouse + installpv + pbefore + pafter + yearsowned + olddownpayment;

% 2. EXPLICITLY BROADCAST ALL VARIABLES TO FULL TENSOR SIZE
% This guarantees that all logical masks will perfectly align with pre-allocated full-size arrays.
zero_tensor = zeros(size(master_tensor), 'like', master_tensor);
installpv      = installpv + zero_tensor;
buyhouse       = buyhouse + zero_tensor;
aprime         = aprime + zero_tensor;
hprime         = hprime + zero_tensor;
a              = a + zero_tensor;
h              = h + zero_tensor;
solarpv        = solarpv + zero_tensor;
pbefore        = pbefore + zero_tensor;
pafter         = pafter + zero_tensor;
yearsowned     = yearsowned + zero_tensor;
olddownpayment = olddownpayment + zero_tensor;
z              = z + zero_tensor;

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
rate_factor = (1+r).^mortgageduration;
pmt_factor = (r .* rate_factor) ./ (rate_factor - 1);

mortgagepayment(paying_mask) = originalmortgage(paying_mask) .* pmt_factor;
outstandingdebt(paying_mask) = originalmortgage(paying_mask) .* (rate_factor - (1+r).^(yearsowned(paying_mask) + 1)) ./ (rate_factor - 1);

% --- 2. Transactions and Costs ---
costofnewhouse = zeros(size(master_tensor), 'like', master_tensor);
move_mask = (hprime ~= h);
costofnewhouse(move_mask) = relevantdownpayment(move_mask) .* pbefore(move_mask) .* pafter(move_mask) .* hprime(move_mask) - outstandingdebt(move_mask);

htc = zeros(size(master_tensor), 'like', master_tensor);
htc(move_mask) = f_htc * pafter(move_mask) .* hprime(move_mask);

pvinstallcost = zeros(size(master_tensor), 'like', master_tensor);
inst_mask = (installpv == 1);

pv_buy0 = inst_mask & (buyhouse == 0);
pvinstallcost(pv_buy0) = Inf;

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
    c = w * kappa_j * z + (1+r)*a - aprime - costofnewhouse - htc - rentalcosts - mortgagepayment - pvinstallcost - (energy_pct_cost * (1 - solarpv/30));
else
    c = pension + (1+r)*a - aprime - costofnewhouse - htc - rentalcosts - mortgagepayment - pvinstallcost - (energy_pct_cost * (1 - solarpv/30));
end

% --- 4. Utility Assembly ---
F = -inf(size(master_tensor), 'like', master_tensor);

valid_c = (c > 0);
F(valid_c) = (((c(valid_c).^(1 - sigma_h)) .* (s(valid_c).^sigma_h)).^(1 - sigma)) ./ (1 - sigma);

% --- 5. Invalid State Bans ---
F((installpv == 1) & (solarpv > 0)) = -Inf;

if agej >= Jr
    F(a < 0) = -Inf;
end

F((hprime == 0) & (buyhouse ~= 0)) = -Inf;
F((hprime ~= 0) & (hprime == h) & (buyhouse ~= 4)) = -Inf;
F((hprime ~= 0) & (hprime ~= h) & ~(buyhouse > 0 & buyhouse < 4)) = -Inf;


end