function F = ElectrifyHousingsemizV_ReturnFn_single(installpv, buyhouse, aprime, hprime, a, h, solarpv, ...
    pbefore, pafter, yearsowned, olddownpayment, z, ...
    w, r, sigma, agej, Jr, pension, kappa_j, sigma_h, f_htc, minhouse, rentprice, houseservices, mortgageduration, pv_pct_cost, energy_pct_cost)

single_1 = single(1);
single_minf = single(-Inf);

% --- 1. House and Mortgage Setup ---
% Pure Arithmetic Masking: No pre-allocation of zeros. MATLAB's JIT will 
% implicitly expand these tensors purely in the GPU compute registers.

relevantdownpayment = (buyhouse == 4) .* olddownpayment + (buyhouse < 4) .* (single(0.2) .* buyhouse);

housevalueatpurchase = (buyhouse == 4) .* (h .* pbefore) + ...
                       (buyhouse > 0 & buyhouse < 4) .* (h .* pbefore .* pafter);

originalmortgage = (buyhouse > 0) .* (single_1 - relevantdownpayment) .* housevalueatpurchase;

paying_mask = (buyhouse > 0) & (yearsowned < 20);
rate_factor = (single_1 + r).^mortgageduration;
pmt_factor  = (r .* rate_factor) ./ (rate_factor - single_1);

mortgagepayment = paying_mask .* originalmortgage .* pmt_factor;
outstandingdebt = paying_mask .* originalmortgage .* ...
                  (rate_factor - (single_1 + r).^(yearsowned + single_1)) ./ (rate_factor - single_1);

% --- 2. Transactions and Costs ---
move_mask = (hprime ~= h);
costofnewhouse = move_mask .* (relevantdownpayment .* pbefore .* pafter .* hprime - outstandingdebt);
htc = move_mask .* (f_htc .* pafter .* hprime);

pv_buy13 = (installpv == 1) & (buyhouse > 0 & buyhouse < 4);
pv_buy4  = (installpv == 1) & (buyhouse == 4);
pvinstallcost = pv_buy13 .* (pv_pct_cost .* h .* pbefore) + ...
                pv_buy4  .* (single(1.1) .* pv_pct_cost .* h .* max(pbefore, pafter));

h_zero = (h == 0);
s = (~h_zero) .* (houseservices .* h) + h_zero .* (single(0.5) .* houseservices .* minhouse);
rentalcosts = h_zero .* rentprice;

% --- 3. Budget Constraint (Consumption) ---
energy_cost = energy_pct_cost .* (single_1 - solarpv ./ single(30));

if agej < Jr
    c = w .* kappa_j .* z + (single_1 + r).*a - aprime - costofnewhouse - htc - rentalcosts - mortgagepayment - pvinstallcost - energy_cost;
else
    c = pension + (single_1 + r).*a - aprime - costofnewhouse - htc - rentalcosts - mortgagepayment - pvinstallcost - energy_cost;
end

% --- 4. Utility Assembly ---
% Preallocate F based strictly on the implicitly expanded size of c
F = -inf(size(c), 'like', c);

valid_c = (c > 0);

% 1. Protect the GPU JIT from evaluating fractional powers of negative numbers
% (We clamp it to a tiny positive number. The invalid states will be overwritten by -Inf anyway).
c_safe = max(c, single(1e-8));

% 2. Calculate utility globally. 
% By omitting the (valid_c) index here, MATLAB's JIT natively broadcasts 
% the smaller 's' tensor against the massive 'c_safe' tensor!
utility = (((c_safe.^(single_1 - sigma_h)) .* (s.^sigma_h)).^(single_1 - sigma)) ./ (single_1 - sigma);

% 3. Slot the valid results into the pre-allocated F tensor
F(valid_c) = utility(valid_c);

% --- 5. Invalid State Bans ---
F((installpv == 1) & (solarpv > 0)) = single_minf;
if agej >= Jr
    F(a < 0) = single_minf;
end
F((hprime == 0) & (buyhouse ~= 0)) = single_minf;
F((hprime ~= 0) & (hprime == h) & (buyhouse ~= 4)) = single_minf;
F((hprime ~= 0) & (hprime ~= h) & ~(buyhouse > 0 & buyhouse < 4)) = single_minf;
F((installpv == 1) & (buyhouse == 0)) = single_minf;

end