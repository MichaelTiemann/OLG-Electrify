function prob = ElectrifyHousing_SemiExoStateFn(pbefore,pafter,yearsowned,downpayment,pbeforeprime,pafterprime,yearsownedprime,downpaymentprime,buyhouse, probhousepricerise, probhousepricefall,pbeforespacing, pafterspacing, maxpbefore, minpbefore,maxpafter, minpafter,mortgageduration)

% 1. Setup Tolerance for floating-point coordinate matching
tol = 1e-4;

% 2. Cast all inputs to double for pristine probability math!
pbefore = double(pbefore);
pafter = double(pafter);
pbeforeprime = double(pbeforeprime);
pafterprime = double(pafterprime);
pnew = pbefore * pafter;
pbeforespacing = double(pbeforespacing);
pafterspacing = double(pafterspacing);

probp1 = 0.0;
probp2 = 0.0;
proby  = 0.0;
probd  = 0.0;

% =========================================================================
% First, the probabilities for pbefore (probp1)
% =========================================================================
if buyhouse == 4
    if abs(pbeforeprime - pbefore) < tol
        probp1 = 1.0;
    end
elseif buyhouse == 1 || buyhouse == 2 || buyhouse == 3
    if pnew >= double(maxpbefore) - tol
        if abs(pbeforeprime - double(maxpbefore)) < tol
            probp1 = 1.0;
        end
    elseif pnew <= double(minpbefore) + tol
        if abs(pbeforeprime - double(minpbefore)) < tol
            probp1 = 1.0;
        end
    else
        dist = abs(pbeforeprime - pnew);
        if dist < pbeforespacing + tol
            probp1 = max(0.0, 1.0 - (dist / pbeforespacing));
        end
    end
elseif buyhouse == 0
    if abs((pbeforeprime - pbefore) - pbeforespacing) < tol
        probp1 = probhousepricerise;
    elseif abs(pbeforeprime - pbefore) < tol
        probp1 = 1.0 - probhousepricerise - probhousepricefall;
    elseif abs((pbefore - pbeforeprime) - pbeforespacing) < tol
        probp1 = probhousepricefall;
    end

    if abs(pbefore - double(maxpbefore)) < tol
        if abs(pbeforeprime - pbefore) < tol
            probp1 = 1.0 - probhousepricefall;
        end
    end
    if abs(pbefore - double(minpbefore)) < tol
        if abs(pbeforeprime - pbefore) < tol
            probp1 = 1.0 - probhousepricerise;
        end
    end
end

% =========================================================================
% Second, the probabilities for pafter (probp2)
% =========================================================================
if buyhouse == 0
    if abs(pafterprime - pafter) < tol
        probp2 = 1.0;
    end
elseif buyhouse == 1 || buyhouse == 2 || buyhouse == 3
    if abs(pafterprime - 1.0) < tol
        probp2 = 1.0;
    end
elseif buyhouse == 4
    if abs((pafterprime - pafter) - pafterspacing) < tol
        probp2 = probhousepricerise;
    elseif abs(pafterprime - pafter) < tol
        probp2 = 1.0 - probhousepricerise - probhousepricefall;
    elseif abs((pafter - pafterprime) - pafterspacing) < tol
        probp2 = probhousepricefall;
    end

    if abs(pafter - double(maxpafter)) < tol
        if abs(pafterprime - pafter) < tol
            probp2 = 1.0 - probhousepricefall;
        end
    end
    if abs(pafter - double(minpafter)) < tol
        if abs(pafterprime - pafter) < tol
            probp2 = 1.0 - probhousepricerise;
        end
    end
end

% =========================================================================
% Third, years owned (proby)
% =========================================================================
if buyhouse == 0
    if yearsownedprime == 0
        proby = 1.0;
    end
elseif buyhouse == 1 || buyhouse == 2 || buyhouse == 3
    if yearsownedprime == 0
        proby = 1.0;
    end
elseif buyhouse == 4
    if yearsownedprime == yearsowned + 1
        proby = 1.0;
    end
    if yearsowned == mortgageduration - 1
        if yearsownedprime == 100
            proby = 1.0;
        end
    end
    if yearsowned == 100
        if yearsownedprime == 100
            proby = 1.0;
        end
    end
end

% =========================================================================
% Fourth, downpayment (probd)
% =========================================================================
if buyhouse == 0
    if abs(double(downpaymentprime) - 0.2) < tol
        probd = 1.0;
    end
elseif buyhouse == 1
    if abs(double(downpaymentprime) - 0.2) < tol
        probd = 1.0;
    end
elseif buyhouse == 2
    if abs(double(downpaymentprime) - 0.4) < tol
        probd = 1.0;
    end
elseif buyhouse == 3
    if abs(double(downpaymentprime) - 0.6) < tol
        probd = 1.0;
    end
elseif buyhouse == 4
    if abs(double(downpaymentprime) - double(downpayment)) < tol
        probd = 1.0;
    end
end

% Calculate final probability strictly in double precision
prob = probp1 * probp2 * proby * probd;

end