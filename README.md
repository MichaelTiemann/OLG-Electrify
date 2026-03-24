# OLG-Electrify -- An Overlapping Generations (OLG) model for electrification

This repository contains VFIToolkit MATLAB models for analyzing the economic evolution of energy transition.
The principal model's starting point is a pair of demonstration files written by Robert Kirkby (author of the VFIToolkit, see https://github.com/vfitoolkit/VFIToolkit-matlab).
One model demonstrates portfolio choice (buying houses vs. risk-free assets) using LifeCycle analysis ([LifeCycleModel35.m](https://github.com/vfitoolkit/IntroToLifeCycleModels/blob/main/Models31to35/LifeCycleModel35.m)).
The other solves a general equilibrium across overlapping generations of households and firms ([OLGModel14.m](https://github.com/vfitoolkit/IntroToOLGModels/blob/main/OLGModel14.m)).
Together, they have all the precursors necessary to evaluate the economics of households and firms deciding to electrify over a period of time.

Lifecycle analysis and overlapping generations allow us to evaluate demographic impacts by age (or age cohort).
Transition path analysis allows us to evaluate economic evolution over time.
Value Function Iteration allows us to use powerful GPUs to solve equilibrium and optimization problems that present as a snapshot in time (across ages) or across time (and across ages).
The VFI Toolkit makes it easy to separate economic cohorts (households, firms, and energy suppliers) so that they can interact in fine-grained ways without creating impractically large state- or action-spaces.

The Reserve Bank of New Zealand ([studying interest rates](https://www.rbnz.govt.nz/news-and-events/news/2025/05/research-investigates-why-the-natural-interest-rate-has-fallen-in-new-zealand-over-recent-decades))
and the UK Government's Office of Budgetary Responsibility ([analyzing pension funding and demographic changes](https://obr.uk/box/the-uk-overlapping-generations-model-uk-olg/))
 now both use Kirkby's VFIToolkit for production modeling of economic forecasts.  Neither of these models incorporate energy in any way.
The Australian Government has also developed an OLG model for their economy ([OLG-A](https://treasury.gov.au/publication/p2023-437296-olga)) that does contemplate elements of energy transition as an economic consideration.  Good on them!

Finally, though our initial examples follow the literature and use Cobb-Douglas production functions (optimizing the marginal contribution of capital and labor), we use an energy-aware production function published by Australian economist [Steve Keen](https://en.wikipedia.org/wiki/Steve_Keen).
Professor Keen argues forcefully that [the classic neoliberal treatment of energy as a rival industry within an economy is wrong](https://profstevekeen.substack.com/p/the-role-of-energy-in-economics), but it can be fixed by moving energy terms from "inside" the economy to astride the economy:

Cobb-Douglas: $`Y = A * K^{1-\alpha} * L^{\alpha}`$

"Wrong" energy-aware Cobb-Douglas: $`Y = A * K^{1-\alpha-\beta} * L^{\alpha} * E^{\beta}`$

Keen's corrected production function: $`Y = (E_k * e_k) * K^{\alpha} * L^{1-{\alpha}}`$

At this point in the model's development, the focus has been on the basic wiring of the model, calibrating parameters only to the extent that the model produces stable equilibria showing economic and programmatic connectivity.  Once we have validated the scope and connectivity of the parameters, we can go about calibrating them and truly analyzing the model.

This project is Open Source Software, and we invite both peer review and participation in that spirit and tradition.

We now present the three major model components: households, firms, the energy sector, and a government that provides abstract benefits and services for the various taxes it levies (payroll, dividends, and capital gains) via the solution to various general equilibrium equations.

## Households

Households in the model make two fundamental "decisions":
* their labor participation rate (from 0% to 100%)
* buying/selling houses (including options to upsize/downsize between several size choices)

They can also make asset allocation choices that cover:
* if they own or buy a house, whether or not to install, retrofit, or upgrade the solar solarPV system (with multiple size choices)
* whether to buy/sell/trade a car, with options for petrol-based or EV
* investing savings into shares of firms (described below), keeping them at a bank (which pays interest on deposits)
* taking out loans collateralized by houses.

Buying and selling houses incur housing transaction costs, and borrowing money exposes a credit wedge between loan and deposit rates.
Households that do not own houses rent.  Renters cannot install their own solarPV systems, but house owners may (with 0, 15kW, 30kW, 45kW and 60kW generation options).
solarPV systems degrade, albeit slowly (1% per year).

Examples from VFIToolkit repositories (linked above, and also even more [here](https://github.com/vfitoolkit/VFItoolkit-matlab-examples)) give rich example
 of many other household questions that can be studied:
* labor participation differences between single male, single female, and married households
* housing purchases and risky asset purchases based on both risk appetite/aversion and exogenous shocks
* the choice to pursue higher education, have children, or both
* investments in personal health, and precautionary saving for medical expenses
* tracking the unpaid labor of childcare and elder care

Alas, while these and many other problems are very interesting, we believe the aspects most salient to energy transition must focus on the elements that are most affected by energy, which come down to powering personal transport, household combustion, and household electricity use.

### Rough parameterization:

Households earn a "standard" wage of 1, which is influenced in the usual way by "work experience" (abbreviated as $`\kappa_j`$), ranging from 0.5 to 2.0,
and also influenced by exogenous $`z`$ and idiosyncratic $`e`$ shocks.  We note that there are three parameterizations to consider for linking this abstract wage to realities in New Zealand:
the living wage (approximately NZD \$55,000/a), the median wage (NZD \$70,000/a), and the average wage (NZD \$80,000).

Houses of size 1, 2, 3, and 4 provide increasingly valuable "housing services" that factor into the consumption function.
Renters achieve 50% of the "housing services" satisfaction of a size 1 house owner).
They are priced at 4x the unit wage based on their size (using the \$80,000 average wage as a base, the range is from NZD \$320,000 to NZD \$1.28 million).
Large houses use more energy in proportion to the square root of their size.

EVs are priced at 50% of the standard wage and petrol cars at 25%.  Cars can be traded/sold for 50% of their value.
Car maintenance costs are priced at 2% of the standard wage.  (We could scale up based on $`\kappa_j`$.)

Housing also incurs annual costs: renters pay 30% of their salary for an all-inclusive deal.
Homeowners spend a maintenance cost of 1% of their house's cost (which amounts to a range of 4% of the standard wage for a house of size 1 and 16% of the standard wage for a house of size 4).

All households are responsible for energy costs, both transportation and housing.  Households without a car budget 20% of the basic wage for public transporation, hiding that energy cost.  But if they have a car, they pay based on the respective energy source.

Renters pay the same energy costs as a house of size 1, and energy costs for larger houses are scaled by the square root of the size of the house.  If the household owns a house, it can install a PV array.  Households are paid the market rate for energy they generate, meaning that sub-sufficient households earn an offset to their energy use, self-sufficient households can zero out their electricity costs, and super-sufficient households actually earn additional money.

Households pay a carbon tax based on obvious carbon expenditures (a petrol vehicle) and on unoffset energy costs based on the exogenous variables *energy_pct_brown* and *carbon_tax*.
*energy_pct_brown* is a function of the energy transition model.  Initially it is 80% and drops to 5% after the full model is finished (100 year horizon).
*carbon_tax* is initially based on a very lightweight ETS pricing model (NDZ \$35/tonne CO2e) and evolves over the course of the model to the true social cost of carbon (NZD \$2450/tonne).

Finally, because survival from age *j* to period *j+1*, some of the population dies within the model.  These result in "accidental bequests" of shares and/or houses.  Following other models, these bequests are valued in the final period by the deceased with a *warmglow_of_bequest*.

## Firms

In the original household/firm model, firms had one decision to make: the size of the dividend they offered.  After discussion, it was decided that trying to use a general equilibrium solver to steer this "decision" was very inefficient.  Instead, the decision is made programmatically.

As mentioned at the outset, we use Steve Keen's modified Cobb-Douglas production model, and based on calibrated $`\alpha`$ parameters from other literature, the firm computes the labor required for a given capital allocation. From there we compute the usual: output, a fraction of energy needed to produce that output, and from that energy amount, a computation of CO2e emissions and thus a cost based on the carbon tax.

In addition to managing its capital stock, firms can decide to purchase industry-scale PVs.
In New Zealand, while most electricity is sourced from renewable energy, energy is just a fraction of energy required.
In fact, only 45% of energy use in New Zealand is from renewable sources.
Using 2023 data, this means that of the 125 TWh of primary energy used by New Zealand for all purposes, 69 TWh is non-renewable.
The model splits this as 20 TWh used by firms and the remainder used by households, the energy service provider, and the government.
With apologies to wind turbines (which are wonderful complements to solarPV arrays), we use PVs as shorthand for net-new renewable resources, and we use recent pricing from industrial-scale installations for initial calibration.

We constrain firms to buy only as much PV resource as required to completely offset their own energy use.
Any non-offset energy use is costed at the prevailing rate of power, and a carbon tax is paid based on total use times the combine factors of *brown_energy_pct* and *carbon_tax* .

After doing the usual tax and expense calculations, the firm is left with the joint operation of paying a dividend and selling shares to raise any necessary monies to pay that dividend as follows (edited and simplified):

```
shares_to_sell=0;
mid_dividend_pp=0.2;
dividend=shares_to_sell+(profit-tau_corp*T)-invest-capitaladjcost;
if dividend<0
    % We will issue new shares and provide a discounted dividend
    low_dividend=0.1;
    % We will completely make up the shortfall (including a paradoxically negative dividend) by selling shares
    shares_to_sell=low_dividend-dividend;
    dividend=low_dividend;
elseif dividend<mid_dividend
    % We will issue new shares and provide a full dividend
    shares_to_sell=mid_dividend-dividend;
    dividend=mid_dividend;
else
    % We don't need to issue shares and can pay rich dividend
end
```

Needless to say, energy and carbon costs have both simple effects on dividend and share allocation, but also strategic impact on total production.

In our model, each time the firm prices its dividend, households can see and act on that price vs everything else the household must decide.

## Energy

The energy sector of the model is mostly a sketch at this point.  It receives revenues from households and firms, and it sees the carbon taxes that households and firms pay.

The revenues received by sector can be translated to energy demand based on per-sector pricing (households pay more than 2x the firm rate).

The energy supplier can also make transition investment decisions, buying renewable generators in time to meet demand forecasts.
In time, it will be the energy sector's transition speed that will govern *brown_energy_pct* that is presently just a model parameter.
When this is done, we can see how accelerating the exposure of the true social cost of carbon further accelerates transition.

## Government benefits, services, and taxes

As mentioned above, we have an invisible government sector.
That is, we don't explicitly model it as an entity unto itself.
But all entities are exposed to taxation, and each entity may receive some explicit or implicit benefit.

As tempting as it would be to put a parlimentary budget negotiation into the model, we have settled on a simpler structure:

* wage taxes on labor pay for the pensions of retirees
* capital gains and dividend taxes pay for government services
* carbon taxes pay for energy transition

## Statistical References and Data

The following are data collected informally for parameterizing the model.  Better sources and methods gladly accepted!

```
%% Basic statistical abstract (NZD)
% NZ GDP: $440B ($80K per capita, $152K per employed worker)
% NZ Wages: $55K living, $70K median, $80K average * 2.9M workers = $232B wages
% NZ Energy:    - 525 PJ/year
%   Oil         - 270 PJ
%   Electricity - 144 PJ
%   Gas         -  58 PJ
%   Biomass     -  40 PJ
%   Coal        -  18 PJ
% NZ Electricity retail: $350/MWh
% NZ HH Energy: 20 kWh/day electricity =>  7 MWh/year =>  $2500/year => 3.5% wages
% NZ HH Energy: 73 kWh/day overall     => 27 MWh/year => $10000/year => 14.0% wages
% NZ Firm Energy retail: $150/MWh
%   Transport    - 200 PJ
%   Industrial   - 160 PJ
%   Commercial   -  55 PJ
%   Ag,Forest,Fish- 30 PJ
%   Total: 445 PJ => 125 TWh => $20B energy costs => 4.5% of 440B GDP
% Energy is 8.6% of Labor costs
% Net capital stocks of NZ $1,329B less $690B real estate = $630B
% K/L = $440B/232B = 1.9
```

### Energy cost and supply allocations

Just as we imagine households able to use their own roofs to generate power (offsetting costs and potentially generating additional income),
we imagine firms able to use their facilities for locally sourcing electricity.  The model contemplates virtually complete electrification (including transportation),
but it allows firms (which may include transport companies) to purchase power from the energy sector rather than being forced to install their own PVs.
Depending on the settings, there may be a role for energy suppliers to generate and sell electricity, or it may be that every energy consumer prefers to be self-reliant.

```
% If Y is full GDP ($440B), then Ek=125 TWh and ek=$440B/125TWh=$3.52 GDP/kWh
% NZ emissions from energy sector: 76.4 Mt CO2e
% 69 TWh to be electrified (56 TWh already renewable); need 46,000 MW generation
% Choosing unit industrial PV to be 200GWh/year generation...
% 200GWh/year = 133MW*1500h/yr = $220M cost @ $1.65M/MW; $220M/$440B = 0.0005 max GDP
% 200GWh PV/year * 1000 MWh/GWh * $150/MWh = Firm PV Energy Cost Offset $30M/PV/year (vs $440B)
% 69 TWh/year to electrify = 345*200GWh/year * $220M/200GWh/year = $75900M cost of 100% conversion
```

At present we allocate firms with the choice of installing 0-100 industry-scale PVs, and the energy sector to install the remaining 245 as a choice between 0-300 industry-scale PVs.
