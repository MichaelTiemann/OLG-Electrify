# Design of an an Energy Transition model using Over-Lapping Generations (OLG)

With the understanding that the VFI Toolkit's Transition Paths are designed to model the evolution of _Prices_ in equilibrium based on exogenous _Parameters_ defining the assumptions of the transition path,
we no longer attempt to _discover_ an optimal transition path; rather, we _analyze the implications_ of the transition path.

## Price Paths

We want to learn the effects of a given transition plan on the following prices

* **w**: the wages paid to labor so that markets clear
* **r**: the interest rate at which capital markets clear
* **dividend**: the dividend rate paid by firms to shareholders
* **pension**: the fraction of income paid to retirees
* **G**: Government/public services as a fraction of total economy
* **housing**: cost of a unit of housing (initially 4x **w**)

## Parameter Paths

MBIE provide a matrix of energy supply and demand that distributes across households and commercial uses by energy type.
We also have rough statistics on the CO2 footprint of various energy / industry uses.  We can project all of the above into:

* **energy use by source**: household and firm usage of green electricity, brown electricity, petrol, gas, coal, and other
* **energy cost by source**: projected prices of green/brown electricity, petrol, gas, coal, and other across the next 50 years (solar becoming cheaper, petrol becoming more dear)
* **energy use by sector**: simplifying as much as possible, sectors include households, agricultural, industrial, commercial, government.
* **cpi**: estimate of the consumer price index (ex energy) over time
* **carbon costs**: whether priced by the ETS or by social cost of carbon (actual cost computed in 2024 as $2450/tonne), we expect this cost to rise over time.  The energy source dictates the carbon burden, which does not change over time.
* **commercial transport electrification**: percentage of electric/zero-emission vehicles to total vehicle fleet
* **$s_j$**: probability of surviving from one lifecycle period to the next
* **${\mu}_j$**: population density of agents by age

We allocate transportation energy use 50/50 between households (with light passenger vehicles) and firms (with commercial vehicles, ships, planes, etc).
We do not (yet) estimate in any detail how actual transportation usage (passenger km, goods tonnes km) may change over time, but do estimate how electrification may change the energy requirements (and GHG footprint) of transportation over time.
We give households the choice of petrol vs. evs at each stage of the observed transition.  The remainder of the transport transition is given by the Parameter Path.

## Open Questions (a growing list)

### Effect of balance of Imports vs. Exports

When we consider the economy of New Zealand, which is very dependent on foreign trade, how does the balance of exports and imports affect our transition modeling?
New Zealand spends $20B/year on energy imports, which must be made up in export goods.
If and when renewable energy such as solar and wind replace $10B of that (or more), what does this mean for households (offering labor, consuming goods, enjoying leisure, etc)?
And what does this mean for firms and the pressure on them to produce and compete?
More specifically, what prices, parameters, and/or equilibrium conditions do we need to add or adjust as we shift up to 5% of our GDP from foreign outflows to domestic uses?

### Equity prices and house ownership in a transition model

The model from which we are building solves a general equilibrium of households and firms as they seek their respective optimums in terms of labor participation, capital investment, wages, dividends, stock price, etc.
While we can calculate the stationary distributions of the agents and their optimal policies based on their respective value functions for a single equilibrium, the model was not designed to operate with history.
One particularly thorny problem is how to manage the value of equity in the firms over time (or the contintued ownership of houses by households).
It is likely that to handle this we need to add a third type of variable to the Transition Path solution: a transitional _experience asset_.

The VFI Toolkit already supports the concept of experience assets (with and without uncertainty), but their meaning is solidly within the context of a stationary distribution and general equilibrium, and not something that lives across transition periods.

Agent-based simulation systems get around this problem by literally instantiating agents that can be endowed with assets.  The VFI Toolkit works at a higher level of abstraction than that.
