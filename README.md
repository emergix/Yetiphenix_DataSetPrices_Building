# YetiPhoenix Autocalls – Instrument Overview and Machine Learning Pricing Goal

## 1. The Instruments: YetiPhoenix Autocalls

A **YetiPhoenix autocall** is a type of structured product combining features of *autocallable notes* and *Phoenix notes*.  
They are usually linked to one or more underlying assets (equities, indices, or baskets) and have the following characteristics:

### Key Features

1. **Autocall (Early Redemption) Feature**
   - On scheduled observation dates, if the underlying(s) price is above a given *autocall barrier*, the product redeems early, paying back the notional plus a coupon.
   - The observation can be based on:
     - The **worst-of** performance among several underlyings.
     - A single underlying's price.

2. **Phoenix Coupon Feature**
   - On each observation date (whether autocall is triggered or not), a coupon may be paid if the underlying is above a *coupon barrier*.
   - Coupons can be **memory**-enabled: if missed previously but conditions are later met, they are paid retroactively.

3. **Protection / Capital at Risk**
   - At maturity, if the product has not been called early:
     - If the underlying is above the *protection barrier* (a.k.a. down-and-in level), 100% of the notional is repaid.
     - If it is below, capital loss occurs, usually in proportion to the underlying's negative performance.

4. **Possible Variations**
   - Basket type: **worst-of** (most common), **best-of**, or weighted average.
   - Payoff legs: additional PDI (put down-and-in) or exotic coupon legs.
   - Observation frequency: quarterly, semiannual, or annual.

---

## 2. Pricing Challenges

Pricing YetiPhoenix autocalls is complex because:
- They are **path-dependent**: autocall and coupon events depend on the trajectory of the underlying(s).
- They often involve **multiple correlated underlyings** (basket worst-of).
- Volatility smiles/skews and stochastic dynamics influence the probabilities of hitting barriers.
- They can require **multi-dimensional integration** or **Monte Carlo** simulation.

Two common modelling approaches:

1. **Local Volatility (VolLoc) Models**
   - Fit to the market smile via Dupire's formula.
   - Deterministic volatility as a function of spot and time.
   - Captures the implied vol surface exactly for vanillas, but may misestimate path-dependent option values.

2. **Stochastic Volatility (VolSto) Models**
   - Example: SABR, Heston.
   - Volatility is itself a random process, adding skew/kurtosis to the path distribution.
   - Often gives more realistic dynamics for barrier and autocall events.

---

## 3. Machine Learning Objective

We want to create a **large dataset** of prices for these products, computed under **both** modelling approaches:
- Price under **VolLoc** (local volatility) model → \( \sigma_{loc} \)
- Price under **VolSto** (stochastic volatility) model → \( \sigma_{sto} \)

Then compute the **spread** between them:

- **Absolute spread:**  
  \[ y = \sigma_{sto} - \sigma_{loc} \]

- **Relative spread:**  
  \[ y = \frac{\sigma_{sto} - \sigma_{loc}}{\sigma_{loc}} \]

### Why?
- The spread is **systematic**: certain market conditions make VolLoc consistently over/under-price relative to VolSto.
- If we can train a **neural network** to predict the spread from market data and product parameters, we can:
  - Quickly adjust a VolLoc price to approximate the more expensive VolSto price.
  - Avoid running slow stochastic simulations for each trade.

---

## 4. Dataset Generation Process

1. **Scenario Sampling**
   - Randomly generate market states: forward prices, interest rates, dividends, correlations, and smile parameters.
   - Randomly generate product parameters: barriers, coupon rates, observation dates, basket composition.

2. **Dual Pricing**
   - For each scenario:
     1. Compute price under **VolLoc**.
     2. Compute price under **VolSto**.

3. **Data Storage**
   - Save to CSV: market features + product features + vol_loc + vol_sto + spread.
   - Ensure positive-definite correlation matrices and realistic market parameters.

4. **Machine Learning Training**
   - Input: market & product features.
   - Target: spread (absolute or relative).
   - Use model to *correct* fast VolLoc prices to match VolSto.

---

## 5. Expected Benefits

- **Speed**: Replace costly stochastic simulations with a fast neural net inference.
- **Accuracy**: Maintain the realism of stochastic-vol prices for path-dependent payoffs.
- **Scalability**: Price large books of structured products quickly in risk systems.

---

*This document explains the context and purpose of the dataset generation using the two programs: one for YetiPhoenix autocalls under local volatility, and one under stochastic volatility.*
