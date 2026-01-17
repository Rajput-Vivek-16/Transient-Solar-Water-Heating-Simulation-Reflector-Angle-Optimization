## Transient Solar Water Heating Simulation & Reflector Angle Optimization (Regression-Based)

### Overview

This project performs a **transient simulation of a solar water heating system** to evaluate the effect of **reflector angle (ϕ)** on system performance. **Polynomial regression is first applied to experimental irradiance and ambient temperature data** to obtain continuous time-dependent inputs. These regressed profiles are then used in a transient energy balance model to predict water temperature and compute the **second figure of merit (F₂)**. The goal is to **identify the reflector angle that maximizes the average F₂**, indicating optimal thermal performance.

---

### Key Objectives

* Apply **regression analysis** to experimental solar irradiance and ambient temperature data
* Generate smooth, continuous functions ( I(t) ) and ( T_a(t) ) from discrete measurements
* Solve the **transient energy balance equation** using ODE integration
* Compute **dynamic F₂ values** and determine the reflector angle yielding maximum average F₂

---

### Methodology

1. **Regression of Experimental Data**

   * Time-series data sampled at 10-minute intervals
   * Second-order polynomial regression (`polyfit`) applied to:

     * Solar irradiance ( I(t) )
     * Ambient temperature ( T_a(t) )
   * Regression ensures smooth input functions for stable transient ODE solving and avoids numerical discontinuities

2. **Transient Thermal Model**

   * Water temperature evolution obtained by solving:
     [
     M C_p \frac{dT_w}{dt} = A_i \left[\tau\alpha I(t) - U_L (T_w - T_a)\right]
     ]
   * Reflector geometry affects the effective absorber area ( A_i ), which is a function of reflector angle ( \phi )

3. **Figure of Merit Calculation**

   * First figure of merit:
     [
     F_1 = \frac{\tau\alpha}{U_L}
     ]
   * Second figure of merit ( F_2 ) computed dynamically between successive time steps using regressed ( I(t) ) and ( T_a(t) )
   * Mean ( F_2 ) evaluated via numerical integration (`trapz`) over the entire operating period

4. **Optimization of Reflector Angle**

   * Reflector angle ( \phi ) varied from **40° to 240°**
   * For each angle:

     * Effective area recalculated
     * Transient temperature predicted
     * Average ( F_2 ) computed
   * Optimal reflector angle identified based on **maximum average F₂**

---

### Inputs

* Experimental transient water temperature data
* Solar irradiance and ambient temperature data (used for regression)
* System parameters:

  * Water mass (M)
  * Specific heat capacity (Cp)
  * Absorber dimensions (a, b)
  * Overall heat loss coefficient (( U_L ))
  * Optical efficiency (( \tau\alpha ))

---

### Outputs

* Regressed profiles of:

  * Solar irradiance ( I(t) )
  * Ambient temperature ( T_a(t) )
* Optimal reflector angle ( \phi_{opt} )
* Maximum average second figure of merit ( F_2 )
* Plot of **Average F₂ vs Reflector Angle**

---

### Results

The code reports:

```text
Maximum avg F2 = <value> at phi = <value>°
```

and generates a plot illustrating the dependence of mean F₂ on reflector angle.

---

### MATLAB Requirements

* MATLAB R2018b or later
* Toolboxes:

  * ODE Suite (`ode45`)
  * Curve fitting utilities (`polyfit`, `polyval`)

---

### Applications

* Regression-based solar thermal performance analysis
* Reflector geometry optimization in concentrating solar collectors
* Transient modeling of solar water heating systems
