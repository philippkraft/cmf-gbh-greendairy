# Gladbacherhof Root-Zone Water Balance Model (`gbh_model`)

A plot-scale root-zone soil water balance and groundwater recharge model implemented with the **Catchment Modelling Framework (CMF)** ([Kraft et al., 2011](https://doi.org/10.1016/j.envsoft.2011.02.006)) and **PyTables/HDF5** ensemble workflows for the **Greendairy** experimental fields at **Gladbacherhof**.

---

## Scientific Context

This code supports the joint manuscript investigating nitrate leaching dynamics in crop rotations of an organic dairy farming system, specifically addressing:

1. Quantification of nitrate leaching from soils under crop rotations in an organic dairy farming system.
2. Comparison of nitrate leaching estimates derived from soil $N_{\min}$ sampling against suction cup measurements combined with model-based soil water balance estimates.
3. Identification of dominant nitrogen sources contributing to nitrate in leachate water using $\delta^{15}\text{N}$ natural abundance.

Within this workflow, the hydrological model simulates root-zone water fluxes and groundwater recharge (drainage leaving the bottom of the soil profile) to be coupled with measured nitrate-N concentration time series across an ensemble of behavioural parameter sets.

---

## Repository Structure

| File | Purpose |
| :--- | :--- |
| `gbh_model.py` | Core CMF root-zone water balance model (`Model`), parameter definitions (`Parameters`), and input data loaders (`load_weather`, `load_soilmoisture`, `load_kc`). |
| `gbh_runner.py` | Parallel Latin Hypercube Sampling (LHS) ensemble runner (`Runner`, `run_lhs`) writing simulation results and summary metrics to HDF5 (`gbh_<code>.h5`) via PyTables. |
| `gbh_plot.py` | Post-processing, GLUE-style thresholding, parameter/flux diagnostics (dotty plots, equifinality, validation scatter, soil moisture & boundary flux plots), and Excel exports. |

---

## Model Overview (`gbh_model.py`)

### Spatial Domain & Boundary Conditions
* **Plot-scale cell:** Single CMF cell (`1000 m²` area) with hourly meteorological and rainfall forcing (`cmf.h`).
* **Soil profile:** Six soil layers down to `1.10 m` depth (`0.1, 0.3, 0.5, 0.7, 0.9, 1.1 m` lower boundaries), plus a `macropores` storage (`self.cell.add_storage('macropores')`). *(Note: Ensure the manuscript text matches the 6 layers implemented in `create_nodes()`.)*
* **Boundary outlets:** Dirichlet boundary outlets for groundwater (`gw` at `z = -2 m`) and surface runoff (`runoff` at `z = 0 m`).

### Hydrological Processes
* **Soil water capacity & initial wetness:** Layer porosity is derived from dead water (`DW`, $> 4.2\ \text{pF}$), plant-available field capacity (`aFC`, $2.5\text{–}4.2\ \text{pF}$), and air capacity (`AC`, $< 2.5\ \text{pF}$, with depth-dependent adjustments in the top two intervals), initialized at `wetness = 0.5`.
* **Evapotranspiration:** Driven by hourly weather inputs (`air_temperature`, `rel_humidity`, `wind_speed`, `global_radiation`, `sunshine_fraction`) and daily interpolated crop coefficients (`Kc`) and rooting depths (`RootDepth80`), using `cmf.PenmanMonteithET` and `cmf.ContentStress`.
* **Surface infiltration & macropores:** Surface water routing via `cmf.ConceptualInfiltration` and `cmf.WaterbalanceFlux`, with optional BROOK90-style infiltration partitioning (`infexp`) into `macropores` and deeper layers via `cmf.PartitionFluxRoute`.
* **Vertical percolation (tipping-bucket / power-law):** Layer-to-layer and bottom drainage connections (`cmf.PowerLawConnection`) parameterized by saturated hydraulic conductivity (`Ksat`), field capacity storage (`Sfc`), and percolation exponent (`percexp`), integrated over hourly steps with `cmf.CVodeDense(project, 1e-9)`.

### Model Parameters (`Parameters`)
* `DW`: Dead water (% of soil, $> 4.2\ \text{pF}$) — uniform around mean `18` (`sd = 2`, $\pm 3\ \text{SD}$).
* `aFC`: Plant-available field capacity (% of soil, $2.5\text{–}4.2\ \text{pF}$) — uniform around mean `14.5` (`sd = 1.2`, $\pm 3\ \text{SD}$).
* `AC`: Air capacity (% of soil, $< 2.5\ \text{pF}$) — uniform around mean `9.6` (`sd = 2.5`, $\pm 3\ \text{SD}$).
* `Ksat`: Saturated hydraulic conductivity ($\text{m d}^{-1}$) — uniform on `[0.05, 1.0]`.
* `percexp`: Exponent for power-law percolation — uniform on `[1.0, 5.0]` (default `2.0`).
* `infexp`: Infiltration exponent for BROOK90 macropore infiltration — uniform on `[0, 0.5]` (default `0`).

---

## Input Data Requirements

The scripts expect local Excel input files in the working directory (adjust paths in code if needed):

1. **`weather-2025.xlsx`** (loaded by `load_weather`):
   * Time-indexed hourly weather table. Column names are parsed before `' at '` (lowercased, spaces replaced by `'_'`, periods removed) to match:
     `air_temperature`, `rel_humidity`, `wind_speed`, `global_radiation`, `sunshine_fraction`, and `rainfall`.
2. **`moisture.xlsx`** (loaded by `load_soilmoisture`):
   * Time-indexed soil water content observations (e.g., TEROS 10 sensors at `40 cm` and `60 cm` depth). Even-indexed columns (`::2`) are aggregated as `40 cm` (`d40cm_mean`, `d40cm_std`) and odd-indexed columns (`1::2`) as `60 cm` (`d60cm_mean`, `d60cm_std`), split around `'2024-05-01'`.
3. **`Kc.xlsx`** (loaded by `load_kc`):
   * Sheet `'Kc_Verlauf'` indexed by `'date'`, containing `'Kc'` and `'RootDepth80'`, resampled daily (`'d'`) and linearly interpolated.

---

## Installation & Dependencies

Install `cmf` and the supporting Python packages:

```bash
pip install cmf numpy pandas scipy matplotlib tables spotpy openpyxl
```
