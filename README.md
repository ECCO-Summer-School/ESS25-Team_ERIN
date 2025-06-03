# st<ins>E</ins><ins>R</ins>ic <ins>I</ins><ins>N</ins>vestigations (ERIN): How well does ECCO capture the steric sea-level signal?
### Collaborators

| Name | Personal goals |
| ------------- | ------------- |
| Erin Robson | To gain experience in accessing and analysing ECCO data products for sea-level research

## Files and folders in project repository

* **`notebooks/`**
<br> notebook contains the core workflow
* **`scripts/`**
<br> functions and relevant .mat files
* `README.md`
<br> project description

## Introduction

One of the key components of sea-level variability is the steric contribution, which arises from changes in ocean volume due to variations in temperature (thermosteric) and salinity (halosteric). Over the 20<sup>th</sup> and early-21<sup>st</sup> century, steric sea-level rates have risen from slightly negative values around 1900, to almost 2 mm yr<sup>-1</sup> in 2021. Noteably, since the 1970s, when examining the individual components that contribute to global mean sea-level (GMSL) variability, the steric component has emerged as the dominant driver of the observed GMSL acceleration (Dangendorf et al., 2024). 

Given its increasing importance, the overarching aim of this project is to **compare the steric sea-level signal derived from the ECCO v4r4 ocean state estimate to the steric signal from the most recent GMSL reconstruction by Dangendorf et al. (2024)**, in order to assess the consistency between model-driven and observation-driven estimates.

### The problem
Studies (e.g., Dangendorf et al., 2021) have highlighted limitations in understanding historical and projected sea-level changes in coastal regions. These limitations stem both from constraints of available observations (e.g., the inability of satellite altimetry to resolve conditions near the coast), and from uncertainties owing to the fact that the impacts of large-scale ocean dynamics are not well constrained. While steric processes are the main driver of sea-level variability in the open ocean, this contribution diminishes at the coast due to shallower open depths. As such, coastal sea-level variability cannot be directly inferred from open-ocean steric changes. Many observation-driven estimates rely on tide-gauge records to provide measurements of relative sea-level at the coast, in an attempt to better constrain coastal variability. As the reconstruction from Dangendorf et al. (2024) utilises tide-gauge data, whereas the ECCO estimate is primarily open-ocean model-based, I wanted to investigate how the inclusion (or exclusion) of coastal data influences the characterisation of the steric sea-level signal. 

## Data and Methods

### Data

* **ECCO version 4 release 4 (v4r4) output via the matlab toolbox gcmfaces** (https://github.com/MITgcm/gcmfaces). The v4r4 datasets provide global coverage from 1992-2017 at a monthly timestep.
* **GMSL reconstruction from Dangendorf et al. (2024**) can be accessed via https://zenodo.org/records/10621070, which provides the sea-level fields for each of the contributing components for both global and regional-scale reconstructions at an annual timestep (1900-2021). 

### Methodology/work flow

- [x] Isolate steric SL signal from both ECCO (SSL = absolute sea level [SSHDYN] - (manometric - IB [OBPNOPAB])) and the reconstruction (SDSL = RSOI - BaryGRD - GIA - IB)
- [x] Regrid ECCO output (360 x 360 x 312) onto the same grid as D24 (74742 x 26)
- [x] Compute linear trend of each signal (1992-2017) 
- [x] Run various comparison analyses (e.g., pattern/temporal correlations, RMSE, trend differences)
- [x] Subdivide grid by region and run area-weighted correlations 
- [x] Nearest neighbour to compare signals at tide-gauge locations

## Project Results
The one-slide summary of this project is available here: [ERIN_ECCO_summaryslide.pdf](https://github.com/user-attachments/files/20576634/ERIN_ECCO_summaryslide.pdf)

##
* Dangendorf, S., Frederikse, T., Chafik, L., Klinck, J.M., Ezer, T. and Hamlington, B.D. (2021). Data-driven reconstruction reveals large-scale ocean circulation control on coastal sea level. Nature Climate Change, 11(6), pp.514–520. doi:https://doi.org/10.1038/s41558-021-01046-1.
* Dangendorf, S., Sun, Q., Wahl, T., Thompson, P., Mitrovica, J.X. and Hamlington, B. (2024). Probabilistic reconstruction of sea-level changes and their causes since 1900. Earth system science data, 16(7), pp.3471–3494. doi:https://doi.org/10.5194/essd-16-3471-2024.

