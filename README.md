# Sample Project

This is an example of how teams can structure their project repositories and format their project README.md file.

When creating a project repository from this template choose "Public" so other participants can follow progress. Add a "topic" to your repository details (click on the gear icon next to the "About" section on the repository page) to help others find your work (e.g. `ecco-hackweek-2024`).


## Files and folders in your project repository

This template provides the following suggested organizaiton structure for the project repository, but each project team is free to organize their repository as they see fit.

* **`contributors/`**
<br> Each team member can create their own folder under contributors, within which they can work on their own scripts, notebooks, and other files. Having a dedicated folder for each person helps to prevent conflicts when merging with the main branch. This is a good place for team members to start off exploring data and methods for the project.
* **`notebooks/`**
<br> Notebooks that are considered delivered results for the project should go in here.
* **`scripts/`**
<br> Code that is shared by the team should go in here (e.g. functions or subroutines). These will be files other than Jupyter Notebooks such as Python scripts (.py).
* `.gitignore`
<br> This file sets the files that will be globally ignored by `git` for the project. (e.g. you may want git to ignore temporary files or large data files, [read more about ignoring files here](https://docs.github.com/en/get-started/getting-started-with-git/ignoring-files))
* `environment.yml`
<br> `conda` environment description needed to run this project.
* `README.md`
<br> Description of the project (see suggested headings below)
* `model-card.md`
<br> Description (following a metadata standard) of any machine learning models used in the project

# stERic INvestigations (ERIN): How well does ECCO capture the steric sea-level signal?

## Introduction

One of the key components of sea-level variability is the steric contribution, which arises from changes in ocean volume due to variations in temperature (thermosteric) and salinity (halosteric). Over the 20th and early-21st century, steric sea-level rates have risen from slightly negative values around 1900, to almost 2 mm yr-1 in 2021. Noteably, since the 1970s, when examining the individual components that contribute to global mean sea-level (GMSL) variability, the steric component has emerged as the dominant driver of the observed GMSL acceleration (Dangendorf et al., 2024). 

Given its increasing importance, the overarching aim of this project was to compare the steric sea-level signal derived from the ECCO v4r4 ocean state estimate to the steric signal from the most recent GMSL reconstruction by Dangendorf et al. (2024), in order to assess the consistency between model-driven and observation-driven estimates of steric variability. 

### Collaborators

| Name | Personal goals |
| ------------- | ------------- |
| Erin Robson | To gain experience in accessing and analysing ECCO data products for sea-level research

### The problem
Studies (e.g., Dangendorf et al., 2021) have highlighted limitations in understanding historical and projected sea-level changes in coastal regions. These limitations stem both from constraints of availible observations (e.g., the inability of satellite altimetry to resolve conditions near the coast), and from uncertainties owing to the fact that the impacts of large-scale ocean dynamics are not well constrained. While steric processes are the main driver of sea-level variability in the open ocean, this contribution diminishes at the coast due to shallower open depths. As such, coastal sea-level variability cannot be directly inferred from open-ocean steric changes. Many observation-driven estimates rely on tide-gauge records to provide measurements of relative sea-level at the coast, in an attempt to better constrain coastal variability. As the reconstruction from Dangendorf et al. (2024) utilises tide-gauge data, whereas the ECCO estimate is primarily open-ocean model-based, I wanted to investigate how the inclusion (or exclusion) of coastal data influences the characterisation of the steric sea-level signal. 

## Data and Methods

### Data

I used the ECCO version 4 release 4 (v4r4) output via the matlab toolbox gcmfaces (https://github.com/MITgcm/gcmfaces). The v4r4 datasets provide global coverage from 1992-2017 at a monthly timestep. The data for the GMSL reconstruction from Dangendorf et al. (2024) was accessed via https://zenodo.org/records/10621070, which provides the sea-level fields for each of the contributing components for both global and regional-scale reconstructions at an annual timestep (1900-2021). 

### Methods/tools

I isolated the steric signal from ECCO (SSL = absolute sea level [SSHDYN] - (manometric - IB [OBPNOPAB])) and the reconstruction (SDSL = RSOI - BaryGRD - GIA - IB), redgridded the ECCO output (360 x 360 x 312) onto the same grid as D24 (74742 x 26), and calculated the linear trend of each signal for 1992-2017 before running various comparison analyses.

Used nearest-neighbour interpolation to find the closest point in the ECCO grid to the tide-gauge locations, in order to compare the steric signals directly on the coast. 

## Project Results
[ERIN_ECCO_summaryslide.pptx](https://github.com/user-attachments/files/20576222/ERIN_ECCO_summaryslide.pptx)

  ## References
Dangendorf, S., Frederikse, T., Chafik, L., Klinck, J.M., Ezer, T. and Hamlington, B.D. (2021). Data-driven reconstruction reveals large-scale ocean circulation control on coastal sea level. Nature Climate Change, 11(6), pp.514–520. doi:https://doi.org/10.1038/s41558-021-01046-1.

Dangendorf, S., Sun, Q., Wahl, T., Thompson, P., Mitrovica, J.X. and Hamlington, B. (2024). Probabilistic reconstruction of sea-level changes and their causes since 1900. Earth system science data, 16(7), pp.3471–3494. doi:https://doi.org/10.5194/essd-16-3471-2024.

