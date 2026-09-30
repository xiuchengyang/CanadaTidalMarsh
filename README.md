# Canadian tidal marsh dynamics, 1988–2023

Explore the [interactive Canada tidal marsh map](https://xiucheng.projects.earthengine.app/view/canadamarsh).

The statistics in this repository are based on our Canadian tidal-marsh maps, produced using [DECODE v2](https://github.com/xiuchengyang/DECODE_v2.0) and available to download from [Google Drive](https://drive.google.com/drive/folders/1TWbaCa5TZByWJ0GrmmRfGs0s8JqU5NVM?usp=sharing). We combined these maps with protected and conserved area (PCA) boundaries and Global Human Modification (GHM) layers to examine marsh extent and change in relation to protection and human pressure. The resulting tables are in `Statistics`, covering individual pixels, PCAs, basins, ecoregions, and both coasts.

The MATLAB scripts in `Code` use these tables to reproduce the analyses for our Canada-wide study and save the figures in `Analysis_Figures`.

## Statistics

### Area-related

Annual tidal-marsh extent and trends by coast, ecoregion, basin, and latitude.

| File | Data |
| --- | --- |
| `Annual_Area_by_Coast_1988_2023.xlsx` | Original and adjusted annual area (ha), with standard deviations of adjusted area, for the Atlantic and Pacific coasts. |
| `Annual_Area_by_Ecoregion_1988_2023.xlsx` | Annual area (ha) by ecoregion, with ecoregion names and abbreviations. |
| `Basin_Area_and_Trends_1988_2023.xlsx` | Annual area (ha), trend rates, confidence intervals, acceleration, and trend categories by basin. |
| `Pixel_Counts_by_Latitude_1988_2023.xlsx` | Annual marsh-pixel counts in 0.1° latitude bands, separately for each coast and both combined. |
| `Spatial_boundaries/Basin_lev06.shp` | Boundaries of the 131 HydroBASINS level-6 basins containing mapped tidal marsh. |
| `Spatial_boundaries/Meow_Ecos_CAMarsh.shp` | Boundaries of the nine marine ecoregions used in the study. |

### GHM-related

Marsh extent and change records linked to Global Human Modification values.

| File | Data |
| --- | --- |
| `Change_Pixels_with_GHM.mat` | Marsh-change pixels with locations, change years, change categories, and matched GHM values. |
| `Marsh_Pixels_with_GHM.mat` | GHM values at marsh pixels for 1990–2020 at five-year intervals, identified by year and coast. |
| `Change_Pixels_with_GHM_by_Basin_1986_2024.csv` | Marsh-change pixels with basin IDs, locations, change years and categories, pixel areas, and matched GHM values and years. |
| `Basin_GHM_Statistics_1990_2020.csv` | Mean and 75th-percentile GHM values and sampled-pixel counts by basin and GHM year. |
| `Annual_Area_by_Basin_1985_2024.csv` | Annual mapped marsh area (km²) by basin; also used for PCA comparisons. |

### PCAs-related

Marsh area overlapping protected and conserved areas, with establishment dates and boundaries.

| File | Data |
| --- | --- |
| `Annual_Area_by_PCA_and_Basin_1985_2024.csv` | Annual marsh area (ha) within each PCA–basin overlap, with coast, ecoregion, and establishment-year attributes. |
| `CPCAD_TidalMarsh_1988_2023/CPCAD_TidalMarsh_1988_2023.shp` | Marsh-overlapping PCA boundaries and attributes for the IUCN I–IV groups used in the category analysis. |

## Software and source data

The code was tested in MATLAB R2026a with the Statistics and Machine Learning Toolbox and Mapping Toolbox. Bundled functions and their licenses are in `Code/Dependencies`.

| Source | Data used and reference |
| --- | --- |
| [Global Human Modification](https://zenodo.org/records/14449495) | GHM v3 temporal layers, 300 m, 1990–2020 at five-year intervals. [Theobald et al. (2025), *Scientific Data*](https://doi.org/10.1038/s41597-025-04892-2). |
| [Canadian Protected and Conserved Areas Database](https://www.canada.ca/en/environment-climate-change/services/national-wildlife-areas/protected-conserved-areas-database.html) | Environment and Climate Change Canada, 2024 release: PCA boundaries, establishment years, and IUCN categories. |
| [HydroBASINS](https://www.hydrosheds.org/products/hydrobasins) | North America, level 6, version 1c. [Lehner and Grill (2013), *Hydrological Processes*](https://doi.org/10.1002/hyp.9740). |
| [Marine Ecoregions of the World](https://databasin.org/datasets/3b6b12e7bcca419990c9081c0af254a2/) | Coastal ecoregion boundaries. [Spalding et al. (2007), *BioScience*](https://doi.org/10.1641/B570707). |

## Citations

**Canadian tidal-marsh study:** Yang, X., McHenry, J., Zhu, Z., Murray, N. J., Li, M., Valenti, V., Park, A., Qiu, S., Kroeger, K. D., O’Connor, M. I., Baum, J. K., and Knox, S. H. *Accelerating tidal marsh loss in Canada reveals the limits of static area-based protection.* Manuscript under review.

**DECODE v2:** Yang, X., Knox, S. H., Qiu, S., Kroeger, K. D., Zhu, Z., Covington, S., and Zhu, Z. (2026). *Tracking US Tidal Marsh Extent and Change with Dense Landsat Time Series and Comparisons with Existing Products.* Remote Sensing of Environment (under review).

## Contact

[Xiucheng Yang](https://xiuchengyang.github.io/).
