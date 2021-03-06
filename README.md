# Automated identification of urban substructure for comparative analysis
by Rohan L. Aras, Nicholas T. Ouellette, and Rishee K. Jain

[![DOI](https://zenodo.org/badge/302712577.svg)](https://zenodo.org/badge/latestdoi/302712577)

This is the code repository for an upcoming paper. Plots from the main manuscript are generated in `manuscript_plots.ipynb`; plots for the supporting information are generated in `s1_synthetic_plots.ipynb`. Supporting code is located underneath `utils`. Data can be found [here](https://drive.google.com/file/d/1dfoinkUviuv-fa8cTw-h3pAHn98fMiDA/view?usp=sharing) as a zip file. To run the code, unzip `data.zip` at the top level of the repository.

## Data Sources
Data is sourced from [City of New York](https://www1.nyc.gov/site/planning/data-maps/open-data/dwn-pluto-mappluto.page), [City of San Francisco](https://data.sfgov.org/Housing-and-Buildings/Land-Use/us3s-fp9q), [City of Houston](https://cohegis.houstontx.gov/cohgispub/rest/services/PD/), and [US Census databases](https://catalog.data.gov/dataset/tiger-line-shapefile-2013-county-harris-county-tx-all-roads-county-based-shapefile). Data for Houston was subset using a custom shapefile tracing the 610 loop in [geojson.io](http://geojson.io/#map=10/29.7793/-95.4135). Please note that the exact ArcGIS endpoint  used to collect data in `data.zip` (`/EconomicDevelopment_wm/MapServer/9`) appears to have been deprecated by the City of Houston Planning Department. A similar endpoint is now available at `/Cadastral_wm/MapServer/0`.
