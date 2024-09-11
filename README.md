# Mt-Info
This module calculates basic mountain parameters such as average mountain length, width, elevation, Moho depth and erosion rate for three modern collisional mountains, the Himalaya-Tibetan plateau, Alps and Zagros. The details to calculate these parameters are provided in the Supplementary Information of Zhu et al. (2024, in review).

 # Datasets
- [Natural Earth’s physical vectors](https://www.naturalearthdata.com/downloads/50m-physical-vectors/50m-physical-labels/) (version 4.1.0)
- [Global topography data from ETOPO1](https://www.ncei.noaa.gov/products/etopo-global-relief-model) at a spatial resolution of one arc-minute (Amante and Eakins, 2009).
- [Global Moho depth data from CRUST1.0](https://igppweb.ucsd.edu/~gabi/crust1.html) at one-degree resolution (Laske et al., 2013).

# Packages
- [geopandas](https://geopandas.org/)
- [shapely](https://shapely.readthedocs.io/en/stable/manual.html) (package for computational geometry)
- [rasterio](https://rasterio.readthedocs.io/en/latest/intro.html) (For accessing the many different kind of raster data files used in the GIS field)

# References
- Amante, C. and B.W. Eakins, 2009. ETOPO1 1 Arc-Minute Global Relief Model: Procedures, Data Sources and Analysis. NOAA Technical Memorandum NESDIS NGDC-24. National Geophysical Data Center, NOAA. doi:10.7289/V5C8276M [2022]
- Laske, G., Masters., G., Ma, Z. and Pasyanos, M., Update on CRUST1.0 - A 1-degree Global Model of Earth's Crust, Geophys. Res. Abstracts, 15, Abstract EGU2013-2658, 2013.
