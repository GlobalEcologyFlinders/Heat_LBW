# Assessing the impacts of heat on low birth weight: evidence from low- and middle-income countries  
<img align="right" src="www/low birth weight transp.png" alt="low birth weight" width="200" style="margin-top: 20px">

**Lead investigator**
- Dr Syeda Hira Fatima, [Global Ecology | <em>Partuyarta Ngadluku Wardli Kuu</em>, Flinders University, Adelaide, Australia](https://globalecologyflinders.com/people/#SHF), *Email: syeda.fatima@flinders.edu.au*  

**Collaborators**
- Professor Corey J. A. Bradshaw, [Global Ecology | <em>Partuyarta Ngadluku Wardli Kuu</em>, Flinders University, Adelaide, Australia](https://globalecologyflinders.com/people/#DIRECTOR), *Email: corey.bradshaw@flinders.edu.com*
- Associate Professor Zohra Lassi, [Research Profile](https://researchers.adelaide.edu.au/profile/zohra.lassi)  
- Professor Zulfiqar Ali Bhutta, [Research Profile](https://www.sickkids.ca/en/staff/b/zulfiqar-bhutta/)  
- Professor Peng Bi, [Research Profile](https://researchers.adelaide.edu.au/profile/peng.bi)  
- Dr Jai Das, [Research Profile](https://www.aku.edu/mcpk/faculty/Pages/profile.aspx?ProfileID=307&Name=Jai++Das)  
- Associate Professor Salima MeherAli, [Research Profile](https://apps.ualberta.ca/directory/person/meherali)  

## Data  
This repository contains the code and documentation for analysing national datasets from [Demographic Health Surveys](https://dhsprogram.com) and [Multiple Indicator Cluster Surveys](https://mics.unicef.org/surveys), merged with geospatial gridded meteorological data from [ERA5-Land](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-monthly-means?tab=overview), to assess the impacts of hot weather conditions on low birth weight in Pakistan. This study is the most comprehensive assessment of heat impacts on low birth weights for low- and middle-income countries.

## Methodological overview  
We used a flexible modelling approach to account for uncertainties and spatio-temporally sparse data by applying distributed-lag nonlinear models within a space-time series study design. We employed generalised linear mixed-effects models with random intercepts to account for provincial-level variation and to explore delayed and nonlinear associations between ambient temperatures and low birth weight.  

Main steps:  
- quadratic functions  applied to address nonlinearity;  
- model uncertainty accounted for by predicting model-averaged outputs;  
- province-level estimates of the heat-related attributable fraction derived using weighted coefficients;  
- for policy implications and targeted interventions, we combined province-level estimates with district-level health and environmental indicators (e.g., live births, child mortality, poverty index, and heat index) to calculate a district-level **heat vulnerability index**;  
- subgroups of women at the population level with higher risk of heat-related low birth weights identified.  
