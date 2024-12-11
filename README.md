# Assessing the impacts of heat on low birth weight: evidence from low- and middle-income countries  
<a href="https://doi.org/10.5281/zenodo.14373384"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.14373384.svg" alt="DOI"></a>
<img align="right" src="www/low birth weight transp.png" alt="low birth weight" width="200" style="margin-top: 20px">

**Lead investigator**
- Dr [Syeda Hira Fatima](https://globalecologyflinders.com/people/#SHF), Global Ecology | <em>Partuyarta Ngadluku Wardli Kuu</em>, Flinders University, Adelaide, Australia, [Email](mailto:syeda.fatima@flinders.edu.au)  

**Collaborators**
- Professor [Corey J. A. Bradshaw](https://globalecologyflinders.com/people/#DIRECTOR), Global Ecology | <em>Partuyarta Ngadluku Wardli Kuu</em>, Flinders University, Adelaide, Australia, [Email](mailto:corey.bradshaw@flinders.edu.com)
- Associate Professor [Zohra Lassi](https://researchers.adelaide.edu.au/profile/zohra.lassi) 
- Professor [Zulfiqar Ali Bhutta](https://www.sickkids.ca/en/staff/b/zulfiqar-bhutta/)
- Professor [Peng Bi](https://researchers.adelaide.edu.au/profile/peng.bi)
- Dr [Jai K. Das](https://www.aku.edu/mcpk/faculty/Pages/profile.aspx?ProfileID=307&Name=Jai++Das)
- Associate Professor [Salima MeherAli](https://apps.ualberta.ca/directory/person/meherali)<br>
<br>
accompanies paper:<br>
<br>
Fatima, SH, CJA Bradshaw, ZA Bhutta, P Bi, JK Das, S MeherAli, ZS Lassi. Disproportionate climate burden of rising temperature on low birth weights in a health-compromised country<br>

### Abstract
Climate change is a major threat to global health, with profound implications for vulnerable populations, especially pregnant women and their newborns. In Pakistan — a country already grappling with socio-economic challenges — the impacts of climate change on maternal health are high. Low birth weight (< 2.5 kg) is a leading contributor to neonatal health problems, heavily burdening Pakistan’s healthcare system and society. Infants born with low birth weight have higher mortality and long-term health impairments, including stunted growth and cognitive deficits — a situation worsened by the escalating effects of climate change. We aimed to evaluate the impact of extreme weather conditions on low birth weight and identify high-risk areas for heat exposure across Pakistan. Additionally, we assessed the projected future risk of low birth weight associated with extreme weather conditions. We used a flexible modelling approach to account for uncertainties and spatiotemporally sparse data by applying distributed-lag nonlinear models within a space-time series study design. We employed generalised linear mixed-effects models with random intercepts to account for  provincial-level variation to explore delayed and nonlinear associations between ambient temperatures and low birth weight. We applied quadratic functions to address nonlinearity, and accounted for model uncertainty by predicting model-averaged outputs. We derived province-level estimates of the heat-related attributed fraction using weighted coefficients. For policy implications and targeted interventions, we combined province-level estimates with district-level health and environmental indicators, such as live births, child mortality, poverty index, and heat index, to calculate a district-level heat vulnerability index. We also identified subgroups of women at the population level who are at higher risk of heat-related low birth weights. Our analysis of national survey data revealed that 28.5% of children in Pakistan (~ 3.81 million from 2008–2017) are born with low birth weight. The cumulative risk of low birth weight due to extreme heat ranges between 1.51 (1.08–2.13; 95% confidence interval) and 2.10 (1.32–3.36) across provinces, highlighting regional disparities. Women of childbearing age in southern Punjab and the northern districts of Baluchistan and Sindh have a higher risk of adverse heat-related effects on the birth weight of their children. The heat-related attributable fraction ranged from 8.2–13.8%, translating to approximately 405,493 cases from 2008 to 2017. Moreover, this fraction is projected to increase further, (1.2–3.0%) with future climate change. This study is the most comprehensive assessment of heat impacts on low birth weights for low-middle-income countries. Our results highlight the vulnerability of this large population to climate change and underscores the need for strong policy support and international cooperation to address the multifaceted challenges posed by climate change on maternal and child health in Pakistan and other low- and middle-income countries.

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

## Scripts
- <code>main_analysis.R</code>: main R code for running models
- <code>LowBirthWeightFigs.R</code>: R code to prepare data for figures in paper

## Required R libraries
<code>data.table</code>, <code>dlnm</code>, <code>dplyr</code>, <code>ggplot2</code>, <code>glmmTMB</code>, <code>here</code>, <code>lubridate</code>, <code>MASS</code>, <code>mice</code>, <code>MuMIn</code>, <code>patchwork</code>, <code>readr</code>, <code>sf</code>, <code>splines</code>, <code>tidyverse</code>, <code>viridis</code>, <code>weathermetrics</code>
<br>
<br>
<p><a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Stacked_RGB_Master.png" alt="Flinders University" width="90" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL" width="85" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://github.com/FutureChildHealth"><img align="bottom-left" src="www/FCHlogo06122024transp.png" alt="Future Child Health" width="90" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://www.sickkids.ca"><img align="bottom-left" src="www/sickkids-logo.webp" alt="Hospital for Sick Children" width="120" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://www.adelaide.edu.au/"><img align="bottom-left" src="www/UAlogo.png" alt="U Adelaide" width="70" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://www.aku.edu/"><img align="bottom-left" src="www/agakhanlogo.png" alt="Aga Khan University" width="70" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://www.ualberta.ca/"><img align="bottom-left" src="www/UAlblogo.png" alt="U Alberta" width="70" style="margin-top: 20px"></a></p>
