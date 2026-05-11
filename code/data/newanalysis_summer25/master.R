# 1. Install required packages
install.packages(c("tidycensus", "tidyverse", "scales", "ggrepel", 
                   "stargazer", "lfe", "viridis", "patchwork", 
                   "gganimate", "sf", "plotly"))

# 2. Set your Census API key
#census_api_key("", install = TRUE)

# 3. Run the main analysis
source("/Users/tommasodesanto/Desktop/Projects/Fertility/Codes/R/newanalysis_summer25/fertility_urban_analysis.R")