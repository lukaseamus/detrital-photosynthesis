# Detrital photosynthesis dataset
This dataset is a compilation of the repsonse of photosynthesis and chlorophyll to detachment in various plants. `Detrital_Photosynthesis.csv` contains all data and `References.pdf` gives a full bibliography of the cited studies. The data are structured into these variables:

-`Reference` Citation
-`DOI` Digital object identifier
-`Group` Non-taxonomic plant group with levels "Terrestrial", "Freshwater", "Seagrass" and "Seaweed"
-`Phylum` `Order` `Family` `Species` Taxonomy according to Plants of the World Online (https://powo.science.kew.org) and AlgaeBase (https://www.algaebase.org)
-`Light` Binary classifier of light availablity during the experiment
-`Water` Binary classifier of water availablity during the experiment
-`Series` Timeseries/experiment within study 
-`Day` Time after detachment/excision given in days
-`Mean` Observation or mean of observations
-`SEM` Standard error of the mean
-`N` Number of observations
-`Response` Response variable with levels "Photosynthesis" and "Chlorophyll"
-`Method` Measurement method
-`Unit` Measurement unit
-`Source` Data source within study

Luka Seamus Wright, 20th March 2025
