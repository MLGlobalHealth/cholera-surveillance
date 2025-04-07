# Cholera Surveillance

In this repo we collect the work on cholera dynamics modelling developed as part of the Wellcome Trust grant "Integrated Surveillance Suite for Targeting Interventions to Cholera Outbreaks".

Currently we have focused on the study of cholera outbreaks occurring in a Uvira, Democratic Republic of the Congo (DRC) in collaboration with the group of Prof. Andrew Azman (Department of Epidemiology, John Hopkins Bloomber School of Public Health and Geneva Centre for Emerging Viral Diseases and Division of Tropical and Humanitarian Medicine). Uvira is in an area of the Eastern DRC characterized by endemic cholera patterns.

The Python and R code is organized in notebooks, in the folder with the same name, as follows:
- data_exploration:
    - eda_uvira_data.ipynb: a short exploratory data analysis of the data. In particular, we highlight potentially significant variation in the effectiveness in testing the suspected cholera cases through time.
- data_processing:
    - climate_data.ipynb: extract climate related indicators to use as covariates in the cholera dynamics models.
    - migration_conflict_disaster_data.ipynb: extract and combine additional indicators concerning Internally Displaced People and the impact of conflicts and natural disasters.
    - combine_data.ipynb: combine all data into a single dataset with daily resolution.
- models: 
    - covariates_model.ipynb: we develop an approach to identify the features of the covariate time series more strongly related to the cholera case load and construct one single combined effect.
    - pomp_single_outbreak.qmd: initial tests using a partially observed markov process (pomp) model to fit a single cholera outbreak. Used only for testing purposes.
    - pomp_multi_outbreak.qmd: the pomp approach is extended to all data and several models with different levels of complexities are implemented. Each model is separately implemented and analysed. In future iterations we will develop a more standardized way to define and test different models. Note that, as of now, all modelling has been focused on all the suspected cholera cases rather than those that were confirmed through some testing mechanism. The main reason for this is the seeming inconsitency noticed in the testing performance through time.
    - pomp_multi_outbreak_vacc.ipynb: some initial attempts to add the effect of vaccination campaigns to the model.  

