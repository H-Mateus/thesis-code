# thesis-code

Scripts from my thesis

## Notes for EMSCI collaborators

The script for modelling of the data is in the "500_patient_modelling" folder and is called "500_modelling.R"

I created some fake data with the same structure as the what we used.
The first few lines of the script do a little cleaning, extract the predictor and model target cols, and then use a function in the "500_functions.R" script to generate the actual models.
If you just want to see the parameter of the models, I copied that part of the function into the script "model_parameters.R"

We built two models for each model target, one using plain linear regression (named "linreg_model" in the script), the other using linear regression with elastic net penalization (named "glmnet_model").
I also tired some logistic regression models whilst initially model building for predicting ASIA conversion, but the models performed poorly so we didn't take them forward.
I tried to remove all the code relevant to the logistic regression experimentation I did, but apologies if I missed anything.

The models themselves are structured in a tibble where each row is a model.
Hopefully the columns of the tibble are self-explanatory, but please let me know if you're not sure of anything or run into any issues.


