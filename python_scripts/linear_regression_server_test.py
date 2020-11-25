#%%
import feather
import numpy as np
import pandas as pd

desired_width = 250
pd.set_option('display.width', desired_width) # Change use more of the space in the console when printing
pd.set_option("display.max_columns", 500) # Change to view more columns with print
pd.set_option('display.max_rows', 10) # Change to view more rows with print

df_imputed = feather.read_dataframe("python_data/500_df_wide_scaled_locf_imputed.feather")
df = feather.read_dataframe("python_data/500_df_wide_scaled_bloods_imputed.feather")
df.set_index('id')
df_imputed.set_index('id')


print(type(df))

print(df.shape)
print(df.columns)
print(df.index)
print(df.head(10))
print(df.info())
test = df.describe()
print(test)
#print(test.iloc.describe[1, 5],[]) # not wokring - attempt to select mean and median columns from .describe
print("End of df checks")

print(df_imputed.info())
#check if there are any NAs left 
print(df.isnull().values.any())
print(df_imputed.isnull().values.any())
print(df.isnull().sum().sum())
print(df_imputed.isnull().sum().sum())

#%%
#This gives a R^2 of 1 for the complete df, which must be wrong, but it gives an R^2 of 0.8 for the imputed dataset

#compare imputed vs non_imputed data with linear regression models
import statsmodels.api as sm

# convert id to numeric to prevent it from being encoded 
df['id'] = pd.to_numeric(df['id'])
df_imputed['id'] = pd.to_numeric(df_imputed['id'])

# Use hot-one-encoding to convert categorical variables to numbers 
df_encoded = pd.get_dummies(df, drop_first = True)
df_imputed_encoded = pd.get_dummies(df_imputed, drop_first = True)
#make a dataframe with only complete cases
df_cc = df_encoded.dropna(how = 'any')

#remove other year 1 columns for predicting motor
year1_cols = "X1_year", "X1year_motor", "X1year_sensor_prick", "X1year_sensor_touch", "year1_SCIM", "by_year_1" # inserting these objects doesn't in [] below doesn't work
year1_motor_cols = "X1_year", "X1year_sensor_prick", "X1year_sensor_touch", "year1_SCIM", "by_year_1" # inserting these objects in the [] below doesn't work

df_cc.drop(["X1_year", "X1year_sensor_prick", "X1year_sensor_touch", "year1_SCIM", "by_year_1"], axis = 1, inplace = True)
df_imputed_encoded.drop(["X1_year", "X1year_sensor_prick", "X1year_sensor_touch", "year1_SCIM", "by_year_1"], axis = 1, inplace = True)

    # The code below uses statsmodels package builds a linear regression model on the data with all NA rows removed, and another on an imputed version of this data - Note: output is not correct 
y = df_cc['X1year_motor'].values #select variable to model
X = df_cc.drop('X1year_motor', axis=1).values #select all other variables 
# print("x1 = ", X[1])
# print("X", X)
# print("y1 = ", y[1])
# print("x2 = ", X[2])
# print("y2 = ", y[2])
# print("last x", X[-1])
# print("last y", y[-1])
lm = sm.OLS(y, X).fit() #build linear model 

#print(lm.summary()) # - this prints a very long output that fills the terminal
print(lm.rsquared_adj)
#print(lm.params) # this shows coeffiecients of the model 

 # Create arrays for the features and the response variable
y = df_imputed_encoded['X1year_motor'].values
X = df_imputed_encoded.drop('X1year_motor', axis=1).values

lm_imputed = sm.OLS(y, X).fit()

#print(lm_imputed.summary()) # - this prints a very long output that fills the terminal 
print(lm_imputed.rsquared_adj)
 #print(lm_imputed.params) 

# # to compare multiple models r-squared easily use the following - Note: the df with complete cases has an R-Squared of 1 - something must be wrong
print(pd.DataFrame({'complete_cases': lm.rsquared_adj, 'locf_imputed': lm_imputed.rsquared_adj}, index=['r_squared_adj']))

print(pd.DataFrame({'complete_cases': lm.params, 'locf_imputed': lm_imputed.params}))

#to create density plots to compare the imputed values to the complete cases
import matplotlib.pyplot as plt
df_cc['discharge_SCIM'].plot(kind='kde', c='red', linewidth=3)
df_imputed_encoded['discharge_SCIM'].plot(kind='kde')
labels = ['Baseline (Complete Case)', 'LOCF Imputation']
plt.legend(labels)
plt.xlabel('discharge_SCIM')
plt.show()


print(df_cc.shape)
print(df_cc.info())
print(df_cc.head(10))


#%%
# This says the R^2 is 1 - something is not right
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_val_score
# X and y for rows with no NAs - modelling 1 year motor score 
y = df_cc['X1year_motor'].values #select variable to model
X = df_cc.drop('X1year_motor', axis=1).values #select all other variables 

# define test and training sets - test_size determins the amount of data removed from testing 
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state=42)

reg_all = LinearRegression()
reg_all.fit(X_train, y_train)
y_pred = reg_all.predict(X_test)

#perform 10-fold cross-validation 
cv_results = cross_val_score(reg_all, X, y, cv = 10)

print(cv_results) #by default this prints R^2
print(np.mean(cv_results))

#to get R-squared
print("R^2: {}".format(reg_all.score(X_test, y_test)))
rmse = np.sqrt(mean_squared_error(y_test, y_pred))
print("Root Mean Squared Error: {}".format(rmse))


#%%
# This gives a negative R^2 of less than -700,000 - something is wrong (R^2 can't be less than -1 (perfect negative correlation))
# sklearn with imputed data 
y = df_imputed_encoded['X1year_motor'].values
X = df_imputed_encoded.drop('X1year_motor', axis=1).values
# define test and training sets - test_size determins the amount of data removed from testing 
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state=42)

reg_all = LinearRegression()
reg_all.fit(X_train, y_train)
y_pred = reg_all.predict(X_test)

#perform 10-fold cross-validation 
cv_results = cross_val_score(reg_all, X, y, cv = 10)

print(cv_results)
print(np.mean(cv_results))

#to get R-squared
print("R^2: {}".format(reg_all.score(X_test, y_test)))
rmse = np.sqrt(mean_squared_error(y_test, y_pred))
print("Root Mean Squared Error: {}".format(rmse))