function [b,bint,r,rint,stats] = MLR_SkySWin(esky,swin,bias)
%{
Estimate multi-linear regession coefficients for sub-canopy longwave
radiation bias with predictors:
- effective emissivity of the sky,
- insolation, and
- interaction term of the above.

Note from Matlab documentation:
When computing statistics, X should include a column of 1s so that the 
model contains a constant term. The F statistic and its p value are 
computed under this assumption, and they are not correct for models without
a constant.
The F statistic is the test statistic of the F-test on the regression 
model, for a significant linear regression relationship between the 
response variable and the predictor variables.
The R2 statistic can be negative for models without a constant, indicating 
that the model is not appropriate for the data.
%}
X = [ones(size(bias)) esky swin swin.*esky];
[b,bint,r,rint,stats] = regress(bias,X);    % removes NaN data
end