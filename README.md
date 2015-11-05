FUNCTIONS FOR GENERATING EMPIRICAL PRIORS OF CLADE AGE USING THE FOSSIL RECORD

Author
Santiago Claramunt
E-mail: claramunt.bio@gmail.com

DETAILS

Given a set of fossil ages for a specific clade, these functions return a probability density for the age of the clade (Claramunt & Cracraft 2015). Specifically, if a set of fossil ages (x) is uniformly distributed, the probability density of the true age t is proportional to the joint probability of x given t, which is 1/t^n for x < t (Wang & Everson 2007, Wang 2010). This formula, used in function PD.age, requires a baseline, which could be the present (baseline=0) or another reference stratum. If a baseline is not provided, the function uses the youngest fossil age as the baseline and n-1 instead of n (Solow 2003). PD.age.unc uses pseudoreplicates to accommodate cases in which fossils ages are not know precisely but assigned to a chronostratigraphic interval; pseudoreplicates of fossil ages are generated by sampling uniformly from the time interval and then a final distribution is computed by averaging across pseudoreplicates. Function PD.age.plot produced a plot that in combination with line plots can be used to find an appropriate parametric probability density function to represent the empirical distribution to be incorporated as clade age priors in Bayesian analysis software (see Examples for details). A wrapper with an alternative parameterization of the log-normal probability density function is provided to facilitate curve fitting by eye.

Suggestions and contributions are welcomed.
