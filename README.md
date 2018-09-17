FUNCTIONS FOR GENERATING EMPIRICAL PRIORS OF CLADE AGE USING THE FOSSIL RECORD

Santiago Claramunt

Department of Natural History, Royal Ontario Museum and Department of Ecology and Evolutionary Biology, University of Toronto, Ontario, Canada.

E-mail: sclaramunt@rom.on.ca


DETAILS

Given a set of fossil ages for a specific clade, these functions model the ucertainty about the age of the clade. The basic method is based on the model of Strauss & Sadler (1989) and assumes a uniform distribution of presicely known fossil ages but there are options for incorporating fossil age uncertatinty (Claramunt & Cracraft 2015) and non-uniform fossil age distributions (Solow 2003). The first version of the functions aimed at generating probability density functions but version 3 implements a more efficient strategy involving the generation of random samples of possible clade ages using quantile functions. These numbers can then be used to generate a histogram and fit probability functions for using as priors in Bayesian time tree estiamtion (Claramunt & Cracraft 2015).

REFERENCES

Claramunt, S., & J. Cracraft 2015. A new time-tree reveals Earth history's imprint on the evolution of modern birds. ScienceAdvances in press.

Solow, A. R. 2003. Estimation of stratigraphic ranges when fossil finds are not randomly distributed. Paleobiology 29(2):181-185.

Strauss D, & P. M. Sadler 1989. Classical confidence intervals and Bayesian probability estimates for ends of local taxon ranges. Mathematica Geology 21:411-427.
