# Applied study

## Code overview

Original public-use data files from the Food Acquisition and Purchase Survey (FoodAPS) are found in the "data'' directory; the file with household variables is "faps_household_puf.csv", the one with individual-level variables (for our case, the primary respondent characteristics) is "faps_individual_puf.csv", and the file containing the survey weights is "faps_hhweights.csv"

First, run "step0-exploration-faps.R" to produce "faps_cleaned.RDS" to be used as an input in "step1-pencomp-analysis-faps.R" where we do the weight trimming step. The final data used for all analyses is "dat_use_faps.RDS"

The code files (step 1 and step 2) are left in their cluster-friendly forms, which means that minor modifications are needed to reproduce results on a local desktop. The postprocessing codes (step 3) are provided to document how the summary statistics and confidence intervals are produced.

Detailed code for running the methods themselves can be found in the "CODE-publish" directory, where they are sourced by the step codes. 

---

FoodAPS (April 2012 - mid January 2013):
 * 	Survey goal: track food acquisition patterns over a 7-day period for non-institutionalized households in the contiguous US.
 *	Note: definition of a FoodAPS "household" is broader than what is used in other surveys. "The FoodAPS household is defined as all persons who live together and
share food and who expect to be present at the sampled address during at least part of the data collection week."

 * 	Survey methods: 
	*	Probability sample
	*	Each household has a "primary respondent," the main food shopper or meal planner in the household.
	*	PR completed two in-person interviews, three telephone interviews. HOusehold members above age 11 scanned barcodes on foods, saved receipts from stores, and recorded the information in their provided food boks.

*	Survey design:
	* Target groups were defined by total reported household income (relative to poverty guideline) and SNAP participation status.
	*	PSU: county or groups of contiguous counties, or MSAs.
	*	MOS determined by number of SNAP households in the PSU, number of non-SNAP households from the three income groups (below 100%, between 100-184%, and 185% and above of poverty guideline)
	* 	SSUs (within boundaries of incorporated cities, towns, villages, boroughs, and unincorporated areas) selected with MOS similar to PSU
	*	Sampling frame for each SSU constructed with indicators for SNAP and non SNAP addresses.
	*	The sampling (and nonreponse) weights are based on the SNAP participation status revised to match administrative data (SNAPNOWHH). Weights were raked based on race (Hispanic, Non-Hispanic White, Non-Hispanic Black, Non-Hispanic other), annual income (<$15K, $15K-49K, $50K+), SNAP participation, household size, number of chidren in the household, and presence of at least one person age 60 or older in the household. Race, income, household size, and number of children were calibrated to 2013 CPS ASEC totals while SNAP receipt and inclusion of a person age 60 or older were obtained from the 2012 ACS.

---
## Main Analysis
*	Z = HH SNAP participation (SNAPNOWHH)
		Race, presumably of the PR (Hispanic, Non-Hispanic White, Non-Hispanic Black, Non-Hispanic other)
		HH Size (HHSIZE)
		HH number of children (HHSIZECHILD)
		HH presence of at least one person age 60 or older.
		HH Annual Income (PCTPOVGUIDEHH_R)
*	X = Household income (<$15K, $15K-49K, $50K+), 
		Race (Non-Hispanic White, Non-Hispanic Black, Hispanic, Other), 
		Education (No HS diploma, HS grad or GED, Some college, College graduate or above)
		Household size 
		Primary store travel time (PRIMSTORETRAVELTIME)
*	A = Household has been food assisted at least once in the past year (Yes/No) =1 if any of the following apply
		If SNAPNOWHH = 1 (n = 1581), then the household is receiving food assistance
		Otherwise, if SNAPEVER = 1 and SNAP12MOS =1 (n = 198), then the household is receiving SNAP food assistance
		USDA FOODS=1 (n=183) (anyone in household receiving USDA foods from local program or distribution site)
		WICHH = 1 (461) (Anyone in household is receiving benefits from WIC)
		MEALFACILITY = 1 (114) (Has anyone in the household received meals at a community center in the past month)
		MEALDELIVERY = 1 (53) (Anyone in household receiving meals at home from community programs)
*	Y = HH adult food security (30 day, Food insecure / secure, ADLTFSCAT)

**Note**: Residential unit = all individuals on the roster for the address. Household = all individuals except guests. Family = PR, all household members related to PR except guests. 

We are interested in the relationship between recent food assistance and household adult food security status, and suspect that the effect of food assistance is confounded by race. 

**Model selection:**

In R notation, the outcome model we use:

$$I(\textrm(food insecure)|\textrm{assist}) \sim race + age + income + educ + hhsize + assist + race*assist + spline(P^{*assist})$$

The treatment propensity model was selected using the BIC criterion:

$$I(\textrm{assist}|X,Z) \sim income + age + WICelig + race + educ + hhsize + age*WIC$$

We will also include the estimate where the spline for sampling weights is included in this model. We include 3 splines for the sampling weights in the treatment propensity model to avoid .


### Methods

We will use the survey-weighted estimator after imputing the potential outcomes in the sample using the following imputation methods:

*	classical PENCOMP
*	WFPBB-PENCOMP without spline in the treatment propensity model.
*	WFPBB-PENCOMP with spline for sampling weights in model for $P(assist|X,Z)$
*	AIPTW (augmented inverse probability of treatment weighting; double-robustness properties)


## Checks:
* Are we seeing treatment heterogeneity, i.e. is there a significant interaction term between Z and A?
	* Yes. If we consider $Z$ as race and $A$ as recent food assistance and fit the survey-weighted model (in R syntax)
	$Y \sim income + age + race + educ + assist + race*assist$
	we find that there is a significant interaction term for Black and food assistance ($\beta_{black:assist} =-0.6, pvalue = 0.01$)
* Is there informative sampling, i.e. is there a significant difference in the weight distributions between levels of Y?
	* Yes. If we model $Y \sim wts$, we find that the coefficient for weights is significant ($p<0.001$).
* Are the variables influencing selection and treatment mechanisms shared? 
	* Yes. Race (Black; $p<0.001$), income (all $p<0.001$), and age (65+; $p<0.001$) are significant in the models $A \sim race + income + race$. All categories are significant at $p<0.01$ except for the Hispanic term for the regression $wts \sim race + income + race$
