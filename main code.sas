%include '/folders/myfolders/sasuser.v94/assignment_macros.sas';
data bchc; 
set '/folders/myfolders/sasuser.v94/assignment/bchc.sas7bdat'; 
keep mortality place race_ethnicity year; 
run;

/* (1) data preparation */
/* We start with removing empty data points from our data */

data bchc;
set bchc;
if cmiss(of mortality race_ethnicity) then delete;
run;

/* We continue to conduct the outlier analysis. First we must check whether our data is normally  */
/* distributed as this analysis relies of normal data */
/* A high number of ties would imply that Anderson-Darling is more reliable than Shapiro-Wilk */

%ties(data=bchc, var=mortality); 

/* Only 5 ties have been found which is 3.38% of the data. We assume this is a low number of ties */
/* ==> we use Shapiro-Wilk */
ods output "Tests for Normality";
proc univariate data=bchc NORMAL;
var mortality;
ods select probplot "Tests for Normality"; 
probplot mortality/Normal(MU=est SIGMA=est);
run;

/* We also create a histogram to visually analyse the normality */
proc univariate data=bchc;
ods select histogram; *suppresses other univariate output;
histogram mortality /normal;
run;

/* Quadratic curve with a significant departure from normality p-value < 0.0001 and.test statistic 0.922236 */
/* The histogram shows that the distribution is skewed to the left */
/* We reject normality */

/* We continue to transform the data to follow a normal distribution by means of the box-cox transformation */
%box_cox(data=bchc, variable=mortality);

/* The positive 1/3 box cox transformation appears the most normal. We test this observation below:  */
ods output "Tests for Normality";
proc univariate data=bchc normal;
var positive_onethird;
ods select probplot "Tests for Normality"; 
probplot positive_onethird/Normal(MU=est SIGMA=est);
run;

/* We obtain a p-value of 0.2342 and a test statistic of 0.988026 ==> fail to reject normality */

/* Now we can perform the outlier analysis */
/* We apply Doornbos */
%Doornbos(data=bchc, variable=positive_onethird, output_name=positive_onethird, iteration=1);

/* No outliers are found so we continue to use Hampel's rule because we want to avoid the masking effect */
%Hampel(data=bchc, variable=positive_onethird, output_name=bchc);

/* Four outliers have been removed (see results) */

/* (2) ANOVA model */
/* we want to apply two-way ANOVA and therefore we have to test its assumptions */
/* 	(i) normality of residuals; ✔*/
/* 	(ii) homogeneity of residual variance across groups; ✔ */
/* 	(iii) normality of random effects; ✔*/
/*  (iv) the two variables have to be in categorical and independent groups. ✔*/

/* the ANOVA model */
ods output "Solution for Random Effects"=EBLUP;
ods output "Covariance Parameter Estimates"=two_ANOVA_place;
proc mixed data=bchc method=type3;
	class race_ethnicity place;
	model mortality = race_ethnicity /solution outp=residuals_data ddfm=satterthwaite;
	random place /solution;
	lsmeans race_ethnicity /diff adjust=tukey cl;
run;

/* (i) */
ods output "Tests for Normality";
proc univariate data=residuals_data NORMAL;
var resid;
ods select probplot "Tests for Normality"; 
probplot resid/Normal(MU=est SIGMA=est);
run;

proc means data=residuals_data skewness kurtosis n; var resid; run;
%normality_SK(skewness=0.2054122, kurtosis=1.0312354, n=144);

/* (ii) */
proc glm data=residuals_data;
class race_ethnicity;
model resid=race_ethnicity;
means race_ethnicity / hovtest=bf;
RUN;

/* (iii) */
ods output "Tests for Normality";
proc univariate data=EBLUP NORMAL;
var estimate;
ods select probplot "Tests for Normality"; 
probplot estimate/Normal(MU=est SIGMA=est);
run;

/* (iv) */
%runscuc(data=residuals_data, var=resid, alpha=0.05);
%runscuc(data=EBLUP, var=estimate, alpha=0.05);

/* .......................................................................... */
/* we proceed to quantify the amount of variability between the cities  */

/* random ANOVA model for σ^2(β~)*/
ods output "Covariance Parameter Estimates"=one_ANOVA_place;
proc mixed data=bchc method=type3;
	class place;
	model mortality=;
	random place;
run;

proc iml;
use one_ANOVA_place; read all var _ALL_; close one_ANOVA_place;
do i = 1 to nrow(CovParm);
	if CovParm[i] = "Place" then one_way_place_var  = Estimate[i];
end;
use two_ANOVA_place; read all var _ALL_; close one_ANOVA_place;
do i = 1 to nrow(CovParm);
	if CovParm[i] = "Place" then two_way_place_var  = Estimate[i];
end;
variability = ((one_way_place_var  - two_way_place_var )/one_way_place_var)*100 ;
print (one_way_place_var || two_way_place_var || variability)[colname={"est. σ^2(β~)", "est. σ^2(β)", "impact percentage"} format=12.2];
run;