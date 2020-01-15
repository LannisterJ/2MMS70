%MACRO runscuc(data=,var=,alpha=);
proc iml;
use &data;
read all var {&var} ;
close &data;

X=&var;
n=nROW(X);
MED=median(X);

XC=X;
DO i=1 to n by 1;
	IF (XC[i] >= MED) then XC[i]=1;
	ELSE XC[i]=0;
END;

n1C=sum(XC);
n2C=n-n1C;

RC=1;
DO i=2 to n by 1;
	if(XC[i] ^= XC[i-1]) then RC=RC+1;
END;

MUC=1+(2*n1C*n2C)/(n1C+n2C);
VARC=2*n1C*n2C*(2*n1C*n2C-n1C-n2C)/((n1C+n2C-1)*(n1C+n2C)**2);

SC=(RC-MUC)/SQRT(VARC);
TC=QUANTILE('NORMAL',&alpha/2);
TCU=QUANTILE('NORMAL',1-&alpha/2);
PC=(1-CDF('NORMAL',abs(SC)))*2;

XUC=REPEAT(0,n-1,1);
TIES=0;
DO i=1 to (n-1) by 1;
	IF (X[i+1] > X[i]) then XUC[i]=1;
	IF (X[i+1] < X[i]) then XUC[i]=0;
	IF (X[i+1] = X[i]) then XUC[i]=XUC[i-1];
	IF (X[i+1] = X[i]) then TIES=TIES+1;
END;

RUC=1;
DO i=2 to (n-1) by 1;
	if(XUC[i] ^= XUC[i-1]) then RUC=RUC+1;
END;

MUUC=(2*(n-TIES)-1)/3;
VARUC=(16*(n-TIES)-29)/90;

SUC=(RUC-MUUC)/SQRT(VARUC);
TUC=QUANTILE('NORMAL',&alpha/2);
TUCU=QUANTILE('NORMAL',1-&alpha/2);
PUC=(1-CDF('NORMAL',abs(SUC)))*2;

c={'runs', 'mu(runs)', 'std(runs)', 'p-value', 'normalised stat.', 'critical value lower', 'critical value upper', 'n-ties'};

PRINT("Median based (conditional) runs test");
PRINT (RC||MUC||sqrt(VARC)||PC||SC||TC||TCU||n)[colname=c];
PRINT("(unconditional) runst test for serial randomness");
PRINT(TIES);
PRINT (RUC||MUUC||sqrt(VARUC)||PUC||SUC||TUC||TUCU||(n-TIES))[colname=c];
quit;
%MEND;

%macro normality_SK(skewness=, kurtosis=, n=); *tests normality with skewness and kurtosis;
DATA approx;
N=&n;
G1=&skewness;
G2=&kurtosis;
b1=(N-2)*G1/(sqrt(N*(N-1)));
b2=G2*((N-2)*(N-3))/((N+1)*(N-1))+3*(N-1)/(N+1);
*JB=N*(b1**2/6+(G2-3)**2/24);
Cn=(3*(N**2+27*N-70)*(N+1)*(N+3))/((N-2)*(N+5)*(N+7)*(N+9));
Wn2=-1+SQRT(2*(Cn-1));
Alphan=SQRT(2/(Wn2-1));
Dn=1/sqrt(log(sqrt(Wn2)));
Bn=sqrt((N+1)*(N+3)/(6*(N-2)))*b1;
Ts=Dn*log(Bn/Alphan+sqrt(1+(Bn/Alphan)**2));
Mun=3*(N-1)/(N+1);
Sigman=sqrt((24*N*(N-2)*(N-3))/((N+3)*(N+5)*(N+1)**2));
Gamma1n=((6*(N**2-5*N+2))/((N+7)*(N+9)))*sqrt(6*(N+3)*(N+5)/(N*(N-2)*(N-3)));
An=6+(8/(Gamma1n))*(2/Gamma1n+sqrt(1+4/(Gamma1n**2)));
Un=(b2-Mun)/Sigman;
Tk=sqrt(9*An/2)*((9*An-2)/(9*An)-((1-2/An)/(1+Un*sqrt(2/(An-4))))**(1/3));
K2=Tk**2+Ts**2;
Ps=2*min(cdf('Normal',Ts,0,1),1-cdf('Normal',Ts,0,1));
Pk=2*min(cdf('Normal',Tk,0,1),1-cdf('Normal',Tk,0,1));
PK2=1-cdf('chisq',K2,2);
*PJB=1-cdf('chisq',JB,2);
run;

proc iml;
use approx;
read all var _ALL_;
close approx;

print (G1 || Ts || Ps)[colname={'skewness approx', 'test statistic', 'p-value'}];
print (G2 || Tk || Pk)[colname={'kurtosis approx', 'test statistic', 'p-value'}];
run;

proc delete data=approx; run;
%mend;

%macro box_cox(data=, variable=);
title &data;
DATA &data;
set &data;
negative_three = (-1/3)*(&variable**-3 -1);
negative_two= (-1/2)*(&variable**-2 -1);
negative_one= (-1)*(&variable**-1 -1);
negative_half= (-2)*(&variable**-(0.5)-1);
negative_onethird = (-3)*(&variable**(-1/3)-1);
zero = log(&variable);
positive_onethird = (3)*(&variable**(1/3) -1);
positive_half= (2)*(&variable**(1/2) -1);
positive_one = &variable -1;
positive_two = (0.5)*(&variable**(2) -1);
positive_three = (1/3)*(&variable**3 -1);
RUN;

proc univariate data=&data;
histogram negative_three /normal;
histogram negative_two /normal;
histogram negative_one /normal;
histogram negative_half /normal;
histogram negative_onethird /normal;
ods select histogram;
histogram zero /normal;
histogram positive_onethird /normal;
histogram positive_half /normal;
histogram positive_one /normal;
histogram positive_two /normal;
histogram positive_three /normal;
run;
%mend;

%macro Doornbos(data=, variable=, output_name=, iteration=);
data L10;
set &data;
run;
proc iml; print (&iteration)[colname={'Holm-Bonferroni (Doornbos) iteration#:'}]; run;
proc means data=L10 mean std n;
var &variable;
output out=L10_sumstat mean=mean median=median std=std n=n;
run;

data L10;
set L10;
if _n_ = 1 then set L10_sumstat;
drop _TYPE_ _FREQ_;
run;

data &output_name;
SET L10;
U = (&variable-MEAN)/STD;
W = SQRT((N*(N-2)*U**2)/((N-1)**2-N*U**2));
DOORNBOS_Y=ABS(W);
CRITER= QUANTILE('T', 1-0.05/(2*N), N-2);
P= MIN(2*MIN(CDF('T', DOORNBOS_Y, N-2),1-CDF('T',DOORNBOS_Y, N-2))*N,1);
RUN;

proc sort data=&output_name;
by P;
run;

proc sql;
create table stripped as 
select * 
from &output_name
where p not in (select min(p) from &output_name where p < 0.05);

create table removed as
select *
from &output_name
where p not in (select p from stripped);
run;

data &output_name;
set stripped;
drop U W DOORNBOS_Y CRITER mean median n std;
run;

proc iml; print('Deleted outliers with Doornbos'); run;
proc print data=removed; run;
proc iml; use removed; read all var{P}; close removed; if nrow(P) = 0 then print("Done: no more outliers");
proc delete data=removed; run;
proc delete data=stripped; run;
proc delete data=L10; run;
proc delete data=L10_sumstat; run;
%mend;

%macro Hampel(data=, variable=, output_name=);
data L10; set &data; run;

proc means data=L10 median n;
var &variable;
output out=L10_sumstat median=median;
run;

data L10;
set L10;
if _n_= 1 then set L10_sumstat;
drop _TYPE_ _FREQ_;
run;

title &data;
proc iml; print('Hampel'); run;

data Hampel;
set L10;
D = abs(&variable-median);
RUN;

proc means data=Hampel median n;
var D;
output out=Hampel_med median=medianD;
run;

data Hampel;
set Hampel;
if _n_= 1 then set Hampel_med;
drop _TYPE_ _FREQ_;
run;

DATA Hampel;
SET Hampel;
Z = ABS(&variable-median)/medianD;
H = (Z>3.5);
RUN;

data removed; 
set hampel;
if H = 0 then delete;
run;

PROC PRINT DATA=removed;RUN;

data &output_name;
set hampel;
if H = 1 then delete;
run;

proc delete data=L10; run;
proc delete data=hampel; run;
proc delete data=L10_sumstat; run;
proc delete data=hampel_med; run;
proc delete data=removed; run;
title;
%mend;

%MACRO TIES(data=,var=); *counts the number of ties of previously unsorted data;
data X;
set &data;
run;

proc sort data=X;
by &var;
run;

proc iml;
use X;
read all var {&var} ;
close X;
Y=&var;
n=nrow(Y);
ties=0;
DO i=1 to (n-1) by 1;
	IF (Y[i+1] = Y[i]) then TIES=TIES+1;
END;
print (ties || (ties/n))[colname={'ties', 'ratio'}];
quit;

proc delete data=X; run;
%MEND;