*************************************************************
* Basic FVA model for Arabidopsis root
*************************************************************
* Lohani Beer
*************************************************************
$INLINECOM /*  */

OPTIONS
    limrow = 10000,
    limcol = 10000,
    optCR = 1E-9,
    optCA = 0.0,
    iterlim = 100000,
    decimals = 8,
    reslim = 100000,
    work = 50000000
;

********* Defining Sets **************************************
SETS
    i   set of metabolites
$include "metabolites.txt"

    j   set of reactions
$include "reactions.txt"

    alias(j,j1);

*************************************************************

*********** Defining Parameters *******************************
PARAMETERS
    S(i,j)              stoichiometric matrix
$include "Sij.txt"


    v_max_std(j)        maximum flux of v(j) for standard condition
$include "UB_STD.txt"

    v_min_std(j)        minimum flux of v(j) for standard condition
$include "LB_STD.txt"

    c(j)
;

**************************************************************

********* Defining Equations ***********************************
EQUATIONS
    objective                objective function
    mass_balance1(i)         steady state mass balance
    lower_bound_std(j)       lower bounds on reactions for standard condition
    upper_bound_std(j)       upper bounds on reactions for standard condition

;

**************************************************************

********* Defining Variables ***********************************
FREE VARIABLES
    v(j)                  reaction flux
    z                    objective value
;

**************************************************************

*************** Defining Model ***********************************
v.lo('Root_Biomass[GS1]') = 0.716;
v.up('Root_Biomass[GS1]') = 0.716;

objective..            z =e= sum(j, c(j) * v(j));

mass_balance1(i)..     sum(j, S(i,j) * v(j)) =e= 0;

lower_bound_std(j)..   v_min_std(j) =l= v(j);
upper_bound_std(j)..   v(j) =l= v_max_std(j);


Model fva_std / objective, mass_balance1, lower_bound_std, upper_bound_std /;

**************************************************************

********** Solving FVA under Standard Conditions *************
FILE RESULTS_STD /results_FVA_STD.txt/;
PUT RESULTS_STD;

PUT "reaction     min_rate    max_rate" /;

LOOP(j1,
    c(j) = 0;
    c(j1) = 1;
    
    PUT j1.tl:0:100;
    
    solve fva_std using lp MINIMIZING z;
    PUT z.l:20:5;
    
    solve fva_std using lp MAXIMIZING z;
    PUT z.l:20:5 /;
);

PUTCLOSE;

**************************************************************

