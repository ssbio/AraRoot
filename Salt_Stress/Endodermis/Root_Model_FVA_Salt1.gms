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

    v_max_salt1h(j)       maximum flux of v(j) for salt condition
$include "UB_Salt1h.txt"

    v_min_salt1h(j)       minimum flux of v(j) for salt condition
$include "LB_Salt1h.txt"

    c(j)
;

**************************************************************

********* Defining Equations ***********************************
EQUATIONS
    objective                objective function
    mass_balance1(i)         steady state mass balance
    lower_bound_salt1h(j)      lower bounds on reactions for salt condition
    upper_bound_salt1h(j)      upper bounds on reactions for salt condition

;

**************************************************************

********* Defining Variables ***********************************
FREE VARIABLES
    v(j)                  reaction flux
    z                    objective value
;

**************************************************************

*************** Defining Model ***********************************
v.lo('Root_Biomass[GS1]') = 0.635;
v.up('Root_Biomass[GS1]') = 0.635;

objective..            z =e= sum(j, c(j) * v(j));

mass_balance1(i)..     sum(j, S(i,j) * v(j)) =e= 0;


lower_bound_salt1h(j)..  v_min_salt1h(j) =l= v(j);
upper_bound_salt1h(j)..  v(j) =l= v_max_salt1h(j);



Model fva_salt1h / objective, mass_balance1, lower_bound_salt1h, upper_bound_salt1h /;


**************************************************************

********** Solving FVA under Salt Conditions 1h ****************
FILE RESULTS_SALT1h /results_FVA_SALT1h.txt/;
PUT RESULTS_SALT1h;

PUT "reaction     min_rate    max_rate" /;

LOOP(j1,
    c(j) = 0;
    c(j1) = 1;
    
    PUT j1.tl:0:100;
    
    solve fva_salt1h using lp MINIMIZING z;
    PUT z.l:20:5;
    
    solve fva_salt1h using lp MAXIMIZING z;
    PUT z.l:20:5 /;
);

PUTCLOSE;
**************************************************************
