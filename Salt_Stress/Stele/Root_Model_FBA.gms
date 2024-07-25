*************************************************************
* Basic FBA model for Arabidopsis root
*************************************************************
* Lohani Beer
*************************************************************
$INLINECOM /*  */

OPTIONS
    decimals = 8,
    lp = cplex
;

********* Defining Sets **************************************
SETS
    i   set of metabolites
$include "metabolites.txt"    

    j   set of reactions
$include "reactions.txt"
;

*************************************************************

*********** Defining Parameters *******************************
PARAMETERS
    S(i,j)              stoichiometric matrix
$include "Sij.txt"

    v_max_salt1h(j)       maximum flux of v(j) for salt condition
$include "UB_Salt1h.txt"

    v_min_salt1h(j)       minimum flux of v(j) for salt condition
$include "LB_Salt1h.txt"

    v_max_salt48h(j)       maximum flux of v(j) for salt condition
$include "UB_Salt48h.txt"

    v_min_salt48h(j)       minimum flux of v(j) for salt condition
$include "LB_Salt48h.txt"

    v_max_std(j)        maximum flux of v(j) for standard condition
$include "UB_STD.txt"

    v_min_std(j)        minimum flux of v(j) for standard condition
$include "LB_STD.txt"
;

**************************************************************

********* Defining Equations ***********************************
EQUATIONS
    objective                objective function
    mass_balance1(i)         steady state mass balance
    lower_bound_std(j)       lower bounds on reactions for standard condition
    upper_bound_std(j)       upper bounds on reactions for standard condition
    lower_bound_salt1h(j)      lower bounds on reactions for salt condition
    upper_bound_salt1h(j)      upper bounds on reactions for salt condition
    lower_bound_salt48h(j)      lower bounds on reactions for salt condition
    upper_bound_salt48h(j)      upper bounds on reactions for salt condition

;

**************************************************************

********* Defining Variables ***********************************
FREE VARIABLES
    v(j)                  reaction flux
    obj                   objective value
;

**************************************************************

*************** Defining Model ***********************************
objective..            obj =e= v('Root_Biomass[GS1]');

mass_balance1(i)..     sum(j, S(i,j) * v(j)) =e= 0;

lower_bound_std(j)..   v_min_std(j) =l= v(j);

upper_bound_std(j)..   v(j) =l= v_max_std(j);

lower_bound_salt1h(j)..  v_min_salt1h(j) =l= v(j);

upper_bound_salt1h(j)..  v(j) =l= v_max_salt1h(j);

lower_bound_salt48h(j)..  v_min_salt48h(j) =l= v(j);

upper_bound_salt48h(j)..  v(j) =l= v_max_salt48h(j);


Model arabidopsis_root_std / objective, mass_balance1, lower_bound_std, upper_bound_std /;
Model arabidopsis_root_salt1h / objective, mass_balance1, lower_bound_salt1h, upper_bound_salt1h /;
Model arabidopsis_root_salt48h / objective, mass_balance1, lower_bound_salt48h, upper_bound_salt48h /;

**************************************************************

********** Solving Model under Standard Conditions ************
Solve arabidopsis_root_std using lp maximizing obj;

FILE RESULTS_STD /results_FBA_STD.txt/;
PUT RESULTS_STD;
PUT "BiomassSTD: " obj.l:0:8 //;

loop(j,
    put j.tl:0:30, "    ", v.l(j):0:8 /;
);

PUTCLOSE;

**************************************************************

********** Solving Model under Salt 1h Conditions ****************
Solve arabidopsis_root_salt1h using lp maximizing obj;

FILE RESULTS_SALT1h /results_FBA_SALT1h.txt/;
PUT RESULTS_SALT1h;
PUT "BiomassSALT1h: " obj.l:0:8 //;

loop(j,
    put j.tl:0:30, "    ", v.l(j):0:8 /;
);

PUTCLOSE;
**************************************************************
********** Solving Model under Salt 48h Conditions ****************
Solve arabidopsis_root_salt48h using lp maximizing obj;

FILE RESULTS_SALT48h /results_FBA_SALT48h.txt/;
PUT RESULTS_SALT48h;
PUT "BiomassSALT48h: " obj.l:0:8 //;

loop(j,
    put j.tl:0:30, "    ", v.l(j):0:8 /;
);

PUTCLOSE;
**************************************************************

