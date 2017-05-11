#!/bin/bash

grep 'Condition number' out | awk '{ print $3 }' > cond.out

lines=`wc -l cond.out | awk '{ print $1 }'`

grep 'Residual PIPS' out | awk '{ print $3 }' > resPIPS.out
grep 'Residual Julia' out | awk '{ print $3 }' > resJulia.out

grep 'Weighed residual PIPS' out | awk '{ print $4 }' > wresPIPS.out
grep 'Weighted residual Julia' out | awk '{ print $4 }' > wresJulia.out

grep 'Error bound PIPS' out | awk '{ print $4 " " $5 }' > errorbound.out

grep 'Diff RHS' out | awk '{ print $3 }' > diffRHS.out

grep 'Cosmin term 1' out | awk '{ print $4 }' > term1.out
grep 'Cosmin term 2' out | awk '{ print $4 }' > term2.out

rm -f count.out

start=1
for i in $(eval echo "{$start..$lines}") ; do echo $i >> count.out ; done

pr -mts count.out cond.out > ./plot/cond.dat
pr -mts count.out resPIPS.out resJulia.out > ./plot/res.dat
pr -mts count.out wresPIPS.out wresJulia.out > ./plot/wres.dat
pr -mts count.out errorbound.out > ./plot/errorbound.dat
pr -mts count.out diffRHS.out > ./plot/rhs.dat
pr -mts count.out term1.out term2.out > ./plot/term.dat

cd plot
gnuplot -c cond.cfg
gnuplot -c errorbound.cfg
gnuplot -c res.cfg
gnuplot -c wres.cfg 
gnuplot -c rhs.cfg 
gnuplot -c term.cfg 

