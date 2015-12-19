set term x11
# set logscale
# set log y
set xrange [0:3.5E3]
set xrange [0:.35E3]
set xrange [0:100.E3]
set yrange [0:0.04]
set xlabel "E_p [keV]"
set ylabel "A [MeV/micron]"
set label "Te=Ti=10 keV ne=1.E25  E_p=100 keV" at 150, 0.045

f(x)=0

plot  \
"main.dat" using 2:3 with lines lt 1 title "A_tot", \
"main.dat" using 2:4 with lines lt 2 title "A_e",   \
"main.dat" using 2:5 with lines lt 3 title "A_I"
#f(x)



#
# create postscript file of plot
#
#set term postscript eps color 
#set output "gr008.eps"
#plot  \
#"main.dat" using 2:3 with lines lt 1 title "A_tot", \
#"main.dat" using 2:4 with lines lt 2 title "A_e",   \
#"main.dat" using 2:5 with lines lt 3 title "A_I"


pause -1
