set terminal x11

mylw = 2
index_check = 24

set xlabel 'Temperature [K]'

set ylabel 'Cp/R'

plot 'prob_iso_octane_therm.txt' index index_check using 1:2 \
     title 'original' \
     with lines lt 1 lw mylw lc rgb "#ff0000",\
     '' index index_check using 1:5 \
     title 'modified' \
     with lines lt 1 lw mylw lc rgb "#0000ff" 


pause -1
