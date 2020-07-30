set terminal x11

mylw = 2

set xlabel 'Temperature [K]'

set ylabel 'S/R'

plot 'thermo_check.stdout' index 0 using 1:4 \
     title 'original' \
     with lines lt 1 lw mylw lc rgb "#ff0000",\
     '' index 1 using 1:4 \
     title 'modified' \
     with lines lt 1 lw mylw lc rgb "#0000ff" 


pause -1
