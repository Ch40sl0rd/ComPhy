# in x.dat stehen: t*dt x[t] dx/dt[t] E_k
set xlab 't'
set yrange [-2:2]
plot 'x.dat' using 1:2 w l t "Position"
replot 'x.dat' using 1:3 w l t "Geschwindigkeit"
replot 'x.dat' using 1:4 w l t "Kinetische Energie"
