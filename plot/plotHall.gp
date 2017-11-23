
set yrange [-10.0:10.0]

plot 'conductivities.dat' u 2:($3*$3/$4) w lp, x, \
     'conductivities.dat' u 2:($3*$3/(2.*$6/3.)) w lp, \
     'conductivities.dat' u 2:($3*$3/((2.*$6/3.)+2*($7)/3.)) w lp
     
     
#plot 'conductivities.dat' u 2:($3*$3/$4) w lp, x, \
#     'conductivities.dat' u 2:($3*$3/(2.*($6-$7)/3.)) w lp, \
#     'conductivities.dat' u 2:($3*$3/((2.*($6-$7)/3.)+4*($8-$9)/3.)) w lp


     
pause -1
