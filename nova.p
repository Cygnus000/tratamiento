set terminal pngcairo size 1200,750 enhanced font 'Verdana,23'
set output 'nova.png'
set title "Modelo de Stepanova"
set xlabel "Tiempo (d)"
set ylabel "concentración (células/litro)"
set grid
plot "nova.dat" w l lt 6 lw 5 title "Carga Tumoral","nova.dat" u 1:3 w l lt 1 lw 5 title "Células T-killer"

