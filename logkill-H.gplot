set terminal pngcairo size 1200,750 enhanced font 'Verdana,23'
set output 'logkill-H.png'
set title "Modelo Efecto Proporcional con Fármaco"
set xlabel "Tiempo (d)"
set ylabel "Concentración (celulas/litro  mg/litro)"
set grid
plot "logkill-H.dat" w l lt 6 lw 5 title "Carga Tumoral", "logkill-H.dat" u 1:3 w l lt 1 lw 5 title "Fármaco"

