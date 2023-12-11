set terminal pngcairo size 1200,750 enhanced font 'Verdana,23'
set output 'logkill.png'
set title "Modelo Efecto Proporcional"
#unset key
set xlabel "Tiempo (d)"
set ylabel "Concentración (celulas/litro  mg/litro)"
set grid
plot "logkill.dat" w l lt 6 lw 5 title "Carga Tumoral", "logkill.dat" u 1:3 w l lt 1 lw 5 title "Fármaco"

