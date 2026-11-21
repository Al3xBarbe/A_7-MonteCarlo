set terminal pngcairo size 1000,600 enhanced font "Arial,12"
set output "Validacion_C0_T0.500.png"

set xlabel "Tiempo (MCS)"
set ylabel "C_0(t) / C_0(0)"
set title "Validación: Relajación de C_0(t) (T = 0.5)"
set xrange [0:50]
set yrange [-0.1:1.1]
set grid

plot "Comparacion_N256_T0.500_MCS500_R200.txt" using 1:2:3 with yerrorbars pt 7 ps 0.7 lc rgb "magenta" title "Metrópolis", \
     "" using 1:4:5 with yerrorbars pt 5 ps 0.7 lc rgb "blue" title "n-fold way", \
     "" using 1:6 with lines lw 2 lc rgb "dark-violet" title "Predicción teórica"

unset output