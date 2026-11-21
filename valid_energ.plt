set terminal pngcairo size 1000,600 enhanced font "Arial,12"
set output "Validacion_Energia_T0.500.png"

set xlabel "Tiempo (MCS)"
set ylabel "Energía por espín e(t)"
set title "Validación: Evolución de la energía (T = 0.5)"
set xrange [0:500]
set yrange [-0.90:-0.68]
set grid

# 'every 5' pinta 1 de cada 5 puntos para que no se apelotonen las barras de error
plot "Comparacion_N256_T0.500_MCS500_R200.txt" every 5 using 1:7:8 with yerrorbars pt 7 ps 0.7 lc rgb "magenta" title "Metrópolis", \
     "" every 5 using 1:9:10 with yerrorbars pt 5 ps 0.7 lc rgb "blue" title "n-fold way", \
     "" using 1:11 with lines dt 2 lw 2 lc rgb "dark-violet" title "Energía de equilibrio"

unset output