reset
set datafile separator ','
set term pdf enhanced lw 2.0
set output 'DosesAnticorpos.pdf'
set xlabel '{/*1.4 tempo (dias)}'
set ylabel '{/*1.4 Anticorpos - mUI/ml}'
set xrange [0:100]
set style data lines
set key right top
set key vertical
set grid
plot 'outputd6.csv' u 1:13 title 'V0 = 0.82 ou 31 UI' lc 'blue', 'outputd5.csv' u 1:13 title 'V0 = 4.16 ou 158 UI' lc 'dark-spring-green', 'outputd4.csv' u 1:13 title 'V0 = 15.47 ou 587 UI' lc 'red', 'outputd3.csv' u 1:13 title 'V0 = 79.39 ou 3013 UI' lc 'cyan', 'outputd2.csv' u 1:13 title 'V0 = 275.28 ou 10447 UI' lc 'dark-magenta', 'outputd1.csv' u 1:13 title 'V0 = 724 ou 27476 UI' lc 'dark-yellow'
