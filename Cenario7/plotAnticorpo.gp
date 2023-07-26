reset
set datafile separator ','
set term pdf enhanced lw 1.0
set output 'anticorpos.pdf'
set xlabel '{/*1.4 tempo (dias)}'
set ylabel '{/*1.4 Concentration}'
set style data lines
set key left top
set key horizontal
set grid
plot 'output.csv' u 1:13 
