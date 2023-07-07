reset
set datafile separator ','
set term pdf enhanced lw 1.0
set output 'viremia.pdf'
set xlabel '{/*1.4 tempo (dias)}'
set ylabel '{/*1.4 Viremia}'
#set xrange [0:0.1]
set style data lines
set key left top
set key horizontal
set grid
plot 'output.csv' u 1:2 
