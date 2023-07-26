reset
set datafile separator ','
set term pdf enhanced lw 1.0
set output 'resultados/anticorposCurvas.pdf'
set xlabel '{/*1.4 tempo (dias)}'
set ylabel '{/*1.4 Concentration}'
set style data lines
set key right top
set key vertical
set grid
plot 'outputCs.csv' u 1:13 lc "blue", 'outputI1.csv' u 1:13 lc "dark-chartreuse", 'outputI2.csv' u 1:13 lc "red",'outputI3.csv' u 1:13 lc "cyan", '../Cenario1/output.csv' u 1:13 lc "violet"
