reset
set datafile separator ','
set term pdf enhanced lw 1.0
set output 'resultados/viremiaCurvas.pdf'
set xlabel '{/*1.4 tempo (dias)}'
set ylabel '{/*1.4 Concentration}'
set style data lines
set key right top
set key vertical
set grid
plot 'outputAz.csv' u 1:2 lc "blue", 'outputDark.csv' u 1:2 lc "dark-chartreuse", 'outputVer.csv' u 1:2 lc "red",'outputCyan.csv' u 1:2 lc "cyan", 'outputMag.csv' u 1:2 lc "violet", 'outputAm.csv' u 1:2 lc "light-goldenrod"
