set term post eps enh color solid
set output 'plot.eps'
set xlabel 'N'
set title 'real/complex pade versus eig matrix exponentials'
set ylabel 'time'
set key left
plot 'test.out' u 1:2 ti 'real eig','' u 1:3 ti 'real pade','' u 1:4 ti 'complex eig','' u 1:5 ti 'complex pade'
set output 'ratio.eps'
set xlabel 'N'
set title 'real/complex pade versus eig matrix exponentials'
set ylabel 'pade/eig'
set key left
plot 'test.out' u 1:($3/$2) ti 'real','' u 1:($5/$4) ti 'complex'
