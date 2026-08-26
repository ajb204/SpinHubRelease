set term post eps enh color solid
set output 'plot.eps'
set xlabel 'N'
set title 'real/complex pade versus eig matrix exponentials'
set ylabel 'time'
set key left
plot 'test.out' u 1:2 ti 'real eig','' u 1:3 ti 'real pade','' u 1:4 ti 'real pad expkit','' u 1:5 ti 'complex eig','' u 1:6 ti 'complex complex pade','' u 1:7 ti 'complex pade expo'
set output 'ratio.eps'
set xlabel 'N'
set title 'real/complex pade versus eig matrix exponentials'
set ylabel 'New/OldPade'
set key left
plot 'test.out' u 1:($3/$2) ti 'real eig','' u 1:($3/$4) ti 'real pade expo','' u 1:($6/$5) ti 'complex eig','' u 1:($6/$7) ti 'complex pade expkit' 
