reset

#set terminal png

min = -30
max =  30

n = 50#floor(sqrt(10000))
width = 1 #(max - min)/n

hist(x,width) =  width*(floor((x-min)/width)+0.5) + min
set boxwidth width
set style fill solid 0.5 # fill style

set size sq
set xlabel "Q" font ",16"

#do for [i=101:600]{
#    set output sprintf("%.3d.png",i) 
p "topchar_hist.dat" i 0 u (hist($1,width)):(1.0) smooth freq w boxes lc rgb "red"
#}
#unset output
