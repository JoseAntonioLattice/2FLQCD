reset

#set terminal png

min = -30
max =  30

n = 50#floor(sqrt(10000))
width = 1 #(max - min)/n

jconf = 85
Lt = 32
Lx = 16
V = Lt*Lx**3
beta = 5.2

set title sprintf("V = %d×%d^3, {/Symbol b = %1.1f}, t/a^2 = %1.2f",Lt,Lx,beta,0.02*jconf) font ",16"


hist(x,width) =  width*(floor((x-min)/width)+0.5) + min
set boxwidth width
set style fill solid 0.5 # fill style

set size sq
set xlabel "Q" font ",16"

#do for [i=101:600]{
#    set output sprintf("%.3d.png",i) 
p "topchar_hist.dat" i jconf u (hist($1,width)):(1.0) smooth freq w boxes lc rgb "red" notitle
#}
#unset output
