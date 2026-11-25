reset

Lt = 32
Lx = 16
V = Lt*Lx**3
beta = 5.2

set key font ",11"
set xlabel "t/a^2" font ",14"
set ylabel "a^4{/Symbol c}_t" font ",14"
set format y "%.1e"
set title sprintf("V = %d×%d^3, {/Symbol b = %1.1f}",Lt,Lx,beta) font ",16"

set arrow from 1.70, graph 0 to 1.70, graph 1 nohead dt 2 lw 1.5

set label "t_0/a^2" at 2,0.00008 center font ",14"

set key left font ",14"
plot 'chi_t.dat' every 5 u 1:2:3 w yerrorbars t 'slab', '' every 5 u 1:4:5 w yerrorbars t '⟨Q^2⟩/V'
