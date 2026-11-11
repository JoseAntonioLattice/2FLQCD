reset


a = 0.25
b = 0.75
Q = 2
jconf = 300
Lt = 32
Lx = 16
V = Lt*Lx**3

f(x) = chi*V*x + (Q**2-chi*V)*x**2


set key font ",12"
set xlabel "x" font ",14"
set ylabel "<q^2>" font ",14"
set title sprintf("Q = %d",Q) font ",18"
fit [a:b] f(x) "WF_slab.dat" i jconf u ($1/Lt):6:7 via chi
p "WF_slab.dat" i jconf u ($1/Lt):6:7 w e notitle, [a:b] f(x) t "{/Symbol c}_tVx + (Q^2-{/Symbol c}_tV)x^2" 
