reset


a = 8./32
b = 16./32
jconf = 150
Lt = 32
Lx = 16
V = Lt*Lx**3

f(x) = chi*V*(1-x)*x


set key font ",12"
set xlabel "x" font ",14"
set ylabel "<q^2 - Q^2x^2>" font ",14"

set object 1 rect from a,0 to b,10   
set object 1 fillstyle solid fillcolor "light-green"
set object 1 behind
fit [a:b] f(x) "WF_slab_new.dat" i jconf u 1:2:3 via chi
p "WF_slab_new.dat" i jconf u 1:2:3 w e notitle, f(x) t "{/Symbol c}_tVx(1-x)" 
