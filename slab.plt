reset


a = 8./32
b = 16./32
jconf = 85
Lt = 32
Lx = 16
V = Lt*Lx**3
beta = 5.2

f(x) = chi*V*(1-x)*x


set key font ",11"
set xlabel "x" font ",14"
set ylabel "⟨q^2 - Q^2x^2⟩" font ",14"
set title sprintf("V = %d×%d^3, {/Symbol b = %1.1f}, t/a^2 = %1.2f ≈ t_0/a^2",Lt,Lx,beta,0.02*jconf) font ",16"

set object 1 rect from a,0 to b,10   
set object 1 fillstyle solid fillcolor "light-green"
set object 1 behind
fit [a:b] f(x) "WF_slab_new.dat" i jconf u 1:2:3 via chi
p "WF_slab_new.dat" i jconf u 1:2:3 w e notitle, [a:b]f(x) lw 2 lc 7 t sprintf("{/Symbol c}_tVx(1-x)"), keyentry t sprintf("{/Symbol c}_t = 1.29(7)×10^{-4}"),\
    [0:a]f(x) lw 2 lc 7 dt 3 notitle, [b:1]f(x) lw 2 lc 7 dt 3 notitle
 #, keyentry t sprintf("{/Symbol c}^2/dof = %1.4f",FIT_STDFIT**2) 
