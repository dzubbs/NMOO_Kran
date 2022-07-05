t1= 0.699
t2= 0.894
tN= 1.592
N=20
N1 = round(N*t1/tN)
N2 = round(N*(t2-t1)/tN)
N3 = round(N*(tN-t2-t1)/tN)
t_gitter = [linspace(0,t1,10)]
