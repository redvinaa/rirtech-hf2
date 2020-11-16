clc
clear
disp('file: calc5.m')


%% 5/a
disp('5/a')

syms w P;
un = 36
T1 = 0.0145
T2 = 1.3825e-4
A = 618.34/un
phi_t = pi/3
w_nom = 4430/(2*pi)

p1 = -1/T1
p2 = -1/T2

s = tf('s');
Wp = A / ( (1+T1*s) * (1+T2*s) )
Hss = ss(Wp)
