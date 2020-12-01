clc
clear
digits 10
disp('file: calc5.m')


%% 5/a
disp('5/a')

syms xi T wn
a2 = 1;
dv_num = .1
T_num = .05
th3_num = .03

% xi
eq = exp(-xi*pi / sqrt(1-xi^2)) - dv_num;
xi = solve(eq, xi);
xi = vpa(xi( vpa(xi)>0 ))

% wn
eq = log(1/th3_num) / (xi*wn) - T_num
wn = solve(eq, wn);
wn = vpa(wn( vpa(wn)>0 ))

beta = wn*xi
wd = wn*sqrt(1-xi^2)

p1 = double(-beta+1i*wd)
p2 = double(-beta-1i*wd)

a1 = -(p1+p2)
a0 = p1*p2
pause

syms k1 k2

T1 = 0.0145;
T2 = 1.3825e-4;

% k1
eq = -1/(T1*T2) - k1 + a0;
k1 = solve(eq, k1);
k1 = double(k1)

% k2
eq = -(T1+T2)/(T1*T2) - k2 + a1;
k2 = solve(eq, k2);
k2 = double(k2)
pause


%% 5/c
disp('5/c')

un = 36;
Psi = 618.34/un;
w_nom = 4430*(2*pi)/60;
w_noload = 5860*(2*pi)/60;

Kx = [k1 k2]
A = [0 1; -1/(T1*T2) -(T1+T2)/(T1*T2)]
B = [0; 1]
C = [Psi/(T1*T2) 0]
D = [0]

A_new = A - B*Kx

sys = ss(A_new, B, C, D);

opt = stepDataOptions;
opt.StepAmplitude = w_noload;

step(sys, opt);grid;title('')
stepinfo(sys)
pause


%% 5/d
disp('5/d')

syms s

[Num,Den] = tfdata(sys,'v')
syms s

Wcl = poly2sym(Num,s)/poly2sym(Den,s)

X = w_noload/s;
Y = Wcl * X;

y_lim = limit(vpa(Y*s), s, 0);
y_lim = double(y_lim)
pause


%% 5/e
disp('5/e')

Kr = -inv(C*inv(A_new)*B)
% Kr = rscale(A, B, C, D, K)
B_new = B*Kr;
sys = ss(A_new, B_new, C, D);
pause


%% 5/f
disp('5/f')

opt = stepDataOptions;
opt.StepAmplitude = w_noload;

step(sys, opt);grid;title('')
stepinfo(sys)
pause


%% 5/g
disp('5/g')

syms s

[Num,Den] = tfdata(sys,'v')
syms s

Wcl = poly2sym(Num,s)/poly2sym(Den,s)

X = w_noload/s;
Y = Wcl * X;

y_lim = limit(vpa(Y*s), s, 0);
y_lim = double(y_lim)
pause


%% 5/h
disp('5/h')

Kaj = inv([A B; C D])*[0;0;1];
Krx = Kaj(1:2)
Kru = Kaj(3)
pause


%% 5/i
disp('5/i')

close

syms x1 x2 U
r = w_noload

num = 1/s * B * (Kru + Kx*Krx)% * r
den = eye(2) - 1/s * (A - B*Kx)

W_sym = C * inv(den) * num;
W_sym = simplify(vpa(W_sym));
pretty(W_sym)

[n,d] = numden(W_sym);
n = sym2poly(n);
d = sym2poly(d);
W = tf(n, d);

step(W*r);grid;title('')
stepinfo(tf(n, d))
pause


%% 5/j
disp('5/j')

syms s
% W * (r/s) * s = W*r
y_lim = limit(W_sym*r, s, 0);
y_lim = double(y_lim)
