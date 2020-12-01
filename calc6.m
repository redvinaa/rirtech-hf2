clc
clear
digits 10
disp('file: calc6.m')


%% 6/a
disp('6/a')

% W = b0 / (a2*s^2 + a1*s + a0)
% T = sqrt(a2/a0)
% xi = a1 / (2*sqrt(a0*a2))
% dv = exp(-xi*pi / sqrt(1-xi^2))

T1 = 0.0145;
T2 = 1.3825e-4;
A = [0 1; -1/(T1*T2) -(T1+T2)/(T1*T2)]
A_new = [0 1 0; -1/(T1*T2) -(T1+T2)/(T1*T2) 0; -1 -1 0]

B_new = [0; 1; 0]

a0 = 400
a1 = 23.65

p2 = roots([1 a1 a0])
p1 = p2(1)
p2 = p2(2)
a0 = 400
a1 = 23.6462

p2_ = roots([1 a1 a0])
p1_ = p2_(1)
p2_ = p2_(2)

p1-p1_
p2-p2_


pause

syms a0 a1 % a2
syms xi T
a2 = 1;
dv_num = .1
T_num = .05

% xi
eq = exp(-xi*pi / sqrt(1-xi^2)) - dv_num;
xi = solve(eq, xi);
xi = vpa(xi( vpa(xi)>0 ))

% a0
eq = sqrt(a2/a0) - T_num;
a0 = solve(eq, a0);
a0 = double(a0)

% a1
eq = xi - a1 / (2*sqrt(a0*a2));
a1 = solve(eq, a1);
a1 = double(a1)


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


%% 6/c
disp('6/c')

un = 36;
Psi = 618.34/un;
w_nom = 4430*(2*pi)/60;
w_noload = 5860*(2*pi)/60;

Kx = [k1 k2]
A = [0 1; -1/(T1*T2) -(T1+T2)/(T1*T2)]
B = [0; 1]
C = [Psi/(T1*T2) 0]
D = [0]

A_new = [A - B*Kx B*KI; -C 0]

sys = ss(A_new, B, C, D);

opt = stepDataOptions;
opt.StepAmplitude = w_noload;

step(sys, opt);grid;title('')
stepinfo(sys)


%% 6/d
disp('6/d')

syms s

[Num,Den] = tfdata(sys,'v')
syms s

Wcl = poly2sym(Num,s)/poly2sym(Den,s)

X = w_noload/s;
Y = Wcl * X;

y_lim = limit(vpa(Y*s), s, 0);
y_lim = double(y_lim)


%% 6/e
disp('6/e')

Kr = -inv(C*inv(A_new)*B)
% Kr = rscale(A, B, C, D, K)
B_new = B*Kr;
sys = ss(A_new, B_new, C, D);


%% 6/f
disp('6/f')

opt = stepDataOptions;
opt.StepAmplitude = w_noload;

step(sys, opt);grid;title('')
stepinfo(sys)


%% 6/g
disp('6/g')

syms s

[Num,Den] = tfdata(sys,'v')
syms s

Wcl = poly2sym(Num,s)/poly2sym(Den,s)

X = w_noload/s;
Y = Wcl * X;

y_lim = limit(vpa(Y*s), s, 0);
y_lim = double(y_lim)


%% 6/h
disp('6/h')

Kaj = inv([A B; C D])*[0;0;1];
Krx = Kaj(1:2)
Kru = Kaj(3)


%% 6/i
disp('6/i')

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


%% 6/j
disp('6/j')

syms s
% W * (r/s) * s = W*r
y_lim = limit(W_sym*r, s, 0);
y_lim = double(y_lim)
