clc
clear
digits 10
disp('file: calc6.m')


%% 6/a
disp('6/a')

% {{{ Előző feladatból
syms xi T wn s

a2 = 1;
dv_num = .1
T_num = .05
th3_num = .03

un = 36;
Psi = 618.34/un;
w_nom = 4430*(2*pi)/60;
w_noload = 5860*(2*pi)/60;

T1 = 0.0145;
T2 = 1.3825e-4;

% xi
eq = exp(-xi*pi / sqrt(1-xi^2)) - dv_num;
xi = solve(eq, xi);
xi = vpa(xi( vpa(xi)>0 ))

% wn
eq = log(1/th3_num) / (xi*wn) - T_num
wn = solve(eq, wn);
wn = double(wn( double(wn)>0 ))

beta = wn*xi
wd = wn*sqrt(1-xi^2)

p1 = double(-beta+1i*wd)
p2 = double(-beta-1i*wd)
%}}}

p3 = 3* real(p1)


%% 6/b
disp('6/b')

kar_pol = (s-p1)*(s-p2)*(s-p3)
coef = coeffs(kar_pol, s)

A_new = double([0 1 0; 0 0 1; -abs(coef(1)) -abs(coef(2)) -abs(coef(3))])
B = [0; 0; 1]

p1_orig = -1/T1
p2_orig = -1/T2

A = [0 1 0; 0 0 1; 0 -abs(p1_orig*p2_orig + p1_orig*p3) -abs(p1_orig+p2_orig)]

Kx = sym('Kx', [1 3])
eq = A_new - (A - B*Kx)
sol = solve(eq)
Kx = [sol.Kx1 sol.Kx2 sol.Kx3]
pause


%% 6/c
disp('6/c')

un = 36;
Psi = 618.34/un;
w_nom = 4430*(2*pi)/60;
w_noload = 5860*(2*pi)/60;

C = [Psi/(T1*T2) 0 0]
D = [0]

sys = ss(A_new, B, C, D);

step(sys*2*pi);grid;title('')
ylabel('Szögelfordulás (rad)')
stepinfo(sys*2*pi)
pause


%% 6/d
disp('6/d')

phi_lim = dcgain(sys) * 2*pi
pause


%% 6/e
disp('6/e')

Kr = 1/dcgain(sys)
B_new = B*Kr;
sys = ss(A_new, B_new, C, D);
pause


%% 6/f
disp('6/f')

step(sys*2*pi);grid;title('')
ylabel('Szögelfordulás (rad)')
stepinfo(sys*2*pi)
pause


%% 6/g
disp('6/g')

phi_lim = dcgain(sys) * 2*pi
pause


%% 6/h
disp('6/h')

Kaj = inv([A B; C D])*[0;0;0;1];
Krx = Kaj(1:3)
Kru = Kaj(4)
pause


%% 6/i
disp('6/i')

close

num = 1/s * B * (Kru + Kx*Krx)
den = eye(3) - 1/s * (A - B*Kx)

W_sym = C * inv(den) * num;
W_sym = simplify(vpa(W_sym));
pretty(W_sym)

[n,d] = numden(W_sym);
n = sym2poly(n);
d = sym2poly(d);
W = tf(n, d);

step(W*2*pi);grid;title('')
ylabel('Szögelfordulás (rad)')
stepinfo(W*2*pi)
pause


%% 6/j
disp('6/j')

y_lim = dcgain(W*2*pi)
