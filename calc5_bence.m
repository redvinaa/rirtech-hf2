clc
clear

disp('file: calc5.m')

%% Kapott paraméterek
th_0=30;  % 'V'
th_1=45;  % '°'
th_2=10;  % '%'
th_3=2;   % '%'
th_4=50;  % 'ms'

U_n=30;
T_1=0.0145;
T_2=1.3719e-4;

psi=645.07/U_n;
phi_t=pi/4;

w_noload=6160/60*(2*pi);

%% 5/A feladat

syms X_i w_n_sym

dv = th_2/100;
T_be = th_4/1000;

% Zeta
eq_X_i = exp(-X_i*pi / sqrt(1-X_i^2)) - dv;
x_i = vpasolve(eq_X_i, X_i);
% Omega
eq_Omega_n = 1/(x_i*w_n_sym)*log(100/th_3) == T_be;
w_n = vpasolve(eq_Omega_n);

Beta = w_n*x_i;
w_d = double(w_n*sqrt(1-x_i^2));

p_1 = double(-Beta+1i*w_d); %PIPA
p_2 = double(-Beta-1i*w_d); %PIPA

%% 5/B feladat

%Ellen: s^2 + 7358 s + 5.044e05
A_hullam=[0         1;
          -p_1*p_2  p_1+p_2];
A = [0            1;
    -1/(T_1*T_2) -(T_1+T_2)/(T_1*T_2)];
B = [0;
     1];
C = [psi/(T_1*T_2)  0];
D = 0;

syms k_1 k_2

K= [0   0;
    k_1 k_2];

eq_A_hullam_AK = A_hullam == (A-K);
solution = solve(eq_A_hullam_AK);
k_1 = double(solution.k_1);
k_2 = double(solution.k_2);

%% 5/C feladat

sys=ss(A_hullam, B,C,D);

step(sys*w_noload);grid;
title('');
xlabel('Idő [s]');
ylabel('Szögsebesség [rad/s]');

%% 5/D feladat

dcgain(sys*w_noload);

%% 5/E feladat

k_r=-inv(C*(A_hullam\B));

%% 5/F feladat

step(sys*w_noload*k_r);grid;
title('5/F');
xlabel('Idő [s]');
ylabel('Szögsebesség [rad/s]');

%% 5/G feladat
dcgain(sys*w_noload*k_r);

%% 5/H feladat

Kaj = inv([A B; C D])*[0;0;1];
Krx = Kaj(1:2);
Kru = Kaj(3);

%% 5/I feladat
syms s

num = 1/s * B * (Kru + [k_1 k_2]*Krx);
den = eye(2) - 1/s * (A - B*[k_1 k_2]);

W_sym = C * inv(den) * num;
W_sym = simplify(vpa(W_sym));
pretty(W_sym)

[n,d] = numden(W_sym);
n = sym2poly(n);
d = sym2poly(d);
W = tf(n, d);

step(W*w_noload);grid;
title('5/I')
xlabel('Idő [s]');
ylabel('Szögsebesség [rad/s]');


%% 5/J feladat

 syms s
 y_lim = limit(W_sym*w_noload, s, 0);
 y_lim = double(y_lim)


%#========
%% 5/J feladat
%dcgain(W_cl*w_noload)




%%% 5/I feladat
% 
% close
% 
% syms x1 x2 U
% r = w_noload
% 
% num = 1/s * B * (Kru + Kx*Krx)% * r
% den = eye(2) - 1/s * (A - B*Kx)
% 
% W_sym = C * inv(den) * num;
% W_sym = simplify(vpa(W_sym));
% pretty(W_sym)
% 
% [n,d] = numden(W_sym);
% n = sym2poly(n);
% d = sym2poly(d);
% W = tf(n, d);
% 
% step(W*r);grid;title('')
% stepinfo(tf(n, d))
% 
% 
% %% 5/J feladat
% 
% syms s
% % W * (r/s) * s = W*r
% y_lim = limit(W_sym*r, s, 0);
% y_lim = double(y_lim)
