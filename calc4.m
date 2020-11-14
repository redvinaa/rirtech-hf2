clc
clear
disp('file: calc4.m')


%% 4/a
disp('4/a')

syms s w P TI TD T1 T2 A;
un = 36;
T1_num = 0.0145;
T2_num = 1.3825e-4;
A_num = 618.34/un;
n = 0.1;
phi_t = pi/3;
w_nom = 4430/(2*pi);

Wc = P * ( 1 + 1/(s*TI) + s*TD/(1+s*n*TD) )
Wp = A/( (1+s*T1) * (1+s*T2) );
Wx = simplify(Wc*Wp / s);
Wcl = Wx/(1+Wx);

pretty(Wx)
latex(Wx)

phi = -pi - atan(T1*w) - atan(T2*w) + atan(TI*w)
phi_num = -pi - atan(T1_num*w) - atan(T2_num*w) + atan(TI*w);
phi_d = diff(phi, w)
phi_d_num = diff(phi_num, w);
wci = solve(phi_d_num, w);

wci = simplify(wci)

return

digits 10
TIi = zeros(4, 1);

% phi_subs = subs(phi, w, wci(1))
% eq = pi - phi_t + phi_subs
% solve(vpa(eq), TI)

% phi_subs = subs(phi, w, wci(2))
% eq = pi - phi_t + phi_subs
% solve(vpa(eq), TI)

% phi_subs = subs(phi, w, wci(3))
% eq = pi - phi_t + phi_subs
% solve(vpa(eq), TI)

% csak a wci(4) megoldás helyes, a többinek vagy nincs megoldása,
% vagy w<0, vagy w nem valós eredményt ad, ami értelmetlen
phi_subs = subs(phi_num, w, wci(4))
eq = pi - phi_t + phi_subs
TI_num = double(solve(vpa(eq), TI))
w_num = double(subs(wci(4), [TI], [TI_num]));
assert( w_num >= 0);       % pipa
assert( imag(w_num) == 0); % pipa

Wx_num = subs(Wx, [T1 T2 TI s], [T1_num T2_num TI_num w_num]);
P_num = solve(abs(Wx_num)-1, P);
P_num = double(P_num)

s = tf('s')
Wc = P_num* (1+s*TI_num)/(s*TI_num)
Wp = A/( (1+s*T1_num) * (1+s*T2_num) );
Wx = minreal(Wc*Wp / s);
Wcl = Wx/(1+Wx);

% margin(Wx);grid


%% 4/b
disp('4/b')

step(2*pi*Wcl)
title('')
grid


%% 4/c
disp('4/c')

syms s
Wc = P_num* (1+s*TI_num)/(s*TI_num)
Wp = A/( (1+s*T1_num) * (1+s*T2_num) );
Wx = simplify(Wc*Wp / s);
Wcl = Wx/(1+Wx);

X = 2*pi/s;
Y = Wcl * X;

theta_inf = limit(Y*s, s, 0);
theta_inf = double(theta_inf)
