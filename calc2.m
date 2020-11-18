clc
clear
disp('file: calc2.m')


%% 2/a
disp('2/a')

syms w P;
un = 36;
T1 = 0.0145;
T2 = 1.3825e-4;
Psi = 618.34/un;
n = 0.1;
phi_t = pi/3;
w_nom = 4430/(2*pi);

phi = - atan(T1*w) - atan(n*T2*w);
eq = phi + pi - phi_t;
wc = solve(eq, w);
wc = double(wc)


Wx = Psi*P/( (1+i*wc*T1) * (1+i*wc*n*T2) );
P = solve(abs(Wx)-1, P);
P = double(P)

s = tf('s');
Wx = Psi*P/( (1+s*T1)* (1+s*n*T2) )
margin(Wx);grid


%% 2/b
disp('2/b')

Wcl = Wx/(1+Wx)
step(Wcl*w_nom);grid;title('')
ylabel('Szögsebesség (rad/s)')
xlabel('Idő')


%% 2/c
disp('2/c')

syms s

Wx = Psi*P/ ( (1+s*T1) * (1+s*n*T2) );
Wcl = Wx/(1+Wx);

X = w_nom/s;
Y = Wcl * X;

w_inf = limit(Y*s, s, 0);
w_inf = double(w_inf)
