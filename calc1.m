clc
clear
disp('file: calc1.m')


%% 1/a
disp('1/a')

syms w P;
un = 36
T1 = 0.0145
T2 = 1.3825e-4
A = 618.34/un
phi_t = pi/3
w_nom = 4430/(2*pi)

phi = -pi/2 - atan(T2*w);
eq = phi + pi - phi_t;
wc = solve(eq, w);
wc = abs(double(wc))


Wx = A*P/T1 * 1/(i*wc*(1+i*wc*T2));
P = solve(abs(Wx)-1, P);
P = double(P)

s = tf('s');
Wx = A*P/T1 * 1/(s*(1+s*T2))
margin(Wx);grid


%% 1/b
disp('1/b')

Wcl = Wx/(1+Wx)
step(Wcl*w_nom);grid;title('')
ylabel('Szögsebesség (rad/s)')
xlabel('Idő')


%% 1/c
disp('1/c')

syms s

Wx = A*P/T1 * 1/(s*(1+s*T2));
Wcl = Wx/(1+Wx);

X = w_nom/s;
Y = Wcl * X;

w_inf = limit(Y*s, s, 0);
w_inf = double(w_inf)
