clc
disp('file: calc3.m')

num_Ra = 11.1;
num_La = 1.52e-3;
num_km = 0.0582;
num_ks = 17.17;
num_ke = 0.05822;
num_Ja = 4.46e-6;
num_wn = 463.91;
num_in = 0.804;
num_un = 36;
num_TI = 0.0145


%% 3/a
disp('3/a')

syms s Ra La km ke Ja TI P;


Wel = 1/(Ra + La*s);
Wme = 1/(Ja*s);


Wx_p = Wel*Wme*km;
Wf_p = ke;
Wp = simplify(Wx_p/(1+Wx_p*Wf_p));

Wc = P*(1 + 1/(TI*s));

Wx = Wc*Wp;
Wo = Wx/(1+Wx);
Wo = simplify(Wo);
disp('Closed loop transfer function')
latex(Wo)

[Wo_num, Wo_den] = numden(Wo);
disp('Karakterisztikus egyenlet')
latex(Wo_den)

num_Wo_den = subs(Wo_den, [Ja, Ra, TI, La, ke, km], [num_Ja, num_Ra, num_TI, num_La, num_ke, num_km])
disp('Poles')
pi = solve(num_Wo_den, s);
pi = vpa(pi)

arr_P = linspace(-20, 40, 1000);
fun_p1 = symfun(pi(1), P);
fun_p2 = symfun(pi(2), P);
fun_p3 = symfun(pi(3), P);
hold on
plot(arr_P, real(fun_p1(arr_P)))
plot(arr_P, real(fun_p2(arr_P)))
plot(arr_P, real(fun_p3(arr_P)))
legend('p1', 'p2', 'p3');grid;
xlim([-1 1]);
ylim([-1e3 1e3]);
hold off
title('')
xlabel('körerősítés')
ylabel('pólus valós része')

%==============================================


% clc
% disp('file: num_calc3')

% s = tf('s');
% P = 4.063
% TI = 0.0145
% parameters


% Wel = 1/(Ra + La*s);
% Wme = 1/(Ja*s);


% %% 3/a
% disp('3/a')


% Wx_p = Wel*Wme*km;
% Wf_p = ke;
% % Wp = minreal(ke * Wx_p/(1+Wx_p*Wf_p));
% Wp = minreal(Wx_p/(1+Wx_p*Wf_p));

% Wc = P*(1 + 1/(TI*s));

% Wx = minreal(Wc*Wp);
% Wo = Wx/(1+Wx);
% Wo = minreal(Wo)

% % bode(Wo, [1, 1e6]);grid;
% % margin(Wx);grid;

% % impulse(Wo*wn/2);grid;
% step(Wo*wn/2);grid;
% title('')
% ylabel('szögsebesség (rad/s)')
