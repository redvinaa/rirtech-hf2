clc
disp('file: calc4.m')

num_Ra = 11.1;
num_La = 1.52e-3;
num_km = 0.0582;
num_ks = 17.17;
num_ke = 0.05822;
num_Ja = 4.46e-6;
num_wn = 463.91;
num_in = 0.804;
num_un = 36;
num_TD = 1.3825e-4;
num_n  = 40.827;


%% 4/a
disp('4/a')

syms s Ra La km ke Ja TD P n;


Wel = 1/(Ra + La*s);
Wme = 1/(Ja*s);


Wo_p = Wel*Wme*km;
Wf_p = ke;
Wp = simplify(Wo_p/(1+Wo_p*Wf_p));

Wc = P*(TD*s + 1)/(n*TD*s + 1);

Wo = Wc*Wp;
Wx = Wo/(1+Wo);
Wx = simplify(Wx);
disp('Closed loop transfer function')
latex(Wx)

[Wx_num, Wx_den] = numden(Wx);
disp('Karakterisztikus egyenlet')
latex(simplify(Wx_den))

num_Wx_den = subs(Wx_den, [Ja, Ra, TD, n, La, ke, km], [num_Ja, num_Ra, num_TD, num_n, num_La, num_ke, num_km])
disp('Poles')
pi = solve(num_Wx_den, s);
pi = vpa(pi)

arr_P = linspace(-2, 1, 1000);
fun_p1 = symfun(pi(1), P);
fun_p2 = symfun(pi(2), P);
fun_p3 = symfun(pi(3), P);
hold on
plot(arr_P, real(fun_p1(arr_P)))
plot(arr_P, real(fun_p2(arr_P)))
plot(arr_P, real(fun_p3(arr_P)))
legend('p1', 'p2', 'p3');grid;
xlim([-2 1]);
ylim([-200 100]);
hold off
title('')
xlabel('körerősítés')
ylabel('pólus valós része')

%==============================================


% clc
% disp('file: num_calc4')

% s = tf('s');
% n = 0.05
% TD = 1.3825e-4
% P = 40.827
% parameters


% Wel = 1/(Ra + La*s);
% Wme = 1/(Ja*s);


% %% 4/a
% disp('4/a')


% Wx_p = Wel*Wme*km;
% Wf_p = ke;
% Wp = minreal(Wx_p/(1+Wx_p*Wf_p));

% Wc = P*(TD*s + 1)/(n*TD*s + 1);

% Wx = minreal(Wc*Wp);
% Wo = Wx/(1+Wx);
% Wo = minreal(Wo)

% % bode(Wo, [1, 1e6]);grid;
% % margin(Wx);grid;

% impulse(Wo*wn/2);grid;
% % step(Wo*wn);grid;
% title('')
% ylabel('szögsebesség (rad/s)')
