clc
clear all
digits 6
disp('file: calc6.m')

num_Ra = 11.1;
num_La = 1.52e-3;
num_km = 0.0582;
num_ks = 17.17;
num_ke = 0.05822;
num_Ja = 4.46e-6;
num_wn = 463.91;
num_in = 0.804;
num_un = 36;
% num_TI = 0.0145;
num_T0 = .005;
num_n = 0.05;
num_TD = 1.3825e-4;
num_P  = 40.827;


syms s Ra La km ke Ja TD n P T0 tau0 wn t;


Wel = 1/(Ra + La*s);
Wme = 1/(Ja*s);

%------------------------------

Uw = wn/(2*s);
Ut = tau0/s * exp(-T0*s);

Wx = Wel*km*Wme
Wf = ke
Wp = simplify(Wx/(1+Wx*Wf));

Wx = -Wme;
Wf = -km*ke*Wel;
Wtw = simplify(Wx/(1+Wx*Wf));

Wc = P*(TD*s + 1)/(n*TD*s + 1);


%{{{ Time responses %}

Yw = (Wp*Wc*Uw + Wtw*Ut)  /  (Wp*Wc + 1);
Yw = simplify(Yw)

Yi = ke*Wel*(Uw - Yw);
Yi = simplify(Yi);

num_Yw = subs(Yw, [Ra, La, km, ke, Ja, TD, n, P, T0, wn], [num_Ra, num_La, num_km, num_ke, num_Ja, num_TD, num_n, num_P, num_T0, num_wn]);
num_Yw = vpa(num_Yw);
yw = ilaplace(num_Yw, s, t);
fun_yw = symfun(yw, [t tau0])

num_Yi = subs(Yi, [Ra, La, km, ke, Ja, TD, n, P, T0, wn], [num_Ra, num_La, num_km, num_ke, num_Ja, num_TD, num_n, num_P, num_T0, num_wn]);
num_Yi = vpa(num_Yi)
yi = ilaplace(num_Yi, s, t)
fun_yi = symfun(yi, [t tau0])
%}}}


%{{{ tau_max %}

% t0 = linspace(num_T0+0.005, num_T0+0.005, 70);
% tau_i = linspace(30, 35, 50);

% imax = 0;
% taumax = 0;
% t0max = 0;
% for i = 1:length(tau_i)
%     for j = 1:length(t0)
%         I = abs(fun_yi(t0(j), tau_i(i)));
%         if I > imax
%             imax = I;
%             taumax = tau_i(i);
%             t0max = t0(j);
%         end
%         if I > num_in
%             disp('found max')
%             disp('imax:')
%             imax
%             disp('taumax')
%             taumax
%             disp('tmax')
%             t0max
%             % return
%             break
%         end
%     end
%     if I > num_in
%         break
%     end
% end

taumax = 32.86

t0 = linspace(0, .01, 1000);
tau_i = taumax * [1/3., 2/3., 1.]
%}}}

%{{{ Plot 5-w %}

hold on
for i = 1:length(tau_i)
	plot(t0, fun_yw(t0, tau_i(i)))
end
xlabel('idő (s)');
ylabel('szögsebesség (rad/s)');
legend('10,96', '21,91', '32,68');
grid;
hold off;
pause;
%}}}

%{{{ Plot 5-i %}

close
hold on
for i = 1:length(tau_i)
	plot(t0, fun_yi(t0, tau_i(i)))
end
xlabel('idő (s)');
ylabel('armatúra áram (A)');
legend('10,96', '21,91', '32,68');
grid;
hold off;
%}}}
