clc
clear

%%

addpath("given_functions")
addpath("Systems/PEECmodel")
addpath("Systems/ButterflyGyro-dim1e5-gyro")

%%
%spiral inductor

[peec.E, rows, cols, entries] = mmread('spiral_inductor_peec.E');

[peec.A, rows, cols, entries] = mmread('spiral_inductor_peec.A');

[peec.B, rows, cols, entries] = mmread('spiral_inductor_peec.B');

%%
E = peec.E;
A = peec.A;
B = peec.B;

%%
m = 3;
tol = 1e-2;

sigma0 = -linspace(1e7,1e10,m);

freq = 10.^[1:0.1:10];
s = i*2*pi()*freq;

[Ahat3, Ehat3, bhat3, chat3, V3] = irka(A, E, B, B, sigma0, tol);
resp_irka_peec3 = bode_from_system(Ahat3,Ehat3, bhat3, chat3,s);

%%


resp_spiral_inductor = bode_from_system(A,E,B,B,s);


%%
m = 4;
sigma0 = -linspace(1e7,1e10,m);
[Ahat4, Ehat4, bhat4, chat4, V4] = irka(A, E, B, B, sigma0, tol);
resp_irka_peec4 = bode_from_system(Ahat4,Ehat4, bhat4, chat4,s);

%%
m = 6;
sigma0 = -linspace(1e7,1e10,m);
[Ahat5, Ehat5, bhat5, chat5, V] = irka(A, E, B, B, sigma0, tol);
resp_irka_peec5 = bode_from_system(Ahat5,Ehat5, bhat5, chat5,s);
%%
figure(31)

loglog(freq, abs(resp_irka_peec3), '-g'); % Amplitude plot
hold on
loglog(freq, abs(resp_irka_peec4), '-r');
loglog(freq, abs(resp_irka_peec5), '-k');
loglog(freq, abs(resp_spiral_inductor), '-b');
title('Amplitude Response');
xlabel('frequency')
ylabel('|H(s)|')
legend('irka3','irka4','irka6','original')
grid on;

figure(32)
semilogx(freq, angle(resp_irka_peec3), '-g'); % Phase plot
hold on
semilogx(freq, angle(resp_irka_peec4), '-r');
semilogx(freq, angle(resp_irka_peec5), '-k');
semilogx(freq, angle(resp_spiral_inductor), '-b');
title('Phase Response');
xlabel('frequency')
ylabel('arg(H(s))')
legend('irka3','irka4','irka6','original')
grid on;

%%
H_diff_inf3 = max(abs(resp_irka_peec3-resp_spiral_inductor))/max(abs(resp_spiral_inductor)) %0.0012
H_diff_inf4 = max(abs(resp_irka_peec4-resp_spiral_inductor))/max(abs(resp_spiral_inductor)) %2.1267e-04
H_diff_inf5 = max(abs(resp_irka_peec5-resp_spiral_inductor))/max(abs(resp_spiral_inductor)) %1.4714e-05

%%
[resis_spiral_peec, induc_spiral_peec] = res_ind(E,A,B);
[resis3, induc3] = res_ind(Ehat3,Ahat3, bhat3);
[resis4, induc4] = res_ind(Ehat4,Ahat4, bhat4);
[resis5, induc5] = res_ind(Ehat5,Ahat5, bhat5);

freq = 10.^[1:0.1:10];

figure(1);
subplot(3,2,1); semilogx(freq,resis_spiral_peec,'-',freq,resis3,'o'); 
title('Resistance'); legend(cat(2,'spiral-peec-',num2str(entries)),cat(2,'irka-reduc-3'));
subplot(3,2,2); semilogx(freq,induc_spiral_peec,'-',freq,induc3,'o'); 
title('Inductance'); legend(cat(2,'spiral-peec-',num2str(entries)),cat(2,'irka-reduc-3'));
subplot(3,2,3); semilogx(freq,resis_spiral_peec,'-',freq,resis4,'o'); 
title('Resistance'); legend(cat(2,'spiral-peec-',num2str(entries)),cat(2,'irka-reduc-4'));
subplot(3,2,4); semilogx(freq,induc_spiral_peec,'-',freq,induc4,'o'); 
title('Inductance'); legend(cat(2,'spiral-peec-',num2str(entries)),cat(2,'irka-reduc-4'));
subplot(3,2,5); semilogx(freq,resis_spiral_peec,'-',freq,resis5,'o'); 
title('Resistance'); legend(cat(2,'spiral-peec-',num2str(entries)),cat(2,'irka-reduc-6'));
subplot(3,2,6); semilogx(freq,induc_spiral_peec,'-',freq,induc5,'o'); 
title('Inductance'); legend(cat(2,'spiral-peec-',num2str(entries)),cat(2,'irka-reduc-6'));