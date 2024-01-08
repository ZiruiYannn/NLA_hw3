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
tol = 1e-2;

i=sqrt(-1);
sigma1 = -1e7;
smin=-1e10;
smax=-1e5;
scount=200;

freq = 10.^[1:0.1:10];
s = i*2*pi()*freq;

[Ahat1, Ehat1, bhat1, chat1, i1] = grka(A, E, B, B, sigma1, smin, smax, scount, tol);
%%
resp_grka_peec1 = bode_from_system(Ahat1,Ehat1, bhat1, chat1,s);

%%
tol = 1e-2;

i=sqrt(-1);
sigma1 = -1e-6;
smin=-1e10;
smax=0;
scount=200;


[Ahat2, Ehat2, bhat2, chat2, i2] = grka(A, E, B, B, sigma1, smin, smax, scount, tol);

%%
resp_grka_peec2 = bode_from_system(Ahat2,Ehat2, bhat2, chat2,s);

%%
tol = 1e-2;

i=sqrt(-1);
sigma1 = -1;
smin=-1e5;
smax=0;
scount=100;


[Ahat3, Ehat3, bhat3, chat3, i3] = grka(A, E, B, B, sigma1, smin, smax, scount, tol);

%%
resp_grka_peec3 = bode_from_system(Ahat3,Ehat3, bhat3, chat3,s);
%%
resp_spiral_inductor = bode_from_system(A,E,B,B,s);

%%
figure(31)

loglog(freq, abs(resp_grka_peec1), '-g'); % Amplitude plot
hold on
loglog(freq, abs(resp_grka_peec2), '-r');
loglog(freq, abs(resp_grka_peec3), '-k');
loglog(freq, abs(resp_spiral_inductor), '-b');
title('Amplitude Response');
xlabel('frequency')
ylabel('|H(s)|')
legend('sigma1=-1e7,i1=229','sigma1=-1e-6,i1=231','sigma1=-1,i1=3','original')
grid on;

figure(32)
semilogx(freq, angle(resp_grka_peec1), '-g'); % Phase plot
hold on
semilogx(freq, angle(resp_grka_peec2), '-r');
semilogx(freq, angle(resp_grka_peec3), '-k');
semilogx(freq, angle(resp_spiral_inductor), '-b');
title('Phase Response');
xlabel('frequency')
ylabel('arg(H(s))')
legend('sigma1=-1e7,i1=229','sigma1=-1e-6,i1=231','sigma1=-1,i1=3','original')
grid on;