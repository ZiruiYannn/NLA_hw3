clc
clear

%%

addpath("given_functions")
addpath("Systems/PEECmodel")
addpath("Systems/ButterflyGyro-dim1e5-gyro")

%%
%gyro
beta = 1e-6;

[gyro.B, rows, cols_B, entries] = mmread('gyro.B');
[gyro.C, cols_C, rows, entries] = mmread('gyro.C');
[gyro.K, rows, cols, entries] = mmread('gyro.K');
[gyro.M, rows, cols, entries] = mmread('gyro.M');

%%
num_out = 2;
M = gyro.M;
K = gyro.K;
D = beta*K;
c = gyro.C(2,:)';
b = gyro.B;

E = [M zeros(size(M)); D eye(size(M,1))];
A = [zeros(size(M)) eye(size(M,1)); -K zeros(size(M))];


%%
m = 4;
tol = 1e-2;

sigma0 = -linspace(1e3,1e7,m);
i = sqrt(-1);
freq = 10.^[1:0.1:4];
s = i*2*pi()*freq;

[Mhat4, Dhat4, Khat4,bhat4, chat4, V4, Ehat4, Ahat4] = qirka(M,D,K,b, c, sigma0, tol);
%%
sys4_irka = mechss(Mhat4,Dhat4,Khat4,bhat4,chat4');
%resp_irka_gyro4 = bode_from_system(Ehat4,Ahat4, [zeros(4,1);bhat4], [chat4;zeros(4,1)],s);

%%
sys = mechss(M,D,K,b,c');
%resp_gyro = bode_from_system(E,A,[zeros(size(M,1),1);b], [c;zeros(size(M,1),1)],s)

%%
m = 6;
tol = 1e-2;

sigma0 = -linspace(1e3,1e7,m);
i = sqrt(-1);
freq = 10.^[1:0.1:4];
s = i*2*pi()*freq;

[Mhat6, Dhat6, Khat6,bhat6, chat6, V6, Ehat6, Ahat6] = qirka(M,D,K,b, c, sigma0, tol);

%%
sys6_irka = mechss(Mhat6,Dhat6,Khat6,bhat6,chat6');
resp_irka_gyro6 = bode_from_system(Ehat6,Ahat6, [zeros(6,1);bhat6], [chat6;zeros(6,1)],s);

%%
m = 2;
tol = 1e-2;

sigma0 = -linspace(1e3,1e7,m);
i = sqrt(-1);
freq = 10.^[1:0.1:4];
s = i*2*pi()*freq;

[Mhat2, Dhat2, Khat2,bhat2, chat2, V2, Ehat2, Ahat2] = qirka(M,D,K,b, c, sigma0, tol);
%%
sys2_irka = mechss(Mhat2,Dhat2,Khat2,bhat2,chat2');
resp_irka_gyro2 = bode_from_system(Ehat2,Ahat2, [zeros(2,1);bhat2], [chat2;zeros(2,1)],s);
%%
bode(sys,sys2_irka,sys4_irka,sys6_irka,2*pi()*10.^[0:0.1:4])

legend('original','order2','order4','order6')

%%
%H_diff_inf4 = max(abs(resp_irka_gyro4-resp_gyro))/max(abs(resp_gyro)) %0.0032
%H_diff_inf6 = max(abs(resp_irka_gyro6-resp_gyro))/max(abs(resp_gyro)) %5.1420e-04
%H_diff_inf2 = max(abs(resp_irka_gyro2-resp_gyro))/max(abs(resp_gyro)) %0.06.3