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

%%
tol = 1e-2;

i=sqrt(-1);
sigma1 = -1;
smin=-1e4;
smax=0;
scount=200;

freq = 10.^[1:0.1:4];
s = i*2*pi()*freq;

[Mhat1, Dhat1, Khat1, bhat1, chat1, i1] = qgrka(M,D,K,b,c, sigma1, smin, smax, scount, tol);

%%
tol = 1e-2;

i=sqrt(-1);
sigma1 = -1e4;
smin=-1e6;
smax=0;
scount=200;


[Mhat2, Dhat2, Khat2, bhat2, chat2, i2] = qgrka(M,D,K,b,c, sigma1, smin, smax, scount, tol);


%%
tol = 1e-2;

i=sqrt(-1);
sigma1 = -1e7;
smin=-1e10;
smax=0;
scount=200;


[Mhat3, Dhat3, Khat3, bhat3, chat3, i3] = qgrka(M,D,K,b,c, sigma1, smin, smax, scount, tol);

%%
sys = mechss(M,D,K,b,c');
sys1 = mechss(Mhat1,Dhat1,Khat1,bhat1,chat1');
sys2 = mechss(Mhat2,Dhat2,Khat2,bhat2,chat2');

%%
bode(sys,sys1,sys2,2*pi()*freq)

legend('original','sigma1 = -1 i1=3','sigma1 = -1e4 i2=23')

