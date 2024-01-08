clc
clear

%%
addpath("given_functions")
addpath("Systems/PEECmodel")

%%
load iss12a.mat

%%
b = B(:,1);
c = C(1,:);

%%
sys = ss(A,b,c,0);


%%
m=160;
tol = 1e-2;

i = sqrt(-1);
sigma0 = [i*linspace(1,50,m/2) -i*linspace(1,50,m/2)];

freq = 10.^[-2:0.1:2];
s = i*2*pi()*freq;

[Ahat1, Ehat1, bhat1, chat1] = irka(A, eye(size(A,1)), b, c', sigma0, tol);
%%
resp_irka1 = bode_from_system(Ahat1,Ehat1, bhat1, chat1,s);

%%
sys1 = dss(Ahat1,bhat1,chat1',0,Ehat1);

%%
m=40;
tol = 1e-2;

i = sqrt(-1);
sigma0 = linspace(1,m*3,m);

freq = 10.^[-2:0.1:2];
s = i*2*pi()*freq;

[Ahat2, Ehat2, bhat2, chat2] = irka(A, eye(size(A,1)), b, c', sigma0, tol);
%%
resp_irka2 = bode_from_system(Ahat2 ,Ehat2 , bhat2 , chat2,s);

%%
sys2 = dss(Ahat2 ,bhat2 ,chat2' ,0,Ehat2);

%%
m=80;
tol = 1e-2;

i = sqrt(-1);
sigma0 = linspace(1,m*3,m);

freq = 10.^[-2:0.1:2];
s = i*2*pi()*freq;

[Ahat3, Ehat3, bhat3, chat3] = irka(A, eye(size(A,1)), b, c', sigma0, tol);
%%
resp_irka3 = bode_from_system(Ahat3 ,Ehat3 , bhat3 , chat3,s);

%%
sys3 = dss(Ahat3 ,bhat3 ,chat3' ,0,Ehat3);

%%
figure(40)
bode(sys,sys1,sys3,sys2,2*pi()*freq)
legend('original','order160','order80','order40')

