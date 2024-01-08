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
tol = 1e-11;

i=sqrt(-1);
sigma1 = -1;
smin=-1e2;
smax=1e-1;
scount=1000;

freq = 10.^[-2:0.1:2];
s = i*2*pi()*freq;

[Ahat1, Ehat1, bhat1, chat1, i1] = grka(A, eye(size(A,1)), b, c', sigma1, smin, smax, scount, tol);

%%


i=sqrt(-1);
sigma1 = -1e-6;
smin=-1e10;
smax=0;
scount=1000;


[Ahat2, Ehat2, bhat2, chat2, i2] = grka(A, eye(size(A,1)), b, c', sigma1, smin, smax, scount, tol);

%%


i=sqrt(-1);
sigma1 = -1e2;
smin=-1e3;
smax=0;
scount=1000;


[Ahat3, Ehat3, bhat3, chat3, i3] = grka(A, eye(size(A,1)), b, c', sigma1, smin, smax, scount, tol);

%%

sys1 = dss(Ahat1,bhat1,chat1',0,Ehat1);
sys2 = dss(Ahat2 ,bhat2 ,chat2' ,0,Ehat2);
sys3 = dss(Ahat3 ,bhat3 ,chat3' ,0,Ehat3);

%%
bode(sys,sys1,sys2,sys3,2*pi()*freq)
legend('original','sigma1 = -1 i1=67','sigma1 = -1e-6 i2=7','sigma1 = -1e2 i3=45')

    
