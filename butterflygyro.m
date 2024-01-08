clc
clear

%%

addpath("given_functions")
addpath("Systems/PEECmodel")
addpath("Systems/ButterflyGyro-dim1e5-gyro")

%%

beta = 1e-6;

[gyro.B, rows, cols_B, entries] = mmread('gyro.B');
[gyro.C, cols_C, rows, entries] = mmread('gyro.C');
[gyro.K, rows, cols, entries] = mmread('gyro.K');
[gyro.M, rows, cols, entries] = mmread('gyro.M');

%%
num_out=2;
beta = 1e-6;

B = gyro.B;
C = gyro.C(num_out,:);
K = gyro.K;
M = gyro.M;


%%

sys = mechss(M,beta*K,K,B,C);

%%
figure(1)
bode(sys,10.^[0:0.1:4])

%%
block_size = size(B,2);
fprintf(1,'doing invK_MB\n');
invK_M = K\M;
invK_B = K\B;
%%
fprintf(1,'calling barnoldi\n');
[X, H, R1] = barnoldi(invK_M, invK_B, 100, 1);

fprintf(1,'computing reduced system\n');

%%

%change order m
m = 2;


Xm= X(:,1:m+1);
Cr = Xm'*M*Xm;
Gr = Xm'*K*Xm;

M_prima2 = Cr(1:m,1:m);
K_prima2 = Gr(1:m,1:m);
B_prima2 = transpose(Xm(:,1:m))*B;
C_prima2 = C*Xm(:,1:m);
%%
sys2 = mechss(M_prima2,beta*K_prima2,K_prima2,B_prima2,C_prima2);


%%

%change order m
m = 3;


Xm= X(:,1:m+1);
Cr = Xm'*M*Xm;
Gr = Xm'*K*Xm;

M_prima3 = Cr(1:m,1:m);
K_prima3 = Gr(1:m,1:m);
B_prima3 = transpose(Xm(:,1:m))*B;
C_prima3 = C*Xm(:,1:m);
%%
sys3 = mechss(M_prima3,beta*K_prima3,K_prima3,B_prima3,C_prima3);

%%

%change order m
m = 4;


Xm= X(:,1:m+1);
Cr = Xm'*M*Xm;
Gr = Xm'*K*Xm;

M_prima4 = Cr(1:m,1:m);
K_prima4 = Gr(1:m,1:m);
B_prima4 = transpose(Xm(:,1:m))*B;
C_prima4 = C*Xm(:,1:m);
%%
sys4 = mechss(M_prima4,beta*K_prima4,K_prima4,B_prima4,C_prima4);
%%
figure(1)
bode(sys,sys2,sys3,sys4,2*pi()*10.^[0:0.1:4])

legend('original','order2','order3','order4');
%%
m = 10;


Xm= X(:,1:m+1);
Cr = Xm'*M*Xm;
Gr = Xm'*K*Xm;

M_prima10 = Cr(1:m,1:m);
K_prima10 = Gr(1:m,1:m);
B_prima10 = transpose(Xm(:,1:m))*B;
C_prima10 = C*Xm(:,1:m);
%%

sys10 = mechss(M_prima10,beta*K_prima10,K_prima10,B_prima10,C_prima10);
%%
%save sys_gyro10 sys10
%%
figure(2)
bode(sys,sys10,2*pi()*10.^[0:0.1:4])
legend('original','order10');