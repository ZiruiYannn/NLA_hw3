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
freq = 10.^[-2:0.1:2];
%%
sys150 = balred(sys,150,FreqIntervals=[10^-2,10^2]);

%%
sys120 = balred(sys,120,FreqIntervals=[10^-2,10^2]);

%%
sys100 = balred(sys,100,FreqIntervals=[10^-2,10^2]);

%%
sys50 = balred(sys,50,FreqIntervals=[10^-2,10^2]);

%%
resp_iss = bode_from_system(sys.A,eye(size(sys.A,1)),sys.B(:,1),sys.C(1,:)',i*2*pi().*freq);
%%
resp_iss50 = bode_from_system(sys50.A,eye(size(sys50.A,1)),sys50.B(:,1),sys50.C(1,:)',i*2*pi().*freq);
%%
resp_iss100 = bode_from_system(sys100.A,eye(size(sys100.A,1)),sys100.B(:,1),sys100.C(1,:)',i*2*pi().*freq);
%%
resp_iss120 = bode_from_system(sys120.A,eye(size(sys120.A,1)),sys120.B(:,1),sys120.C(1,:)',i*2*pi().*freq);
%%
resp_iss150 = bode_from_system(sys150.A,eye(size(sys150.A,1)),sys150.B(:,1),sys150.C(1,:)',i*2*pi().*freq);
%%
figure(1)
loglog(freq, abs(resp_iss50), '*r');
hold on
loglog(freq, abs(resp_iss100), '*g'); % Amplitude plot
loglog(freq, abs(resp_iss150), '*k');
loglog(freq, abs(resp_iss), '-b');
title('Amplitude Response');
xlabel('frequency')
ylabel('|H(s)|')
legend('order50','order100','order150','original')
grid on;

figure(2)
semilogx(freq, angle(resp_iss50), '*r'); 
hold on
semilogx(freq, angle(resp_iss100), '*g'); % Phase plot
semilogx(freq, angle(resp_iss150), '*k');
semilogx(freq, angle(resp_iss), '-b'); % Phase plot
title('Phase Response');
xlabel('frequency')
ylabel('arg(H(s))')
legend('order50','order100','order150','original')
grid on;
%%
sys_full = ss(A,B,C,0);
sys150_full = balred(sys_full,150,FreqIntervals=[10^-2,10^2]);
resp_iss150_full = bode_from_system(sys150_full.A,eye(size(sys150_full.A,1)),sys150_full.B(:,1),sys150_full.C(1,:)',i*2*pi().*freq);
%%
sys260_full = balred(sys_full,260,FreqIntervals=[10^-2,10^2]);
resp_iss260_full = bode_from_system(sys260_full.A,eye(size(sys260_full.A,1)),sys260_full.B(:,1),sys260_full.C(1,:)',i*2*pi().*freq);
%%
figure(11)
loglog(freq, abs(resp_iss150), '*k');
hold on
loglog(freq, abs(resp_iss150_full), '*r');
loglog(freq, abs(resp_iss260_full), '*g');
loglog(freq, abs(resp_iss), '-b');
title('Amplitude Response');
xlabel('frequency')
ylabel('|H(s)|')
legend('input1-output1 with order 150','all inputs and outputs with order 150','all inputs and outputs with order 260','original')
grid on;

figure(12)
semilogx(freq, angle(resp_iss150), '*k'); 
hold on
semilogx(freq, angle(resp_iss150_full), '*r');
semilogx(freq, angle(resp_iss260_full), '*g');
semilogx(freq, angle(resp_iss), '-b'); % Phase plot
title('Phase Response');
xlabel('frequency')
ylabel('arg(H(s))')
legend('input1-output1 with order 150','all inputs and outputs with order 150','all inputs and outputs with order 260','original')
grid on;

%%
save iss_sys.mat