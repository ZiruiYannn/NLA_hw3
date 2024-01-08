clc
clear

%%

addpath("given_functions")
addpath("Systems/PEECmodel")

%%

% m is PRIMA reduction order. 
m = 50;

run("plot_spiral.m")

%%
%bode_from_system(A,E,b,c,s)

resp_spiral_inductor = zeros(1,size(freq,2));

i = sqrt(-1);

resp_spiral_inductor = bode_from_system(spiral_inductor_peec.A,spiral_inductor_peec.E, ...
        spiral_inductor_peec.B,spiral_inductor_peec.B,i*2*pi().*freq);


figure(10)
subplot(1, 2, 1);
loglog(freq, abs(resp_spiral_inductor), '-b'); % Amplitude plot
title('Amplitude Response');
xlabel('frequency')
ylabel('|H(s)|')
grid on;

subplot(1, 2, 2);
semilogx(freq, angle(resp_spiral_inductor), '-b'); % Phase plot
title('Phase Response');
xlabel('frequency')
ylabel('arg(H(s))')
grid on;

%%

i = sqrt(-1);

resp_spiral_inductor = zeros(1,size(freq,2));
resp_spiral_inductor = bode_from_system(spiral_inductor_peec.A,spiral_inductor_peec.E, ...
        spiral_inductor_peec.B,spiral_inductor_peec.B,i*2*pi().*freq);
%%
m = 50;
run("plot_spiral.m")
resp_spiral_inductor50 = zeros(1,size(freq,2));
resp_spiral_inductor50 = bode_from_system(-R_prima,L_prima, ...
        B_prima,B_prima,i*2*pi().*freq);
%%
m = 10;
run("plot_spiral.m")
resp_spiral_inductor10 = zeros(1,size(freq,2));
resp_spiral_inductor10 = bode_from_system(-R_prima,L_prima, ...
        B_prima,B_prima,i*2*pi().*freq);
%%
m = 3;
run("plot_spiral.m")
resp_spiral_inductor3 = zeros(1,size(freq,2));
resp_spiral_inductor3 = bode_from_system(-R_prima,L_prima, ...
        B_prima,B_prima,i*2*pi().*freq);


%%
figure(12)

loglog(freq, abs(resp_spiral_inductor50), '-g'); % Amplitude plot
hold on
loglog(freq, abs(resp_spiral_inductor10), '-r');
loglog(freq, abs(resp_spiral_inductor3), '-y');
loglog(freq, abs(resp_spiral_inductor), '-b');
title('Amplitude Response');
xlabel('frequency')
ylabel('|H(s)|')
legend('order50','order10','order3','original')
grid on;

figure(13)
semilogx(freq, angle(resp_spiral_inductor50), '-g'); % Phase plot
hold on
semilogx(freq, angle(resp_spiral_inductor10), '-r'); % Phase plot
semilogx(freq, angle(resp_spiral_inductor3), '-y'); % Phase plot
semilogx(freq, angle(resp_spiral_inductor), '-b');
title('Phase Response');
xlabel('frequency')
ylabel('arg(H(s))')
legend('order50','order10','order3','original')
grid on;

%%
sys_peec3 = dss(-R_prima,B_prima,B_prima',0,L_prima)
%%
pole(sys_peec3)
%-3.1469 e8
%-0.1403 e8
%-1.0908 e8

%   1.0e+09 *

%   -5.8959
 %  -0.5734
  % -0.3959
  % -0.2803
  % -0.0140
   %-0.1749
   %-0.1369
   %-0.0846
   %-0.1020
   %-0.1112