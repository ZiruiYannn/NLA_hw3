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
%dpa
freq = 2.^[1:1:21]

lamda_peecs = zeros(size(freq));
x_peecs = zeros(rows,size(freq,2));
y_peecs = zeros(rows,size(freq,2));

i =sqrt(-1);

%%
tol = 1e-5;

for j = 1:size(freq,2)
    [lambda_peecs(j), x_peecs(:,j), y_peecs(:,j)] = dpa(peec.E, peec.A, peec.B, peec.B, i*2*pi()*freq(j),tol);
    freq(j)
end

%%
figure(30)
plot3(freq, real(lambda_peecs), imag(lambda_peecs) ,'*b'); 
set(gca,'XScale','log')
xlabel('The imag part of initial points')
ylabel('The real part of lambda')
zlabel('The imag part of lambda')
grid on

%%
%sadpa
s0 = -logspace(7,12,64);
nwanted = 50;
options = struct("nwanted",nwanted,"tol",1e-5, "displ",1,"strategy",'LR',"kmin",1,"kmax",15,"maxrestarts",100,...
    "f_ax",'N',"f_ex",'N',"f_semax",'N',"f_semax_s",'N',"use_lu",0,"use_lu_w_amd",0,"dpa_bordered",0,"yEx_scaling",0, ...
    "rqitol",1e-4,"turbo_deflation",0);
[poles, residues, rightev, leftev, nr_solves, ress] = sadpa(peec.A, peec.E, peec.B, peec.B, 0, s0, options);

%%
%save sadpa_peec50.mat poles residues rightev leftev nr_solves
%%
m=10;

freq = 10.^[1:0.1:10];
s = i*2*pi()*freq;

n = length(s);
resp_sadpa_peec10 = zeros(n,1);
for j=1:n
    temp=0;
    for k=1:m
            temp=temp+ residues(k)/(s(j)-poles(k));
    end
    resp_sadpa_peec10(j,:) = temp;
end
%%
m=30;

freq = 10.^[1:0.1:10];
s = i*2*pi()*freq;

n = length(s);
resp_sadpa_peec30 = zeros(n,1);
for j=1:n
    temp=0;
    for k=1:m
            temp=temp+ residues(k)/(s(j)-poles(k));
    end
    resp_sadpa_peec30(j,:) = temp;
end

%%
m=50;

freq = 10.^[1:0.1:10];
s = i*2*pi()*freq;

n = length(s);
resp_sadpa_peec50 = zeros(n,1);
for j=1:n
    temp=0;
    for k=1:m
            temp=temp+ residues(k)/(s(j)-poles(k));
    end
    resp_sadpa_peec50(j,:) = temp;
end
%%
freq = 10.^[1:0.1:10];
s = i*2*pi()*freq;

resp_spiral_inductor = bode_from_system(peec.A,peec.E, ...
        peec.B,peec.B,s);

%%
figure(31)

loglog(freq, abs(resp_sadpa_peec10), '-g'); % Amplitude plot
hold on
loglog(freq, abs(resp_sadpa_peec30), '-r');
loglog(freq, abs(resp_sadpa_peec50), '-k');
loglog(freq, abs(resp_spiral_inductor), '-b');
title('Amplitude Response');
xlabel('frequency')
ylabel('|H(s)|')
legend('order10','order30','order50','original')
grid on;

figure(32)
semilogx(freq, angle(resp_sadpa_peec10), '-g'); % Phase plot
hold on
semilogx(freq, angle(resp_sadpa_peec30), '-r');
semilogx(freq, angle(resp_sadpa_peec50), '-k');
semilogx(freq, angle(resp_spiral_inductor), '-b');
title('Phase Response');
xlabel('frequency')
ylabel('arg(H(s))')
legend('order10','order30','order50','original')
grid on;


%%
H_diff_inf = zeros(1,50);
for m=1:50    
    n = length(s);
    resp_sadpa_peec_temp = zeros(n,1);
    for j=1:n
        temp=0;
        for k=1:m
                temp=temp+ residues(k)/(s(j)-poles(k));
        end
        resp_sadpa_peec_temp(j,:) = temp;
    end
    H_diff_inf(m) = max(abs(resp_sadpa_peec_temp-resp_spiral_inductor));
end

H_diff_inf = H_diff_inf/max(abs(resp_spiral_inductor));
%%

figure(33)
plot(H_diff_inf)
title('$\frac{|H(s)-\hat{H(s)}|_{\infty}}{|H(s)|_{\infty}}$ for reduced models of using different number of dominant poles','interpreter','latex');
xlabel('number of poles')
ylabel('$\frac{|H(s)-\hat{H(s)}|_{\infty}}{|H(s)|_{\infty}}$','interpreter','latex')

%%


%%




