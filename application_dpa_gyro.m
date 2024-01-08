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
%dpa
freq = 2.^[1:1:21];

lamda_gyros = zeros(size(freq));
x_gyros = zeros(rows,size(freq,2));
y_gyros = zeros(rows,size(freq,2));

i =sqrt(-1);

%%
%tol = 1e-2;
%sigma0 = -linspace(1e3,1e7,16);

%parfor j = 1:size(sigma0,2)
%    [lamda_gyros(j), x_gyros(:,j), y_gyros(:,j)] = qdpa(M, D, K, b, c, sigma0(j), tol);
%    sigma0(j)
%end

%%
%figure(30)
%plot3(freq, real(lambda_gyros), imag(lambda_gyros) ,'*b'); 
%set(gca,'XScale','log')
%xlabel('The imag part of initial points')
%ylabel('The real part of lambda')
%zlabel('The imag part of lambda')
%grid on

%%
%sadpa
i=sqrt(-1);
s0 = [-logspace(3,7,16) (-1.2028 + 0.9947i) (-0.0148 + 0.1715i) (-0.0041 + 0.0902i) (-0.0025 + 0.0708i) (-0.0010 - 0.0441i)...
    (-0.0007 - 0.0382i) (-0.0004 + 0.0277i) (-0.0001 - 0.0106i) (-0.0001 + 0.0158i) ]*1e6;
nwanted = 10;
options = struct("nwanted",nwanted,"tol",1e-3, "displ",1,"strategy",'LR',"kmin",1,"kmax",15,"maxrestarts",15,...
    "f_ax",'N',"f_ex",'N',"f_semax",'N',"f_semax_s",'N',"use_lu",0,"use_lu_w_amd",0,"dpa_bordered",0,"yEx_scaling",0, ...
    "rqitol",1e-2,"turbo_deflation",1);
[poles, residues, rightev, leftev, nr_solves] = saqdpa(K, D, M,b, c, 0, s0, options);



%%
%save sadpa_gyro10.mat poles residues rightev leftev nr_solves
%%
%load sadpa_gyro10.mat
load sys_gyro10.mat

%%
i=sqrt(-1);
s = i*2*pi()*10.^[0:0.2:4];

resp_gyro=bode_from_system([zeros(10) eye(10);-sys10.K zeros(10)], [sys10.M zeros(10);sys10.C eye(10)],[zeros(10,1);sys10.B],[sys10.F';zeros(10,1)],s);
%%
%systemp=dss([zeros(10) eye(10);-sys10.K zeros(10)],[zeros(10,1);sys10.B],[sys10.F zeros(1,10)],0,[sys10.M zeros(10);sys10.C eye(10)]);
%%
H_diff_inf = zeros(1,10);
n = length(s);
resp_saqdpa_gyro_temp = zeros(n,10);
for m=1:10      
    for j=1:n
        temp=0;
        for k=1:m
                temp=temp+ residues(k)/(s(j)-poles(k));
        end
        resp_saqdpa_gyro_temp(j,m) = temp;
    end
    H_diff_inf(m) = max(abs(resp_saqdpa_gyro_temp(:,m)-resp_gyro));
end

H_diff_inf = H_diff_inf/max(abs(resp_gyro));
%%

%figure(33)
%plot(H_diff_inf)
%title('$\frac{|H(s)-\hat{H(s)}|_{\infty}}{|H(s)|_{\infty}}$ for reduced models of using different number of dominant poles','interpreter','latex');
%xlabel('number of poles')
%ylabel('$\frac{|H(s)-\hat{H(s)}|_{\infty}}{|H(s)|_{\infty}}$','interpreter','latex')

%%
freq = 10.^[0:0.2:4];
figure(34)
semilogx(freq,abs(resp_saqdpa_gyro_temp)')
title('Amplitude Response');
xlabel('frequency')
ylabel('|H(s)|')
figure(35)
semilogx(freq,abs(resp_gyro))
title('Amplitude Response');
xlabel('frequency')
ylabel('|H(s)|')