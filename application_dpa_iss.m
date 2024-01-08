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
resp_iss = bode_from_system(sys.A,eye(size(sys.A,1)),sys.B(:,1),sys.C(1,:)',i*2*pi().*freq);

%%
%sadpa
s0 = zeros(1,300);
i = sqrt(-1);
for j=1:10
    for k=1:30
        s0(30*(j-1)+k) = -0.02*j + i *(j+k);
    end
end
nwanted = 200;
options = struct("nwanted",nwanted,"tol",1e-5, "displ",1,"strategy",'LR',"kmin",1,"kmax",15,"maxrestarts",100,...
    "f_ax",'N',"f_ex",'N',"f_semax",'N',"f_semax_s",'N',"use_lu",0,"use_lu_w_amd",0,"dpa_bordered",0,"yEx_scaling",0, ...
    "rqitol",1e-4,"turbo_deflation",1);
[poles, residues, rightev, leftev, nr_solves, ress] = sadpa(A, eye(size(A,1)), b, c', 0, s0, options);

%%
s= i*2*pi().*freq;
n = length(s);
H_diff_inf = zeros(1,200);
for m=1:200    
    resp_sadpa_iss_temp = zeros(n,1);
    for j=1:n
        temp=0;
        for k=1:m
                temp=temp+ residues(k)/(s(j)-poles(k));
        end
        resp_sadpa_iss_temp(j,:) = temp;
    end
    H_diff_inf(m) = max(abs(resp_sadpa_iss_temp-resp_iss));
end

H_diff_inf = H_diff_inf/max(abs(resp_iss));

%%

figure(33)
plot(H_diff_inf)
title('$\frac{|H(s)-\hat{H(s)}|_{\infty}}{|H(s)|_{\infty}}$ for reduced models of using different number of dominant poles','interpreter','latex');
xlabel('number of poles')
ylabel('$\frac{|H(s)-\hat{H(s)}|_{\infty}}{|H(s)|_{\infty}}$','interpreter','latex')

%%
resp_sadpa_iss120 = zeros(n,1);
for j=1:n
    temp=0;
    for k=1:120
        temp=temp+ residues(k)/(s(j)-poles(k));
    end
    resp_sadpa_iss120(j,:) = temp;
end

resp_sadpa_iss150 = zeros(n,1);
for j=1:n
    temp=0;
    for k=1:150
        temp=temp+ residues(k)/(s(j)-poles(k));
    end
    resp_sadpa_iss150(j,:) = temp;
end

figure(31)

loglog(freq, abs(resp_sadpa_iss120), '-g'); % Amplitude plot
hold on
loglog(freq, abs(resp_sadpa_iss150), '-r');
loglog(freq, abs(resp_iss), '-b');
title('Amplitude Response');
xlabel('frequency')
ylabel('|H(s)|')
legend('order120','order150','original')
grid on;

figure(32)
semilogx(freq, angle(resp_sadpa_iss120), '-g'); % Phase plot
hold on
semilogx(freq, angle(resp_sadpa_iss150), '-r');
semilogx(freq, angle(resp_iss), '-b');
title('Phase Response');
xlabel('frequency')
ylabel('arg(H(s))')
legend('order120','order150','original')
grid on;