%%% Control de vibracions
% Amb un sistema MKC absorbir vibracions d'un sistema MK
% J.A.Mayugo, UdG, 2012

clear all; close all;

%% Sistema inicial MK
m1 = 150;   % massa [kg]
k1 = 2.65e6;% rigidesa [N/m] 

A0 = [ 0     1;
    -k1/m1  0];
B0 = [ 0 ; 1/m1]; % input una força

C0 = [ 1  0 ];   % defineix output posicio x
D0 = 0;

MK_ss = ss(A0, B0, C0, D0);

%% Sistema MK + MKC de control
m2 = 25;    % massa del sitema MKC de control [kg]

% càcul control òptim
mu = m2/m1
q_opt = 1 /(1+mu)
zeta_opt=sqrt(3*mu/(8*(1+mu)))
  
k2= q_opt^2*k1*mu
c = zeta_opt*2*sqrt(k2*m2)

% espai d'estat
A = [   0       1       0       0;
   -(k1+k2)/m1 -c/m1    k2/m1   c/m1;
        0       0       0       1;
        k2/m2   c/m2   -k2/m2  -c/m2];
B = [ 0 ; 1/m1 ; 0 ; 0 ];

C = [ 1 0 0 0 ];
D = 0;

MK_MKS_ss = ss(A, B, C, D);

H0=figure;
%% Resposta freqüencial (H(W))
h = bodeplot(MK_ss*k1,MK_MKS_ss*k1,[50:0.01:200]);
setoptions(h,'MagUnits','abs','MagScale','linear','FreqScale','linear')
setoptions(h,'YLim',{[0 5],[-180 0]});
saveas(H0,'Ex_bode_abs.eps')

% %% NATURAL FREQUENCY
omega_n0 = abs(pole(MK_ss))
omega_n  = abs(pole(MK_MKS_ss))
% %% RESONANT FREQUENCY
[omega_r_c] = abs(imag(pole(MK_MKS_ss)))
% %% DAMPING RATIO
[Wn,zeta_c] = damp(MK_MKS_ss)

%% Respostes transitòries
t=0:0.001:.4;
F_0=1000;       
omega = omega_n0(1)*1.10  % exitació a un 10% de la W_natural
f = F_0*sin(omega.*t);

[Y0,t1]=lsim(A0,B0,C0,D0,f,t);
[Y ,t1]=lsim(A, B, C, D,f,t);

H1=figure;hold on
plot(t,Y0*1e3,'--r')  % resposta amb MK inicial
plot(t,Y*1e3 ,'-.m')  % resposta amb MKC de control
plot(t,f/k1*1e3)    % si no hi haguès efectes dinàmics
legend('resposta inicial','resposta controlada','f(t)/k1')
xlabel('temps (s)');ylabel('desplaçament (mm)')
hold off
saveas(H1,'Ex_tran_abs.eps')