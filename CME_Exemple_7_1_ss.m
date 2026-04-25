%%% Obtenir l'eqüació 1 gdl en espai d'estat
% Obtenir el diagrama de Bode
% J.A.Mayugo, UdG, 2012

clear all; close all;

m = 0.02; % massa [kg]
c = 1.25; % esmorteïment [Ns/m]
k = 5500; % rigidesa [N/m] 

A = [ 0     1;
    -k/m -c/m];
B = [ 0 ; 1/m ]; % input una força

C = [ k  0 ];   % defineix output posicio x
D = 0;

model = ss(A, B, C, D)

%% Resposta freqüencial (H(W))

omega_n   = abs(pole(model))   % Natural frequency
[zeta]    = abs(real(pole(model)))./abs(pole(model))  % Damping Ratio
[omega_r] = abs(imag(pole(model)))   % Resonant frequency

% FUNCIÓ AMB L'ORDRE BODE de MatLab (solució numèrica)
H1=figure;
h = bodeplot(model,{400,800});
setoptions(h,'MagUnits','abs','FreqScale','linear')
saveas(H1,'Ex_bode.eps')

%% TRANSFER FUNCTION to SS conversion 
%[num, den] = ss2tf(A, B, C, D)
%tutorial_tf = tf(model)
%[A,B,C,D] = tf2ss(num,den)
%tutorial_H=tutorial_tf*k
%h2 = bodeplot(tutorial_H,{400,1000});
%setoptions(h2,'MagUnits','abs')