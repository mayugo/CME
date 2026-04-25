%%% Obtenir l'eqüació 1 gdl en funció de transferència
% Obtenir el diagrama de Bode
% J.A.Mayugo, UdG, 2012

clear all; close all;

m = 0.02; % massa [kg]
c = 1.25; % esmorteïment [Ns/m]
k = 5500; % rigidesa [N/m] 

model = tf([1],[m c k])

%% Respostes transitòries
%u = 4;x0 = [ -0.08; 0.002];
%initial(model, x0)
%step(u * model)
%impulse(u * model)

%% Resposta freqüencial (H(W))

% %% NATURAL FREQUENCY
% The natural frequency in rad/s is equal to the module of complex poles.
omega_n   = sqrt(k/m)   % know solution for 1 dof system
omega_n_c = abs(pole(model)) % computed solution

% %% DAMPING RATIO
% The damping ratio for the complex eigenvalue can also be manually calcu-
% lated as the ratio of the real part of the pole to the natural frequency.
zeta     = c/(2*sqrt(k*m))  % know solution for 1 dof system
[zeta_c] = abs(real(pole(model)))./abs(pole(model))  % computed solution

% %% RESONANT FREQUENCY  
% The resonant frequency in rad/sec is equal to the size of the imaginary 
% part of the complex poles which can be found by examination, or by 
% taking the absolute value of the imaginary part of the pole.
omega_r     = omega_n*sqrt(1-zeta^2)  % know solution for 1 dof system
[omega_r_c] = abs(imag(pole(model)))  % computed solution

% The natural frequency and the damping ratio can be found using MATLAB's 
% damp function.
[Wn,Z] = damp(model)

% A) FUNCIÓ AMB l'equació (solució analítica)
omg =400:0.1:800;
Omega=omg/omega_n;
H_omg=1./sqrt((1-Omega.^2).^2+(2*zeta.*Omega).^2);
Theta_omg=atan2((2*zeta.*Omega),(1-Omega.^2))

H0=figure;
subplot(2,1,1); plot(omg,H_omg);ylabel('H');
axis([400 800 0 10])  
subplot(2,1,2); plot(omg,Theta_omg*180/pi);
xlabel('Freq (rad/s)');ylabel('Desfase (graus)');
axis([400 800 0 180])
saveas(H0,'Ex_H.eps')

% B) FUNCIÓ AMB L'ORDRE BODE de MatLab (solució numèrica)
H1=figure;
h = bodeplot(model*k,{400,800});
setoptions(h,'MagUnits','abs','FreqScale','linear')
saveas(H1,'Ex_bode.eps')