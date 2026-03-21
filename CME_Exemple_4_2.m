%% RESPOSTA TRANSITORIA EN CONDICIONS GENERALS NO NUL·LES
% J.A.Mayugo, UdG, 2015
clear, close all

%% Dades exercici
m = 108     % kg
c = 78      % Ns/m
k = 2925    % N/m

m_eq=(m+1/12*m)/9
c_eq=c
k_eq=k/9

x0 = [1e-3 ];           % CI desplaçament inicial

t_ = [0:0.01:1.6];      % temps, s

%% MÈTODE 1 - Descomposició en fraccions parcials
syms s t real
eqn = s^2*m_eq+c_eq*s+k_eq == 0;
sols = solve(eqn,s)

[residus,pols,kk] = residue([x0*m_eq c_eq*x0],[m_eq c_eq k_eq]);
FF1=residus(1)/(s-pols(1))
FF2=residus(2)/(s-pols(2))

ff=ilaplace(FF1+FF2)
pretty(ff)

ff_=eval(subs(ff,t,t_));
figure;
plot(t_,ff_);
xlabel('t [s]');ylabel('x [m]')
title('Resposta amb Condicions Inicials no nul·les')

%% MÈTODE 2 - Emprant funció Inversa de TF de Matlab (ilaplace)
syms s t real
F = (x0*m_eq*s+c_eq*x0)/(s^2*m_eq+c_eq*s+k_eq);
f=ilaplace(F, t);
pretty(f)

limit(s*F,s,0)      % valor final de x(t)
limit(s*F,s,Inf)    % valor inicial de x(t)

f_=eval(subs(f,t,t_));
figure;plot(t_,f_);
xlabel('t [s]');ylabel('x [m]')
title('Resposta amb Condicions Inicials no nul·les')

figure;plot(t_,f_/3*k);
xlabel('t [s]');ylabel('F [N]')
title('Força molla')

%% MÊTODE 3 - Funció TF de Matlab en espai s

model = tf([x0*m_eq x0*c_eq],[m_eq c_eq k_eq]);

figure;
impulse(model);
title('Resposta amb Condicions Inicials no nul·les');xlabel('t [s]');ylabel('x [m]')