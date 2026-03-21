%%% Resolucio equacio diferencial de segon ordre
% m x'' + c x' + k x = F(t), on F(t) es una funcio escalo
% J.A.Mayugo, UdG, 2012

clear all; close all;syms s t real

m = 1;          % massa equivalent (kg)
c = [2,4,5];    % esmorteiment viscos (Nm/s)
k = 4;          % rigidesa (Nm)
 
F=1/s  % si la funcio escalo es F(t) = 1 N; llavors F(s)=1/s)

x_=[]; v_=[];
t_=[0:0.01:10];  % simulaci? 10 segons
for i=[1:3]
    X=F/(m*s^2+c(i)*s+k); 
    x = ilaplace(X);
    pretty(X)
    pretty(x)
    
    X_dot=s*X;                  % velocitat, deriva posicio
    x_dot = ilaplace(X_dot);

    x_(:,i)=subs(x,t,t_);       % posicio x(t)
    v_(:,i)=subs(x_dot,t,t_);   % velocitat v(t)
    fr_(:,i)=x_(:,i)*k;         % força ressort f_r(t)
    fc_(:,i)=v_(:,i)*c(i);         % força esmorteïdor f_r(t)
end

h=figure;
plot(t_,x_);
legend(['c = ' num2str(c(1)) ' Ns/m'],...
       ['c = ' num2str(c(2)) ' Ns/m'],...
       ['c = ' num2str(c(3)) ' Ns/m']);
xlabel('t [s]');ylabel('x [m]')
%saveas(h,'./figures/Exemple_MatLab_1plot_x.png')

h=figure;
plot(t_,v_);
legend(['c = ' num2str(c(1)) ' Ns/m'],...
       ['c = ' num2str(c(2)) ' Ns/m'],...
       ['c = ' num2str(c(3)) ' Ns/m']);
xlabel('t [s]');ylabel('v [m/s]');
%saveas(h,'./figures/Exemple_MatLab_1plot_v.png')

h=figure;
plot(t_,fr_);
legend(['c = ' num2str(c(1)) ' Ns/m'],...
       ['c = ' num2str(c(2)) ' Ns/m'],...
       ['c = ' num2str(c(3)) ' Ns/m']);
xlabel('t [s]');ylabel('força ressort [N]')
%saveas(h,'./figures/Exemple_MatLab_1plot_fr.png')

h=figure;
plot(t_,fc_);
legend(['c = ' num2str(c(1)) ' Ns/m'],...
       ['c = ' num2str(c(2)) ' Ns/m'],...
       ['c = ' num2str(c(3)) ' Ns/m']);
xlabel('t [s]');ylabel('força esmorteïdor [N]')
%saveas(h,'./figures/Exemple_MatLab_1plot_fc.png')

h=figure;
plot(t_,fr_ + fc_);
legend(['c = ' num2str(c(1)) ' Ns/m'],...
       ['c = ' num2str(c(2)) ' Ns/m'],...
       ['c = ' num2str(c(3)) ' Ns/m']);
xlabel('t [s]');ylabel('força a la base [N]')
%saveas(h,'./figures/Exemple_MatLab_1plot_fc.png')