function kran_5()
close all
global x0 xf umax umin g

scenario = 1;
%# Modellparameter
K_x = 0.99;
l_s = 1.0;
g   = 9.81;
v_soll = 0.4;
v_min = -10.0;
a_max  = 0.33;
a_min = 0.0;
phi_0 = 0;
dphi_0 = 0;
phi_tf = 0;
dphi_tf = 0;
u_max = a_max/K_x;
umax = u_max;
switch scenario
    case 0
        t1= 0.6898200664133287;
        t2= 0.8888066240325581;
        tf= 1.5918924609538354;
        umin = -u_max;
        results_X = load('kran.dat');
        N1 = (size(results_X,1)-1)/3;
        x10 = 0;%results_X(1,3);
        x20 = 0;%results_X(1,4);
        x30 = 0;%results_X(1,5);
        results_dual = load('kran_dual.dat');
        p10 = results_dual(1,2);
        p20 = results_dual(1,3);
        p30 = results_dual(1,4);        
    case 1
        t1= 0.5879036776475854;
        t2= 0.9798394627459757;
        tf= 1.6033736663115965;
        umin = 0.00000;
        results_X = load('kran_a_min_0.dat');
        N1 = (size(results_X,1)-1)/3;
        x10 = results_X(1,3);
        x20 = results_X(1,4);
        x30 = results_X(1,5);
        results_dual = load('kran_a_min_0_dual.dat');
        p10 = results_dual(1,2);
        p20 = results_dual(1,3);
        p30 = results_dual(1,4);
end
% N1 = 100;
x0 = [x10 x20 x30]';
xf = [v_soll, 0,0]';
% solinit = bvpinit([linspace(0, 1.0, N1) linspace(1.0, 2.0, N1), ...
%     linspace(2.0, 3.0, N1)],[0 0 0 0 0 0], [t1; t2;tf]);
solinit = bvpinit([linspace(0, 1.0, N1) linspace(1.0, 2.0, N1), ...
     linspace(2.0, 3.0, N1)],[x10;x20;x30; p10; p20; p30], [t1; t2;tf]);
sol = bvp4c(@ode, @bc, solinit, bvpset('stats', 'on', 'RelTol', 1e-12, ...
    'AbsTol', 1e-10));
t1 = sol.parameters(1);
t2 = sol.parameters(2);
tf = sol.parameters(3);
N2 = 300;
t = [linspace(0, 1-eps(1), N2) linspace(1+eps(1), 2-eps(1), N2)...
     linspace(2+eps(1), 3, N2)]';
y = deval(sol, t);
xs = deval(sol, 1-eps(1)); % Knick
u = zeros(length(t),1);
for i = 1:length(u)
    if t(i) < 1
        u(i) = umax;
    elseif t(i) < 2
        u(i) = umin;
    else
        u(i) = umax;
    end
end
t = [t(t<1.0)*t1; t1+(t(t>1.0& t <= 2.0)-1.0)*(t2-t1);...
     t2+(t(t>2.0)-2.0)*(tf-t2)]; % Entnormierung
% H = 1+p(:,1).*(-0.5*x(:, 1)+x(:,2))+p(:,2).*(-x(:, 2)+u);
x = y(1:3, :)';    
p = y(4:6, :)';
H = 1+p(:,1).*u+p(:,2).*x(:,3)-p(:,3).*(x(:,2)*g+u)/l_s;
figure(1)
tiledlayout(3,2)
ax1 = nexttile;
plot(t, x(:,1))
ax2 = nexttile;
plot(t, x(:,2))
ax3 = nexttile;
plot(t, x(:,3))
ax4 = nexttile;
plot(t, p(:,1))
hold on
plot(results_dual(:,1), results_dual(:,2))
ax5 = nexttile;
plot(t, p(:,2))
hold on
plot(results_dual(:,1), results_dual(:,3))
ax6 = nexttile;
plot(t, p(:,3))
linkaxes([ax1, ax2, ax3,ax4,ax5,ax6], 'x')
hold on
plot(results_dual(:,1), results_dual(:,4))
figure(2)
tiledlayout(2,1)
bx1 = nexttile;
plot(t,u);
bx2 =nexttile;
plot(t,H);
linkaxes([bx1, bx2],'x')

    function yd = ode(t, y, int_nr, par)
        % kanonisches Dgl.-System
        x = y(1:3);
        p = y(4:6);
        t1 = par(1);
        t2 = par(2);
        tf = par(3);
        if int_nr == 1          % 1. Teilintervall
            u = umax;
            dt_dtnorm = t1;     % t = tnorm*ts;
        elseif int_nr == 2                   % 2. Teilintervall
            u = umin;
            dt_dtnorm = t2-t1;  % t = ts+(tnorm-1.0)*(tf-ts);
        else
            u = umax;
            dt_dtnorm = tf -t2;
        end
        yd = [u; x(3); -g/l_s*x(2)-u/l_s;...
            0; g/l_s*p(3); -p(2)]*dt_dtnorm;
    end
    function res = bc(y_0, y_f, par)
        % Residuen Rand- und Transversalitaetsbedingungen
        x_0 = y_0(1:3, 1);              % x(0)
        x_f = y_f(1:3, 3);              % x(tf)
        p_f = y_f(4:6, 3);              % p(tf)
        p_s2 = y_f(4:6, 2);             % p(t2)        
        p_s1 = y_f(4:6, 1);             % p(t1)
        u_f = umax;                     % u(tf)
        H_f = 1+p_f(1)*u_f + p_f(2)*x_f(3) -p_f(3)*(g*x_f(3) + u_f)/l_s; % H(tf)
        res = [x_0-x0; ...              % x(0)  = x0
            x_f-xf; ...                 % x(tf) = xf
            H_f; ...                    % H(tf) = 0
            y_f(:, 1)-y_0(:, 2); ...    % Stetigkeit x, p in ts
            y_f(:, 2)-y_0(:, 3); ...    % -""-
            p_s1(1)-p_s1(3)/l_s; ...    % Schaltbedingung
            p_s2(1)-p_s2(3)/l_s];       % -""-            
    end


end