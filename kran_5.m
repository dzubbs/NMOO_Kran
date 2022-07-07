function kran_5()
global x0 xf umax umin g

t1= 0.697;
t2= .893;
tf= 1.592;
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
u_min = eps(10);
umax = u_max;
umin = -u_max;
N = 200;
x0 = [0 0 0]';
xf = [v_soll, 0,0]';
solinit = bvpinit([linspace(0, 1.0, 20) linspace(1.0, 2.0, 20), ...
    linspace(2.0, 3.0, 20)],[0 0 0 0 0 0], [t1; t2;tf]);
sol = bvp4c(@ode, @bc, solinit, bvpset('stats', 'on', 'RelTol', 1e-8));
t1 = sol.parameters(1);
t2 = sol.parameters(2);
tf = sol.parameters(3);
t = [linspace(0, 1-eps(1), 50) linspace(1+eps(1), 2-eps(1), 50)...
     linspace(2+eps(1), 3, 50)]';
y = deval(sol, t);
xs = deval(sol, 1-eps(1)); % Knick
t = [t(t<1.0)*t1; t1+(t(t>1.0& t <= 2.0)-1.0)*(t2-t1);...
     t2+(t(t>2.0)-2.0)*(tf-t2)]; % Entnormierung
% H = 1+p(:,1).*(-0.5*x(:, 1)+x(:,2))+p(:,2).*(-x(:, 2)+u);
x = y(1:3, :)';    
p = y(4:6, :)';

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
        p_s2 = y_f(4:6, 2);              % p(tf)        
        p_s1 = y_f(4:6, 1);              % p(t1)
        u_f = umax;                     % u(tf)
        H_f = 1+p_f(1)*u_f + p_f(2)*x_f(3) -p_f(3)*(g*x_f(3) + u_f)/l_s; % H(tf)
        res = [x_0-x0; ...              % x(0)  = x0
            x_f-xf; ...                 % x(tf) = xf
            H_f; ...                    % H(tf) = 0
            y_f(:, 1)-y_0(:, 2); ...    % Stetigkeit x, p in ts
            y_f(:, 2)-y_0(:, 3); ...
            p_s1(1)-p_s1(3)/l_s; ...
            p_s2(1)-p_s2(3)/l_s];                    % Schaltbedingung
    end


end