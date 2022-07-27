function kran_5()
% Miniprojekt 5 Zeitoptimale Umsteuerung Brückenkran Numerische Methoden 
%   der Optimierung und Optimalen Steuerung
% Exercise 5, indirect collocation
% Authors: Laura Kleckner and Marcel Dzubba
% Summer term 2022
close all
global x0 xf umax umin g l_s t1 t2 tf  % variables used in subfunctions --> globar
%%%%%%%%%%%%% chose case manually!!! %%%%%%%%%%%%%%%%%
choose_problem = 0; % 0: a_min = - a_max; 1: a_min = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Modellparameter
K_x = 0.99;
l_s = 1.0;
g   = 9.81;
v_soll = 0.4;
a_max  = 0.33;
umax = a_max/K_x;

switch choose_problem
    case 0 % a_min = - a_max
        t1= 0.6898200664133287; % switching time 1 (from julia)
        t2= 0.8888066240325581; % switching time 2 
        tf= 1.5918924609538354; % end time
        umin = -umax; 
        % load data from direct collocation for initializing 
        results_X = load('kran.dat');       
        N1 = (size(results_X,1)-1)/3; % same dimension as direct collocation
        results_dual = load('kran_dual.dat');   
    case 1 % a_min = 0
        t1= 0.5879036776475854;
        t2= 0.9798394627459757;
        tf= 1.6033736663115965;
        umin = 0;
        results_X = load('kran_a_min_0.dat');
        N1 = (size(results_X,1)-1)/3;
        results_dual = load('kran_a_min_0_dual.dat');
end
x0 = [0 0 0]'; % initial values 
xf = [v_soll, 0,0]'; % terminal values

% generate suitable starting values for bvp4c
solinit = bvpinit([linspace(0, 1.0, N1) linspace(1.0, 2.0, N1), ...
     linspace(2.0, 3.0, N1)],[0;0;0; 0; 0; 0],[t1; t2;tf]);
% use values from directe collocation 
solinit.y =  [results_X(1:end-1,3)';results_X(1:end-1,4)';results_X(1:end-1,5)';...
      results_dual(:,2)'; results_dual(:,3)'; results_dual(:,4)'];
% determine optimal solution 
sol = bvp4c(@ode, @bc, solinit, bvpset('stats', 'on', 'RelTol', 1e-13, ...
    'AbsTol', 5e-9));
% visualization 
t1 = sol.parameters(1) % switching time 1
t2 = sol.parameters(2) % switching time 2
tf = sol.parameters(3) % terminal time 
N2 = 200;              % 3*N1 points for visualization
t = [linspace(0, 1-eps(2), N2) linspace(1+eps(2), 2-eps(2), N2)...
     linspace(2+eps(2), 3, N2)]'; %generate normalizion grid  
y = deval(sol, t);      % get optimal solution

% build u vector according to interval
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
     t2+(t(t>2.0)-2.0)*(tf-t2)]; % denormalize time vector
x = y(1:3, :)';     % state vector
p = y(4:6, :)';     % co-state vector
H = 1+p(:,1).*u+p(:,2).*x(:,3)-p(:,3).*(x(:,2)*g+u)/l_s; % hamilton fct
%% plots
figure(1)
tiledlayout(2,2)
ax1 = nexttile; hold on
plot(t, x(:,1))
plot(results_X(:,1),results_X(:,3))
xlabel('time [s]')
ylabel('x2 [m/s]')
ax2 = nexttile; hold on
plot(t, x(:,2))
plot(results_X(:,1),results_X(:,4))
xlabel('time [s]')
ylabel('x2 [rad]')
ax3 = nexttile; hold on
plot(t, x(:,3))
plot(results_X(:,1),results_X(:,5))
xlabel('time [s]')
ylabel('x3 [rad/s]')
ax4 = nexttile;
plot(t,u);hold on
plot(results_X(1:end-1,1),results_X(1:end-1,2))
linkaxes([ax1, ax2, ax3,ax4], 'x')
leg = legend('Indirect Collocation', 'Direct Collocation');
leg.Layout.Tile = 'south';
figure(2)
tiledlayout(2,1)
plot(t,H);
xlabel('time [s]')
ylabel('Hamilton function [-]')
    function yd = ode(t, y, int_nr, par)
        % (co-) states
        x = y(1:3);
        p = y(4:6);
        t1 = par(1);
        t2 = par(2);
        tf = par(3);
        if int_nr == 1          % 1. interval
            u = umax;
            dt_dtnorm = t1;     % t = tnorm*ts;
        elseif int_nr == 2      % 2. interval
            u = umin;
            dt_dtnorm = t2-t1;  
        else
            u = umax;           % 3. interval
            dt_dtnorm = tf -t2; 
        end
        % [xdot;pdot]-Vector
        yd = [u; x(3); -g/l_s*x(2)-u/l_s;...
            0; g/l_s*p(3); -p(2)]*dt_dtnorm;
    end
    function res = bc(y_0, y_f, par)
        % multibound and transersality condition 
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
            y_f(:, 1)-y_0(:, 2); ...    % continuity x, p in ts
            y_f(:, 2)-y_0(:, 3); ...    % -""-
            p_s1(1)-p_s1(3)/l_s; ...    % switching condition 
            p_s2(1)-p_s2(3)/l_s];       % -""-            
    end
end