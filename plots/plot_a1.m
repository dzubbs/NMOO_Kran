clear; close all
standardoptions
addpath('..\')
results_X = load('kran.dat');       
results_X_amin0 = load('kran_a_min_0.dat');

figure(1)
tiledlayout(2,2)
ax1 = nexttile; 
plot(results_X(:,1),results_X(:,3), 'Color',plot_settings.color.UStuttDarkBlue)
hold on
plot(results_X_amin0(:,1),results_X_amin0(:,3),'Color',plot_settings.color.UStuttDarkOrange)
legend('$u_\mathrm{min} = - u_\mathrm{max}$','$u_\mathrm{min} = 0$', 'Direct Collocation', 'Interpreter','latex', 'Location','southeast')
xlim([0,1.6])
xlabel('Zeit $[\si{\s}]$', 'Interpreter','latex')
ylabel('$x \ [\si{\m\per\s}]$', 'Interpreter','latex')
ax2 = nexttile; 
plot(results_X(:,1),results_X(:,4),'Color',plot_settings.color.UStuttDarkBlue)
hold on
plot(results_X_amin0(:,1),results_X_amin0(:,4),'Color',plot_settings.color.UStuttDarkOrange)
xlim([0,1.6])

xlabel('Zeit $[\si{\s}]$', 'Interpreter','latex')
ylabel('$\varphi \ [\si{\radian}]$', 'Interpreter','latex')
ax3 = nexttile; 
plot(results_X(:,1),results_X(:,5),'Color',plot_settings.color.UStuttDarkBlue)
hold on
plot(results_X_amin0(:,1),results_X_amin0(:,5),'Color',plot_settings.color.UStuttDarkOrange)
xlim([0,1.6])

xlabel('Zeit $[\si{\s}]$', 'Interpreter','latex')
ylabel('$\dot{\varphi} \ [\si{\radian\per\s}]$', 'Interpreter','latex')
ax4 = nexttile;
plot(results_X(1:end-1,1),results_X(1:end-1,2),'Color',plot_settings.color.UStuttDarkBlue)
hold on
plot(results_X_amin0(1:end-1,1),results_X_amin0(1:end-1,2),'Color',plot_settings.color.UStuttDarkOrange)
xlim([0,1.6])

xlabel('Zeit $[\si{\s}]$', 'Interpreter','latex')
ylabel('$u \ [\si{\m\per\square\s}]$', 'Interpreter','latex')
% linkaxes([ax1, ax2, ax3,ax4], 'x')
cleanfigure();
name='plot_a1.tikz';
matlab2tikz(name,'height', '\doubledfigureheight', 'width', '\figurewidth',...
     'extraaxisoptions','ylabel absolute,scaled ticks = false, yticklabel style={/pgf/number format/.cd,fixed,precision=2}, xlabel absolute,scaled ticks = false, xticklabel style={/pgf/number format/.cd,fixed,precision=2}')

% 