clear; close all
standardoptions
load('states_and_costates.mat');
t_mamax = t;
p1_mamax = p(:,1);
p2_mamax = p(:,2);
p3_mamax = p(:,3);
H_mamax  = H;
load('states_and_costates_amin_0.mat')
t_0 = t;
p1_0 = p(:,1);
p2_0 = p(:,2);
p3_0 = p(:,3);
H_0  = H;

tiledlayout(2,2)
nexttile
plot(t_mamax, p1_mamax, 'Color',plot_settings.color.UStuttDarkBlue); hold on
plot(t_0, p1_0,'Color',plot_settings.color.UStuttDarkOrange);
legend('$u_\mathrm{min} = - u_\mathrm{max}$','$u_\mathrm{min} = 0$', 'Interpreter','latex', 'Location','west')

xlim([0,1.6])
ylim([-1.64, -1.49])
xlabel('Zeit $[\si{\s}]$', 'Interpreter','latex')
ylabel('$p_1 \ [-]$', 'Interpreter','latex')
nexttile
plot(t_mamax, p2_mamax, 'Color',plot_settings.color.UStuttDarkBlue); hold on
plot(t_0, p2_0,'Color',plot_settings.color.UStuttDarkOrange);
xlim([0,1.6])
ylim([-6, 6])

xlabel('Zeit $[\si{\s}]$', 'Interpreter','latex')
ylabel('$p_2 \ [-]$', 'Interpreter','latex')
nexttile
plot(t_mamax, p3_mamax, 'Color',plot_settings.color.UStuttDarkBlue); hold on
plot(t_0, p3_0,'Color',plot_settings.color.UStuttDarkOrange);
xlim([0,1.6])
xlabel('Zeit $[\si{\s}]$', 'Interpreter','latex')
ylabel('$p_3 \ [-]$', 'Interpreter','latex')
nexttile
plot(t_mamax, H_mamax*1e10, 'Color',plot_settings.color.UStuttDarkBlue); hold on
plot(t_0, H_0*1e10,'Color',plot_settings.color.UStuttDarkOrange);
xlim([0,1.6])
xlabel('Zeit $[\si{\s}]$', 'Interpreter','latex')
ylabel('$H \cdot 10^{10}   \ [-]$', 'Interpreter','latex')
cleanfigure();
name='plot_a2.tikz';
matlab2tikz(name,'height', '\doubledfigureheight', 'width', '\figurewidth',...
     'extraaxisoptions','ylabel absolute,scaled ticks = false, yticklabel style={/pgf/number format/.cd,fixed,precision=2}, xlabel absolute,scaled ticks = false, xticklabel style={/pgf/number format/.cd,fixed,precision=10}')

