close all

%  load DSGE_NL_IRANO1_simul_results.mat 
  dynare DSGE_NL_IRANO1_simul
close all
fig = figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% IRFs To GOVERNMENT EXPENDITURE SHOCK
figure (1)
subplot(2,3,1)
plot(oo_.irfs.GDP_epscg*100/oo_.steady_state(39), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(GDP_epscg*100/oo_.stochastic_steady_state(39), 'r--', 'LineWidth', 2)
yyaxis left
title('GDP')
subplot(2,3,2)
plot(oo_.irfs.yd_epscg*100/oo_.steady_state(22), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(yd_epscg*100/oo_.stochastic_steady_state(22), 'r--', 'LineWidth', 2)
yyaxis left
title('y_d')
subplot(2,3,3)
plot(oo_.irfs.l_epscg*100/oo_.steady_state(27), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(l_epscg*100/oo_.stochastic_steady_state(27), 'r--', 'LineWidth', 2)
yyaxis left
title('L')
subplot(2,3,4)
plot(oo_.irfs.PI_epscg*100/oo_.steady_state(8), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(PI_epscg*100/oo_.stochastic_steady_state(8), 'r--', 'LineWidth', 2)
yyaxis left
title('\pi')
subplot(2,3,5)
plot(oo_.irfs.c_epscg*100/oo_.steady_state(2), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(c_epscg*100/oo_.stochastic_steady_state(2), 'r--', 'LineWidth', 2)
yyaxis left
title('c')
subplot(2,3,6)
plot(oo_.irfs.x_epscg*100/oo_.steady_state(10), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(x_epscg*100/oo_.stochastic_steady_state(10), 'r--', 'LineWidth', 2)
yyaxis left
title('I')
legend({'وضعیت‌پایدار‌غیرتصادفی(محور‌چپ)', 'وضعیت‌پایدارتصادفی(محورراست)' }, 'FontSize', 14, 'Location','SouthOutside','Orientation','horizontal')
suptitle('IRFs To GOVERNMENT EXPENDITURE SHOCK')
saveas(gcf,'CG.pdf')
saveas(gcf,'CG.eps')

% IRFs To MONETARY POLICY SHOCK
fig = figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
figure (2)
subplot(2,3,1)
plot(oo_.irfs.GDP_epsm*100/oo_.steady_state(39), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(GDP_epsm*100/oo_.stochastic_steady_state(39), 'r--', 'LineWidth', 2)
yyaxis left
title('GDP')
subplot(2,3,2)
plot(oo_.irfs.yd_epsm*100/oo_.steady_state(22), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(yd_epsm*100/oo_.stochastic_steady_state(22), 'r--', 'LineWidth', 2)
yyaxis left
title('y_d')
subplot(2,3,3)
plot(oo_.irfs.l_epsm*100/oo_.steady_state(27), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(l_epsm*100/oo_.stochastic_steady_state(27), 'r--', 'LineWidth', 2)
yyaxis left
title('L')
subplot(2,3,4)
plot(oo_.irfs.PI_epsm*100/oo_.steady_state(8), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(PI_epsm*100/oo_.stochastic_steady_state(8), 'r--', 'LineWidth', 2)
yyaxis left
title('\pi')
subplot(2,3,5)
plot(oo_.irfs.c_epsm*100/oo_.steady_state(2), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(c_epsm*100/oo_.stochastic_steady_state(2), 'r--', 'LineWidth', 2)
yyaxis left
title('c')
subplot(2,3,6)
plot(oo_.irfs.x_epsm*100/oo_.steady_state(10), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(x_epsm*100/oo_.stochastic_steady_state(10), 'r--', 'LineWidth', 2)
yyaxis left
title('I')
legend({'وضعیت‌پایدار‌غیرتصادفی(محور‌چپ)', 'وضعیت‌پایدارتصادفی(محورراست)' }, 'FontSize', 14, 'Location','SouthOutside','Orientation','horizontal') 
suptitle('IRFs To MONETARY POLICY SHOCK')
saveas(gcf,'M_dot.pdf')
saveas(gcf,'M_dot.eps')


% IRFs To TECHNOLOGY SHOCK
fig = figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
figure (3)
subplot(2,3,1)
plot(-oo_.irfs.GDP_epsA*100/oo_.steady_state(39), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(-GDP_epsA*100/oo_.stochastic_steady_state(39), 'r--', 'LineWidth', 2)
yyaxis left
title('GDP')
subplot(2,3,2)
plot(-oo_.irfs.yd_epsA*100/oo_.steady_state(22), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(-yd_epsA*100/oo_.stochastic_steady_state(22), 'r--', 'LineWidth', 2)
yyaxis left
title('y_d')
subplot(2,3,3)
plot(-oo_.irfs.l_epsA*100/oo_.steady_state(27), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(-l_epsA*100/oo_.stochastic_steady_state(27), 'r--', 'LineWidth', 2)
yyaxis left
title('L')
subplot(2,3,4)
plot(-oo_.irfs.PI_epsA*100/oo_.steady_state(8), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(-PI_epsA*100/oo_.stochastic_steady_state(8), 'r--', 'LineWidth', 2)
yyaxis left
title('\pi')
subplot(2,3,5)
plot(-oo_.irfs.c_epsA*100/oo_.steady_state(2), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(-c_epsA*100/oo_.stochastic_steady_state(2), 'r--', 'LineWidth', 2)
yyaxis left
title('c')
subplot(2,3,6)
plot(-oo_.irfs.x_epsA*100/oo_.steady_state(10), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(-x_epsA*100/oo_.stochastic_steady_state(10), 'r--', 'LineWidth', 2)
yyaxis left
title('I')
legend({'وضعیت‌پایدار‌غیرتصادفی(محور‌چپ)', 'وضعیت‌پایدارتصادفی(محورراست)' }, 'FontSize', 14, 'Location','SouthOutside','Orientation','horizontal') 
suptitle('IRFs To TECHNOLOGY SHOCK')
saveas(gcf,'A.pdf')
saveas(gcf,'A.eps')

fig = figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% IRFs To GOVERNMENT INVESTMENT SHOCK
figure (4)
subplot(2,3,1)
plot(oo_.irfs.GDP_epsig*100/oo_.steady_state(39), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(GDP_epsig*100/oo_.stochastic_steady_state(39), 'r--', 'LineWidth', 2)
yyaxis left
title('GDP')
subplot(2,3,2)
plot(oo_.irfs.yd_epsig*100/oo_.steady_state(22), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(yd_epsig*100/oo_.stochastic_steady_state(22), 'r--', 'LineWidth', 2)
yyaxis left
title('y_d')
subplot(2,3,3)
plot(oo_.irfs.l_epsig*100/oo_.steady_state(27), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(l_epsig*100/oo_.stochastic_steady_state(27), 'r--', 'LineWidth', 2)
yyaxis left
title('L')
subplot(2,3,4)
plot(oo_.irfs.PI_epsig*100/oo_.steady_state(8), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(PI_epsig*100/oo_.stochastic_steady_state(8), 'r--', 'LineWidth', 2)
yyaxis left
title('\pi')
subplot(2,3,5)
plot(oo_.irfs.c_epsig*100/oo_.steady_state(2), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(c_epsig*100/oo_.stochastic_steady_state(2), 'r--', 'LineWidth', 2)
yyaxis left
title('c')
subplot(2,3,6)
plot(oo_.irfs.x_epsig*100/oo_.steady_state(10), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(x_epsig*100/oo_.stochastic_steady_state(10), 'r--', 'LineWidth', 2)
yyaxis left
title('I')
legend({'وضعیت‌پایدار‌غیرتصادفی(محور‌چپ)', 'وضعیت‌پایدارتصادفی(محورراست)' }, 'FontSize', 14, 'Location','SouthOutside','Orientation','horizontal') 
suptitle('IRFs To GOVERNMENT INVESTMENT SHOCK')
saveas(gcf,'IG.pdf')
saveas(gcf,'IG.eps')

fig = figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% IRFs To GOVERNMENT INVESTMENT SHOCK
figure (5)
subplot(2,3,1)
plot(oo_.irfs.GDP_ue2*100/oo_.steady_state(39), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(GDP_ue2*100/oo_.stochastic_steady_state(39), 'r--', 'LineWidth', 2)
yyaxis left
title('GDP')
subplot(2,3,2)
plot(oo_.irfs.yd_ue2*100/oo_.steady_state(22), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(yd_ue2*100/oo_.stochastic_steady_state(22), 'r--', 'LineWidth', 2)
yyaxis left
title('y_d')
subplot(2,3,3)
plot(oo_.irfs.l_ue2*100/oo_.steady_state(27), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(l_ue2*100/oo_.stochastic_steady_state(27), 'r--', 'LineWidth', 2)
yyaxis left
title('L')
subplot(2,3,4)
plot(oo_.irfs.PI_ue2*100/oo_.steady_state(8), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(PI_ue2*100/oo_.stochastic_steady_state(8), 'r--', 'LineWidth', 2)
yyaxis left
title('\pi')
subplot(2,3,5)
plot(oo_.irfs.c_ue2*100/oo_.steady_state(2), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(c_ue2*100/oo_.stochastic_steady_state(2), 'r--', 'LineWidth', 2)
yyaxis left
title('c')
subplot(2,3,6)
plot(oo_.irfs.x_ue2*100/oo_.steady_state(10), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(x_ue2*100/oo_.stochastic_steady_state(10), 'r--', 'LineWidth', 2)
yyaxis left
title('I')
legend({'وضعیت‌پایدار‌غیرتصادفی(محور‌چپ)', 'وضعیت‌پایدارتصادفی(محورراست)' }, 'FontSize', 14, 'Location','SouthOutside','Orientation','horizontal') 
suptitle('IRFs To OIL SHOCK')
saveas(gcf,'OIL.pdf')
saveas(gcf,'OIL.eps')






fig = figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% IRFs To GOVERNMENT INVESTMENT SHOCK
figure (6)
subplot(2,3,1)
plot(oo_.irfs.GDP_epsoil*100/oo_.steady_state(39), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(GDP_epsoil*100/oo_.stochastic_steady_state(39), 'r--', 'LineWidth', 2)
yyaxis left
title('GDP')
subplot(2,3,2)
plot(oo_.irfs.yd_epsoil*100/oo_.steady_state(22), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(yd_epsoil*100/oo_.stochastic_steady_state(22), 'r--', 'LineWidth', 2)
yyaxis left
title('y_d')
subplot(2,3,3)
plot(oo_.irfs.l_epsoil*100/oo_.steady_state(27), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(l_epsoil*100/oo_.stochastic_steady_state(27), 'r--', 'LineWidth', 2)
yyaxis left
title('L')
subplot(2,3,4)
plot(oo_.irfs.PI_epsoil*100/oo_.steady_state(8), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(PI_epsoil*100/oo_.stochastic_steady_state(8), 'r--', 'LineWidth', 2)
yyaxis left
title('\pi')
subplot(2,3,5)
plot(oo_.irfs.c_epsoil*100/oo_.steady_state(2), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(c_epsoil*100/oo_.stochastic_steady_state(2), 'r--', 'LineWidth', 2)
yyaxis left
title('c')
subplot(2,3,6)
plot(oo_.irfs.x_epsoil*100/oo_.steady_state(10), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(x_epsoil*100/oo_.stochastic_steady_state(10), 'r--', 'LineWidth', 2)
yyaxis left
title('I')
legend({'وضعیت‌پایدار‌غیرتصادفی(محور‌چپ)', 'وضعیت‌پایدارتصادفی(محورراست)' }, 'FontSize', 14, 'Location','SouthOutside','Orientation','horizontal') 
suptitle('IRFs To OIL PRICE SHOCK')
saveas(gcf,'OIL_Prc.pdf')
saveas(gcf,'OIL_Prc.eps')

fig = figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% IRFs To GOVERNMENT INVESTMENT SHOCK
figure (7)
subplot(2,3,1)
plot(oo_.irfs.GDP_eps_sigoil*100/oo_.steady_state(39), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(GDP_eps_sigoil*100/oo_.stochastic_steady_state(39), 'r--', 'LineWidth', 2)
yyaxis left
title('GDP')
subplot(2,3,2)
plot(oo_.irfs.yd_eps_sigoil*100/oo_.steady_state(22), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(yd_eps_sigoil*100/oo_.stochastic_steady_state(22), 'r--', 'LineWidth', 2)
yyaxis left
title('y_d')
subplot(2,3,3)
plot(oo_.irfs.l_eps_sigoil*100/oo_.steady_state(27), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(l_eps_sigoil*100/oo_.stochastic_steady_state(27), 'r--', 'LineWidth', 2)
yyaxis left
title('L')
subplot(2,3,4)
plot(oo_.irfs.PI_eps_sigoil*100/oo_.steady_state(8), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(PI_eps_sigoil*100/oo_.stochastic_steady_state(8), 'r--', 'LineWidth', 2)
yyaxis left
title('\pi')
subplot(2,3,5)
plot(oo_.irfs.c_eps_sigoil*100/oo_.steady_state(2), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(c_eps_sigoil*100/oo_.stochastic_steady_state(2), 'r--', 'LineWidth', 2)
yyaxis left
title('c')
subplot(2,3,6)
plot(oo_.irfs.x_eps_sigoil*100/oo_.steady_state(10), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(x_eps_sigoil*100/oo_.stochastic_steady_state(10), 'r--', 'LineWidth', 2)
yyaxis left
title('I')
legend({'وضعیت‌پایدار‌غیرتصادفی(محور‌چپ)', 'وضعیت‌پایدارتصادفی(محورراست)' }, 'FontSize', 14, 'Location','SouthOutside','Orientation','horizontal') 
suptitle('IRFs To OIL PRICE UNCERTAINTY SHOCK')
saveas(gcf,'OIL_Pr_UNC.pdf')
saveas(gcf,'OIL_Pr_UNC.eps')


fig = figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% IRFs To GOVERNMENT INVESTMENT SHOCK
figure (8)
subplot(2,3,1)
plot(oo_.irfs.GDP_eps_sigeps2*100/oo_.steady_state(39), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(GDP_eps_sigeps2*100/oo_.stochastic_steady_state(39), 'r--', 'LineWidth', 2)
yyaxis left
title('GDP')
subplot(2,3,2)
plot(oo_.irfs.yd_eps_sigeps2*100/oo_.steady_state(22), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(yd_eps_sigeps2*100/oo_.stochastic_steady_state(22), 'r--', 'LineWidth', 2)
yyaxis left
title('y_d')
subplot(2,3,3)
plot(oo_.irfs.l_eps_sigeps2*100/oo_.steady_state(27), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(l_eps_sigeps2*100/oo_.stochastic_steady_state(27), 'r--', 'LineWidth', 2)
yyaxis left
title('L')
subplot(2,3,4)
plot(oo_.irfs.PI_eps_sigeps2*100/oo_.steady_state(8), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(PI_eps_sigeps2*100/oo_.stochastic_steady_state(8), 'r--', 'LineWidth', 2)
yyaxis left
title('\pi')
subplot(2,3,5)
plot(oo_.irfs.c_eps_sigeps2*100/oo_.steady_state(2), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(c_eps_sigeps2*100/oo_.stochastic_steady_state(2), 'r--', 'LineWidth', 2)
yyaxis left
title('c')
subplot(2,3,6)
plot(oo_.irfs.x_eps_sigeps2*100/oo_.steady_state(10), 'k', 'LineWidth', 2)
yyaxis right
hold on
plot(x_eps_sigeps2*100/oo_.stochastic_steady_state(10), 'r--', 'LineWidth', 2)
yyaxis left
title('I')
legend({'وضعیت‌پایدار‌غیرتصادفی(محور‌چپ)', 'وضعیت‌پایدارتصادفی(محورراست)' }, 'FontSize', 14, 'Location','SouthOutside','Orientation','horizontal') 
suptitle('IRFs To OIL UNCERTAINTY SHOCK')
saveas(gcf,'OIL_UNC.pdf')
saveas(gcf,'OIL_UNC.eps')
