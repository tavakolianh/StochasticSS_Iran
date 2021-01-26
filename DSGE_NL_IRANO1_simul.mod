
var 

    d       ${d}$                         //preference shock
    c       $\tilde{c}$                   //consumption
    mu_z 	${\mu}^{z}$                  //trend growth rate of the economy (from neutral and investment specific technology)
    mu_I    ${\mu}^{I}$                //growth rate of investment-specific technology growth
    mu_A    ${A}$                      //growth rate of neutral technology
    lambda  $\tilde{\lambda}$            //Lagrange multiplier
    R       ${R}$                      //Nominal Interest rate
    PI      $\tilde{\Pi}$               //Inflation
    r       $\tilde{r}$                  //rental rate of capital
    x       $\tilde{x}$                  //investment
    u       $\tilde{u}$                  //capacity utilization
    q       $\tilde{q}$                 //Tobin's marginal q
    f       $\tilde{f}$               //variable for recursive formulation of wage setting
    ld      ${l}^{d}$                   //aggregate labor demand
    ldy      ${l}^{dy}$                   //aggregate labor demand in non-oil
    w       $\tilde{w}$               //real wage
    wstar   $\tilde{w*}$           //optimal real wage
    PIstarw $\tilde{\Pi_{w*}}$        //optimal wage inflation
    PIstar  $\tilde{\Pi *}$       //optimal price inflation
    g1      $\tilde{g}^{1}$         //variable 1 for recursive formulation of price setting
    g2      $\tilde{g}^{2}$          //variable 2 for recursive formulation of price setting
    yd      $\tilde{y}^{d}$            //aggregate output
    mc      $\tilde{mc}$              //marginal costs
    k       $\tilde{k}$               //capital
    vp      ${v}^{p}$                 //price dispersion term
    vw      ${v}^{w}$               //wage dispersion term
    l       ${l}$                   //aggregate labor bundle
    phi     ${\varphi}$             //labor disutility shock
    F       $\mathcal{F}$             //firm profits
    m       $m$                     // money
    m_dot   $\dot{m}$                  // modey Growth 
    cg       $\tilde{cg}$             //government expenditures
    ig       $\tilde{ig}$             //government expenditures
    kg       $\tilde{kg}$             //government expenditures
    yobs    $\Delta\log y$             // output data
    cobs    $\Delta\log c$             //consumption data
    xobs    $\Delta\log x$            //investment data
    cgobs   $\Delta\log cg$            //government expenditures data
    GDP Y_oil P_oil  K_oil No OR W_oil eps2 igobs poilobs yoilobs wobs sigoil sigeps2; 
varexo 
    epsd       $\varepsilon^{d}$ 
    epsphi     $\varepsilon^{\varphi}$  
    epsmu_I    $\varepsilon^{I}$  
    epsA       $\varepsilon^{A}$  
    epscg      $\varepsilon^{cg}$  
    epsig      $\varepsilon^{ig}$  
    epsm       $\varepsilon^{m}$
    epsoil ue2 eps_sigoil eps_sigeps2;

predetermined_variables k;

parameters 
           h            $h$                  //consumption habits
           betta        $\beta$              //discount factor
           gammma1      $\gamma_{1}$         //capital utilization, linear term
           gammma2      $\gamma_{2}$        //capital utilization, quadratic term
           delta        $\delta$             //depreciation rate
           kappa        $\kappa$             //capital adjustment costs parameter
           eta          $\eta$              //elasticity of substitution between labor varieties
           epsilon      $\epsilon$          //elasticity of substitution between goods varieties
           varpsi       $\psi$            //labor disutility parameter
           gammma       $\gamma$            //inverse Frisch elasticity
           nu_m         $\nu_{m}$            //money utility parameter
           chiw         $\chi_{w}$            //wage indexation parameter
           chi          $\chi$               //price indexation
           thetap       $\theta_{p}$        //Calvo parameter prices
           thetaw       $\theta_{w}$        //Calvo parameter wages
           alppha       $\alpha$            //capital share
           Rbar         $\bar{R}$            //steady state interest rate
           PIbar        $\bar{\Pi}$       //steady state inflation
           mdotbar     $\bar{\dot{dc}}$   //steady state money growth
           gammmam     $\gamma_{dc}$       //money growth smoothing coefficient Taylor rule
           gammmaPI     $\gamma_{\Pi}$    //feedback coefficient to inflation monetary policy rule
           gammmay      $\gamma_{y}$         //feedback coefficient to output growth deviation in monetary policy rule
           Phi          $\phi$              //firms fixed costs
           rhod         $\rho_{d}$          //autocorrelation preference shock
           rhophi       $\rho_{d}$           //autocorrelation labor disutility shock
           rhocg        $\rho_{cg}$          //autocorrelation government expenditure shock
           rhoig        $\rho_{ig}$          //autocorrelation government expenditure shock
           cgss         $\bar{cg}$          // steady state government
           igss         $\bar{ig}$          // steady state government
           deltag       $\delta_g$          // steady state government
           psig         $\psi_g$          // steady state government
           Lambdamu  	$\Lambda_{\mu}$    //steady state growth rate of investmentment-specific technology
           LambdaA      $\Lambda_{A}$      //steady state neutral technology growth 
           Lambdax      $\Lambda_{x}$       //steady state growth rate of investment
           LambdaYd     $\Lambda_{y^{d}}$    //steady state growth rate of output
           sigma_d      $\sigma_{d}$        //standard deviation preference shock
           sigma_phi    $\sigma_{\phi}$    //standard deviation labor disutility shock
           sigma_mu     $\sigma_{\mu}$    //standard deviation investment-specific technology
           sigma_A      $\sigma_{A}$        //standard deviation neutral technology
           sigma_cg     $\sigma_{cg}$        //standard deviation government expenditure
           sigma_ig     $\sigma_{ig}$        //standard deviation government expenditure
           sigma_m      $\sigma_{m}$    //standard deviation preference shock
           sigma_oil      $\sigma_{m}$    //standard deviation preference shock
rho_oil
alpha_oil gamma_o delta_o P_oil_ss  rho_eps2 rho_sigoil sigoil_ss rho_sigeps2 sigeps2_ss;

    delta     = 0.06; 
    epsilon   = 10;
    eta       = 10;
    Phi       = 0;
    gammma2   = 0.001;
    betta     = 0.97;
    h         = 0.97;
    varpsi    = 1.92;
    gammma    = 1.17;
    kappa     = 9.51;
    alppha    = 0.45;
    thetap    = 0.32;
    chi       = 0.43;
    thetaw    = 0.48;
    chiw      = 0.42;
    gammmam   = 0.2;
    gammmay   = -0.5;
    gammmaPI  = -1.5;
    PIbar     = (1 + 0.178)^(1/4);
    rhod      = 0.12;
    rhophi    = 0.93;
    sigma_A   = -4.97;
    sigma_d   = -5.51;
    sigma_phi = -5.36;
    sigma_mu  = -5.43;
    sigma_m   = -5.85;
    sigma_cg   = -2.85;
    sigma_ig   = -2.85;
    sigma_oil   = -5.85;
    Lambdamu  = 1.4e-3;
    LambdaA   = 4.8e-3;
    nu_m      = 1.5;
    rhocg     = 0.5;
    rhoig     = 0.8;
    cgss      = 0.05;
    igss      = 0.2;
    psig      = 0.65;
    deltag    = 0.048;
rho_oil = 0.20;
alpha_oil = 0.03;
P_oil_ss =1;
rho_eps2 = 0.5;
gamma_o =  0.7634;
delta_o = 0.0047;
rho_sigoil =0.7;
sigoil_ss =1;
rho_sigeps2 =0.7;
sigeps2_ss =1;

model; 

// 1. FOC consumption

d*(c-h*c(-1)*mu_z^(-1))^(-1)-h*betta*d(+1)*(c(+1)*mu_z(+1)-h*c)^(-1)=lambda;

// 2. Euler equation

lambda=betta*lambda(+1)*mu_z(+1)^(-1)/PI(+1)*R;

// 3. FOC capital utilization

r=gammma1+gammma2*(u-1);

// 4. FOC capital

q=betta*lambda(+1)/lambda*mu_z(+1)^(-1)*mu_I(+1)^(-1)*((1-delta)*q(+1)+r(+1)*u(+1)-(gammma1*(u(+1)-1)+gammma2/2*(u(+1)-1)^2));

// 5. FOC investment

1=q*(1-(kappa/2*(x/x(-1)*mu_z-Lambdax)^2)-(kappa*(x/x(-1)*mu_z-Lambdax)*x/x(-1)*mu_z))
  +betta*q(+1)*lambda(+1)/lambda*mu_z(+1)^(-1)*kappa*(x(+1)/x*mu_z(+1)-Lambdax)*(x(+1)/x*mu_z(+1))^2;

// 6-7. Wage setting

f=(eta-1)/eta*wstar^(1-eta)*lambda*w^eta*ld+betta*thetaw*(PI^chiw/PI(+1))^(1-eta)*(wstar(+1)/wstar*mu_z(+1))^(eta-1)*f(+1);

f=varpsi*d*phi*PIstarw^(-eta*(1+gammma))*ld^(1+gammma)+betta*thetaw*(PI^chiw/PI(+1))^(-eta*(1+gammma))*(wstar(+1)/wstar*mu_z(+1))^(eta*(1+gammma))*f(+1);

// 8. demand for money

nu_m/m=((R-1)/R)*lambda;

// 9-11. firm's price setting

g1=lambda*mc*yd+betta*thetap*(PI^chi/PI(+1))^(-epsilon)*g1(+1);

g2=lambda*PIstar*yd+betta*thetap*(PI^chi/PI(+1))^(1-epsilon)*PIstar/PIstar(+1)*g2(+1);

epsilon*g1=(epsilon-1)*g2;

// 12-13. optimal inputs

u*k/ldy=alppha/(1-alppha)*w/r*mu_z*mu_I;

mc=(1/(1-alppha))^(1-alppha)*(1/alppha)^alppha*w^(1-alppha)*r^alppha*(kg(-1)*Y_oil)^(-psig);

// 14. law of motion wages

1=thetaw*(PI(-1)^chiw/PI)^(1-eta)*(w(-1)/w*mu_z^(-1))^(1-eta)+(1-thetaw)*PIstarw^(1-eta);

// 15. law of motion prices
1=thetap*(PI(-1)^chi/PI)^(1-epsilon)+(1-thetap)*PIstar^(1-epsilon);

// 16-17. Monetary Policy

m_dot = m*PI/m(-1);

m_dot/mdotbar=(m_dot(-1)/mdotbar)^gammmam*((PI/PIbar)^gammmaPI*((GDP/GDP(-1)*mu_z)/exp(LambdaYd))^gammmay)^(1-gammmam)*exp(epsm);

// 18. government

kg = (1-deltag)*kg(-1)+ig;

log(ig) = (1-rhoig)*log(igss)+rhoig*(log(ig(-1)))+epsig;

log(cg) = (1-rhocg)*log(cgss)+rhocg*(log(cg(-1)))+epscg;

//19-20. Market clearing

GDP=yd+Y_oil;

GDP=c+x+cg+ig+mu_z^(-1)*mu_I^(-1)*(gammma1*(u-1)+gammma2/2*(u-1)^2)*k;

yd=(mu_A*mu_z^(-1)*(u*k)^alppha*ldy^(1-alppha)*(kg(-1)*Y_oil)^psig-Phi)/vp;

//21-23. Price and wage dispersion terms

l=vw*ld; 

vp=thetap*(PI(-1)^chi/PI)^(-epsilon)*vp(-1)+(1-thetap)*PIstar^(-epsilon);

vw=thetaw*(w(-1)/w*mu_z^(-1)*PI(-1)^chiw/PI)^(-eta)*vw(-1)+(1-thetaw)*(PIstarw)^(-eta);

// 24. Law of motion for capital

k(+1)*mu_z*mu_I-(1-delta)*k-mu_z*mu_I*(1-kappa/2*(x/x(-1)*mu_z-Lambdax)^2)*x=0;

// 25. Profits

F=yd-1/(1-alppha)*w*ldy;

// 26. definition optimal wage inflation

PIstarw=wstar/w;

// exogenous processes

// 27. Preference Shock

log(d)=rhod*log(d(-1))+epsd;

//28. Labor disutility Shock

log(phi)=rhophi*log(phi(-1))+epsphi;

//29. Investment specific technology

log(mu_I)=Lambdamu+epsmu_I;

//30. Neutral technology

log(mu_A)=LambdaA+epsA; 

//31. Defininition composite technology

mu_z=mu_A^(1/(1-alppha))*mu_I^(alppha/(1-alppha));

yobs = mu_z*yd/yd(-1);

xobs = mu_I*x/x(-1);

cobs = mu_z*c/c(-1);
wobs = w*PI/w(-1);

cgobs = mu_z*cg/cg(-1);
igobs = mu_z*ig/ig(-1);
poilobs = P_oil/P_oil_ss;
yoilobs = Y_oil/Y_oil(-1);

P_oil =  (1-rho_oil)* P_oil_ss + rho_oil * P_oil(-1) + sigoil*epsoil;  // 41 P_oil 

Y_oil = exp(eps2)*K_oil(-1)^(gamma_o)*No^(1-gamma_o);   // 42 Y_oil

OR = (1-alpha_oil) * P_oil * Y_oil - W_oil * No;   // 43 Y_oil

No = (1-alpha_oil) * (1-gamma_o) * P_oil * Y_oil/W_oil;   // 44 Y_oil

K_oil  = (1-delta_o) * K_oil(-1)  + alpha_oil * P_oil * Y_oil; // 45 K_oil

W_oil = w;   // 46 W_oil

log(sigoil)= (1-rho_sigoil)* log(sigoil_ss) + rho_sigoil * log(sigoil(-1)) + eps_sigoil;

log(sigeps2)= (1-rho_sigeps2)* log(sigeps2_ss) + rho_sigeps2 * log(sigeps2(-1)) + eps_sigeps2;

ld = ldy+No;

eps2= rho_eps2 * eps2(-1) + sigeps2*ue2;
end;

shocks;
    var     epsd;    stderr exp(sigma_d);
    var     epsphi;  stderr exp(sigma_phi);
    var     epsmu_I; stderr exp(sigma_mu);
    var     epsA;    stderr exp(sigma_A);
    var     epsm;    stderr exp(sigma_m);
    var     epscg;   stderr exp(sigma_cg);
    var     epsig;   stderr exp(sigma_ig);
    var     epsoil;   stderr exp(sigma_oil);
    var     eps_sigoil;   stderr 0.01;
    var     eps_sigeps2;   stderr 0.01;
    var     ue2;   stderr exp(sigma_m);
end;






write_latex_dynamic_model;
write_latex_parameter_table;
write_latex_definitions;
collect_latex_files;

steady;
check;

//-----------------------------------7. PARAMETERS ESTIMATION-----------------------------

estimated_params;

          // PARAM NAME, INITVAL, LB, UB, PRIOR_SHAPE, PRIOR_P1, PRIOR_P2, PRIOR_P3, PRIOR_P4, JSCALE
          // PRIOR_SHAPE: BETA_PDF, GAMMA_PDF, NORMAL_PDF, INV_GAMMA_PDF

          h,              GAMMA_PDF,        0.97,  0.001;
          betta,          BETA_PDF,        0.97,  0.001;
          kappa,          GAMMA_PDF,       9.51,   0.2;
          eta,            GAMMA_PDF,        10,  0.1;
          varpsi,         GAMMA_PDF,       8.92,     0.1;
          gammma,         GAMMA_PDF,       1.17,   0.02;
          nu_m,           GAMMA_PDF,       1.5,   0.2;
          gammma2,          BETA_PDF,        0.001,  0.00001;
          chi,            BETA_PDF,        0.43,  0.02;
          //chiw,           BETA_PDF,        0.62,  0.05;
          thetap,         BETA_PDF,        0.32,  0.005;
          epsilon,         GAMMA_PDF,        10,  0.1;
          //thetaw,         BETA_PDF,        0.48,  0.01;
          alppha,         BETA_PDF,        0.55,  0.002;
          gammmam,        BETA_PDF,        0.2,   0.01;
          gammmaPI,       NORMAL_PDF,     -1.5,  0.02;
          gammmay,        NORMAL_PDF,     -0.5,  0.01;
          rhod,           BETA_PDF,        0.12,  0.01;
          rhophi,         BETA_PDF,        0.93,   0.005;
          //rhocg,          BETA_PDF,        0.8,   0.01;
          rhoig,          BETA_PDF,        0.8,   0.01;
          psig,          BETA_PDF,        0.05,   0.001;
          rho_oil,          BETA_PDF,        0.2,   0.001, 0, 1;
          alpha_oil,          BETA_PDF,        0.01,   0.001;
          rho_eps2,          BETA_PDF,        0.5,   0.01, 0, 1;
          Lambdamu,       GAMMA_PDF,       3.4e-3,  0.0001;
          LambdaA,        GAMMA_PDF,       1.8e-3,  0.0001;
          stderr epsd,    INV_GAMMA_PDF,   0.01,  inf;
          stderr epsphi,  INV_GAMMA_PDF,   0.01,  inf;
          stderr epsmu_I, INV_GAMMA_PDF,   0.01,  inf;
          stderr epsA,    INV_GAMMA_PDF,   0.01,  inf;
          stderr epscg,   INV_GAMMA_PDF,   0.01,  inf;
          stderr epsig,   INV_GAMMA_PDF,   0.01,  inf;
          stderr epsm,    INV_GAMMA_PDF,   0.01,  inf;
          stderr epsoil,    INV_GAMMA_PDF,   0.01,  inf;
          stderr ue2,    INV_GAMMA_PDF,   0.01,  inf;
          
         // dsge_prior_weight,  uniform_pdf, , , 0, 2;

end;

//----------------------------------------------------------------------------------------

//varobs yobs cobs cgobs xobs PI;
varobs  yobs yoilobs cobs PI m_dot xobs  igobs poilobs;//  wobs  poilobs ;// ; //GDPobs
//-----------------------------------8. ESTIMATION----------------------------------------

//dynare_sensitivity;//(lik_init=3);//, pvalue_ks=0.05, pvalue_corr=0.05);


//identification(advanced=0,max_dim_cova_group=3,prior_mc=0, lik_init=4);

/*
estimation
           (datafile=data_nl_dsge, filtered_vars, smoother, order=1,
            mh_drop = 0.5, mh_jscale =0.3, mh_replic=500000, mh_nblocks=1, 
            bayesian_irf, mode_check, irf=10,conditional_variance_decomposition =[1:4],  
            // mode_file = DSGE_VAR_DUAL_mode
            plot_priors =0,
            mode_compute=5,
            optim=('Hessian', 5, 'MaxIter', 20000, 'TolFun', 1e-24,'SaveFiles', 0),
            //mode_compute=6, 
            //optim=('nclimb-mh', 100000, 'ncov-mh', 10000, 'nscale-mh', 100000,'NumberOfMh', 1),
            graph_format=(pdf), silent_optimizer,
            lik_init=3, geweke_interval =[0.2 0.5]) yobs cobs xobs w PI m_dot;
*/

//-----------------------------------9. BVAR ESTIMATION------------------------------------

/*
bvar_density(datafile = data_BVAR, first_obs = 40, bvar_prior_tau=10, 
bvar_prior_decay=0.5, bvar_prior_lambda=10, bvar_prior_mu=4, bvar_prior_omega=0,
bvar_prior_train = 30) 8;
bvar_forecast(forecast = 4, bvar_replic = 10000, nobs =62) 8;
options_.irf=30;
options_.bvar.conf_sig=0.5;
bvar_irf(2,'Cholesky');

*/

//stoch_simul( order=2,irf=60) yd c R PI;

xx=load('DSGE_NL_IRANO1_results.mat');
M_.params=xx.M_.params;

stoch_simul(order=3,pruning,nograph,irf=30);



irflength           = 30;
irfburninlength     = 20000;

oo_.stochastic_steady_state = sss(oo_.dr, irfburninlength, options_.order);

if M_.dynare_version == '4.5.7'
    for ivariable = 1:M_.endo_nbr
    assignin('base',[deblank(M_.endo_names(ivariable,:)) '_sss'],...
                     oo_.stochastic_steady_state(ivariable,:)');                     
end

for ishock = 1:M_.exo_nbr    

    variables_irfsss  =  irfsss(oo_.dr,M_.Sigma_e(:,ishock),irflength,irfburninlength,options_.order)';

    for ivariable = 1:M_.endo_nbr
        assignin('base',[deblank(M_.endo_names(ivariable,:)) '_' deblank(M_.exo_names(ishock,:))],...
                          variables_irfsss(:,ivariable));                         
    end

end

else
    for ivariable = 1:M_.endo_nbr
    assignin('base',[char(M_.endo_names(ivariable,:)) '_sss'],...
                     oo_.stochastic_steady_state(ivariable,:)');                     
end

for ishock = 1:M_.exo_nbr    

    variables_irfsss  =  irfsss(oo_.dr,M_.Sigma_e(:,ishock),irflength,irfburninlength,options_.order)';

    for ivariable = 1:M_.endo_nbr
        assignin('base',[char(M_.endo_names(ivariable,:)) '_' char(M_.exo_names(ishock,:))],...
                          variables_irfsss(:,ivariable));                         
    end

end

end

for kkk=1:length(M_.endo_names)
    M_.endo_names(kkk,2) = {oo_.stochastic_steady_state(kkk)};
    M_.endo_names(kkk,3) = {oo_.steady_state(kkk)};
end




for kkk=1:length(M_.endo_names)
    sss_GDP(kkk,1)=oo_.stochastic_steady_state(kkk,1)'./oo_.stochastic_steady_state(39);
    sss_GDP(kkk,2)=oo_.steady_state(kkk,1)'./oo_.steady_state(39);
end


gammma1=(oo_.stochastic_steady_state(3,1))*oo_.stochastic_steady_state(4,1)/betta-(1-delta);
(oo_.stochastic_steady_state(10,1)+(oo_.stochastic_steady_state(3,1))^(-1)*oo_.stochastic_steady_state(4,1)^(-1)*(gammma1*(oo_.stochastic_steady_state(11,1)-1)+gammma2/2*(oo_.stochastic_steady_state(11,1)-1)^2)*oo_.stochastic_steady_state(24,1))/oo_.stochastic_steady_state(39,1)


