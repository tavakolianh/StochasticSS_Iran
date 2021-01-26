function [ys,params,check] = DSGE_NL_IRANO1_simul_steadystate(ys,exo,M_,options_)
% function [ys,params,check] = NK_baseline_steadystate(ys,exo,M_,options_)
% computes the steady state for the NK_baseline.mod and uses a numerical
% solver to do so
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%   - M_        [structure] Dynare model structure
%   - options   [structure] Dynare options structure
%
% Output: 
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - params    [vector] vector of parameter values
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 of not (allows to impose restrictions on parameters)

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = M_.param_names{ii};
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;

        
%% Enter model equations here

options=optimset(); % set options for numerical solver

% the steady state computation follows FVRR (2006), section 4.1
PI=PIbar;
u=1;
q=1;
d=1;
phi=1;
zeta=1;
mdotbar=PIbar;
m_dot=mdotbar;
LambdaYd= (LambdaA+alppha*Lambdamu)/(1-alppha);
mu_z=exp(LambdaYd);
mu_I=exp(Lambdamu);
mu_A=exp(LambdaA);
cg=cgss;
ig=igss;
kg=ig/deltag;
NNN= 0.1;
%set the parameter Lambdax
Lambdax=mu_z;

%set the parameter gammma1
gammma1=mu_z*mu_I/betta-(1-delta);
if gammma1<0 % parameter violates restriction; Preventing this cannot be implemented via prior restriction as it is a composite of different parameters and the valid prior region has unknown form
    check=1; %set failure indicator
    return; %return without updating steady states
end


r=1*gammma1;
R=1+(PI*mu_z/betta-1);

%set Rbar
Rbar=R;

PIstar=((1-thetap*PI^(-(1-epsilon)*(1-chi)))/(1-thetap))^(1/(1-epsilon));
PIstarw=((1-thetaw*PI^(-(1-chiw)*(1-eta))*mu_z^(-(1-eta)))/(1-thetaw))^(1/(1-eta));
Y_oil= 0.7538;
mc=(epsilon-1)/epsilon*(1-betta*thetap*PI^((1-chi)*epsilon))/(1-betta*thetap*PI^(-(1-epsilon)*(1-chi)))*PIstar;
w=(1-alppha)*(mc*((Y_oil*kg)^psig)*(alppha/r)^alppha)^(1/(1-alppha));
wstar=w*PIstarw;
vp=(1-thetap)/(1-thetap*PI^((1-chi)*epsilon))*PIstar^(-epsilon);
vw=(1-thetaw)/(1-thetaw*PI^((1-chiw)*eta)*mu_z^eta)*PIstarw^(-eta);
tempvaromega=alppha/(1-alppha)*w/r*mu_z*mu_I;



[ld,fval,exitflag]=fzero(@(ld)(1-betta*thetaw*mu_z^(eta-1)*PI^(-(1-chiw)*(1-eta)))/(1-betta*thetaw*mu_z^(eta*(1+gammma))*PI^(eta*(1-chiw)*(1+gammma)))...
-(eta-1)/eta*wstar/(varpsi*PIstarw^(-eta*gammma)*ld^gammma)*((1-h*mu_z^(-1))^(-1)-betta*h*(mu_z-h)^(-1))*...
((mu_A*mu_z^(-1)*vp^(-1)*((Y_oil*kg)^psig)*tempvaromega^alppha-tempvaromega*(1-(1-delta)*(mu_z*mu_I)^(-1)))*(1-NNN)*ld-vp^(-1)*Phi+Y_oil-cg-ig)^(-1),0.317,options);
if exitflag <1
    %indicate the SS computation was not sucessful; this would also be detected by Dynare
    %setting the indicator here shows how to use this functionality to
    %filter out parameter draws
    check=1; %set failure indicator
    return; %return without updating steady states
end


l=vw*ld;
ldy = (1-NNN)*ld;
No = NNN * ld;
k=tempvaromega*ldy;
x=(1-(1-delta)*(mu_z*mu_I)^(-1))*k;
yd=(mu_A/mu_z*k^alppha*ldy^(1-alppha)*(Y_oil*kg)^psig-Phi)/vp;
W_oil  = w ;%%%35
P_oil = 1 ;
QoilperGDP= 0.25;
% Y_oil = QoilperYH * yd/P_oil ;
c=(mu_A*mu_z^(-1)*vp^(-1)*((Y_oil*kg)^psig)*tempvaromega^alppha-tempvaromega*(1-(1-delta)*(mu_z*mu_I)^(-1)))*ldy+Y_oil-vp^(-1)*Phi-cg-ig;
lambda=(1-h*betta*mu_z^(-1))*(1-h/mu_z)^(-1)/c;
F=yd-1/(1-alppha)*w*ldy;
f=(eta-1)/eta*wstar*PIstarw^(-eta)*lambda*ld/(1-betta*thetaw*mu_z^(eta-1)*PI^(-(1-chiw)*(1-eta)));
f2=varpsi*d*phi*PIstarw^(-eta*(1+gammma))*ld^(1+gammma)/(1-betta*thetaw*(PI^chiw/PI)^(-eta*(1+gammma))*(wstar/wstar*mu_z)^(eta*(1+gammma)));
GDP=c+x+cg+ig+mu_z^(-1)*mu_I^(-1)*(gammma1*(u-1)+gammma2/2*(u-1)^2)*k;
g1=lambda*mc*yd/(1-betta*thetap*PI^((1-chi)*epsilon));
g2=epsilon/(epsilon-1)*g1;
m=nu_m*Rbar/((Rbar-1)*lambda);
yobs=mu_z;
cobs=mu_z;
cgobs=mu_z;
igobs=mu_z;
xobs=mu_I;
poilobs = 1;
yoilobs=1;
wobs=PI;
eps2 =0;
eps_sigoil=0;
sigoil = sigoil_ss;
eps_sigeps2=0;
sigeps2 = sigeps2_ss;
gamma_o = fsolve(@(gamma_o)  No - ((1-alpha_oil) * (1-gamma_o) * P_oil  * Y_oil /W_oil), 0.76, options);%32
K_oil= fsolve(@(K_oil)    Y_oil  - exp(eps2)*K_oil  ^(gamma_o)*No ^(1-gamma_o), 0.004, options); % 42,31
delta_o  = alpha_oil * P_oil  * Y_oil /K_oil; 

OR = Y_oil * ((1 - alpha_oil) *gamma_o * P_oil)   ;


%% end own model equations

params=NaN(NumberOfParameters,1);
for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = M_.endo_names{ii};
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end
