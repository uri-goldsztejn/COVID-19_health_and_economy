function [x] = diff_system(t,y,theta,kappa,theta_s, kappa_s,new_patients,new_deaths)
% System of differential equations, which implements the compartmental
% model of our work.


% enum variables
S = 1; IsoS = 2; E=3; IsoE=4; I=5; H=6; IsoI=7; D=8; R=9; IsoR=10; S_S=11; IsoS_S=12; E_S=13;
IsoE_S=14; I_S=15; H_S=16; IsoI_S = 17; D_S=18; R_S=19; IsoR_S=20;
Y_S_L = 21; Y_S_H = 22; Y_G_L = 23; Y_G_H = 24; C_H = 25; C_L = 26; C_D = 27;


N = 330e6; % total population
x = zeros(27,1); %output vector (derivatives of states);

% dS/dt
beta = 1; % probabilty of encounter with asymptomatic patient
epsilon = beta/10; % probabilty of encounter with symptomatic patient
v = [y(H) + y(H_S); y(D) + y(D_S); new_patients; new_deaths; sum(y(21:24)) - sum(y(25:27))]/N; %v vector as described in the Appendix

tS1 = beta*(y(I) + y(I_S))*y(S)/N;
tS2 = epsilon*(y(H) + y(H_S))*y(S)/N;
x(S) = -tS1 -tS2  + kappa*y(IsoS) - theta'*v*y(S);

%dIso^S/dt
x(IsoS) =  -kappa*y(IsoS) + theta'*v*y(S);

%dE/dt
gamma_rate = 1/6.4; %1/days untill becoming contagious
tE1 = (gamma_rate*y(E));
x(E) = tS1 + tS2 + kappa*y(IsoE) - theta'*v*y(E) - tE1;

%dIso^E/dt
tIsoE1 =  (gamma_rate*y(IsoE));
x(IsoE) =  -kappa*y(IsoE) + theta'*v*y(E) - tIsoE1;

%dI^S/dt
p = 0.03; % probability of becoming symptomatic
chi = 1/7; %1/days to transition from asymptomatic to symptomatic
rho = 1/14; %1/days to recover after showing symptoms
delta = 1/7;

max_HFR = 0.015; %initial HFR
min_HFR = 0.005; % final HFR

% The HFR decreases over time thanks to medical advances
if t<365
    HFR = (365-t)/365*(max_HFR - min_HFR) + min_HFR ; % hospitalization fatality rate without health care saturation
else
    HFR = min_HFR;
end
m = 1; %slope of HFR as a funtion of patients when the system saturates

if (y(H) + y(H_S) < N*(0.08/100)) % health care system H not strained or saturated
    HFR_adjusted = HFR; 
elseif (y(H)+ y(H_S) < N*(0.22/100)) % health care system strained 
    HFR_adjusted = HFR*1.5;
else    % health care system saturated
    HFR_adjusted  = HFR*2;
end

tH1 =  y(H);
x(H) =   -delta*HFR_adjusted*tH1 - rho*(1-HFR_adjusted)*y(H) + p*chi*y(I) + p*chi*y(IsoI);

%dD/dt
x(D) = delta*HFR_adjusted*tH1;

%dI/dt
mu = 1/14; %1/days to recover without showing symptoms
tI1 = (mu*y(I));
x(I) =  tE1  + kappa*y(IsoI) - theta'*v*y(I) - (1-p)*tI1 - p*chi*y(I);

%dIso^I/dt
tIsoI1 =  mu*y(IsoI);
x(IsoI) =  -kappa*y(IsoI) + theta'*v*y(I) + tIsoE1 - (1-p)*tIsoI1 - p*chi*y(IsoI);

%dR/dt
x(R) =  (1-p)*tI1 +  rho*(1-HFR_adjusted)*tH1 + kappa*y(IsoR) - theta'*v*y(R);%(1-p)*mu*y(I)+

%dIso^R/dt
x(IsoR) = -kappa*y(IsoR) + theta'*v*y(R) + (1-p)*tIsoI1;

%% Equations for the senior population

%dS/dt
x(S_S) = -beta*(y(I) + y(I_S))*y(S_S)/N - epsilon*(y(H) + y(H_S))*y(S_S)/N + kappa_s*y(IsoS_S) - theta_s'*v*y(S_S);

%dIso^S/dt
x(IsoS_S) = -kappa_s*y(IsoS_S) + theta_s'*v*y(S_S);

%dE/dt
gamma_rate = 1/6.4; %1/days till becoming contagious
x(E_S) = beta*(y(I) + y(I_S))*y(S_S)/N + epsilon*(y(H) + y(H_S))*y(S_S)/N + kappa_s*y(IsoE_S) - theta_s'*v*y(E_S) - gamma_rate*y(E_S);

%dIso^E/dt
x(IsoE_S)=  -kappa_s*y(IsoE_S) + theta_s'*v*y(E_S) - gamma_rate*y(IsoE_S);

%dI/dt
p = 0.15; % probability of becoming sy]mptomatic being senior
chi = 1/7*1.3; %1/days to transition from asymptomatic to symptomatic
rho = 1/14; %1/days to recover after showing symptoms

max_HFR_S = 0.12; %initial HFR
min_HFR_S = 0.09; % final HFR

if t<365
    HFR_s = (365-t)/365*(max_HFR_S - min_HFR_S) + min_HFR_S ; % hospitalization fatality rate without health care saturation
else
    HFR_s = min_HFR_S;
end

HFR_s = (365-t)/365*(max_HFR_S - min_HFR_S) + min_HFR_S ; % hospitalization fatality rate without health care  saturation
% HFR_s = 0.12;
m = 1; %slope of HFR as a funtion of patients when the system saturates

if (y(H) + y(H_S)< N*(0.08/100))% health system H not strained or saturated
    HFR_adjusted_S = HFR_s; 
elseif (y(H) + y(H_S)< N*(0.22/100)) % health care system strained
    HFR_adjusted_S = HFR_s*1.5;
%     disp('Hospital capacity strained')
else    % health care system saturated
    HFR_adjusted_S  = HFR_s*2;

end

x(H_S) =  - delta*HFR_adjusted_S*y(H_S) - rho*(1-HFR_adjusted_S)*y(H_S) + p*chi*y(I_S) + p*chi*y(IsoI_S);

%dD/dt
x(D_S) = delta*HFR_adjusted_S*y(H_S);

%dI/dt
mu = 1/14; %1/days to recover without showing symptoms
x(I_S) = gamma_rate*y(E_S)  + kappa_s*y(IsoI_S) - theta_s'*v*y(I_S) - mu*(1-p)*y(I_S) - p*chi*y(I_S) ;%*(1-p) - chi*y(I_S) 

% tE1  + kappa*y(IsoI) - theta'*media_info*y(I) - (1-p)*tI1 - p*chi*y(I);

%dIso^I/dt
x(IsoI_S) = -kappa_s*y(IsoI_S) + theta_s'*v*y(I_S) +  gamma_rate*y(IsoE_S) - (1-p)*mu*y(IsoI_S) - p*chi*y(IsoI_S);

%dR/dt
x(R_S) = (1-p)*mu*y(I_S) + rho*(1-HFR_adjusted_S)*y(H_S) + kappa_s*y(IsoR_S) - theta_s'*v*y(R_S);

%dIso^R/dt
x(IsoR_S) = -kappa_s*y(IsoR_S) + theta_s'*v*y(R_S) + (1-p)*mu*y(IsoI_S);

%% Economic model

% Four groups: S_L (seniors with lower ability to work in an isolated state), S_H (seniors with higher ability to work in an isolated state), 
% G_L (non-seniors (general population) with lower ability to work in an isolated state), G_H (non-seniors (general population) with higher ability to work in an isolated state)

% Public policy will change as a function of Y/N, and that H the only 
% part of the economic model that H explicitly feeding into the public 
% policy and disease spread. 

% Productivity depends on current portion of (S + E + I), (H + D), (R + IsoR)
% and on (S_S + E_S + I_S), (H_S + D_S), (R_S + IsoR_S), and on theta

% Change in productivity depends in change in proportion of the above
% groups
POP_S_L = .2354;
POP_S_H = .0953;
POP_G_L = .464;
POP_G_H = .2053;

% 25.4% of S_L are in labor force, 36.8% of S_H, 69.3% of G_L, 88.6% of G_H
LAB_S_L = .254;
LAB_S_H = .368;
LAB_G_L = .693;
LAB_G_H = .886;

%% Conditional on being in the labor force:
% Average income for S_L H 37950.38
% Average income for S_H H 76255.1
% Average income for G_L H 33919.6
% Average income for G_H H 79094.35
INC_S_L = 37950.38;
INC_S_H = 76255.1;
INC_G_L = 33919.6;
INC_G_H = 79094.35;

phi = 0.7; % productivity penalty from lockdown for high type, function of theta, will vectorize
iota = 0.1; % boost in productivity when recovered for high type, function of theta, will vectorize

psi = 0.5; % productivity penalty from lockdown for low type, function of theta, will vectorize
xi = 0.2; % boost in productivity when recovered for low type, function of theta, will vectorize

% t_h = 0.4*(t_temp/segments); % will be a function of t
t_h=0.0;
% t_l = 0.3*(t_temp/segments);  % will be a function of t
t_l=0.0;

% What should cost be as a fraction of productivity?

c_h = 18000; % cost of treating higher-risk 
c_l = 12000; % cost of treating lower risk
c_d = 25; % based on a expected interest rate of 4%


% Ultimately, theta only factors into the economic part of the model through
% its impact on the value of phi, delta, psi, and gamma


%dY_S_L/dt
x(Y_S_L) = ((x(IsoS_S) + x(IsoE_S) + x(IsoI_S)) * (psi * (1-t_h)) ...
+ (x(S_S) + x(E_S) + x(I_S)) * (1-t_h) + (x(H_S) + x(D_S))*0 ...
+ x(IsoR_S) * ((psi+xi)* (1-t_h)) + ...
+ x(R_S)*(1-t_h))*(POP_S_L)*(LAB_S_L)*(INC_S_L);

%dY_S_H/dt
x(Y_S_H) = ((x(IsoS_S) + x(IsoE_S) + x(IsoI_S))* (phi * (1-t_l)) ...
+ (x(S_S) + x(E_S) + x(I_S)) *(1-t_l) + (x(H_S) + x(D_S))*0 ...
+ x(IsoR_S)  * ((phi + iota) * (1-t_l)) ... 
+ x(R_S) *(1-t_l))*(POP_S_H)*(LAB_S_H)*(INC_S_H);

%dY_G_L/dt
x(Y_G_L) = ((x(IsoS) + x(IsoE) + x(IsoI)) * (psi * (1-t_h)) ...
+ (x(S) + x(E) + x(I)) * (1- t_h) + (x(H) + x(D))*0 ...
+ x(IsoR)  * ((psi+xi) * (1-t_h)) ...
+ x(R) *(1-t_h))*(POP_G_L)*(LAB_G_L)*(INC_G_L);

%dY_G_H/dt
x(Y_G_H) = ((x(IsoS) + x(IsoE) + x(IsoI))  * (phi * (1-t_l)) ...
+ (x(S) + x(E) + x(I)) *(1-t_l) + (x(H) + x(D))*0 ... 
+ x(IsoR) * ((phi + iota) * (1-t_l)) ... 
+ x(R)*(1-t_l))*(POP_G_H)*(LAB_G_H)*(INC_G_H); 

%dc_h/dt
x(C_H)= c_h*(x(H_S));

%dc_l/dt
x(C_L) = c_l*(x(H));

%dc_d/dt
x(C_D) = c_d*(x(D)*((POP_G_L)/(POP_G_L + POP_G_H))*(LAB_G_L)*(INC_G_L) + ((POP_G_H)/(POP_G_L + POP_G_H))*(LAB_G_H)*(INC_G_H) ... 
+x(D_S)*((POP_S_L)/(POP_S_L + POP_S_H))*(LAB_S_L)*(INC_S_L) + ((POP_S_H)/(POP_S_L + POP_S_H))*(LAB_S_H)*(INC_S_H));

end