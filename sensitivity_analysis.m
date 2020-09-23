%% This code accompanies the work in "Public policy and economic dynamics of COVID-19 spread: a mathematical modeling study"
% By Uri Goldsztejn, David Schwartzman, and Arye Nehorai. 2020

% If you find this code useful in your research, please consider citing.


% This code is used for the sensitivity analysis
% The model is run several times, with various parameter values

%% initialization
clear all;  clc ;close all;

% initial conditions (individuals in each state)
S = 1; IsoS = 2; E=3; IsoE=4; I=5; H=6; IsoI=7; D=8; R=9; IsoR=10; S_S=11; IsoS_S=12; E_S=13;
IsoE_S=14; I_S=15; H_S=16; IsoI_S = 17; D_S=18; R_S=19; IsoR_S=20;
Y_S_L = 21; Y_S_H = 22; Y_G_L = 23; Y_G_H = 24; C_H = 25; C_L = 26; C_D = 27;

N = 330e6; % total population
y0 = zeros(27,1); %initial values

y0_seed = [220499400,0,385569,0,170325,2068,0,14,42621,0,108604180,0,189907,0,80479,6224,0,368,18817,0,247091512391.000,291213254477.000,2411496327298.00,3180920390834.00,112093778,24848801,94687527];

y0(S) = y0_seed(S)*0.15; y0(IsoS) = y0_seed(S)*0.85;
y0(E) =  y0_seed(E)*0.15; y0(IsoE) = y0_seed(E)*0.85;
y0(I) =  y0_seed(I)*0.15; y0(IsoI) = y0_seed(I)*0.85;
y0(R) =  y0_seed(R)*0.15; y0(IsoR) = y0_seed(R)*0.85;
y0(D) = y0_seed(D);
y0(H) = y0_seed(H);

y0(S_S) = y0_seed(S_S)*0.05; y0(IsoS_S) = y0_seed(S_S)*0.95;
y0(E_S) =  y0_seed(E_S)*0.05; y0(IsoE_S) = y0_seed(E_S)*0.95;
y0(I_S) =  y0_seed(I_S)*0.05; y0(IsoI_S) = y0_seed(I_S)*0.95;
y0(R_S) =  y0_seed(R_S)*0.05; y0(IsoR_S) = y0_seed(R_S)*0.95;
y0(D_S) = y0_seed(D_S);
y0(H_S) = y0_seed(H_S);
%%


y0 = y0/sum(y0)*N; % percent of the population


% From 2018 American Community Survey data:
POP_S_L = .2354;
POP_S_H = .0953;
POP_G_L = .464;
POP_G_H = .2053;

POP_H = POP_S_L + POP_S_H;
POP_L = POP_G_L + POP_G_H;

% 25.4% of S_L are in labor force, 36.8% of S_H, 69.3% of G_L, 88.6% of G_H
LAB_S_L = .254;
LAB_S_H = .368;
LAB_G_L = .693;
LAB_G_H = .886;

%% Conditional on being in the labor force:
% Average income for S_L is 37950.38
% Average income for S_H is 76255.1
% Average income for G_L is 33919.6
% Average income for G_H is 79094.35
INC_S_L = 37950.38;
INC_S_H = 76255.1;
INC_G_L = 33919.6;
INC_G_H = 79094.35;

phi = 0.7; % productivity penalty from lockdown for high type (0.7)
iota = 0.1; % boost in productivity when recovered for high type (0.1)

psi = 0.5; % productivity penalty from lockdown for low type *(0.5)
xi = 0.2; % boost in productivity when recovered for low type (0.2)


y0(Y_S_L) = (y0(S_S) + y0(E_S) + y0(I_S)) * (POP_S_L)*(LAB_S_L)*(INC_S_L) ... 
+ (y0(IsoS_S) + y0(IsoE_S) + y0(IsoI_S))*psi*(POP_S_L)*(LAB_S_L)*(INC_S_L) ...
+ (y0(IsoR_S))*(POP_S_L)*(LAB_S_L)*(INC_S_L)*(psi + xi) + y0(R_S)*(POP_S_L)*(LAB_S_L)*(INC_S_L); 

y0(Y_S_H) = (y0(S_S) + y0(E_S) + y0(I_S)) * (POP_S_H)*(LAB_S_H)*(INC_S_H) ... 
+ (y0(IsoS_S) + y0(IsoE_S) + y0(IsoI_S))*phi*(POP_S_H)*(LAB_S_H)*(INC_S_H) ...
+ (y0(IsoR_S))* (POP_S_H)*(LAB_S_H)*(INC_S_H) * (phi + iota) + y0(R_S) * (POP_S_H)*(LAB_S_H)*(INC_S_H);

y0(Y_G_L) = (y0(S) + y0(E) + y0(I)) * (POP_G_L)*(LAB_G_L)*(INC_G_L) ... 
+ (y0(IsoS) + y0(IsoE_S) + y0(IsoI))*psi*(POP_G_L)*(LAB_G_L)*(INC_G_L) ...
+ (y0(IsoR))* (POP_G_L)*(LAB_G_L)*(INC_G_L) * (psi + xi) + y0(R) * (POP_G_L)*(LAB_G_L)*(INC_G_L);

y0(Y_G_H) = (y0(S) + y0(E) + y0(I)) * (POP_G_H)*(LAB_G_H)*(INC_G_H) ... 
+ (y0(IsoS) + y0(IsoE_S) + y0(IsoI))*phi*(POP_G_H)*(LAB_G_H)*(INC_G_H) ...
+ (y0(IsoR))* (POP_G_H)*(LAB_G_H)*(INC_G_H) * (phi + iota) + y0(R) * (POP_G_H)*(LAB_G_H)*(INC_G_H);

y0(C_H) = 0; y0(C_L) = 0; y0(C_D) = 0;

% Baseline pre-quarantine economy size:
econ_previrus = (((POP_S_L*LAB_S_L*INC_S_L + POP_S_H*LAB_S_H*INC_S_H)*(sum(y0(11:20))/sum(y0(1:20))) ...
+ (POP_G_L*LAB_G_L*INC_G_L + POP_G_H*LAB_G_H*INC_G_H)*(sum(y0(1:10))/sum(y0(1:20)))) - y0(C_H) -y0(C_L) - y0(C_D))*N;

y0_econ = y0(Y_S_L) + y0(Y_S_H) + y0(Y_G_L) + y0(Y_G_H) - y0(C_H) - y0(C_L) - y0(C_D);


%% Evaluate strictness of isolation

theta = [1;2;3;5;0];
theta_s = [1;2;3;5;0];
kappa = 0.001;
kappa_s = 0.00001;
% 
theta = theta/sum(theta);% normalization
theta_s = theta/sum(theta_s);
rate_into_quarantine = 1; % higher number reflects a more strict enforcement of qurantine
theta = theta*1;
theta_s = theta_s*1;



y_econ_init = econ_previrus/N;

%% solve the system 

% the system is solved for 7 days at a time. 

weeks = 76;
t_end = round(weeks*7);

segments = floor(t_end/7);
y = zeros(t_end,length(y0));
Y_econ = zeros(t_end,1);

factors = [100 50 10 5 2 1 0 1e-3]; 

outputs_q = zeros(length(factors),5);
for s = 1:length(factors)
    
    infected = 0;
    if factors(s) == 0
        theta_tmp = theta *factors(s);
        theta_s_tmp = theta_s *factors(s);

        kappa_tmp = kappa*factors(s);
        kappa_s_tmp = kappa_s*factors(s);

        
    else
        theta_tmp = theta /factors(s);
        theta_s_tmp = theta_s /factors(s);

        kappa_tmp = kappa*factors(s);
        kappa_s_tmp = kappa_s*factors(s);
    end
    for i = 1 : segments

        curr_tspan = (i-1)*8+1:i*8; %current 14 days being evaluated
        if i == 1
            [t_temp,y_temp] = ode45(@(t,y) diff_system_sensitivity(t,y,theta_tmp,kappa_tmp,theta_s_tmp, kappa_s_tmp,0,0,1,0.1), curr_tspan ,y0);
            y_temp = floor(y_temp);
            y(curr_tspan,:) = y_temp;

            Y_econ(curr_tspan) = (sum(y_temp(:,21:24),2) - sum(y_temp(:,25:27),2)) ./N;
        else
    shift = -i +1;
       
       current_time = curr_tspan+shift;
            new_patients = sum(diff(y_temp(:,H)) + diff(y_temp(:,H_S)));
            new_deaths = sum(diff(y_temp(:,D)) + diff(y_temp(:,D_S)));

            [t_temp,y_temp] = ode45(@(t,y) diff_system_sensitivity(t,y,theta_tmp,kappa_tmp,theta_s_tmp, kappa_s_tmp,new_patients,new_deaths,1,0.1), curr_tspan,y_temp(end,:));
            y_temp = floor(y_temp);
            y(curr_tspan+shift,:) = y_temp;
            Y_econ(curr_tspan+shift) = (sum(y_temp(:,21:24),2) - sum(y_temp(:,25:27),2)) ./N;

        end
        
    end

y = y(1:end-1,:);
Y_econ = Y_econ(1:end-1,:);


hospitalized = [y(:,H) + y(:,H_S)];


p = 0.03; % probability of becoming symptomatic
chi = 1/7; 

p_S = 0.15; % probability of becoming sy]mptomatic being senior
chi_S = 1/7*1.3; %1/days to transition from asymptomatic to symptomatic

infected = sum(p*chi*y(:,I) + p*chi*y(:,IsoI) + p_S*chi_S*y(:,I_S) + p_S*chi_S*y(:,IsoI_S));

outputs_q(s,:) = [factors(s) ,infected, max(hospitalized), max(y(t_end,D) + y(t_end,D_S)),(Y_econ(end) - y_econ_init)/y_econ_init*100];
end


%% Evaluate relative importance to the economy

theta(:,1) = [0;0;0;0;0];

theta(:,2) = [0;0;0;0;0];
theta(:,3) = [0;0;0;0;0];

theta(:,4) = [0;0;0;0;0];


theta_s = theta;
kappa = 0.00;
kappa_s = 0.000;
theta(end,4) = -1e-8;
theta(end,3) = -1e-9;
theta(end,2) = -1e-10;
theta(end,1) = 0;

y_econ_init = econ_previrus/N;


weeks = 76;
t_end = round(weeks*7);

segments = floor(t_end/7);
y = zeros(t_end,length(y0));
Y_econ = zeros(t_end,1);


outputs_rel_theta = zeros(4,5);
for s = 1:4
    
    infected = 0;
    theta_tmp = theta(:,s);
    theta_s_tmp = theta_s(:,s);
    kappa_tmp = kappa;
    kappa_s_tmp = kappa_s;

    for i = 1 : segments

        curr_tspan = (i-1)*8+1:i*8;
        if i == 1
            [t_temp,y_temp] = ode45(@(t,y) diff_system_sensitivity(t,y,theta_tmp,kappa_tmp,theta_s_tmp, kappa_s_tmp,0,0,1,0.1), curr_tspan ,y0);
            y_temp = floor(y_temp);
            y(curr_tspan,:) = y_temp;

            Y_econ(curr_tspan) = (sum(y_temp(:,21:24),2) - sum(y_temp(:,25:27),2)) ./N;
        else
            shift = -i +1;
       
            current_time = curr_tspan+shift;
       
            new_patients = sum(diff(y_temp(:,H)) + diff(y_temp(:,H_S)));
            new_deaths = sum(diff(y_temp(:,D)) + diff(y_temp(:,D_S)));

            [t_temp,y_temp] = ode45(@(t,y) diff_system_sensitivity(t,y,theta_tmp,kappa_tmp,theta_s_tmp, kappa_s_tmp,new_patients,new_deaths,1,0.1), curr_tspan,y_temp(end,:));
            y_temp = floor(y_temp);
            y(curr_tspan+shift,:) = y_temp;
            Y_econ(curr_tspan+shift) = (sum(y_temp(:,21:24),2) - sum(y_temp(:,25:27),2)) ./N;

        end
        
    end


y = y(1:end-1,:);
Y_econ = Y_econ(1:end-1,:);
hospitalized = [y(:,H) + y(:,H_S)];

p = 0.03; % probability of becoming symptomatic
chi = 1/7; 
p_S = 0.15; % probability of becoming sy]mptomatic being senior
chi_S = 1/7*1.3; %1/days to transition from asymptomatic to symptomatic

infected = sum(p*chi*y(:,I) + p*chi*y(:,IsoI) + p_S*chi_S*y(:,I_S) + p_S*chi_S*y(:,IsoI_S));

outputs_rel_theta(s,:) = [factors(s) ,infected, max(hospitalized), max(y(t_end,D) + y(t_end,D_S)),(Y_econ(end) - y_econ_init)/y_econ_init*100];
end

%%
theta = [1;2;3;5;0];
theta_s = [1;2;3;5;0];
kappa = 0.00;
kappa_s = 0.000;

theta = theta/sum(theta);% normalization
theta_s = theta/sum(theta_s);
rate_into_quarantine = 1; % higher number reflects a more strict enforcement of qurantine

theta = theta*00;
theta_s = theta_s*000;



%% Evaluate different values of beta

kappa = 0.00;
kappa_s = 0.0000;

theta = theta*0;
theta_s = theta_s*0;

betas = 2:-0.1:0;


outputs_beta = zeros(length(betas),5);
for s = 1:length(betas)
   
for i = 1 : segments

    curr_tspan = (i-1)*8+1:i*8; 
    if i == 1
        [t_temp,y_temp] = ode45(@(t,y) diff_system_sensitivity(t,y,theta,kappa,theta_s, kappa_s,0,0,betas(s),0.1), curr_tspan ,y0);
        y_temp = floor(y_temp);
        y(curr_tspan,:) = y_temp;
        
        Y_econ(curr_tspan) = (sum(y_temp(:,21:24),2) - sum(y_temp(:,25:27),2)) ./N;
    else
        shift = -i +1;
       
       current_time = curr_tspan+shift;
        
        new_patients = sum(diff(y_temp(:,H)) + diff(y_temp(:,H_S)));
        new_deaths = sum(diff(y_temp(:,D)) + diff(y_temp(:,D_S)));
               
        [t_temp,y_temp] = ode45(@(t,y) diff_system_sensitivity(t,y,theta,kappa,theta_s, kappa_s,new_patients,new_deaths,betas(s),0.1), curr_tspan,y_temp(end,:));
        y_temp = floor(y_temp);
        y(curr_tspan+shift,:) = y_temp;
        Y_econ(curr_tspan+shift) = (sum(y_temp(:,21:24),2) - sum(y_temp(:,25:27),2))./N;
   
       
    end
end



y = y(1:end-1,:);
Y_econ = Y_econ(1:end-1,:);

infectedPeople = [y(:,I) + y(:,H) + y(:,E) + y(:,IsoI) + y(:,IsoE) + y(:,I_S) + y(:,H_S) + y(:,E_S) + y(:,IsoI_S) + y(:,IsoE_S)];
hospitalized = [y(:,H) + y(:,H_S)];

p = 0.03; % probability of becoming symptomatic
chi = 1/7; 
p_S = 0.15; % probability of becoming sy]mptomatic being senior
chi_S = 1/7*1.3; %1/days to transition from asymptomatic to symptomatic

infected = sum(p*chi*y(:,I) + p*chi*y(:,IsoI) + p_S*chi_S*y(:,I_S) + p_S*chi_S*y(:,IsoI_S));
outputs_beta(s,:) = [betas(s), infected, max(hospitalized), y(t_end,D) + y(t_end,D_S), (Y_econ(end) - y_econ_init)/y_econ_init*100];
end

%% Evaluate differnet values of epsilon
ep = 2:-0.1:0.3;

outputs_epsilon = zeros(length(ep),5);
for s = 1:length(ep)
   
for i = 1 : segments

    curr_tspan = (i-1)*8+1:i*8; %current 14 days being evaluated
    if i == 1
        [t_temp,y_temp] = ode45(@(t,y) diff_system_sensitivity(t,y,theta,kappa,theta_s, kappa_s,0,0,1,ep(s)), curr_tspan ,y0);
        y_temp = floor(y_temp);
        y(curr_tspan,:) = y_temp;
        
        Y_econ(curr_tspan) = (sum(y_temp(:,21:24),2) - sum(y_temp(:,25:27),2)) ./N;%(N - y_temp(:,D) - y_temp(:,D_S));
    else
        shift = -i +1;
       
       current_time = curr_tspan+shift;
        new_patients = sum(diff(y_temp(:,H)) + diff(y_temp(:,H_S)));
        new_deaths = sum(diff(y_temp(:,D)) + diff(y_temp(:,D_S)));
               
        [t_temp,y_temp] = ode45(@(t,y) diff_system_sensitivity(t,y,theta,kappa,theta_s, kappa_s,new_patients,new_deaths,1,ep(s)), curr_tspan,y_temp(end,:));
        y_temp = floor(y_temp);
        y(curr_tspan+shift,:) = y_temp;
        Y_econ(curr_tspan+shift) = (sum(y_temp(:,21:24),2) - sum(y_temp(:,25:27),2))./N;%(N - y_temp(:,D) - y_temp(:,D_S));
   
       
    end
end

y = y(1:end-1,:);
Y_econ = Y_econ(1:end-1,:);
infectedPeople = [y(:,I) + y(:,H) + y(:,E) + y(:,IsoI) + y(:,IsoE) + y(:,I_S) + y(:,H_S) + y(:,E_S) + y(:,IsoI_S) + y(:,IsoE_S)];
hospitalized = [y(:,H) + y(:,H_S)];

p = 0.03; % probability of becoming symptomatic
chi = 1/7; 
p_S = 0.15; % probability of becoming sy]mptomatic being senior
chi_S = 1/7*1.3; %1/days to transition from asymptomatic to symptomatic

infected = sum(p*chi*y(:,I) + p*chi*y(:,IsoI) + p_S*chi_S*y(:,I_S) + p_S*chi_S*y(:,IsoI_S));
outputs_epsilon(s,:) = [ep(s),infected,  max(hospitalized), y(t_end,D) + y(t_end,D_S), (Y_econ(end) - y_econ_init)/y_econ_init*100];
end

%% plot

y1 = (Y_econ(1)-y_econ_init)/y_econ_init*100;

lineWidth = 3;
line_dotted_width = 2;
set(0, 'DefaultAxesFontName', 'Arial');
figure('DefaultAxesFontSize',11)

set(gcf,'Color','white');

subplot(3,3,1)
plot(linspace(0,1,length(outputs_beta(:,2))),outputs_beta(:,3),'k','LineWidth',lineWidth);hold on
plot(linspace(0,1,length(outputs_epsilon(:,2))),outputs_epsilon(:,3),'--k','LineWidth',lineWidth);hold on
xline(0.5,':b','LineWidth',line_dotted_width)
ylabel('Hospitalizations')

subplot(3,3,2)
plot(linspace(0,1,length(outputs_beta(:,2))),outputs_beta(:,4),'k','LineWidth',lineWidth);hold on
plot(linspace(0,1,length(outputs_epsilon(:,2))),outputs_epsilon(:,4),'--k','LineWidth',lineWidth);hold on
ylabel('Deaths');
xline(0.5,':b','LineWidth',line_dotted_width)

subplot(3,3,3)
plot(linspace(0,1,length(outputs_beta(:,2))),outputs_beta(:,5),'k','LineWidth',lineWidth);hold on
plot(linspace(0,1,length(outputs_epsilon(:,2))),outputs_epsilon(:,5),'--k','LineWidth',lineWidth);hold on
ylabel('\Delta Y');
xline(0.5,':b','LineWidth',line_dotted_width)
legend('\beta','\epsilon')


subplot(3,3,4)
plot(linspace(0,1,length(outputs_q(:,2))),outputs_q(:,3),'k','LineWidth',lineWidth);hold on
ylabel('Hospitalizations')
xline(0.85,':b','LineWidth',line_dotted_width)

subplot(3,3,5)
plot(linspace(0,1,length(outputs_q(:,2))),outputs_q(:,4),'k','LineWidth',lineWidth);hold on
ylabel('Deaths')
xline(0.85,':b','LineWidth',line_dotted_width)

subplot(3,3,6)
plot(linspace(0,1,length(outputs_q(:,2))),outputs_q(:,5),'k','LineWidth',lineWidth);hold on
ylabel('\Delta Y');
xline(0.85,':b','LineWidth',line_dotted_width)



subplot(3,3,7)
plot(linspace(0,1,length(outputs_rel_theta(:,2-1))),outputs_rel_theta(1:end,3),'k','LineWidth',lineWidth);hold on
ylabel('Hospitalizations')
xline(0.02,':b','LineWidth',line_dotted_width)

subplot(3,3,8)
plot(linspace(0,1,length(outputs_rel_theta(:,2-1))),outputs_rel_theta(1:end,4),'k','LineWidth',lineWidth);hold on
ylabel('Deaths')
xline(0.02,':b','LineWidth',line_dotted_width)

subplot(3,3,9)
plot(linspace(0,1,length(outputs_rel_theta(:,2))),outputs_rel_theta(1:end,5),'k','LineWidth',lineWidth);hold on
ylabel('\Delta Y');
xline(0.02,':b','LineWidth',line_dotted_width)

% make figure nice

for i = 1:3
    subplot(3,3,i)
    ax = gca;
    set(gca,'box','off')
    xlabel("Control of contagiousness")
    ax.FontUnits = 'normalized';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'off';
    ax.TickLength = [0.05 0.1];
    set(gca, 'FontName', 'Arial')
end

for i = 4:6
    subplot(3,3,i)
    ax = gca;
    set(gca,'box','off')
    xlabel("Strictness of social isolation")
    ax.FontUnits = 'normalized';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'off';
    ax.TickLength = [0.05 0.1];
    set(gca, 'FontName', 'Arial')
end


for i = 7:9
    subplot(3,3,i)
    ax = gca;
    set(gca,'box','off')
    xlabel("Economic weight in policy")
    ax.FontUnits = 'normalized';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'off';
    ax.TickLength = [0.05 0.1];
    set(gca, 'FontName', 'Arial')
end

annotation('textbox', [0.09, 0.92, 0.05, 0.06], 'String', "A",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.37, 0.92, 0.05, 0.06], 'String', "B",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.66, 0.92, 0.05, 0.06], 'String', "C",'FontName','Arial','FontSize',13,'LineStyle','none')

annotation('textbox', [0.09, 0.62, 0.05, 0.06], 'String', "D",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.37, 0.62, 0.05, 0.06], 'String', "E",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.66, 0.62, 0.05, 0.06], 'String', "F",'FontName','Arial','FontSize',13,'LineStyle','none')

annotation('textbox', [0.09, 0.32, 0.05, 0.06], 'String', "G",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.37, 0.32, 0.05, 0.06], 'String', "H",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.66, 0.32, 0.05, 0.06], 'String', "I",'FontName','Arial','FontSize',13,'LineStyle','none')

