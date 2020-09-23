%% This code accompanies the work in "Public policy and economic dynamics of COVID-19 spread: a mathematical modeling study"
% By Uri Goldsztejn, David Schwartzman, and Arye Nehorai. 2020

% If you find this code useful in your research, please consider citing.


% This file is used to simulate the gradual release of the non-senior
% population
%% initialization
clear all; clc ;close all;

% initial conditions (individuals in each state)
S = 1; IsoS = 2; E=3; IsoE=4; I=5; H=6; IsoI=7; D=8; R=9; IsoR=10; S_S=11; IsoS_S=12; E_S=13;
IsoE_S=14; I_S=15; H_S=16; IsoI_S = 17; D_S=18; R_S=19; IsoR_S=20;
Y_S_L = 21; Y_S_H = 22; Y_G_L = 23; Y_G_H = 24; C_H = 25; C_L = 26; C_D = 27;


N = 330e6; % total population
y0 = zeros(27,1); %initial values

economy_init = 4.0803e12;
%% Initial conditions
% initial conditions for the epidemiological model
% here we definte the fraction of the population in each compartment

% y0(S) = 0.1; y0(IsoS) = 0.6; y0(E) =  2.4e-4; y0(IsoE) = 9.6e-4; y0(H) = 2e-4; y0(I) = 2.4e-4; y0(IsoI) = 9.6e-04; y0(D) = 3.5e-06;
% y0(R) = 7e-07; y0(IsoR) = 5.8e-6;
% 
% y0(S_S) = 0.01; y0(IsoS_S) = 0.21; y0(E_S) = 1.4e-4; y0(IsoE_S) = 2.4e-4; y0(H_S) = 2e-4; y0(I_S) = 2e-4; y0(IsoI_S) = 2e-4; y0(D_S) = 3.5e-06;
% y0(R_S) = 7e-07; y0(IsoR_S) = 5.8e-6;
% 


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

y0 = y0/sum(y0)*N; % percent of the population

% From 2018 American Community Survey data:
% 23.5 percent is 60+ with a lower ability to work in an isolated state (S_L), 9.5% is 60+ with a higher ability to work in an isolated state (S_H),
% 46.4% is under 60 (general population or non-seniors) with a lower ability to work in an isolated state (G_L), 20.5% is under 60 (general population or non-seniors) with a higher ability to work in an isolated state (G_H)
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


theta = [1;2;3;5;0];
theta_s = [1;2;3;5;0];
kappa = 0.00;
kappa_s = 0.000;
% 
theta = theta/sum(theta);% normalization
theta_s = theta/sum(theta_s);
rate_into_quarantine = 1; % higher number reflects a more strict enforcement of qurantine

theta = theta*0;
theta_s = theta_s*0;



%% solve the system 
% the system is solve for 14 days at a time. Then the number of new deaths and
% new patients is counted.

% 76 weeks of simulation
weeks = 76;
t_end = round(weeks*7);

segments = floor(t_end/7);
y = zeros(t_end,length(y0));
Y_econ = zeros(t_end,1);

reopen_flag = 0;
for i = 1 : segments

    curr_tspan = (i-1)*8+1:i*8; %current 14 days being evaluated
    if i == 1
        [t_temp,y_temp] = ode45(@(t,y) diff_system(t,y,theta,kappa,theta_s, kappa_s,0,0), curr_tspan ,y0);
        y_temp = floor(y_temp);
        y(curr_tspan,:) = y_temp;
        
        Y_econ(curr_tspan) = (sum(y_temp(:,21:24),2) - sum(y_temp(:,25:27),2)) ./N;%(N - y_temp(:,D) - y_temp(:,D_S));
    else
        shift = -i +1;
        

        current_time = curr_tspan+shift;
          if (y_temp(end,I)+y_temp(end,I_S))<2e-3*N & (y_temp(2,I)+y_temp(2,I_S)) -(y_temp(1,I)+y_temp(1,I_S)) <0 & reopen_flag ==0
                 release_time = current_time(1);  
                 reopen_flag = 1;
    
            theta = [1;2;3;5;0];
            theta_s = [1;2;3;5;0];
            theta = theta/sum(theta);% normalization
            theta_s = theta/sum(theta_s);
            theta = theta*10;
            theta_s = theta_s*0;
       
            
        kappa = 0.001;
        kappa_s = 0.000;

          end
          
        new_patients = sum(diff(y_temp(:,H)) + diff(y_temp(:,H_S)));
        new_deaths = sum(diff(y_temp(:,D)) + diff(y_temp(:,D_S)));
               
        [t_temp,y_temp] = ode45(@(t,y) diff_system(t,y,theta,kappa,theta_s, kappa_s,new_patients,new_deaths), curr_tspan,y_temp(end,:));
        y_temp = floor(y_temp);
        y(curr_tspan+shift,:) = y_temp;
        Y_econ(curr_tspan+shift) = (sum(y_temp(:,21:24),2) - sum(y_temp(:,25:27),2))./N;%(N - y_temp(:,D) - y_temp(:,D_S));
    
     end
end

y = y(1:end-1,:);
Y_econ = Y_econ(1:end-1,:);


%% poorly implemented

theta = [1;2;3;5;0];
theta_s = [1;2;3;5;0];
kappa = 0.00;
kappa_s = 0.000;
% 
theta = theta/sum(theta);% normalization
theta_s = theta/sum(theta_s);
rate_into_quarantine = 1; % higher number reflects a more strict enforcement of qurantine

theta = theta*0;
theta_s = theta_s*0;


%% solve the system 


% 18  months of simulation
weeks = 76;
t_end = round(weeks*7);

segments = floor(t_end/7);
y_poorly = zeros(t_end,length(y0));
Y_econ_poorly = zeros(t_end,1);

reopen_flag = 0;

for i = 1 : segments

    curr_tspan = (i-1)*8+1:i*8; %current 14 days being evaluated
    if i == 1
        [t_temp,y_temp] = ode45(@(t,y) diff_system(t,y,theta,kappa,theta_s, kappa_s,0,0), curr_tspan ,y0);
        y_temp = floor(y_temp);
        y_poorly(curr_tspan,:) = y_temp;
        
        Y_econ_poorly(curr_tspan) = (sum(y_temp(:,21:24),2) - sum(y_temp(:,25:27),2)) ./N;
    else
        shift = -i +1;
        
    
    current_time = curr_tspan+shift;
     if (y_temp(end,I)+y_temp(end,I_S))<2e-3*N & (y_temp(2,I)+y_temp(2,I_S)) -(y_temp(1,I)+y_temp(1,I_S)) <0 & reopen_flag ==0
                 release_time = current_time(1);  
                 reopen_flag = 1;
    
            theta = [1;2;3;5;0];
            theta_s = [1;2;3;5;0];
            theta = theta/sum(theta);% normalization
            theta_s = theta/sum(theta_s);
            theta = theta*10;
            theta_s = theta_s*10;
       
            
        kappa = 0.001;
        kappa_s = 0.001;
        
     end
        new_patients = sum(diff(y_temp(:,H)) + diff(y_temp(:,H_S)));
        new_deaths = sum(diff(y_temp(:,D)) + diff(y_temp(:,D_S)));
               
        [t_temp,y_temp] = ode45(@(t,y) diff_system(t,y,theta,kappa,theta_s, kappa_s,new_patients,new_deaths), curr_tspan,y_temp(end,:));
        y_temp = floor(y_temp);
        y_poorly(curr_tspan+shift,:) = y_temp;
        Y_econ_poorly(curr_tspan+shift) = (sum(y_temp(:,21:24),2) - sum(y_temp(:,25:27),2))./N;
    
    end
end

y_poorly = y_poorly(1:end-1,:);
Y_econ_poorly = Y_econ_poorly(1:end-1,:);


%% Summary of results at the end of the simulation
y_econ_init = econ_previrus/N;

fprintf("Summary of results at the end of the simulated pedriod:\n\n")
fprintf("Deaths of non-senior population: %d\n",max(y(:,D)))
fprintf("Deaths of senior population: %d\n",max(y(:,D_S)))
fprintf("Maximum hospitalized people at any given time: %d\n",max(y(:,H) + y(:,H_S)))
fprintf("Net productivity change (dY): %3.1f%%\n",(Y_econ(end) -y_econ_init)/y_econ_init*100 )


%% plot for paper 

strain_level = N*(0.08/100);
saturation_level = N*(0.22/100);

lineWidth = 4;
line_dotted_width = 2;

set(0, 'DefaultAxesFontName', 'Arial');
figure('DefaultAxesFontSize',11)
set(gcf,'Color','white');

subplot(2,4,1)
plot(1:t_end,y(:,H) + y(:,I),'k','LineWidth',lineWidth);hold on
ylabel('Non-seniors infected');  
plot(1:t_end,y_poorly(:,H) + y_poorly(:,I),'LineWidth',lineWidth)
xline(release_time,':b','LineWidth',line_dotted_width)


subplot(2,4,2)
plot(1:t_end,y(:,D),'k','LineWidth',lineWidth);hold on
ylabel('Non-senior deaths');
plot(1:t_end,y_poorly(:,D),'LineWidth',lineWidth)
xline(release_time,':b','LineWidth',line_dotted_width)
legend('Extreme restrictions on seniors','Weak restrictions on seniors')

subplot(2,4,3)
plot(1:t_end,y(:,H_S) + y(:,H),'k' ,'LineWidth',lineWidth);hold on
plot(1:t_end,y_poorly(:,H_S) +y_poorly(:,H) ,'LineWidth',lineWidth)
ylabel('Total hospitalizations')
xline(release_time,':b','LineWidth',line_dotted_width)
plot([1, t_end],[1,1]*strain_level,':', 'Color',[0.7,0.7,0.7],'LineWidth',line_dotted_width );
plot([1, t_end],[1,1]*saturation_level,':k','LineWidth',line_dotted_width  );


subplot(2,4,4)

ax = gca;
ax.Visible = 'off';
subplot(2,4,5)
plot(1:t_end,y(:,H_S) + y(:,I_S),'k','LineWidth',lineWidth);hold on
ylabel('Seniors infected');  
plot(1:t_end,y_poorly(:,H_S) + y_poorly(:,I_S),'LineWidth',lineWidth)
xline(release_time,':b','LineWidth',line_dotted_width)


subplot(2,4,6)
plot(1:t_end,y(:,D_S),'k','LineWidth',lineWidth);hold on
ylabel('Senor deaths');
plot(1:t_end,y_poorly(:,D_S),'LineWidth',lineWidth)
xline(release_time,':b','LineWidth',line_dotted_width)

subplot(2,4,7)
%y_econ_init = 2.0425e+04;
y_econ_init = econ_previrus/N;
plot(1:t_end,(Y_econ - y_econ_init)/y_econ_init*100 ,'k','LineWidth',lineWidth);hold on
ylabel('Y [%]')
plot(1:t_end,(Y_econ_poorly - y_econ_init)/y_econ_init*100,'LineWidth',lineWidth)
xline(release_time,':b','LineWidth',line_dotted_width)

subplot(2,4,8)

ax = gca;
ax.Visible = 'off';

for i = 1:8
    subplot(2,4,i)
    ax = gca;
    set(gca,'box','off')
    xlabel("Time [days]")
    ax.FontUnits = 'normalized';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'off';
    ax.XTick = [0 100 200 300 400 500 600];
    ax.TickLength = [0.03 0.1];
    xlim([0 532])
    set(gca, 'FontName', 'Arial')
end

annotation('textbox', [0.08, 0.92, 0.05, 0.06], 'String', "A",'FontName','Arial','FontSize',12,'LineStyle','none')
annotation('textbox', [0.29, 0.92, 0.05, 0.06], 'String', "B",'FontName','Arial','FontSize',12,'LineStyle','none')
annotation('textbox', [0.5, 0.92, 0.05, 0.06], 'String', "C",'FontName','Arial','FontSize',12,'LineStyle','none')
annotation('textbox', [0.08, 0.45, 0.05, 0.06], 'String', "D",'FontName','Arial','FontSize',12,'LineStyle','none')
annotation('textbox', [0.29, 0.45, 0.05, 0.06], 'String', "E",'FontName','Arial','FontSize',12,'LineStyle','none')
annotation('textbox', [0.5, 0.45, 0.05, 0.06], 'String', "F",'FontName','Arial','FontSize',12,'LineStyle','none')
