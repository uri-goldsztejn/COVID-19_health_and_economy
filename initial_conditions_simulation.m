%% This code accompanies the work in "Public policy and economic dynamics of COVID-19 spread: a mathematical modeling study"
% By Uri Goldsztejn, David Schwartzman, and Arye Nehorai. 2020

% If you find this code useful in your research, please consider citing.



%% initialization
clear all;clc ;close all;

%% Constants used in the simulation

% this are enum variables
S = 1; IsoS = 2; E=3; IsoE=4; I=5; H=6; IsoI=7; D=8; R=9; IsoR=10; S_S=11; IsoS_S=12; E_S=13;
IsoE_S=14; I_S=15; H_S=16; IsoI_S = 17; D_S=18; R_S=19; IsoR_S=20;
Y_S_L = 21; Y_S_H = 22; Y_G_L = 23; Y_G_H = 24; C_H = 25; C_L = 26; C_D = 27;



N = 330e6; %total population
y0 = zeros(27,1); %This variable will store the initial values

economy_init = 4.0803e12;
%% Initial seed
% We start our simulation with only a few infected individuals.
% Then we use the values in y to define the intial conditions for our
% scenarios.
% We used the values of day 35 as initial conditions for our scenarios.

y0(S) = 0.67; y0(I)= 1e-07;

y0(S_S) = 0.33;

y0 = y0/sum(y0)*N; % percent of the population

% Initial conditions for the economic model
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

% Conditional on being in the labor force:
% Average income for S_L  37950.38
% Average income for S_H  76255.1
% Average income for G_L  33919.6
% Average income for G_H  79094.35
INC_S_L = 37950.38;
INC_S_H = 76255.1;
INC_G_L = 33919.6;
INC_G_H = 79094.35;

% Average earnings for member of each group are POP_n * LAB_n * INC_n

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

y0_econ/econ_previrus;


%% Define transition rates into and out of quarantine

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
% the system is solved for 7 days at a time. 

% About 18  months of simulation
weeks = 76; % at the end the simulation will be shorter due to the overlaps
t_end = round(weeks*7);


segments = floor(t_end/7); % each segment corresponds to a week in the simulation
y = zeros(t_end,length(y0));
Y_econ = zeros(t_end,1);


% The last values of each segment are used as initial conditions for the
% next segment. Therefore, we calculate 8 days at a time and the override
% the last day.

for i = 1 : segments

    curr_tspan = (i-1)*8+1:i*8; %current 8 days being evaluated
    if i == 1
        [t_temp,y_temp] = ode45(@(t,y) diff_system(t,y,theta,kappa,theta_s, kappa_s,0,0), curr_tspan ,y0);
        y_temp = floor(y_temp);
        y(curr_tspan,:) = y_temp;
        
        Y_econ(curr_tspan) = (sum(y_temp(:,21:24),2) - sum(y_temp(:,25:27),2)) ./N;
    else
        shift = -i +1; % Shift in time due to overriding the last day in each segment.
        
        % we keep track of the new patients and deaths
        new_patients = sum(diff(y_temp(:,H)) + diff(y_temp(:,H_S)));
        new_deaths = sum(diff(y_temp(:,D)) + diff(y_temp(:,D_S)));
               
        % Here we solve the system of equations
        [t_temp,y_temp] = ode45(@(t,y) diff_system(t,y,theta,kappa,theta_s, kappa_s,new_patients,new_deaths), curr_tspan,y_temp(end,:));
        y_temp = floor(y_temp);
        y(curr_tspan+shift,:) = y_temp;
        Y_econ(curr_tspan+shift) = (sum(y_temp(:,21:24),2) - sum(y_temp(:,25:27),2))./N;
   
    end
end

y = y(1:end-1,:);
Y_econ = Y_econ(1:end-1,:);

%% plot the number of individuals in each compartment (Supplementary fig 1)
lineWidth = 5;

figure
set(gcf,'Color','white');

subplot(4,4,1)
plot(1:t_end,y(:,S),'k','LineWidth',lineWidth)
ylabel('Suceptible');  

subplot(4,4,2)
plot(1:t_end,y(:,E),'k','LineWidth',lineWidth)
ylabel('Exposed');  

subplot(4,4,3)
plot(1:t_end,y(:,H),'LineWidth',lineWidth);hold on
ylabel('Infected');  
plot(1:t_end,y(:,I),'LineWidth',lineWidth)
ylabel('Infected');  
legend('Hospitalized','Not hospitalized','FontSize',7) 

subplot(4,4,4)
plot(1:t_end,y(:,R),'k','LineWidth',lineWidth);hold on;
ylabel('Recovered');
yyaxis right
plot(1:t_end,y(:,D),'LineWidth',lineWidth)
ylabel('Dead');

subplot(4,4,5)
plot(1:t_end,y(:,IsoS),'k','LineWidth',lineWidth)
ylabel('Isolated suceptible');  

subplot(4,4,6)
plot(1:t_end,y(:,IsoE),'k','LineWidth',lineWidth)
ylabel('Isolated exposed');  

subplot(4,4,7)
plot(1:t_end,y(:,IsoI),'k','LineWidth',lineWidth)
ylabel('Isolated infected');  

subplot(4,4,8)
plot(1:t_end,y(:,IsoR),'k','LineWidth',lineWidth)
ylabel('Isolated recovered');  

subplot(4,4,9)
plot(1:t_end,y(:,S_S),'k','LineWidth',lineWidth)
ylabel('Suceptible_S');  

subplot(4,4,10)
 
plot(1:t_end,y(:,E_S),'k','LineWidth',lineWidth)
ylabel('Exposed_S');  

subplot(4,4,11)
plot(1:t_end,y(:,H_S),'LineWidth',lineWidth);hold on
ylabel('Infected_S');  
plot(1:t_end,y(:,I_S),'LineWidth',lineWidth)
ylabel('Infected_S');  
legend('Hospitalized','Not hospitalized','FontSize',7) 

subplot(4,4,12)
plot(1:t_end,y(:,R_S),'k','LineWidth',lineWidth);hold on;
ylabel('Recovered_S');
yyaxis right
plot(1:t_end,y(:,D_S),'LineWidth',lineWidth)
ylabel('Dead_S');

subplot(4,4,13)
plot(1:t_end,y(:,IsoS_S),'k','LineWidth',lineWidth)
ylabel('Isolated suceptible_S');  

subplot(4,4,14)
plot(1:t_end,y(:,IsoE_S),'k','LineWidth',lineWidth)
ylabel('Isolated exposed_S');  

subplot(4,4,15)
plot(1:t_end,y(:,IsoI_S),'k','LineWidth',lineWidth)
ylabel('Isolated infected_S');  

subplot(4,4,16)
plot(1:t_end,y(:,IsoR_S),'k','LineWidth',lineWidth)
ylabel('Isolated recovered_S');  

for i = 1:16
    subplot(4,4,i)
    ax = gca;
    set(gca,'box','off')
    xlabel("Time [days]")
    ax.FontUnits = 'normalized';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'off';
    ax.XTick = [0 250 500 ];
    ax.TickLength = [0.05 0.1];
     if i ==1 || i==3 ||i == 5
 ytickformat('%1.0f')
     else
        ytickformat('%1.0f')
     end
end


annotation('textbox', [0.08, 0.9, 0.05, 0.06], 'String', "A",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.28, 0.9, 0.05, 0.06], 'String', "B",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.49, 0.9, 0.05, 0.06], 'String', "C",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.7, 0.9, 0.05, 0.06], 'String', "D",'FontName','Arial','FontSize',13,'LineStyle','none')

annotation('textbox', [0.08, 0.68, 0.05, 0.06], 'String', "E",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.28, 0.68, 0.05, 0.06], 'String', "F",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.49, 0.68, 0.05, 0.06], 'String', "G",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.7, 0.68, 0.05, 0.06], 'String', "H",'FontName','Arial','FontSize',13,'LineStyle','none')

annotation('textbox', [0.08, 0.45, 0.05, 0.06], 'String', "I",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.28, 0.45, 0.05, 0.06], 'String', "J",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.49, 0.45, 0.05, 0.06], 'String', "K",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.7, 0.45, 0.05, 0.06], 'String', "L",'FontName','Arial','FontSize',13,'LineStyle','none')

annotation('textbox', [0.08, 0.25, 0.05, 0.06], 'String', "M",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.28, 0.25, 0.05, 0.06], 'String', "N",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.49, 0.25, 0.05, 0.06], 'String', "O",'FontName','Arial','FontSize',13,'LineStyle','none')
annotation('textbox', [0.7, 0.25, 0.05, 0.06], 'String', "P",'FontName','Arial','FontSize',13,'LineStyle','none')
   
%% plot the population evolution in the baseline scenario (fig 2 in the manuscript)

strain_level = N*(0.08/100);
saturation_level = N*(0.22/100);

lineWidth = 4;
line_dotted_width = 2;
% figure
set(0, 'DefaultAxesFontName', 'Arial');
figure('DefaultAxesFontSize',11)
set(gcf,'Color','white');

subplot(2,3,1)
plot(1:t_end,y(:,I),'k','LineWidth',lineWidth)
ylabel('Non-seniors infected');  


subplot(2,3,2)
plot(1:t_end,y(:,D),'k','LineWidth',lineWidth)
ylabel('Non-senior deaths');

subplot(2,3,3)
plot(1:t_end,y(:,H_S) +y(:,H),'k' ,'LineWidth',lineWidth);hold on
plot([1, t_end],[1,1]*strain_level,':', 'Color',[0.7,0.7,0.7],'LineWidth',line_dotted_width );
plot([1, t_end],[1,1]*saturation_level,':k','LineWidth',line_dotted_width  );
ylabel('Total hospitalizations')


subplot(2,3,4)
plot(1:t_end,y(:,I_S),'k','LineWidth',lineWidth)

ylabel('Seniors infected ');  


subplot(2,3,5)
plot(1:t_end,y(:,D_S),'k','LineWidth',lineWidth)
ylabel('Senior deaths');

% lineWidth = 4;

y_econ_init = 2.0425e+04;
subplot(2,3,6)
plot(1:t_end,(Y_econ - y_econ_init)/y_econ_init *100,'k','LineWidth',lineWidth);
ylabel('Y [%]')


for i = 1:6
    subplot(2,3,i)
    ax = gca;
    set(gca,'box','off')
    xlabel("Time [days]")
    ax.FontUnits = 'normalized';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'off';
    ax.XTick = [0 100 200 300 400 500 600 700 800 900 1000];
    ax.TickLength = [0.03 0.1];
    set(gca, 'FontName', 'Arial')
end
annotation('textbox', [0.06, 0.92, 0.05, 0.06], 'String', "A",'FontName','Arial','FontSize',12,'LineStyle','none')
annotation('textbox', [0.35, 0.92, 0.05, 0.06], 'String', "B",'FontName','Arial','FontSize',12,'LineStyle','none')
annotation('textbox', [0.64, 0.92, 0.05, 0.06], 'String', "C",'FontName','Arial','FontSize',12,'LineStyle','none')
annotation('textbox', [0.06, 0.45, 0.05, 0.06], 'String', "D",'FontName','Arial','FontSize',12,'LineStyle','none')
annotation('textbox', [0.35, 0.45, 0.05, 0.06], 'String', "E",'FontName','Arial','FontSize',12,'LineStyle','none')
annotation('textbox', [0.64, 0.45, 0.05, 0.06], 'String', "F",'FontName','Arial','FontSize',12,'LineStyle','none')

