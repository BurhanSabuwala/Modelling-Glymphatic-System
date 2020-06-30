% crk
% codes for the model that studies the combination of basic fluid mechanics
% and role of apoe/lrp1 on beta amyloid levels in the brain. the isf is
% modeled as well as the perivascular spaces.

clear all

% the model here describes a mm^3 space of brain within the ca1 region of
% the hippocampus. all cell numbers have been dervied from the literature

% simulation time step = 24 hrs x 365 days x 50 years
t = 24*365*50;
x = zeros(t,1);

% Studying the effect of Caffeine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cell # definition
sleep_hours = 7; % Average number of hours the person sleeps
N = zeros(t,1);
M = zeros(t,1);
EC = zeros(t,1);

N(1,1) = 5275;  % neurons
M(1,1) = 4500;  % microglia
EC(1,1) = 22050; %brain endothelial cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition of the matrices for the chemical species

Ccaf = zeros(t,1);  % conc. of Caffeine in serum
C1 = zeros(t,1);    % conc. of Ab40 in brain isf
C2 = zeros(t,1);    % conc of Ab40 in perivascular space
C3 = zeros(t,1);    % conc. of Ab40 deposited in vessels as CAA
C9 = zeros(t,1);

C4 = zeros(t,1);    % conc. of Ab42 in brain isf
C5 = zeros(t,1);    % conc. of Ab42 in PVS
C6 = zeros(t,1);    % conc. of Ab42 deposited in vessels as CAA
C10 = zeros(t,1);

LRPEC = zeros(t,1); % LRP density on brain ECs
ApoE = zeros(t,1); % ApoE conc in ISF/PVS (assumption)

Fi40 = zeros(t,1);    % Flow of Ab40 between ISF and PVS
Fi42 = zeros(t,1);    % Flow of Ab42 between ISF and PVS

%stiffness = zeros(t,1);

% initializing conditions
HR = 60;        % heart rate per minute
HR = 72;        % heart rate per minute
%for n = 1:t
%     stiffness = 1*exp(3.8*(n-1));  % relative vessel stiffness, increases with aging
%end
stiffness = 1;
velocity = 0.14*(HR/60)*(1/stiffness)/24;      % mm per day, velocity of bulk ISF flow (from Abbott 2004)
v_factor = 0.1;            % Average awake/sleep factor
Ccaf(1,1) = 0;

C1(1,1) = 6.26e7;          % Ab40 in brain
C2(1,1) = 6.26e6;          % Ab40 in PVS
C3(1,1) = 0;               % Ab40 deposited in vessel
C9(1,1) = 0;

C4(1,1) = 9.74e6;          % Ab42 in brain
C5(1,1) = 9.74e5;          % Ab42 in PVS
C6(1,1) = 0;               % Ab42 deposited in vessel
C10(1,1) = 0;

LRPEC(1,1) = 9.3e11;       % # of LRP1 receptors per endothelial cell
ApoE(1,1) = 1;             % not used in this modeling yet

Fi40(1,1) = velocity*C1(1,1);           % mass flux out of brain of Ab40 (#,hr)
Fi42(1,1) = velocity*C4(1,1);           % mass flux out of brain of Ab42 (#,hr)

% rate constants

d = 0.000027/24;   % neuronal death rate per hour assuming non-Apoe4 carrier
k1 = 0.000093/24;  % additional neuronal death attributed to elevated Ab
%k1 = 0.000050/24;
k2 = 0.0000027/24; % conversion of ECs to senesnce (basically death rate) without additional injury to ECs

k3 = (1.6e4)/24;     % generation rate of Ab40 per neuron per hr
k3_factor = 10;     % kn_factor = kn_wake/kn_sleep

k8 = (3e3)/24;       % generation rate of Ab42 per neuron per hr
k8_factor = 10;     % kn_factor = kn_wake/kn_sleep

kd_caf = log(2)/5.7;  % Half life of Caffeiene is 5.7 - 6 hours - Institute of Medicine (US) Committee on Military Nutrition Research Study
kd_pgp = log(2)/26.7; % Half life of P-Glycoprotein is 26.7 hours  -  Katayama et al.

k4 = 1.5;         % increased production factor for Ab gen when stimulated

k5 = (1.4e3)/24;    % uptake rate of Ab per hr per microglial cell
k11 = (2e2)/24;

k6 = 0.4/24;       % per hr; fraction of Ab removed by LRP1 from PVS
%k6 = 0.12;      % per hr, fraction of Ab removed by LRP1 from PVS when ApoE4 allele present (Bachmeir 2012)

%k7 = 0.004/24;     % deposition rate of Ab40 per hr in brain vessel
k7 = 0.006/24;     % deposition rate of Ab40 per hr in brain vessel
k9 = 0.001/24;     % deposition rate of Ab42 per hr in brain isf
k7a= 0.0015/24;     % deposition rate of Ab40 per hr in brain isf when elevated Ab
k7a2 = 0.001/24;
%k9a = 0.002/24;    % deposition rate of Ab42 per hr in brain vessel
k9a = 0.003/24;    % deposition rate of Ab42 per hr in brain vessel

k10 = (2.7e-5)/24;   % normal loss of LRP1 with aging, fraction/hr

k12 = 0.000073/24; % loss of effective microglia over time (become about 80% dystrophic by 60 yo)

kpgp_increment = 1.45; % Increase in Flow rates  From Qosa et al results for Mouse

% equations

m = length(C1);

for n = 1:m
    
    if rem(n,24)<16
        velocity_wake = velocity*24*v_factor/((24-sleep_hours)*v_factor +sleep_hours);
        velocity_now = velocity_wake;
        if n == 1:175200 
            Fi40(n+1,1) = 2*velocity_now*C1(n,1);           % mass flux out of brain of Ab40 (#,day)
        else
            Fi40(n+1,1) = velocity_now*C1(n,1);           % mass flux out of brain of Ab40 (#,day)
        end
        Fi42(n+1,1) = velocity_now*C4(n,1);      % mass flow out of brain ISF to PVS
    else
        velocity_sleep = velocity*24/((24-sleep_hours)*v_factor+sleep_hours);
        velocity_now = velocity_sleep;
        Fi42(n+1,1) = velocity_now*C4(n,1);      % mass flow out of brain ISF to PVS
        if n == 1:175200
            Fi40(n+1,1) = 2*velocity_now*C1(n,1);           % mass flux out of brain of Ab40 (#,day)
        else
            Fi40(n+1,1) = velocity_now*C1(n,1);           % mass flux out of brain of Ab40 (#,day)
        end
    end
    
    if Ccaf(n,1) >= 3                                    % Arbitary Cutoff
        Fi40(n,1) = kpgp_increment * Fi40(n,1);
    end
    
    x(n+1,1) = n;
    if C4(n,1)>= 5*C4(1,1) && N(n,1)>=0
        N(n+1,1) = (1-(k1+d))*N(n,1);
    elseif N(n,1)>=0
        N(n+1,1) = (1-d)*N(n,1);
    else
        N(n+1,1) = 0;
    end
    
    M(n+1,1) = M(n,1)*(1-k12);
    
    if EC(n,1) >=0
        EC(n+1,1) = (1-k2)*EC(n,1);         % change in endothelial cells without addition of loss 2/2 Ab
    else
        EC(n+1,1) = 0;
    end
    
    ratio = EC(n,1)/EC(1,1);            % fraction of ECs remaining
    LRPEC(n+1,1) = (1-k10)*LRPEC(n,1); %*ratio;   %(EC(n,1)./(EC(1,1)));  % change in LRP1 # wrt aging
        
    %Caffeine dyanmics
    if Ccaf(n,1)>0 
        Ccaf(n+1,1) =  Ccaf(n,1) - kd_caf*Ccaf(n,1);  % Caffeine degradation
    else
        Ccaf(n,1) = 0;
    end
    if rand()>0.9 && (rem(n,24)==0 | rem(n,24)==8)
        Ccaf(n+1,1) = Ccaf(n+1,1) + 10;   % Coffee Intake
    end    
    
    % brain isf ab40
    if C1(n,1) >=0 && C1(n,1)>=5*C1(1,1)
        if rem(n,24)<16
            k3_wake = 24*k3_factor*k3/((24-sleep_hours)*k3_factor + sleep_hours); % generation rate of Ab40 per neuron per hour during wake
            C1(n+1,1) = C1(n,1) + k4*k3_wake*N(n,1) - k5*M(n,1) - Fi40(n,1) - k7a2*C1(n,1);
            C9(n+1,1) = C9(n,1) + k7a2*C1(n,1);     %Ab40 deposition in brain
        else
            k3_sleep = 24*k3/((24-sleep_hours)*k3_factor + sleep_hours);       % generation rate of Ab40 per neuron per hour during sleep
            C1(n+1,1) = C1(n,1) + k4*k3_sleep*N(n,1) - k5*M(n,1) - Fi40(n,1) - k7a2*C1(n,1);
            C9(n+1,1) = C9(n,1) + k7a2*C1(n,1);     %Ab40 deposition in brain
        end 
    elseif C1(n,1) >=0
        if rem(n,24)<16
            k3_wake = 24*k3_factor*k3/((24-sleep_hours)*k3_factor+sleep_hours);       % generation rate of Ab40 per cell per hour
            C1(n+1,1) = C1(n,1) + k3_wake*N(n,1) - k5*M(n,1) - Fi40(n,1) - k7a2*C1(n,1);
            C9(n+1,1) = C9(n,1) + k7a2*C1(n,1);     %Ab40 deposition in brain
        else
            k3_sleep = 24*k3/((24-sleep_hours)*k3_factor+sleep_hours);       % generation rate of Ab40 per cell per hour
            C1(n+1,1) = C1(n,1) + k3_sleep*N(n,1) - k5*M(n,1) - Fi40(n,1) - k7a2*C1(n,1);
            C9(n+1,1) = C9(n,1) + k7a2*C1(n,1);     %Ab40 deposition in brain
        end
    else
        C1(n+1,1) = 0;
    end
    
    if C2(n,1) >=0
        C2(n+1,1) = C2(n,1) + Fi40(n,1) - velocity_now.*C2(n,1) - k6.*((LRPEC(n,1))/LRPEC(1,1)).*C2(n,1) - k7*C2(n,1);
    else
        C2(n+1,1) = 0;
    end
    
    if C3(n,1) >=0
        C3(n+1,1) = C3(n,1) + k7*C2(n,1);
    else     
        C3(n+1,1) = 0;
    end
     
    Fi42(n+1,1) = velocity_now*C4(n,1);      % mass flow out of brain ISF to PVS
    
    if Ccaf(n,1) >= 2
        Fi42(n,1) = kpgp_increment * Fi42(n,1);
    end
    
    if C4(n,1) >=0 && C4(n,1)>=5*C4(1,1)
        if rem(n,24)<16
            k8_wake = 24*k8_factor*k8/((24-sleep_hours)*k8_factor+sleep_hours);
            C4(n+1,1) = C4(n,1) + k4*k8_wake*1.01*N(n,1) - k11*M(n,1) - Fi42(n,1);
        else
            k8_night = 24*k8/((24-sleep_hours)*k8_factor+sleep_hours);
            C4(n+1,1) = C4(n,1) + k4*k8_night*1.01*N(n,1) - k11*M(n,1) - Fi42(n,1);
        end
    elseif C4(n,1) >=0
        if rem(n,24)<16
            k8_wake = 24*k8_factor*k3/((24-sleep_hours)*k8_factor+sleep_hours);
            C4(n+1,1) = C4(n,1) + k8_wake*1.01*N(n,1) - k11*M(n,1) - Fi42(n,1);
        else
            k8_night = 24*k8/((24-sleep_hours)*k8_factor+sleep_hours);
            C4(n+1,1) = C4(n,1) + k8_night*1.01*N(n,1) - k11*M(n,1) - Fi42(n,1);
        end
    else
        C4(n+1,1) = 0;
    end
    
     C10(n+1,1) = C10(n,1) + k9*C4(n,1);     %Ab42 deposition in brain
    
    if C5(n,1) >=0 
        C5(n+1,1) = C5(n,1) + Fi42(n,1) - velocity_now.*C5(n,1) - 0.5*k6.*((LRPEC(n,1))/LRPEC(1,1)).*C5(n,1) - k9*C5(n,1);
    else
        C5(n+1,1) = 0;
    end
    
    if C6(n,1) >=0
        C6(n+1,1) = C6(n,1) + k9a*C5(n,1);     
    else
        C6(n+1,1) = 0;
    end
    
end

figure(1)
plot(x, C9, 'k', x, C3, 'g', x, C10, 'b', x, C6, 'm', 'Linewidth', 2);
legend('Ab40_I_S_F', 'Ab40_V_e_s_s_e_l', 'Ab42_I_S_F', 'Ab42_V_e_s_s_e_l')
xlabel('Days'), ylabel('Ab #')
hold on

figure(2)
plot(x, C1, 'k', x, C2, 'g', x, C4, 'b', x, C5, 'm', 'Linewidth', 2);
legend('Ab40_b_r_a_i_n', 'Ab40_I_S_F','Ab42_b_r_a_i_n', 'Ab42_I_S_F')
xlabel('Hours'), ylabel('Ab #')
hold on
% subplot(2,2,1), plot(C9(1:t), 'r'), xlabel('days'), ylabel('[Ab40]'), title('Ab40 deposition in brain ISF'), hold on
% subplot(2,2,3), plot(C3(1:t,1), 'g'), xlabel('days'), ylabel('[Ab40]'), title('Ab40 deposition in cerebral vasculature'), hold on
% subplot(2,2,2), plot(C6(1:t,1), 'm'), xlabel('days'), ylabel('[Ab42]'), title('Ab42 deposition in cerebral vasculature'), hold on
% subplot(2,2,4), plot(C10(1:t), 'b'), xlabel('days'), ylabel('[Ab42]'), title('Ab42 deposition in brain ISF'), hold on
% % subplot(3,2,4), plot(C5(1:t), 'b'), xlabel('days'), ylabel('[Ab42]'), title('Ab42 in perivascular space'), hold on
%subplot(3,2,6), plot(C6(1:t), 'r'), xlabel('days'), ylabel('[Ab42]'), title('Ab42 deposition in cerebral vasculature'), hold on

figure(3)
plot(x,N, 'r', x, M, 'g', x, EC, 'b', 'Linewidth', 1);
legend('Neuron #', 'Microglia #', 'EC #')
xlabel('Hours'), ylabel('Cell #')
hold on

figure(4)
plot(x,C1./C4, 'r', x, C2./C5, 'g', x, C3./C6, 'b', 'Linewidth', 1);
legend('Brain Parenchyma Ratio', ' PVS ratio', 'blood vessel ratio')
xlabel('Hours'), ylabel('Cell #')
hold on

% subplot(1,4,1), plot(N(1:t), 'r'), xlabel('days'), ylabel('Neuron #'), hold on
% subplot(1,4,2), plot(M(1:t), 'k'), xlabel('days'), ylabel('Microglia #'), hold on
% subplot(1,4,3), plot(EC(1:t), 'b'), xlabel('days'), ylabel('Brain endothelial cell #'), hold on
% subplot(1,4,4), plot(LRPEC(1:t), 'g'), xlabel('days'), ylabel('LRP #'), hold on

figure (5)
plot(x, LRPEC, 'Linewidth', 2)
xlabel('Hours'), ylabel('Total LRP #')
hold on

% figure(4)
% plot(x,C1, '+')
% hold on

%output final values of important chemical species
% C9(t,1) % Ab40 brain
% C3(t,1) %Ab40 vessel
% C10(t,1) %Ab42 brain
% C6(t,1) %Ab42 vessel

% LRPEC(t,1)
% EC(t,1)
