function [tsleep_vec,C,x] = sleep_irregular(t,sleep_reg)
% codes for the model that studies the combination of basic fluid mechanics
% and role of apoe/lrp1 on beta amyloid levels in the brain. the isf is
% modeled as well as the perivascular spaces.

% clear all

% the model here describes a mm^3 space of brain within the ca1 region of
% the hippocampus. all cell numbers have been dervied from the literature

% simulation time step = 24 hrs x 365 days x 50 years
% t = 365*24*50;
% t = 24*365;
x = zeros(t,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cell # definition

N = zeros(t,1);
M = zeros(t,1);
EC = zeros(t,1);

N(1,1) = 5275;  % neurons
M(1,1) = 4500;  % microglia
EC(1,1) = 22050; %brain endothelial cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% definition of the matrices for the chemical species

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

%% initializing conditions
HR = 60;        % heart rate per minute
HR = 72;        % heart rate per minute
%for n = 1:t
%     stiffness = 1*exp(3.8*(n-1));  % relative vessel stiffness, increases with aging
%end
stiffness = 1;
velocity = 0.14*(HR/60)*(1/stiffness)/24;      % mm per day, velocity of bulk ISF flow (from Abbott 2004)
v_factor = 0.01;           % V_wake/V_sleep
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


%% rate constants

d = 0.000027/24;   % neuronal death rate per hour assuming non-Apoe4 carrier
k1 = 0.000093/24;  % additional neuronal death attributed to elevated Ab
%k1 = 0.000050;
k2 = 0.0000027/24; % conversion of ECs to senesnce (basically death rate) without additional injury to ECs

k3 = 1.6e4/24;     % generation rate of Ab40 per neuron per hr
k3_factor = 10;     % kn_factor = kn_wake/kn_sleep

k8 = 3e3/24;       % generation rate of Ab42 per neuron per hr
k8_factor = 10;    % kn_factor = kn_wake/kn_sleep
k4 = 1.5;         % increased production factor for Ab gen when stimulated

k5 = 1.4e3/24;    % uptake rate of Ab per hr per microglial cell
k11 = 2e2/24;

k6 = 0.4/24;       % per hr; fraction of Ab removed by LRP1 from PVS
%k6 = 0.12;      % per hr, fraction of Ab removed by LRP1 from PVS when ApoE4 allele present (Bachmeir 2012)

k7 = 0.004/24;     % deposition rate of Ab40 per hr in brain vessel
k9 = 0.001/24;     % deposition rate of Ab42 per hr in brain isf
k7a= 0.001/24;     % deposition rate of Ab40 per hr in brain isf when elevated Ab
k7a2 = 0.001/24;
k9a = 0.002/24;    % deposition rate of Ab42 per hr in brain vessel

k10 = (2.7e-5)/24;   % normal loss of LRP1 with aging, fraction/hr

k12 = 0.000073/24; % loss of effective microglia over time (become about 80% dystrophic by 60 yo)

%% equations

m = length(C1);
tsleep_reg = 8;
tsleep_act = 8;
tsleep_vec = [tsleep_act];
%%
for n = 1:m
    if sleep_reg == 0
        if rem(n,24) == 0        
            tsleep_act = randi(8);
            tsleep_vec(n) = tsleep_act;
        else
            if n>1
                tsleep_act = tsleep_vec(n-1);
                tsleep_vec(n) = tsleep_act;
            end
        end
    else
            tsleep_act = tsleep_reg;
    end
    if rem(n,24)<(24-tsleep_act)
        velocity_wake = velocity*24*v_factor/((24-tsleep_reg)*v_factor +tsleep_reg);
        velocity_now = velocity_wake;
        if n == 1:175200
            Fi40(n+1,1) = 2*velocity_now*C1(n,1);           % mass flux out of brain of Ab40 (#,day)
        else
            Fi40(n+1,1) = velocity_now*C1(n,1);           % mass flux out of brain of Ab40 (#,day)
        end

    else
        velocity_sleep = velocity*24/((24-tsleep_reg)*v_factor +tsleep_reg);
        velocity_now = velocity_sleep;
        Fi42(n+1,1) = velocity_now*C4(n,1);      % mass flow out of brain ISF to PVS
        if n == 1:175200
            Fi40(n+1,1) = 2*velocity_now*C1(n,1);           % mass flux out of brain of Ab40 (#,day)
        else
            Fi40(n+1,1) = velocity_now*C1(n,1);           % mass flux out of brain of Ab40 (#,day)
        end
        Fi42(n+1,1) = velocity_now*C4(n,1);      % mass flow out of brain ISF to PVS
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
        
    % brain isf ab40
    if C1(n,1) >=0 && C1(n,1)>=5*C1(1,1)
        if rem(n,24)<(24-tsleep_act)
            k3_wake = 24*k3_factor*k3/((24-tsleep_reg)*k3_factor + tsleep_reg); % generation rate of Ab40 per neuron per hour during wake
            C1(n+1,1) = C1(n,1) + k4*k3_wake*N(n,1) - k5*M(n,1) - Fi40(n,1) - k7a2*C1(n,1);
            C9(n+1,1) = C9(n,1) + k7a2*C1(n,1);     %Ab40 deposition in brain
        else
            k3_sleep = 24*k3/((24-tsleep_reg)*k3_factor + tsleep_reg);       % generation rate of Ab40 per neuron per hour during sleep
            C1(n+1,1) = C1(n,1) + k4*k3_sleep*N(n,1) - k5*M(n,1) - Fi40(n,1) - k7a2*C1(n,1);
            C9(n+1,1) = C9(n,1) + k7a2*C1(n,1);     %Ab40 deposition in brain
        end 
    elseif C1(n,1) >=0
        if rem(n,24)<(24-tsleep_act)
            k3_wake = 24*k3_factor*k3/((24-tsleep_reg)*k3_factor + tsleep_reg);     % generation rate of Ab40 per cell per hour
            C1(n+1,1) = C1(n,1) + k3_wake*N(n,1) - k5*M(n,1) - Fi40(n,1) - k7a2*C1(n,1);
            C9(n+1,1) = C9(n,1) + k7a2*C1(n,1);     %Ab40 deposition in brain
        else
            k3_sleep = 24*k3/((24-tsleep_reg)*k3_factor + tsleep_reg);      % generation rate of Ab40 per cell per hour
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
    
    if C4(n,1) >=0 && C4(n,1)>=5*C4(1,1)
        if rem(n,24)<(24-tsleep_act)
             k8_wake = 24*k8_factor*k8/((24-tsleep_reg)*k8_factor + tsleep_reg); % generation rate of Ab42 per neuron per hour during wake
            C4(n+1,1) = C4(n,1) + k4*k8_wake*N(n,1) - k11*M(n,1) - Fi42(n,1);
            C10(n+1,1) = C10(n,1) + k9*C4(n,1);     %Ab42 deposition in brain
        else
            k8_sleep = 24*k8/((24-tsleep_reg)*k8_factor + tsleep_reg);       % generation rate of Ab40 per neuron per hour during sleep
            C4(n+1,1) = C4(n,1) + k4*k8_sleep*N(n,1) - k11*M(n,1) - Fi42(n,1);
            C10(n+1,1) = C10(n,1) + k9*C4(n,1);     %Ab42 deposition in brain
        end 
    elseif C4(n,1) >=0
        if rem(n,24)<(24-tsleep_act)
             k8_wake = 24*k8_factor*k8/((24-tsleep_reg)*k8_factor + tsleep_reg);       % generation rate of Ab42 per cell per hour
            C4(n+1,1) = C4(n,1) + k8_wake*N(n,1) - k11*M(n,1) - Fi42(n,1);
            C10(n+1,1) = C10(n,1) + k9*C4(n,1);     %Ab42 deposition in brain
        else
            k8_sleep = 24*k8/((24-tsleep_reg)*k8_factor + tsleep_reg);       % generation rate of Ab42 per cell per hour
            C4(n+1,1) = C4(n,1) + k8_sleep*N(n,1) - k11*M(n,1) - Fi42(n,1);
            C10(n+1,1) = C10(n,1) + k9*C10(n,1);     %Ab42 deposition in brain
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
C = [C1 C2 C3 C4 C5 C6 C9 C10];
end
