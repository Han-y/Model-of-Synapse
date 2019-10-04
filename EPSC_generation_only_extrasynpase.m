 hold off;
%%% Parameters can be changed
glu_No = 10000;
N_ampar_PerCluster=100 ; %number of total AMPARs
runs = 80;

period = 12; % P4, P8, P12, P18, P30
ratio = 0.45; % GluA2 / (GluA2 + GluA4)
              % P4 : 0.75;  P8 : 0.7; P12: 0.45
              % P18: 0.35;  P30: 0.25
SpillOver = 0;
NND_No = 0; % number of synapse @NND
      
%%%% Calculation of mEPSC %%%%
timestepSize = 10 * 2000;% 10 ms * 2000 timeStepSize per microsecond
baseline = 2 * 2000; % 2 ms
dt = 0.5; % 0.5us, time step
range = 1500; % the range of fitting, 
              % it should be changed sometimes, because the fitting
              % equation can't be got when the curve larger than zero

%%%% membrane attribute
Vm = -70; %Resting membrane potential in mV, Gupta R, Reneaux M, Karmeshu, 2016
Cm=0.9;   %Membrane capacitance in nanofarad, Gupta R, Reneaux M, Karmeshu, 2016
gm=0.025; %Membrane leak conductance in nanoSieman, Gupta R, Reneaux M, Karmeshu, 2016


%%%% AMPAR attribute
Pr = 1; % release probability
Vampa = 0; % reversal potential

V = zeros(timestepSize,1);
V(1) = -70; %mV, initial volt

V_a2 = zeros(timestepSize,1);
V_a2(1) = -70; %mV, initial volt

V_a4 = zeros(timestepSize,1);
V_a4(1) = -70; %mV, initial volt

open_state_total_a2 = zeros(44000,10);
open_state_total_a4 = zeros(44000,10);
open_state_total_a2_SpillOver = zeros(44000,10);
open_state_total_a4_SpillOver = zeros(44000,10);
I_total_a2 = zeros(20000,runs);% Store the traces of GluA1 in each loop
I_total_a4 = zeros(20000,runs);% Store the traces of GluA4 in each loop
for i = 1: runs
    %%%%  slow GluAs  %%%%
    if ratio == 0
        I_total_a2 = zeros(20000,runs);
        continue
    else
        name = ['GluA2_P',num2str(period),'_',num2str(i) '.mat'];
        load(name);
        open_state_total_a2(:,i) = amparStates(:,4);
        open_state_a2 = open_state_total_a2(:,i);

        g_a2 = 0.021;% nanoSieman, Robert A. and Howe J., 2003
        Po_a2 = open_state_a2./(N_ampar_PerCluster*ratio); % open probability on PSD
        
        if SpillOver == 1
            name_SpillOver = ['GluA2_P',num2str(period),'_',num2str(i) '_SpillOver.mat'];
            load(name_SpillOver);
            open_state_total_a2_SpillOver(:,i) = amparStates(:,4);
            open_state_a2_SpillOver = open_state_total_a2_SpillOver(:,i);
            Po_a2 = Po_a2 + NND_No * open_state_a2_SpillOver./(N_ampar_PerCluster*ratio); % open probability on PSD
            open_state_a2 = open_state_a2 + NND_No * open_state_a2_SpillOver;
        end

        for k = 1:(timestepSize-1)
            V_a2(k+1)=V_a2(k)+dt*(1/Cm)*((-1*gm*(V_a2(k)-Vm))+(-1.*g_a2.*Po_a2(k)*(V_a2(k)-Vampa))); %mV, Gupta R, Reneaux M, Karmeshu, 2016
        end
          
        I_2_a2 = zeros(timestepSize ,1);
        for k = 1:timestepSize
            I_2_a2(k) = g_a2 .* open_state_a2(k) * (V_a2(k)-Vampa);%Ampa excitatory current in picoAmpere.
        end            
    
        I_total_a2(:,i) = I_2_a2; 
    end
 

    %%%%  fast GluAs  %%%%
    if ratio ==1
        I_total_a4 = zeros(20000,runs);
        continue
    else
        name = ['GluA4_P',num2str(period),'_',num2str(i+80), '.mat'];
        load(name);
        open_state_total_a4(:,i) = amparStates(:,4);
        open_state_a4 = open_state_total_a4(:,i);

        g_a4 = 0.024;% nanoSieman, Robert A. and Howe J., 2003
        Po_a4 = open_state_a4./(N_ampar_PerCluster*(1-ratio)); % open probability on PSD
        
        if SpillOver == 1
            name_SpillOver = ['GluA4_P',num2str(period),'_',num2str(i+80) '_SpillOver.mat'];
            load(name_SpillOver);
            open_state_total_a4_SpillOver(:,i) = amparStates(:,4);
            open_state_a4_SpillOver = open_state_total_a4_SpillOver(:,i);
            Po_a4 = Po_a4 + NND_No * open_state_a4_SpillOver./(N_ampar_PerCluster*ratio); % open probability on PSD
            open_state_a4 = open_state_a4 + NND_No * open_state_a4_SpillOver;

        end

        for k = 1:(timestepSize-1)
            V_a4(k+1)=V_a4(k)+dt*(1/Cm)*((-1*gm*(V_a4(k)-Vm))+(-1.*g_a4.*Po_a4(k)*(V_a4(k)-Vampa))); %mV, Gupta R, Reneaux M, Karmeshu, 2016
        end
        

        I_2_a4 = zeros(timestepSize ,1);
        for k = 1:timestepSize
            I_2_a4(k) = g_a4 .* open_state_a4(k) * (V_a4(k)-Vampa);%Ampa excitatory current in picoAmpere.
        end
        I_total_a4(:,i) = I_2_a4;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Calculate the 10-90% rise time and decay time %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmp = timestepSize - baseline; 
        X = (1:tmp)'; 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GluA2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I_2_a2 = mean(I_total_a2,2);
        I_fit1 = I_2_a2;
        I_fit1(1:baseline) = []; % Remove the '0' part

        I_max_a2 = min(I_fit1); % Actually it's the peak value of current              
        index_max_a2 = find(I_fit1 == I_max_a2); 

        % extract rise phase
        X_rise_a2 = X; X_rise_a2((index_max_a2+1):tmp) = [];
        I_rise_a2 = I_fit1; I_rise_a2((index_max_a2+1):tmp) = [];
        % extract decay phase
        X_decay_a2 = X; X_decay_a2(1:(index_max_a2-1)) = [];
        I_decay_a2 = I_fit1; I_decay_a2(1:(index_max_a2-1)) = [];


        f_decay_a2 = fit(X_decay_a2(1:range),I_decay_a2(1:range),'exp2','Robust','LAR','Algorithm','Levenberg-Marquardt'); % double exponential fit, for decay time                                                   
        I_37_a2 = 0.37 * I_max_a2; % 37% of the current peak
        syms x ;
        fun1_a2 = f_decay_a2.a * exp(f_decay_a2.b * x) + f_decay_a2.c * exp(f_decay_a2.d * x) == I_37_a2;  % V = V0 * exp(-t/tau) = a * exp(b * x)
        decayTime_a2 = double(vpasolve(fun1_a2,x,1000) - index_max_a2)/2000;

        f_rise_a2 = fit(X_rise_a2(:),I_rise_a2,'exp2','Robust','LAR','Algorithm','Levenberg-Marquardt'); % double exponential fit, for decay time                                                   
        I_10_a2 = 0.1 * I_max_a2; % 10% of the current peak
        I_90_a2 = 0.9 * I_max_a2; % 90% of the current peak
        syms x ;
        fun10_a2 = f_rise_a2.a * exp(f_rise_a2.b * x) + f_rise_a2.c * exp(f_rise_a2.d * x) == I_10_a2;  % V = V0 * exp(-t/tau) = a * exp(b * x)
        fun90_a2 = f_rise_a2.a * exp(f_rise_a2.b * x) + f_rise_a2.c * exp(f_rise_a2.d * x) == I_90_a2;
        riseTime10_a2 = double(vpasolve(fun10_a2,x));
        riseTime90_a2 = double(vpasolve(fun90_a2,x));

        riseTime_a2 = double(riseTime90_a2-riseTime10_a2)/2000;

        plot(I_total_a2,'Color',[0.65,0.65,0.65]); hold on;
        plot(I_2_a2,'r');hold on;

        decay_Text_a2 = ['GluA2 decay time = ', num2str(decayTime_a2)];
        text(12000,-6,decay_Text_a2);
        rise_Text_a2 = ['GluA2 rise time = ', num2str(riseTime_a2)];
        text(12000,-7,rise_Text_a2);
        amplitude_a2 = ['GluA2 amplitude = ', num2str(I_max_a2)];
        text(12000,-8,amplitude_a2);

        xlim([0 timestepSize]);
        set(gca,'XTick',[0 4000 8000 12000 16000 20000],'XTickLabel',[0 2 4 6 8 10]);

        title('Slow-gating GluA2 ');
        xlabel('time (ms)');ylabel('mEPSC (pA)');
        box off;
        
        %%%%%%%%%%%%%%%%%%%% GluA4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        I_2_a4 = mean(I_total_a4,2);
        I_fit_a4 = I_2_a4;
        I_fit_a4(1:baseline) = []; % Remove the '0' part    

        I_max_a4 = min(I_fit_a4);% Actually it's the peak value of current
        index_max_a4 = find(I_fit_a4 == I_max_a4);

%        extract rise phase
        X_rise_a4 = X; X_rise_a4((index_max_a4 + 1) : tmp) = [];
        I_rise_a4 = I_fit_a4; I_rise_a4 ((index_max_a4 + 1) : tmp) = [];
%        extract decay phase
        X_decay_a4 = X; X_decay_a4(1:(index_max_a4-1)) = [];
        I_decay_a4 = I_fit_a4; I_decay_a4(1:(index_max_a4-1)) = [];                


        f_decay_a4 = fit(X_decay_a4(1:range),I_decay_a4(1:range),'exp2','Robust','LAR','Algorithm','Levenberg-Marquardt'); % double exponential fit, for decay time                                                   
        I_37_a4 = 0.37 * I_max_a4; % 37% of the current peak
        syms x ;
        fun1_a4 = f_decay_a4.a * exp(f_decay_a4.b * x) + f_decay_a4.c * exp(f_decay_a4.d * x) == I_37_a4;  % V = V0 * exp(-t/tau) = a * exp(b * x)
        decayTime_a4 = double(vpasolve(fun1_a4,x,1000) - index_max_a4)/2000;

        f_rise_a4 = fit(X_rise_a4(:),I_rise_a4,'exp2','Robust','LAR','Algorithm','Levenberg-Marquardt'); % double exponential fit, for decay time                                                   
        I_10_a4 = 0.1 * I_max_a4; % 10% of the current peak
        I_90_a4 = 0.9 * I_max_a4; % 90% of the current peak
        syms x ;
        fun10_a4 = f_rise_a4.a * exp(f_rise_a4.b * x) + f_rise_a4.c * exp(f_rise_a4.d * x) == I_10_a4;  % V = V0 * exp(-t/tau) = a * exp(b * x)
        fun90_a4 = f_rise_a4.a * exp(f_rise_a4.b * x) + f_rise_a4.c * exp(f_rise_a4.d * x) == I_90_a4;
        riseTime10_a4 = double(vpasolve(fun10_a4,x));
        riseTime90_a4 = double(vpasolve(fun90_a4,x));

        riseTime_a4 = double(riseTime90_a4 - riseTime10_a4)/2000; 

        plot(I_total_a4,'Color',[0.65,0.65,0.65]); hold on;
        plot(I_2_a4,':r');hold on;     

        decay_Text_a4 = ['GluA4 decay time = ', num2str(decayTime_a4)];
        text(12000,-10,decay_Text_a4);        
        rise_Text_a4 = ['GluA4 rise time = ', num2str(riseTime_a4)];
        text(12000,-11,rise_Text_a4);
        amplitude_a4 = ['GluA4 amplitude = ', num2str(I_max_a4)];
        text(12000,-12,amplitude_a4);

        xlim([0 timestepSize]);
        set(gca,'XTick',[0 4000 8000 12000 16000 20000],'XTickLabel',[0 2 4 6 8 10]);

        title('Fast-gating GluA4 ');
        xlabel('time (ms)');ylabel('mEPSC (pA)');
        box off;
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%% Combine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmp = timestepSize - baseline; 
        X = (1:tmp)'; 
        
        if ratio == 0
            I_2_a2 = zeros(20000,runs);
        elseif ratio ==1
            I_2_a4 = zeros(20000,runs);
        end
        
        I_total = I_2_a2 + I_2_a4;
        I_fit = I_total;
        I_fit(1:baseline) = []; % Remove the '0' part    

        I_max = min(I_fit);% Actually it's the peak value of current
        index_max = find(I_fit == I_max);

        % extract rise phase
        X_rise = X; X_rise((index_max + 1) : tmp) = [];
        I_rise = I_fit; I_rise((index_max + 1) : tmp) = [];
       % extract decay phase
        X_decay = X; X_decay(1:(index_max-1)) = [];
        I_decay = I_fit; I_decay(1:(index_max-1)) = [];                

        
        f_decay = fit(X_decay(1:range),I_decay(1:range),'exp2','Robust','LAR','Algorithm','Levenberg-Marquardt'); % double exponential fit, for decay time                                                   
        I_37 = 0.37 * I_max; % 37% of the current peak
        syms x ;
        fun = f_decay.a * exp(f_decay.b * x) + f_decay.c * exp(f_decay.d * x) == I_37;  % V = V0 * exp(-t/tau) = a * exp(b * x)
        decayTime = double(vpasolve(fun,x,1000) - index_max)/2000;

        f_rise = fit(X_rise(:),I_rise,'exp2','Robust','LAR','Algorithm','Levenberg-Marquardt'); % double exponential fit, for decay time                                                   
        I_10 = 0.1 * I_max; % 10% of the current peak
        I_90 = 0.9 * I_max; % 90% of the current peak
        syms x ;
        fun10 = f_rise.a * exp(f_rise.b * x) + f_rise.c * exp(f_rise.d * x) == I_10;  % V = V0 * exp(-t/tau) = a * exp(b * x)
        fun90 = f_rise.a * exp(f_rise.b * x) + f_rise.c * exp(f_rise.d * x) == I_90;
        riseTime10 = double(vpasolve(fun10,x));
        riseTime90 = double(vpasolve(fun90,x));

        riseTime = double(riseTime90 - riseTime10)/2000; 

        plot(I_total,'k'); hold off; 
        
        decay_Text = ['GluA1&GluA4 decay time = ', num2str(decayTime)];
        text(12000,-14,decay_Text);
        rise_Text = ['GluA1&GluA4 rise time = ', num2str(riseTime)];
        text(12000,-15,rise_Text);
        amplitude = ['Combined amplitude = ', num2str(I_max)];
        text(12000,-16,amplitude);
                
        xlim([0 timestepSize]);
        set(gca,'XTick',[0 4000 8000 12000 16000 20000],'XTickLabel',[0 2 4 6 8 10]);

        title('Combined ')
        xlabel('time (ms)');ylabel('mEPSC (pA)');
        box off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beep;