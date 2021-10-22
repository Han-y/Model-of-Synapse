clear; 

runs = 160; % number of runs
stage =4;
ratio = 0.75; % ratio = GluA2 / (GluA2+GluA4)
              % P4 : 0.75;  P8 : 0.7; P12: 0.45
              % P16: 0.35;  P30: 0.25
factor = 1.0;
              

amparNo = 100; % number of AMPARs
duration = 0.2; % ms; duration of vesicle pore open 
gluPerVesicle = 6000;
D_glu = 0.4 * 1E6; %400000,0.2 (mu m)^2/ms = 400000 nm^2/ms
transpDensity= 2500;

SpillOver = 0; % SpillOver = 0: there is no spillover
               % SpillOver = 1: include spillover with one, two, three, or
               % four neighboring synapses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STROM data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% P8 ctrl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% A1 %%%%%
Nano_Volume_A1 = 2.206; %nm^2
Nano_Density_A1 = 3.577;
Synaptic_Volume_A1 = 15.47;
Synaptic_Density_A1 = 2.399;
ratio_Num_SynNano_STORM_A1 = (Synaptic_Volume_A1 * Synaptic_Density_A1) / (Nano_Volume_A1 * Nano_Density_A1);

%%% used data in simulation
ratio_EM_STORM = 132/(sqrt(30.35)); % 132 (nm) got from EM data
                                         % 30.35 is the synaptic volume of GluA4 under STORM at P16
Synaptic_radius_A1 = ratio_EM_STORM * sqrt(Synaptic_Volume_A1); 
Nano_radius_A1 = ratio_EM_STORM * sqrt(Nano_Volume_A1); % nm


%%%%% A4 %%%%%
Nano_Volume_A4 = 2.105; %nm^2
Nano_Density_A4 = 7.7;
Synaptic_Volume_A4 = 22.93;
Synaptic_Density_A4 = 3.792;
ratio_Num_SynNano_STORM_A4 = (Synaptic_Volume_A4 * Synaptic_Density_A4) / (Nano_Volume_A4 * Nano_Density_A4);
%%% used data in simulation
Synaptic_radius_A4 = ratio_EM_STORM * sqrt(Synaptic_Volume_A4); % 132 nm, got from EM data
Nano_radius_A4 = ratio_EM_STORM * sqrt(Nano_Volume_A4); % nm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% P8 ko %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% A1 %%%%%
% Nano_Volume_A1 = 1.763; %nm^2
% Nano_Density_A1 = 4.328;
% Synaptic_Volume_A1 = 14.69;
% Synaptic_Density_A1 = 2.915;
% ratio_Num_SynNano_STORM_A1 = (Synaptic_Volume_A1 * Synaptic_Density_A1) / (Nano_Volume_A1 * Nano_Density_A1);
% 
% %%% used data in simulation
% ratio_EM_STORM = 132/(sqrt(30.35)); % 132 (nm) got from EM data
%                                          % 30.35 is the synaptic volume of GluA4 under STORM at P16
% Synaptic_radius_A1 = ratio_EM_STORM * sqrt(Synaptic_Volume_A1); 
% Nano_radius_A1 = ratio_EM_STORM * sqrt(Nano_Volume_A1); % nm
% 
% 
% %%%%% A4 %%%%%
% Nano_Volume_A4 = 2.055; %nm^2
% Nano_Density_A4 = 8.024;
% Synaptic_Volume_A4 = 24.05;
% Synaptic_Density_A4 = 3.709;
% ratio_Num_SynNano_STORM_A4 = (Synaptic_Volume_A4 * Synaptic_Density_A4) / (Nano_Volume_A4 * Nano_Density_A4);
% %%% used data in simulation
% Synaptic_radius_A4 = ratio_EM_STORM * sqrt(Synaptic_Volume_A4); % 132 nm, got from EM data
% Nano_radius_A4 = ratio_EM_STORM * sqrt(Nano_Volume_A4); % nm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% P12 ctrl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% A1 %%%%%
% Nano_Volume_A1 = 2.911; %nm^2
% Nano_Density_A1 = 4.0;
% Synaptic_Volume_A1 = 33.64;
% Synaptic_Density_A1 = 1.988;
% ratio_Num_SynNano_STORM_A1 = (Synaptic_Volume_A1 * Synaptic_Density_A1) / (Nano_Volume_A1 * Nano_Density_A1);
% 
% %%%%% used data in simulation
% ratio_EM_STORM = 132/(sqrt(30.35)); % 132 (nm) got from EM data
%                                       %  30.35 is the synaptic volume of GluA4 under STORM at P16
% Synaptic_radius_A1 = ratio_EM_STORM * sqrt(Synaptic_Volume_A1); 
% Nano_radius_A1 = ratio_EM_STORM * sqrt(Nano_Volume_A1); % nm
% 
% 
% %%%% A4 %%%%%
% Nano_Volume_A4 = 1.992; %nm^2
% Nano_Density_A4 = 4.723;
% Synaptic_Volume_A4 = 26.41;
% Synaptic_Density_A4 = 2.435;
% ratio_Num_SynNano_STORM_A4 = (Synaptic_Volume_A4 * Synaptic_Density_A4) / (Nano_Volume_A4 * Nano_Density_A4);
% %%%% used data in simulation
% Synaptic_radius_A4 = ratio_EM_STORM * sqrt(Synaptic_Volume_A4); % 132 nm, got from EM data
% Nano_radius_A4 = ratio_EM_STORM * sqrt(Nano_Volume_A4); % nm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% P12 KO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% A1 %%%%%
% Nano_Volume_A1 = 4.627; %nm^2
% Nano_Density_A1 = 1.833;
% Synaptic_Volume_A1 = 37.27;
% Synaptic_Density_A1 = 1.225;
% ratio_Num_SynNano_STORM_A1 = (Synaptic_Volume_A1 * Synaptic_Density_A1) / (Nano_Volume_A1 * Nano_Density_A1);
% 
% %%% used data in simulation
% ratio_EM_STORM = 132/(sqrt(30.35)); % 132 (nm) got from EM data
%                                          % 30.35 is the synaptic volume of GluA4 under STORM at P16
% Synaptic_radius_A1 = ratio_EM_STORM * sqrt(Synaptic_Volume_A1); 
% Nano_radius_A1 = ratio_EM_STORM * sqrt(Nano_Volume_A1); % nm
% 
% 
% %%%%% A4 %%%%%
% Nano_Volume_A4 = 3.139; %nm^2
% Nano_Density_A4 = 4.359;
% Synaptic_Volume_A4 = 40.11;
% Synaptic_Density_A4 = 2.16;
% ratio_Num_SynNano_STORM_A4 = (Synaptic_Volume_A4 * Synaptic_Density_A4) / (Nano_Volume_A4 * Nano_Density_A4);
% %%% used data in simulation
% Synaptic_radius_A4 = ratio_EM_STORM * sqrt(Synaptic_Volume_A4); % 132 nm, got from EM data
% Nano_radius_A4 = ratio_EM_STORM * sqrt(Nano_Volume_A4); % nm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% P16 ctrl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% A1 %%%%%
% Nano_Volume_A1 = 2.486; %nm^2
% Nano_Density_A1 = 2.805;
% Synaptic_Volume_A1 = 30.35;
% Synaptic_Density_A1 = 1.803;
% ratio_Num_SynNano_STORM_A1 = (Synaptic_Volume_A1 * Synaptic_Density_A1) / (Nano_Volume_A1 * Nano_Density_A1);
% 
% %%% used data in simulation
% ratio_EM_STORM = 132/(sqrt(30.35)); % 132 (nm) got from EM data
%                                          % 30.35 is the synaptic volume of GluA4 under STORM at P16
% Synaptic_radius_A1 = ratio_EM_STORM * sqrt(Synaptic_Volume_A1); 
% Nano_radius_A1 = ratio_EM_STORM * sqrt(Nano_Volume_A1)/1; % nm
% 
% 
% %%%%% A4 %%%%%
% Nano_Volume_A4 = 4.029; %nm^2
% Nano_Density_A4 = 2.983;
% Synaptic_Volume_A4 = 45.43;
% Synaptic_Density_A4 = 1.581;
% ratio_Num_SynNano_STORM_A4 = (Synaptic_Volume_A4 * Synaptic_Density_A4) / (Nano_Volume_A4 * Nano_Density_A4);
% %%% used data in simulation
% Synaptic_radius_A4 = ratio_EM_STORM * sqrt(Synaptic_Volume_A4); % 132 nm, got from EM data
% Nano_radius_A4 = ratio_EM_STORM * sqrt(Nano_Volume_A4)/1; % nm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% P16 KO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% A1 %%%%%
% Nano_Volume_A1 = 3.177; %nm^2
% Nano_Density_A1 = 2.423;
% Synaptic_Volume_A1 = 32.74;
% Synaptic_Density_A1 = 1.532;
% ratio_Num_SynNano_STORM_A1 = (Synaptic_Volume_A1 * Synaptic_Density_A1) / (Nano_Volume_A1 * Nano_Density_A1);
% 
% %%%% used data in simulation
% ratio_EM_STORM = 132/(sqrt(30.35)); % 132 (nm) got from EM data
%                                          30.35 is the synaptic volume of GluA4 under STORM at P16
% Synaptic_radius_A1 = ratio_EM_STORM * sqrt(Synaptic_Volume_A1); 
% Nano_radius_A1 = ratio_EM_STORM * sqrt(Nano_Volume_A1); % nm
% 
% 
% %%%% A4 %%%%%
% Nano_Volume_A4 = 3.931; %nm^2
% Nano_Density_A4 = 3.242;
% Synaptic_Volume_A4 = 43.28;
% Synaptic_Density_A4 = 1.688;
% ratio_Num_SynNano_STORM_A4 = (Synaptic_Volume_A4 * Synaptic_Density_A4) / (Nano_Volume_A4 * Nano_Density_A4);
% %%%% used data in simulation
% Synaptic_radius_A4 = ratio_EM_STORM * sqrt(Synaptic_Volume_A4); % 132 nm, got from EM data
% Nano_radius_A4 = ratio_EM_STORM * sqrt(Nano_Volume_A4); % nm


release_zone = (Synaptic_radius_A1 + Synaptic_radius_A4)/2; %nm; can be changed under different stages


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%      Run Simulations      %%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for j = 1:runs
    fileNo = j;
    period = ['P',num2str(stage)]; 
    synapse_sim(gluPerVesicle, amparNo,ratio,fileNo,period,transpDensity,D_glu,factor,...
                Nano_radius_A1, Synaptic_radius_A1, ratio_Num_SynNano_STORM_A1,...
                Nano_radius_A4, Synaptic_radius_A4, ratio_Num_SynNano_STORM_A4,...                
                duration, release_zone, SpillOver);
end

beep;