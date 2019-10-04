runs = 1; % number of runs
ratio = 0.45; % ratio = GluA2 / (GluA2+GluA4)
              % P4: 0.25; P8: 0.3; P12: 0.45
              % P18: 0.65; P30: 0.75
PSD_factor = 2; % can be changed under different stages
PSD = 164; release_zone = 132; %nm; can be changed under different stages
duration = 0.2; % ms; duration of vesicle release
amparNo = 100; % number of AMPARs 

SpillOver = 0; % SpillOver = 0: there is no spillover
               % SpillOver = 1: include spillover with one, two, three, or
               % four neighboring synapses

load('Vesicle_Distributions_P12'); % get the location of releaes site in each run

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%      Run Simulations      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:runs
    fileNo = j;
    period = ['P',num2str(12)];
    synapse_sim(10000,amparNo,ratio,fileNo,period, random_number_P12(j), PSD_factor,duration, release_zone, PSD, SpillOver);
end

beep;