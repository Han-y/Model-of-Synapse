clear; 

runs = 160; % number of runs
stage =12;
ratio = 0.45; % ratio = GluA2 / (GluA2+GluA4)
              % P4 : 0.75;  P8 : 0.7; P12: 0.45
              % P16: 0.35;  P30: 0.25
              
 group = 'Ctrl';
%  group = 'KO';

factor = 1;

D_glu = 0.4 * 1E6;% glu diffusion constant,
cleftHeight = 28;
transpDensity= 5000;

duration = 0.2; % ms; duration of vesicle pore open 
gluPerVesicle = 8000;
release_zone = 60; %nm; can be changed under different stages

SpillOver = 0; % SpillOver = 0: there is no spillover
               % SpillOver = 1: include spillover with one, two, three, or
               % four neighboring synapses
                             
              
if strcmp(group, 'Ctrl')
    
    amparNo_a1 = 450; % /10 number of AMPARs
    amparNo_a4 = 550;
    
elseif strcmp(group, 'KO')
    
    amparNo_a1 = 423; % 93.8 AMPARs
    amparNo_a4 = 515;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%      Run Simulations      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for j = 1:runs
    fileNo = j;
    period = ['P',num2str(stage)]; 
    synapse_sim_on_PSD_STROM(gluPerVesicle, amparNo_a1, amparNo_a4, ratio, factor,...
                             fileNo,period,group,...
                             D_glu,transpDensity,cleftHeight,...                               
                             duration, release_zone, SpillOver);
end

beep;