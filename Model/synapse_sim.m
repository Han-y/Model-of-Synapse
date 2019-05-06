%
%   simulation
%

function synapse_sim(gluPerVesicle,amparNo, ratio,fileNo,period, random_number,PSD_factor,duration,release_zone,PSD)


%---- globals -------------------------------------------------------------
global  matlabFileInfo  %#ok

global  amparPeakDistr  %#ok
global  amparTags       %#ok
global  amparTagsVar    %#ok
global  amparStates     %#ok
global  amparMean       %#ok
global  amparVar        %#ok
global  amparVarS       %#ok
global  gluStates       %#ok
global  gluDistrES      %#ok

global  amparState_F        %#ok
global  amparState_B        %#ok
global  amparState_D        %#ok
global  amparState_O        %#ok

%---- random seeds --------------------------------------------------------
randn('state', 0);
rand('seed', 0);

%---- parameters ----------------------------------------------------------
clear('S'); % DONT REMOVE THIS LINE

%   S.savePath              = '';
    S.workspaceFilePrefix   = 'sim';    
    S.simulation            = 'new';
        
    
%---- statistics ----------------------------------------------------------
    S.runs                  = 1;
%---- pulses --------------------------------------------------------------
    S.waitBeforeFirstPulse  = 2; % in ms
    S.pulseNo               = 1;
    S.pulseToPulseTime      = 20;% ms, 1500;
    S.gluPerVesicle         = gluPerVesicle;
    S.vesiclesPerPulse      = 1;
    S.vesicleJitter         = 0;
    S.vesFillingFraction    = [1 1];
    S.vesFusionDuration     = duration; % ms
    
    % Kinetic schemes for AMPARs
    % Form: [x1 x2 x3 x4], where this vector contains the fractions of AMPARs
    %       with
    %  x1 = JS scheme
    %  x2 = MN scheme without TARP ligation
    %  x3, x4 = MN scheme with ligations of type 1 and 2, cf. to
    %           the Milstein, Nicoll paper on AMPARs and TARPs
    
%---- time steps ----------------------------------------------------------
    S.timeStepsPerMS        = round(1000*2); %=0.5 microseconds
%---- cleft geometry control constants ------------------------------------
    S.R_Cleft               = 1128; % radius, nm                      % 400/2
    S.R_PSD                 = PSD * PSD_factor*1 ; % radius, nm. Data is from Xue Lei % 200/2
    S.R_ES                  = S.R_Cleft + 40;   % distance between pre/post-synaptic cylinders and glial sheath
    S.R_Reservoir           = 126 * PSD_factor;
    
    S.R_ReleaseZone         = 0 + (release_zone - 0) * random_number; % 0;%S.R_Cleft% radius of release zone, 0=center 

    S.cleftHeight           = 28; % In Calyx of Held
    S.esHeight              = 10000;
    S.esBinSize             = 50; % in nm, for the distribution of glu in the ES
%---- AMPAR control constants ---------------------------------------------
    S.D_ampar_PSD           = 00;       % AMPAR diffusion constant in nm^2/ms, on the PSD % Choquet  0.2 (mu m)^2/s = 200 nm^2/ms
    S.D_ampar_outside       = 00;       % inside the reservoir
    
    a = 0.1;                            % outside/inside density ratio
    p = 1.0;    
    S.P_Reflect_Inside      = 1-p*a; % 1-p*a for a <= 0.1, otherwise 1-p   % prob to be reflected hitting the PSD boundary from inside the PSD
    S.P_Reflect_Outside     = 1-p;   % 1-p/a for a >= 0.1, otherwise 1-p    

  
    S.R_BindingToAMPAR      = 10;   %in nm, binding radius, 2.5nm yields good agreement with analytically manageable situations
%---- glu control constants -----------------------------------------------
    S.D_glu                 = 300000;   % glu diffusion constant, Franks, Sejnowski 0.2 (mu m)^2/ms = 200000 nm^2/ms, OTHER: 0.6  (mu m)^2/ms, 0.76 (mu m)^2/ms
%---- transporter control constants ---------------------------------------
    S.D_transp              = 0;        % NOT implemented, set to 0
    S.transpDensity         = 2500 *1/1000^2; % in number per nm^2
    S.transpBindingTime     = 5.5;      % in ms, IGNORED
    S.transpPumpingTime     = 38;       % in ms, IGNORED

    S.absorbAtGlia          = 0;        % set to 1 if glu should be absorbed upon hitting the glial sheath
    S.absorbAtCleftBD       = 0;        % set to 1 if glu should be absorbed upon leaving the cleft

    S.R_BindingToTransp     = 20;       %in nm, binding radius

    % rates for the kinetic model of a transporter, see Franks, Stevens, Sejnowsky
    S.transpRate_assoc      = 1.8  *10^7 *1/1000;
    S.transpRate_dissoc     = 1.8  *10^2 *1/1000;
    S.transpRate_pump       = 1.8  *10^2 *1/1000;
    S.transpRate_getReady   = 2.57 *10^1 *1/1000;
%--------------------------------------------------------------------------  
  
    for i = 1:2
        if i == 1 % GluA2
            % AMPAR number on the PSD at time 0
            S.amparNo_PSD           = round(amparNo * ratio)*1 ;
            if S.amparNo_PSD == 0
                continue
            end
            S.amparNo_outside       = round(amparNo * ratio)*0;
            S.kineticsFrac_PSD      = [1 0 0 0]; % for AMPARs initially on the PSD
            S.kineticsFrac_Outside  = [1 0 0 0];
            S.subunit = 'GluA2';
            S.factor =1.1;
            S.fileName = fileNo;
               
            S.period = period;
            synapse_fun(S);     % run the simulation with the parameters specified in S
        else % GluA4
            % AMPAR number on the PSD at time 0
            S.amparNo_PSD           = round(amparNo * (1 - ratio))*1 ; 
            if S.amparNo_PSD == 0
                continue
            end
            S.amparNo_outside       = round(amparNo * (1 - ratio))*0;
            S.kineticsFrac_PSD      = [0 1 0 0]; % for AMPARs initially on the PSD
            S.kineticsFrac_Outside  = [0 1 0 0];
            S.subunit = 'GluA4';
            S.factor =1;
            S.fileName = fileNo;  
            S.period = period;
            synapse_fun(S);     % run the simulation with the parameters specified in S
        end
    end