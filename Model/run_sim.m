%%%% ratio = GluA2 / (GluA2+GluA4)
runs = 160;
ratio = 0.75; % P4
amparNo = 100;
PSD_factor = 1;
PSD = 132; % nm
release_zone = 132; %nm

duration = 0.2; %ms

load('Vesicle_Distributions_P04');


for j = 1:runs
    fileNo = j;
    period = ['P',num2str(04)];
    synapse_sim(10000,amparNo,ratio,fileNo,period, random_number_P4(j), PSD_factor,duration, release_zone, PSD);
end

beep;