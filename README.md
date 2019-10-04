# Model-of-Synapse
This is the README for the models associated with the paper:

Ying Han, Lulu Y.Chen, Thomas C.Sudhof, Bo Zhang
Neuroligins maintain clustering of AMPA receptors at a central synapse. (in submission)

To make the code concise, we used GluA2 represents slow-GluAs and GluA4 represents fast-GluAs. 
There are 15 files, the order of operation should be: 
(1) make.m
    To get the .mexw64 file that can be invoked within MATLAB.
(2) random_release_zone.m
    To get the random distribution of release sites. The default value of random sites is 160. Note that, the random release sites in 160 runs we used are shown in "Vesicle_distributions_P04~30.mat".
(3) run_sim.m
    This file can change the parameters of simulation and 
    get the operation result including states of AMPAR, states of glutamates and so on.
    It needs to take hours if you run 160 times and the running time depends on your device.
(4) EPSC_generation.m
    This file can get mEPSCs of fast-/slow-GluAs and their combined traces. This file can also get mEPSCs evoked by spillover.
(5) EPSC_generation_only_extrasynapse.m 
    This file can get mEPSCs of a single AMPAR in extrasynaptic region. It randomly chooses GluA2 or GluA4 in extrasynaptic region in these 160 runs.

The remaining four files don't need to be run.
(6) printLoopStateInfo.m
    This file name is literal. It prints information of loop state.
(7) synapse_C.c
    Run to make compile. Note that the parameters in this file should be the same as in synapse_fun.m.
(8) synapse_fun.m
    The simulation structure.
(9) synapse_sim.m
    Including all the parameters that can be changed in simulation.
 
 If you have any problems, please don't hesitate to contact me: hany1900@pku.edu.cn.
