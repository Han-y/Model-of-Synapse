# Model-of-Synapse
This is the README for the models associated with the paper:

Ying Han, Ran Cao, Liming Qin, Lulu Y.Chen, Ai-Hui Tang, Thomas C. Südhof, Bo Zhang
Neuroligin-3 confines AMPA-receptors into nanoclusters, thereby controlling synaptic strength at the calyx of Held synapses. (in submission)

To make the code concise, we used GluA1 represents slow-GluAs and GluA4 represents fast-GluAs, respectively. 
There are two folders containing seven files under each folder. The program in each folder is based on the STORM data of GluA1/GluA4 or PSD-95, respectively.

The order of operation should be: 
(1) make.m 
    To get the .mexw64 file that can be invoked within MATLAB. This will result in three new files: 
    synapse_absorb_at_cleftbd_C.mexw64, synapse_absorb_at_glia_C.mexw64, and synapse_C.c.
(2) run_sim.m
    This file can change the parameters of simulation and 
    get the operation result including states of AMPARs, states of glutamates and so on.
    It needs to take hours if you run 160 times and the running time depends on your device.
(3) (Optional) synapse_sim.m
    Including all the parameters that can be changed in simulation.
(4) EPSC_generation.m
    This file can get currents of fast- and slow-GluAs and their combined traces. 

The remaining four files don't need to be run.
(5) printLoopStateInfo.m
    This file name is literal. It prints information of loop state.
(6) synapse_C.c
    Run to make compile. Note that the parameters in this file should be the same as in synapse_fun.m.
(7) synapse_fun.m
    The simulation structure.
