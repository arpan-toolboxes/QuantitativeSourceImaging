Put all the files in path. All codes are based on fieldtrip toolbox, therefore add fieldtrip to path as well.

1) Timelock_Example: This simulates mixture of gaussian with different SNR values and drifts. It also computes sources based on the covariance of the timelocked response using eLORETA, MNE. and LCMV.
2) Freqency_Example: This simulates mixture of sinusoidal with different SNR values and power ratios. This also computes sources based on the cross spectral density of the frequency response using eLORETA and DICS.

Functions:
1) fun_simulation: function to simulate sinusoidal or gaussian based on the arguments entered in the data structure 'fun'.
2) SNR_Gol: SNR computed using equation given in Goldenholz2009.
3) fun_SL: The function which calls the source localization commands and interpolates it to the T1- MRI.
4) source_zscore: Computes the z scores of all voxels in each hemisphere. Further statistics employed.
5) coords2surf2: Function to view the computed sources on brain surface.

Workspace variables:
1) Headmodel_MRI: headmodel and mri on Colin27 model. 
2) Grid_Elec: Includes the leadfield computed from the headmodel with the electrode configuration based on 64 channel cap.
3) data_base_TL: baseline resting state data acquired empirically, epoched according to the length of the simulated signal.
4) data_base_freq: baseline resting state data acquired empirically, epoched according to the length of the simulated signal.
5) Template.gii: Used by coords2surf2 to import surface of the brain. 