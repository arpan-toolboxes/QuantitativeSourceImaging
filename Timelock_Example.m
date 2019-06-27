%% Example of localizing gaussian pulse as per distributed dipole condition
clear
clc
% Simulation of EEG data
load('Forward_Var_Colin.mat');
load('data_base_TL');

scale=[95 410 620 885 1000]; % yield different SNR
lag=[0] %[15 30 45]; insert time lags between the signals in each hemisphere, or else random jitters applied

for ii=1:length(scale)
    fun=[];
    fun.grid=grid; %leadfield
    fun.headmodel=vol; %headmodel 
    fun.elec=elec; % electrode structure
    fun.model='dist'; % or 'single' or 'two'
    fun.signal='gauss';
    fun.noise_perc= 0;
    fun.noise_perc_sensor = 0;
    fun.dip_pos=[-60, -28 ,-6; 64, -24 ,6];
    fun.scale=scale(ii);
    fun.lag=lag;
    [data_stim]=fun_simulation(fun);

    % Adding the signal to baseline
    Sim_Signal{ii}=data_stim;
    Sim_Signal{ii}.trial{1}=Sim_Signal{ii}.trial{1}(:,1:1000)+data_base_avg.trial{1};
    Sim_Signal{ii}.time{1}=Sim_Signal{ii}.time{1}(1,1:1000);
    
    ii
    
    % Obtaining SNR
    SNR(ii)=SNR_Gol(Sim_Signal{ii}.trial{1}',data_base_avg.trial{1}'); 
end
SNR
%% Source Localization
for ii=1:5
    fun=[];
    fun.grid=grid;
    fun.headmodel=vol;
    fun.elec=elec;
    fun.analysis='timelock';
    fun.simdata=Sim_Signal{ii};
    fun.basedata=data_base_avg;
    % fun.freq=[40];
    % fun.taper=2;
    fun.time = [375 425];
    fun.mri=mri;
    fun.lambda=1;
    fun.pow='contrast';
    
    [source_eLOR_sim, source_LCMV_sim, source_MNE_sim] = fun_SL(fun);
    
    % Z-scores dim{SNRs}{various z threshold}
    
    z_thr=0.9999;
    [Z_LCMV{ii},Z_LCMV_P{ii}] = source_zscore(source_LCMV_sim,'dist',z_thr);
    [Z_eLOR{ii},Z_eLOR_P{ii}] = source_zscore(source_eLOR_sim,'dist',z_thr);
    [Z_MNE{ii},Z_MNE_P{ii}]  = source_zscore(source_MNE_sim,'dist',z_thr);
    
    dip_pos=[-60, -28 ,-6; 64, -24 ,6];
    
    % plotting on surface
    figure;
    coords2surf2(Z_eLOR{ii}{1},Z_LCMV{ii}{1},dip_pos,Z_MNE{ii}{1}); %(red, green, blue, pink)
end