function [source_elor_sim, source_lcmv_sim, source_mne_sim, source_dics_sim, topo] = fun_SL(fun)

% Source localization. Required fieldtrip toolbox.
% Inputs in [fun]:
% 1. grid, 2. headmodel, 3. elec
% 4. simdata = simulus
% 5. basedata = baseline
% 6. analysis = 'frequency' or 'timelock'
% 7. freq = [freq] frequency if analysis is 'frequency'
% 8. taper = for a frequency band
% 9. time = [time1 time2] time range if analysis is 'timelock'
% 10. fun.mri = MRI data str for interpolation
% 11. fun.pow = 'individual' or 'contrast'
% Output (in sequence):
% 1. source_eLOR_sim = eLORETA simulus source
% 2. source_eLOR_base = eLORETA baseline source
% 3. source2_sim = DICS/LCMV simulus source (based on analysis)
% 4. source2_base = DICS/LCMV baseline source

if isfield(fun,{'grid', 'elec', 'headmodel', 'analysis', 'mri', 'pow', 'simdata', 'basedata'}) ==1
    disp('All arguments present')
else
    
    error('Error: Put in all arguments. Check comments.');
end

source_elor_sim=[];
source_lcmv_sim=[];
source_mne_sim=[];
source_dics_sim=[];

if strcmp(fun.analysis,'frequency')
    
    alldata=ft_appenddata([],fun.simdata, fun.basedata);
    
    cfg = [];
    cfg.method    = 'mtmfft';
    cfg.output    = 'powandcsd';
    cfg.pad       = 'maxperlen';
    cfg.tapsmofrq = fun.taper;
    cfg.taper     = 'dpss';
    cfg.foilim   = fun.freq;
    
    freq_sim = ft_freqanalysis(cfg, fun.simdata);
    freq_base= ft_freqanalysis(cfg, fun.basedata);
    freq_all = ft_freqanalysis(cfg, alldata);  %for common filter
    
    %% Method : DICS
    fun.pow='contrast';
    if strcmp(fun.pow,'contrast')==1
        cfg=[];
        cfg.method='dics';
        cfg.frequency    = 40;
        cfg.grid         =  ft_convert_units(fun.grid,'mm');
        cfg.vol  =  ft_convert_units(fun.headmodel,'mm');
        cfg.dics.keepfilter='yes';
        cfg.grid.unit       = 'mm';
        cfg.elec =  ft_convert_units(fun.elec,'mm');
        cfg.elec.unit='mm';
        cfg.senstype =  'EEG';
        cfg.dics.lambda=fun.lambda;
        dics_all=ft_sourceanalysis(cfg, freq_all);
        
        cfg.dics.keepfilter='no';
        cfg.grid.filter   =  dics_all.avg.filter;
        
        dics_sim=ft_sourceanalysis(cfg, freq_sim);
        dics_base=ft_sourceanalysis(cfg, freq_base);
        
        dics_sim_pow=dics_sim;
        dics_sim_pow.avg.pow=(dics_sim.avg.pow-dics_base.avg.pow)./dics_base.avg.pow;
        
        cfg=[];
        cfg.parameter='pow';
        source_dics_sim=ft_sourceinterpolate(cfg, dics_sim_pow, fun.mri);
        
        
    elseif strcmp(fun.pow,'individual')==1
        cfg=[];
        cfg.method='dics';
        cfg.dics.keepfilter='no';
        cfg.grid         =  ft_convert_units(fun.grid,'mm');
        cfg.grid.unit       = 'mm';
        cfg.headmodel  =  ft_convert_units(fun.headmodel,'mm');
        cfg.elec =  ft_convert_units(fun.elec,'mm');
        cfg.elec.unit='mm';
        cfg.senstype =  'EEG';
        cfg.frequency    = 40;
        cfg.lambda=fun.lambda;
        cfg.projectnoise='yes';
        
        dics_sim=ft_sourceanalysis(cfg, freq_sim);
        dics_base=ft_sourceanalysis(cfg, freq_base);
        
        dics_sim_pow=dics_sim;
        dics_sim_pow.avg.pow=(dics_sim.avg.pow);
        
        cfg=[];
        cfg.parameter='pow';
        source_dics_sim=ft_sourceinterpolate(cfg,dics_sim_pow,fun.mri);
    end
    
    %% Method : eLORETA in frequency domain
    
    if strcmp(fun.pow,'contrast')==1
        cfg=[];
        cfg.method='eloreta';
        cfg.eloreta.keepfilter='yes';
        cfg.grid         =  ft_convert_units(fun.grid,'mm');
        cfg.grid.unit       = 'mm';
        cfg.eloreta.lambda=fun.lambda;
        cfg.vol  =  ft_convert_units(fun.headmodel,'mm');
        cfg.elec =  ft_convert_units(fun.elec,'mm');
        cfg.elec.unit='mm';
        cfg.senstype =  'EEG';
        cfg.frequency    = 40;
        
        elor_all=ft_sourceanalysis(cfg, freq_all);
        
        cfg.eloreta.keepfilter='no';
        cfg.grid.filter   =  elor_all.avg.filter;
        
        elor_sim=ft_sourceanalysis(cfg, freq_sim);
        elor_base=ft_sourceanalysis(cfg, freq_base);
        
        elor_sim_pow=elor_sim;
        elor_sim_pow.avg.pow=(elor_sim.avg.pow-elor_base.avg.pow);
        
        cfg=[];
        cfg.parameter='pow';
        source_elor_sim=ft_sourceinterpolate(cfg, elor_sim_pow, fun.mri);
        
    elseif strcmp(fun.pow,'individual')==1
        
        cfg=[];
        cfg.method='eloreta';
        cfg.eloreta.keepfilter='yes';
        cfg.grid         =  ft_convert_units(fun.grid,'mm');
        cfg.grid.unit       = 'mm';
        cfg.eloreta.lambda=fun.lambda;
        cfg.vol  =  ft_convert_units(fun.headmodel,'mm');
        cfg.elec =  ft_convert_units(fun.elec,'mm');
        cfg.elec.unit='mm';
        cfg.senstype =  'EEG';
        cfg.frequency    = 40;
        
        elor_sim=ft_sourceanalysis(cfg, freq_sim);
        
        elor_sim_pow=elor_sim;
        elor_sim_pow.avg.pow=(elor_sim.avg.pow);
        
        cfg=[];
        cfg.parameter='pow';
        source_elor_sim=ft_sourceinterpolate(cfg, elor_sim_pow, fun.mri);
        
        
        
    end
    
    %Topoplot
    cfg=[];
    cfg.elec=fun.elec;
    layout=ft_prepare_layout(cfg);
    
    freq_diff = freq_sim;
    freq_diff.powspctrm= freq_sim.powspctrm-freq_base.powspctrm;
    
    cfg=[];
    cfg.marker='labels';
    cfg.colorbar='EastOutside';
    cfg.markerfontsize=9;
    cfg.parameter='powspctrm';
    cfg.layout=layout;
    topo=figure,ft_topoplotER(cfg,freq_diff);
    
    
elseif strcmp(fun.analysis,'timelock')
    
    fun.pow='contrast';
    cfg=[];
    cfg.begsample    = fun.time(1,1);
    cfg.endsample    = fun.time(1,2);
    data_sim= ft_redefinetrial(cfg,fun.simdata);
    data_base= ft_redefinetrial(cfg,fun.basedata);
    
    alldata=ft_appenddata([],data_sim, data_base);
    
    cfg = [];
    cfg.covariance = 'yes';
    timelock_sim = ft_timelockanalysis(cfg, data_sim);
    timelock_base = ft_timelockanalysis(cfg, data_base);
    
    timelock_all=ft_timelockanalysis(cfg,alldata);
    
    cfg = [];
    cfg.elec =  ft_convert_units(fun.elec,'mm');
    cfg.elec.unit='mm';
    cfg.vol  =  ft_convert_units(fun.headmodel,'mm');
    cfg.grid    = fun.grid;
    cfg.method  = 'lcmv';
    cfg.lcmv.keepfilter= 'yes';
    cfg.senstype =  'EEG';
    cfg.lcmv.projectnoise='yes';
    cfg.lcmv.lambda=fun.lambda;
    lcmv_all = ft_sourceanalysis(cfg, timelock_all);
    
    cfg.gridfilter=lcmv_all.avg.filter;
    lcmv_sim = ft_sourceanalysis(cfg, timelock_sim);
    lcmv_base = ft_sourceanalysis(cfg, timelock_base);
    
    if strcmp(fun.pow,'individual')
        
        cfg=[];
        cfg.parameter='pow';
        source_lcmv_sim=ft_sourceinterpolate(cfg, lcmv_sim, fun.mri);
        
    elseif strcmp(fun.pow,'contrast')
        lcmv_sim_pow=lcmv_sim;
        lcmv_sim_pow.avg.pow= (lcmv_sim.avg.pow - lcmv_base.avg.pow) ./ lcmv_base.avg.pow;
        
        cfg=[];
        cfg.parameter='pow';
        source_lcmv_sim=ft_sourceinterpolate(cfg, lcmv_sim_pow, fun.mri);
    end
    
    %% Method : eLORETA in time domain
    
    fun.pow='contrast';
    if strcmp(fun.pow,'contrast');
        cfg=[];
        cfg.method='eloreta';
        cfg.eloreta.keepfilter='yes';
        cfg.grid         =  ft_convert_units(fun.grid,'mm');
        cfg.grid.unit       = 'mm';
        cfg.vol  =  ft_convert_units(fun.headmodel,'mm');
        cfg.elec =  ft_convert_units(fun.elec,'mm');
        cfg.elec.unit='mm';
        cfg.senstype =  'EEG';
        cfg.projectnoise='yes';
        cfg.eloreta.lambda=fun.lambda;
        source_all = ft_sourceanalysis(cfg, timelock_all);
        
        cfg.eloreta.keepfilter='yes';
        cfg.grid.filter  = source_all.avg.filter;
        elor_sim = ft_sourceanalysis(cfg, timelock_sim);
        elor_base = ft_sourceanalysis(cfg, timelock_base);
        
        elor_sim_pow=elor_sim;
        elor_sim_pow.avg.pow=(elor_sim.avg.pow)-elor_base.avg.pow;
        
        cfg=[];
        cfg.parameter='pow';
        source_elor_sim=ft_sourceinterpolate(cfg, elor_sim_pow, fun.mri);
        
    elseif strcmp(fun.pow,'individual')
        cfg=[];
        cfg.method='eloreta';
        cfg.eloreta.keepfilter='no';
        cfg.grid         =  ft_convert_units(fun.grid,'mm');
        cfg.grid.unit       = 'mm';
        cfg.eloreta.lambda=fun.lambda;
        cfg.vol  =  ft_convert_units(fun.headmodel,'mm');
        cfg.elec =  ft_convert_units(fun.elec,'mm');
        cfg.elec.unit='mm';
        cfg.senstype =  'EEG';
        cfg.projectnoise='yes';
        cfg.parameter='pow';
        elor_sim = ft_sourceanalysis(cfg, timelock_sim);
        
        source_elor_sim=ft_sourceinterpolate(cfg, elor_sim, fun.mri);
    end
    
    
    
    %% MNE
    
    fun.pow='contrast';
    cfg=[];
    cfg.begsample    = 70;
    cfg.endsample    = 71;
    data_sim= ft_redefinetrial(cfg,fun.simdata);
    data_base= ft_redefinetrial(cfg,fun.basedata); % comment if high computation available
    
    cfg = [];
    cfg.covariance = 'yes';
    timelock_sim = ft_timelockanalysis(cfg, data_sim);
    timelock_base = ft_timelockanalysis(cfg, data_base);
    
    cfg = [];
    cfg.elec =  ft_convert_units(fun.elec,'mm');
    cfg.elec.unit='mm';
    cfg.vol  =  ft_convert_units(fun.headmodel,'mm');
    cfg.grid    = fun.grid;
    cfg.method  = 'mne';
    cfg.mne.keepleadfield = 'yes';
    cfg.senstype =  'EEG';
    cfg.mne.lambda='1';
    cfg.mne.projectnoise='yes';
    mne_sim = ft_sourceanalysis(cfg, timelock_sim);
    mne_base = ft_sourceanalysis(cfg, timelock_base);
    
    if strcmp(fun.pow,'individual')
        
        cfg=[];
        cfg.parameter='pow';
        source_mne_sim=ft_sourceinterpolate(cfg, mne_sim, fun.mri);
        
    elseif strcmp(fun.pow,'contrast')
        mne_sim_pow=mne_sim;
        mne_sim_pow.avg.pow= (mne_sim.avg.pow-mne_base.avg.pow);
        
        cfg=[];
        cfg.parameter='pow';
        source_mne_sim=ft_sourceinterpolate(cfg, mne_sim_pow, fun.mri);
    end
    
    %Topoplot
    cfg=[];
    cfg.elec=fun.elec;
    layout=ft_prepare_layout(cfg);
    
    tl_diff = timelock_sim;
    tl_diff.avg= timelock_sim.avg-timelock_base.avg;
    
    cfg=[];
    cfg.marker='labels';
    cfg.colorbar='EastOutside';
    cfg.markerfontsize=9;
    cfg.parameter='avg';
    cfg.layout=layout;
    topo=figure,ft_topoplotER(cfg,tl_diff);
end