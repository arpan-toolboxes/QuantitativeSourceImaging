function [data_stim dipole_no dipole_mom data_base fig_stim fig_base]=fun_simulation(fun)

% Dipole simulation for 3 models: single, two point and distributed dipole
% models. Requires fieldtrip toolbox.
% Inputs in fun.fieldname:
% Requires 1. grid (leadfield), 2. headmodel and 3. elec (electrode configuration)
% To be specified:
% 4. signal (gauss or sin),
% 5. model (single, two, dist)
% 6. noise_perc (at source level)
% 7. noise_perc_sensor (at sensor level)
% 8. dip_pos (dipole locations for both stimulus and baseline condition, 1
% position each, for single model, 2 positions each for two and dist models)
% Position according to MNI coordinate system.
% Moment assumed to ba radially aligned.
% Output: SNR (in dB)
% data_stim (Stimulus simulated signal)
% data_base (Baseline random noise simulated signal).


if isfield(fun,{'grid', 'elec', 'headmodel', 'dip_pos', 'signal', 'model'}) ==1
    disp('All arguments present')
else
    
    error('Error: Put in all arguments. Check comments');
end



%% Set Dipole Location and moments according to the model
% fun.grid.pos=fun.grid.pos(fun.grid.inside,:);
if strcmp(fun.model,'single'); % single dipole condition
    
    % Stimulus
    dipole_stim_left=find(sqrt((fun.grid.pos(:,1) -fun.dip_pos(1,1)*ones(size(fun.grid.pos(:,2)))).^2+...
        (fun.grid.pos(:,2) -fun.dip_pos(1,2)*ones(size(fun.grid.pos(:,2)))).^2+...
        (fun.grid.pos(:,3) -fun.dip_pos(1,3)*ones(size(fun.grid.pos(:,2)))).^2)<=3)';
    
    dipole_loc_stim=[dipole_stim_left]';
    
    dip_mom_left=[-1 0 0];
    for xx=1:length(dipole_stim_left)-1
        dip_mom_left = cat(2,dip_mom_left,[-1 0 0]);
    end
    
    dip_mom_stim=[dip_mom_left];
    
elseif strcmp(fun.model,'two'); % two-point dipole condition
    
    % Stimulus
    dipole_stim_left=find(sqrt((fun.grid.pos(:,1) -fun.dip_pos(1,1)*ones(size(fun.grid.pos(:,2)))).^2+...
        (fun.grid.pos(:,2) -fun.dip_pos(1,2)*ones(size(fun.grid.pos(:,2)))).^2+...
        (fun.grid.pos(:,3) -fun.dip_pos(1,3)*ones(size(fun.grid.pos(:,2)))).^2)<=3)';
    
    dipole_stim_right=find(sqrt((fun.grid.pos(:,1) -fun.dip_pos(2,1)*ones(size(fun.grid.pos(:,2)))).^2+...
        (fun.grid.pos(:,2) -fun.dip_pos(2,2)*ones(size(fun.grid.pos(:,2)))).^2+...
        (fun.grid.pos(:,3) -fun.dip_pos(2,3)*ones(size(fun.grid.pos(:,2)))).^2)<=3)';
    
    dipole_loc_stim=[dipole_stim_left dipole_stim_right]';
    
    dip_mom_left=[-1 0 0];
    for xx=1:length(dipole_stim_left)-1
        dip_mom_left = cat(2,dip_mom_left,[-1 0 0]);
    end
    
    dip_mom_right=[1 0 0];
    for xx=1:length(dipole_stim_right)-1
        dip_mom_right = cat(2,dip_mom_right,[1 0 0]);
    end
    
    dip_mom_stim=[dip_mom_left dip_mom_right];
    
elseif strcmp(fun.model,'dist'); % distributed model: 12mm around the dipoles
    
    % Stimulus: Locations
    dipole_stim_left=find(sqrt((fun.grid.pos(:,1) -fun.dip_pos(1,1)*ones(size(fun.grid.pos(:,2)))).^2+...
        (fun.grid.pos(:,2) -fun.dip_pos(1,2)*ones(size(fun.grid.pos(:,2)))).^2+...
        (fun.grid.pos(:,3) -fun.dip_pos(1,3)*ones(size(fun.grid.pos(:,2)))).^2)<=12)';
    
    dipole_stim_right=find(sqrt((fun.grid.pos(:,1) -fun.dip_pos(2,1)*ones(size(fun.grid.pos(:,2)))).^2+...
        (fun.grid.pos(:,2) -fun.dip_pos(2,2)*ones(size(fun.grid.pos(:,2)))).^2+...
        (fun.grid.pos(:,3) -fun.dip_pos(2,3)*ones(size(fun.grid.pos(:,2)))).^2)<=12)';
    
    dipole_loc_stim=[dipole_stim_left dipole_stim_right]';
    
    dip_mom_left=[-1 0 0];
    for xx=1:length(dipole_stim_left)-1
        dip_mom_left = cat(2,dip_mom_left,[-1 0 0]);
    end
    
    dip_mom_right=[1 0 0];
    for xx=1:length(dipole_stim_right)-1
        dip_mom_right = cat(2,dip_mom_right,[1 0 0]);
    end
    
    dip_mom_stim=[dip_mom_left dip_mom_right];
    
end
%% Stimulus Simulation
nlength=5000; % Length of the the trial for sinusoidal signal
signal_stim=[];
time = (1:nlength)/1000;

fs=1000; %sampling frequency
sigma=0.05;
ts=0:1/fs:1; %time base

variance=sigma^2;
x1=1/(sqrt(2*pi*variance))*(exp(-(ts-0.4).^2/(2*variance)));
x2=-1/(sqrt(2*pi*variance))*(exp(-(ts-0.55).^2/(2*variance)));

x_gauss = x1+x2;
x_gauss=x_gauss/max(x_gauss);

if strcmp(fun.signal,'sin')==1
    
    if strcmp(fun.model,'single')==1
        for ll=1:size(dipole_loc_stim,1 )
            signal_stim(ll,:) = (1-fun.noise_perc_l)*fun.scale*sin(40*time*2*pi)+(fun.noise_perc_l+0.00001)*randn(1,nlength);
        end
        
    else
        
        qq=1;
        for ll=1:size(dipole_stim_left,2 )
            signal_stim(qq,:) = (1-fun.noise_perc_l)*fun.scale*sin(40*time*2*pi)+(fun.noise_perc_l+0.00001)*randn(1,nlength);
            qq=qq+1;
        end
        
        for ll=1:size(dipole_stim_right,2 )
            signal_stim(qq,:) = (1-fun.noise_perc_r)*fun.scale*sin(40*time*(2*pi) - fun.phase) + (fun.noise_perc_r+0.00001)*randn(1,nlength);
            qq=qq+1;
        end
        
    end
    
elseif strcmp(fun.signal,'gauss')==1
   
    if strcmp(fun.model,'single')==1
        for ll=1:size( dipole_loc_stim,1 )
            signal_stim(ll,:) = (1-fun.noise_perc)*fun.scale*x_gauss+(fun.noise_perc+0.00001)*randn(1,length(x_gauss));
        end
        
    else
        
%         x1=1/(sqrt(2*pi*variance))*(exp(-(ts-0.4).^2/(2*variance)));
%         x2=-1/(sqrt(2*pi*variance))*(exp(-(ts-0.55).^2/(2*variance)));
%         x_gauss = x1+x2;
%         x_gauss=x_gauss/max(x_gauss);
%         
%         qq=1;
%         for ll=1:length( dipole_stim_left )
%             
%             signal_stim(qq,:) = (1-fun.noise_perc)*fun.scale*x_gauss+(fun.noise_perc+0.00001)*randn(1,length(x_gauss));
%             qq=qq+1;
%         end
%          
%         x1=1/(sqrt(2*pi*variance))*(exp(-(ts-(0.4+fun.lag)).^2/(2*variance)));
%         x2=-1/(sqrt(2*pi*variance))*(exp(-(ts-(0.55+fun.lag)).^2/(2*variance)));
%         x_gauss = x1+x2;
%         x_gauss=x_gauss/max(x_gauss);
%         
%         for ll=1:length( dipole_stim_right )
%             signal_stim(qq,:) = (1-fun.noise_perc)*fun.scale*x_gauss+(fun.noise_perc+0.00001)*randn(1,length(x_gauss));
%             qq=qq+1;
%         end
        
        for ll=1:size( dipole_loc_stim,1 )
            r = (0.6-0.2).*rand(1,1) + 0.2;
            x1=1/(sqrt(2*pi*variance))*(exp(-(ts-r).^2/(2*variance)));
            r = (0.75-0.55).*rand(1,1) + 0.55;
            x2=-1/(sqrt(2*pi*variance))*(exp(-(ts-r).^2/(2*variance)));
            
            x_gauss = x1+x2;
            x_gauss=x_gauss/max(x_gauss);
            signal_stim(ll,:) = (1-fun.noise_perc)*fun.scale*x_gauss+(fun.noise_perc+0.00001)*randn(1,length(x_gauss));
        end
    end
    
else
    disp('Error: Put sin or gauss for signal');
end

cfg      = [];
cfg.elec = fun.elec;
cfg.headmodel = fun.headmodel;
cfg.dip.pos = fun.grid.pos(dipole_loc_stim,:);
dipole_no = length(cfg.dip.pos);
cfg.fsample = 1000; % Hz
cfg.dip.mom = [dip_mom_stim]';
cfg.dip.signal = {[signal_stim]}'; % 1 trial only

raw_stim = ft_dipolesimulation(cfg);

raw_stim.trial{1} = raw_stim.trial{1} + fun.noise_perc_sensor*randn(size(raw_stim.trial{1}));

cfg=[];
cfg.channel=fun.elec.label;
% cfg.reref='yes';
% cfg.refchannel='all'
% cfg.method='avg';
data_stim=ft_preprocessing(cfg, raw_stim);

% fig_stim=figure,plot(data_stim.time{1}, data_stim.trial{1});

