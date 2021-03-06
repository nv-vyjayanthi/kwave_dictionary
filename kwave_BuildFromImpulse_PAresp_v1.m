% master build of dictionary
% spacing will be set at 1e-6m
clear;
addpath('~/Documents/MATLAB/k-wave-toolbox-version-1.2.1/k-Wave/');

SaveImpulse = 'no';
GenerateImpulse = 'load';
Dir = '~/Documents/MATLAB/Programming/BlindSpeckleIllumPA/kwave_dict_ImpulseResp/';
SaveFileName = [Dir '/ImpulseResponse_3.mat'];

%% create computational grid and get impusle response

if strcmp(GenerateImpulse,'yes')
    
    Nx = 32;
    dx = 1e-6;
    Ny = Nx; dy = dx;
    % Nz = Nx; dz = dx;
    Nz = 128; dz = dx;

    kgrid = kWaveGrid(Nx,dx,Ny,dy,Nz,dz);

    % define the properties of the propagation medium
    medium.sound_speed = 1500;	% [m/s]

    % create the time array
    t_end = 2.5e-6;
    kgrid.t_array = makeTime(kgrid, medium.sound_speed, [], t_end);
    % kgrid.makeTime(medium.sound_speed);

    % define grid sensor
    n_sensors = 1;
    det_pos = [Nx/2+1 Ny/2+1 10];
    sensor.mask = zeros([Nx,Ny,Nz]);
    sensor.mask(det_pos(1),det_pos(2),det_pos(3)) = 1; % pos in z of detector a constant

    % input arguments
    input_args = {'DataCast', 'single', 'CartInterp', 'nearest', 'DisplayMask', 'off', ...
        'PlotSim', boolean(0), 'PMLInside',boolean(0)};

    % create initial pressure distribution 
    impulse_pos = [Nx/2+1 Ny/2+1 Nz/2+1];
    p01 = zeros([Nx Ny Nz]);
    p01(impulse_pos(1),impulse_pos(2),impulse_pos(3)) = 1; % this will be const.
    source.p0 = p01;
    sensor_data_impulse = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

    if strcmp(SaveImpulse,'yes')
        save(SaveFileName);
    end
    
elseif strcmp(GenerateImpulse,'load')   
    tic;
    load(SaveFileName);
    disp(['Impulse response loaded in ' num2str(toc) 's']);
    
end
%% define object to be built

obj_mag = 1; 

% alter this for x-z objects
no_pts = 6;
obj_pos_x = (Nx/2+1)*ones([1 no_pts]);
obj_pos_y = (Ny/2+1)*ones([1 no_pts]);
obj_pos_z = [36 58 73 85 94 100];
% obj_pos_z = [65 65+6];

% alter this for x-y objects

% alter this for 3d objects

p00_new = zeros([Nx Ny Nz]);
p00_new(obj_pos_x, obj_pos_y, obj_pos_z) = obj_mag;

%% calculate shifting and scaling constants

time_shifts = zeros([1 no_pts]);
scaling_consts = zeros([no_pts kgrid.Nt]);
shifted_waveforms = zeros([no_pts n_sensors kgrid.Nt]);
scaling_factors = zeros([no_pts kgrid.Nt]);
sensor_data_pts = zeros([no_pts n_sensors kgrid.Nt]);

vst_array = medium.sound_speed*kgrid.t_array;
r0 = sqrt((impulse_pos(1)-det_pos(1))^2+(impulse_pos(2)-det_pos(2))^2+(impulse_pos(3)-det_pos(3))^2)*kgrid.dz;

for i_det = 1:n_sensors
    for i_pt = 1:no_pts
        
        % calc time shifts
        r_0pt = sqrt((impulse_pos(1)-obj_pos_x(i_pt))^2+(impulse_pos(2)-obj_pos_y(i_pt))^2+(impulse_pos(3)-obj_pos_z(i_pt))^2)*kgrid.dz;
        r_deltaT = sign(obj_pos_z(i_pt)-impulse_pos(3))*(r_0pt/medium.sound_speed);
        time_shifts(i_pt) = ceil((r_deltaT/kgrid.dt));
        
        % shift the waveform
        shifted_data = circshift(sensor_data_impulse(i_det,:),time_shifts(i_pt));
        shifted_waveforms(i_pt,i_det,:) = shifted_data;
        
        % calc scaling factor
        r_det_pt = sqrt((obj_pos_x(i_pt)-det_pos(1))^2+(obj_pos_y(i_pt)-det_pos(2))^2+(obj_pos_z(i_pt)-det_pos(3))^2)*kgrid.dz;
        scaling_factors(i_pt,:) = (1-vst_array/r_det_pt)./(1-((vst_array-r_deltaT*medium.sound_speed)/r0));
        
    end
end

scaling_factors(find(scaling_factors==Inf)) = 1;
scaling_factors(find(isnan(scaling_factors))) = 1;

% superposition to get final waveform
for i_det = 1:n_sensors
    for i_pt = 1:no_pts
        shifted_waveform = squeeze(shifted_waveforms(i_pt,i_det,:)).';
        scaling_const = scaling_factors(i_pt,:);
        sensor_data_pts(i_pt,i_det,:) = scaling_const.*shifted_waveform.*obj_mag;
    end
end

sensor_data_obj = reshape(sum(sensor_data_pts,1),[n_sensors kgrid.Nt]);

%% generate speckle patterns

N_speckle = 200;
kspeckle = 12*dx; 
l = 1/kspeckle; % length of phase patterns
Nlx = round(Nx*dx*l); % [grid points]; this size is defined in the Fourier domain
Nly = round(Ny*dy*l); % [grid points]; this size is defined in the Fourier domain
Nlz = round(Nz*dz*l); % [grid points]; this size is defined in the Fourier domain

speckle_patterns = zeros([Nx Ny Nz N_speckle]);
phase_patterns = zeros([Nx Ny N_speckle]);
p00_new_speckle = zeros([Nx Ny Nz N_speckle]);
sensor_data_speckle_pts = zeros([no_pts n_sensors kgrid.Nt N_speckle]);
speckle_patterns_long = zeros([10*Nz N_speckle]);

% generting non-uniform speckle in x-y-z
for i_speckle = 1:N_speckle

    SM_3d_long = Generate_3D_Speckle(Nx, Ny, 10*Nz, Nlx, Nly, 10*Nlz);
    SM_3d = SM_3d_long(:,:,round(size(SM_3d_long,3)/2-Nz/2)+1:round(size(SM_3d_long,3)/2+Nz/2),:);
    speckle_patterns(:,:,:,i_speckle) = SM_3d;  
    speckle_patterns_long(:,i_speckle) = SM_3d_long(1,1,:); % take only 1 instance in x,y to measure size  
    p00_new_speckle(:,:,:,i_speckle) = p00_new.*SM_3d;
    
end

[speckle_FWHM,speckle_cov] = SpeckleSize1d(speckle_patterns_long);

% figure;
% plot(speckle_cov); title({'Speckle Covariance function',['Speckle Width ' num2str(speckle_FWHM*kgrid.dz/1e-6) '\mum']})

for i_speckle = 1:N_speckle
    for i_pt = 1:no_pts
            sensor_data_speckle_pts(i_pt,:,:,i_speckle) = sensor_data_pts(i_pt,:,:)*p00_new_speckle(obj_pos_x(i_pt),obj_pos_y(i_pt),obj_pos_z(i_pt),i_speckle);
    end
end

sensor_data_speckle_obj = reshape(sum(sensor_data_speckle_pts,1),[n_sensors kgrid.Nt N_speckle]);

%% add low pass filter reponse 

% keeping everything in MHz generate filter coeff
fs = 1/kgrid.dt/1e6;
nyquist_freq = fs/2;
cutoff_freq = 30;
wn = cutoff_freq/nyquist_freq;
[n_coeff,d_coeff] = butter(2,wn,'low');

% cut-off data; zero pad it; filter it
t_cutoff = 800;
zero_pad = 400;
impulse_resp = double(sensor_data_impulse(1:t_cutoff));
impulse_resp = [zeros([1 zero_pad]) impulse_resp];
impulse_resp_filt = filtfilt(n_coeff,d_coeff,impulse_resp);

obj_resp = sensor_data_obj(1:t_cutoff);
obj_resp = [zeros([1 zero_pad]) obj_resp];
obj_resp_filt = filtfilt(n_coeff,d_coeff,obj_resp);

speckle_resp = squeeze(sensor_data_speckle_obj(1,1:t_cutoff,:));
speckle_resp = [zeros([zero_pad N_speckle]); speckle_resp];

speckle_resp_filt = zeros(size(speckle_resp));
for i_speckle = 1:N_speckle
    speckle_resp_filt(:,i_speckle) = filtfilt(n_coeff,d_coeff, speckle_resp(:,i_speckle));
end

% create new time array
t_array = kgrid.t_array(1:t_cutoff);
t_array = [kgrid.dt.*(-zero_pad:1:-1) t_array];

% add noise to the low pass filtered signals
noise_level = 10;
obj_resp_filt = obj_resp_filt.*(1 + (noise_level/100)/sqrt(N_speckle)*randn(size(obj_resp_filt)));
speckle_resp_filt = speckle_resp_filt.*(1 + (noise_level/100)*randn(size(speckle_resp_filt)));

% % no intragrating filtered responses
% impulse_resp_int = impulse_resp_filt;
% obj_resp_int = obj_resp_filt;
% speckle_resp_int = speckle_resp_filt;

% integrate filtered responses
impulse_resp_int = cumtrapz(t_array,impulse_resp_filt);
obj_resp_int = cumtrapz(t_array,obj_resp_filt);
speckle_resp_int = cumtrapz(t_array,speckle_resp_filt,1);
% keeping only +ve values from integrated response
impulse_resp_int(impulse_resp_int<0) = 0;
obj_resp_int(obj_resp_int<0) = 0;
speckle_resp_int(speckle_resp_int<0) = 0;

% % taking absolute value of integrated sgls
% impulse_resp_int = abs(impulse_resp_int);
% obj_resp_int = abs(obj_resp_int);
% speckle_resp_int = abs(speckle_resp_int);

% offset integrated values so that everything is +ve 
% impulse_resp_int = impulse_resp_int + abs(min(impulse_resp_int));
% obj_resp_int = obj_resp_int + abs(min(obj_resp_int));
% speckle_resp_int = speckle_resp_int + abs(min(speckle_resp_int));

%% visualize data

figure;
subplot(2,2,1); plot(t_array/1e-6,impulse_resp); xlabel('t \mus'); title('Impulse response'); xlim([t_array(1)/1e-6 t_array(end)/1e-6])
subplot(2,2,2); plot(kgrid.z_vec/1e-6,squeeze(p00_new(Nx/2+1,Ny/2+1,:))); xlabel('z \mum');title('Object')
subplot(2,2,3); plot(t_array/1e-6,obj_resp); xlabel('t \mus');title('PA resp unif'); xlim([t_array(1)/1e-6 t_array(end)/1e-6])
subplot(2,2,4); plot(t_array/1e-6,mean(speckle_resp,2)); xlabel('t \mus');
                title(['Mean speckle resp N=' num2str(N_speckle)]); xlim([t_array(1)/1e-6 t_array(end)/1e-6])
suptitle('Un-Filtered Responses')

% figure;
% subplot(2,2,1); plot(t_array/1e-6,impulse_resp_filt); xlabel('t \mus'); title('Impulse response'); xlim([t_array(1)/1e-6 t_array(end)/1e-6])
% subplot(2,2,2); plot(kgrid.z_vec/1e-6,squeeze(p00_new(Nx/2+1,Ny/2+1,:))); xlabel('z \mum');title('Object')
% subplot(2,2,3); plot(t_array/1e-6,obj_resp_filt); xlabel('t \mus');title('PA resp unif'); xlim([t_array(1)/1e-6 t_array(end)/1e-6])
% subplot(2,2,4); plot(t_array/1e-6,mean(speckle_resp_filt,2)); xlabel('t \mus');
%                 title(['Mean speckle resp N=' num2str(N_speckle)]); xlim([t_array(1)/1e-6 t_array(end)/1e-6])
% suptitle('Filtered Responses')
                
% figure;
% subplot(1,3,1); plot(t_array/1e-9,impulse_resp_int); title('Impulse response'); 
% subplot(1,3,2); plot(t_array/1e-9,obj_resp_int); title('PA resp unif'); 
% subplot(1,3,3); plot(t_array/1e-9,mean(speckle_resp_int,2)); title('Speckle response mean'); 
% suptitle('Filtered&Integrated Responses')

figure;
plot(speckle_cov); title({'Speckle Covariance function',['Speckle Width ' num2str(speckle_FWHM*kgrid.dz/1e-6) '\mum'], ...
    ['Speckle Width ' num2str(speckle_FWHM*kgrid.dz/1e-6/1.5) 'ns']})



