clear all; close all
%% Start-up eidors 
run 'C:\Users\erutkovs\OneDrive - University College London\PhD sEIT\EIT general\eidors-v3.11-ng\eidors\startup.m'; % start-up eidors 
%% Load Sciospec EIT data

% Without perturbation
fpath='C:\Users\erutkovs\OneDrive - University College London\PhD sEIT\EIT general\Sciospec comparison\tank_EIT_03062024_v0.3.28\20240603 15.50.14 no pert\setup\';

files=dir([fpath 'setup_*.eit']);
files={files.name};
NrOfFramesRef = length(files);

% Manually define skip pattern
NSkip=0;

for iFiles=1:NrOfFramesRef
    iFiles
    filename=files{iFiles};
    FrameAllRef(iFiles)=fnc_read_SciospecData_v2(fullfile(fpath,filename));
    VoltageRef_16chan = FrameAllRef(iFiles).Voltages.voltage(:,(1:16));
    VoltageRef_temp(:,iFiles) = reshape(VoltageRef_16chan,[],1); % reshape into a single column vector
end

for k = 1:NrOfFramesRef
    VoltageRef(:,:,k)=reshape(VoltageRef_temp(:,k),16,16);
end

% With perturbation
fpath='C:\Users\erutkovs\OneDrive - University College London\PhD sEIT\EIT general\Sciospec comparison\tank_EIT_03062024_v0.3.28\20240603 15.58.26 with pert\setup\';
files=dir([fpath 'setup_*.eit']);
files={files.name};
NrOfFramesPert = length(files);

for iFiles=1:NrOfFramesPert
    iFiles
    filename=files{iFiles};
    FrameAllPert(iFiles)=fnc_read_SciospecData_v2(fullfile(fpath,filename));
    VoltagePert_16chan = FrameAllPert(iFiles).Voltages.voltage(:,(1:16));
    VoltagePert_temp(:,iFiles) = reshape(VoltagePert_16chan,[],1); % reshape into a single column vector
end

for k = 1:NrOfFramesPert
    VoltagePert(:,:,k) = reshape(VoltagePert_temp(:,k),16,16);
end

% Get amplitude of current
disp(['Injected Current amplitude : ' FrameAllRef(1).Amplitude])
amplitude = str2double(regexp(FrameAllRef(1).Amplitude,['\d' '.' '\d*'],'Match')); % Amps

%% Convert Sciospec data to EIT data
% Get number of channels
NChannel=length(FrameAllRef(1).MeasChannels);

% Data without perturbation
for k = 1:NrOfFramesRef
    V = VoltageRef(:,:,k);
    Veit = func_ConvertSciospecToEIT(V',NChannel,NSkip,false);
    vv_wout(:,k) = Veit;
end

% Data with perturbation
for k = 1:NrOfFramesPert
    V=VoltagePert(:,:,k);
    Veit = func_ConvertSciospecToEIT(V',NChannel,NSkip,false);
    vv_with(:,k) = Veit;
end

% Get real parts of complex voltages only
vv_wout = real(vv_wout);
vv_with = real(vv_with);
%% Calculate forward model and sensitivity matrix
fmdl = mk_common_model('j2C', 16); % refined mesh
fmdl = fmdl.fwd_model; 
% Adjust for the injection and measurement patterns
fmdl.stimulation = mk_stim_patterns(16,1,'{ad}', '{ad}', {'no_meas_current'},amplitude);
img_sim = mk_image(fmdl, 1); 
% Not required:
loc = [0.543, 0.089, 0]; % define location of the perturbation
c = find_element_centres(img_sim.fwd_model);
idx = find_perturbation_indices(c,1, 0.4, loc); % Paste the location of interest

imgh = img_sim;
imgi = img_sim; % reserve for perturbed model

vv = imgi.elem_data(idx{1});
imgi.elem_data(idx{1}) = vv.*0.01; % voltage change around the perturbation

v1 = fwd_solve(imgh); % unperturbed 
v2 = fwd_solve(imgi); % perturbed
dv = v2.meas - v1.meas; % voltage diff for the difference EIT

J = calc_jacobian(img_sim); % Calculate Jacobian

%% Average measurements for diff eit
vv_wout_mean = mean(vv_wout,2); % average each measurement channel over time -> difference EIT
vv_with_mean = mean(vv_with,2); 
%% Flip signs according to simulated forward model voltages
vv_with_mean = flip_sign(vv_with_mean, v2.meas);
vv_wout_mean = flip_sign(vv_wout_mean, v1.meas);

vv_diff = vv_with_mean - vv_wout_mean; % array of difference
%% Reconstruct measurements

[S,X] = eit_recon_tik0((vv_diff'),J,logspace(-25,2,250)); % Tikhonov's EIT reconstruction using the array of difference
imgh.elem_data = S; % Assign reconstructed conductivities to the image elements
figure(1);
show_slices(imgh); 
show_fem(imgh,[0 1.02])
title("Reconstructed image")
% Threshold
figure(2);
S2 = S; S2(S > min(S)/2) = 0;
img2 = imgh;
img2.elem_data = S2;
show_slices(img2)
show_fem(img2,[0 1.02])
title("Thresholded reconstrucred image")

%% Show location of the perturbation and simulate reconstruction
[S,X] = eit_recon_tik0((dv'),J,logspace(-25,2,250)); % reconstruct simulation
figure(3);
subplot(1,2,1)
show_fem(imgi,[0 1.02]);
title("Location of simulated perturbation")
subplot(1,2,2);
imgi.elem_data = S;
show_slices(imgi);
title("Reconstructed simlated perturbation")
