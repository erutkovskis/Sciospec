clear all; close all
%% Start-up eidors 
run 'C:\Users\erutkovs\OneDrive - University College London\PhD sEIT\EIT general\eidors-v3.11-ng\eidors\startup.m'; % start-up eidors 
%% Load Sciospec EIT data

% Reference
fpath='C:\Users\erutkovs\OneDrive - University College London\PhD sEIT\EIT general\Sciospec comparison\tank_EIT_03062024_v0.3.28\20240603 15.50.14 no pert\setup\';


for i=1:190
    fname=['setup_' sprintf('%05d',i) '.eit'];
    FrameAll(i)=fnc_read_SciospecData_v2(fullfile(fpath,fname));
    VoltageRef_16chan = FrameAll(i).Voltages.voltage(:,(1:16));
    VoltageRef_temp(:,i) = reshape(VoltageRef_16chan,[],1); % reshape into a single column vector
    %VoltageRef_temp(:,i)=FrameAll(i).Voltages.voltage(:,(1:16));
end


% Data
fpath='C:\Users\erutkovs\OneDrive - University College London\PhD sEIT\EIT general\Sciospec comparison\tank_EIT_03062024_v0.3.28\20240603 15.58.26 with pert\setup\';

for i=1:302
    fname=['setup_' sprintf('%05d',i) '.eit'];
    FrameAll(i)=fnc_read_SciospecData_v2(fullfile(fpath,fname));
    VoltageAnoMany_16chan = FrameAll(i).Voltages.voltage(:,(1:16));
    VoltageAnoMany_temp(:,i) = reshape(VoltageAnoMany_16chan,[],1);

    %VoltageAnoMany_temp(:,i)=FrameAll(i).Voltages.voltage(:);
end

disp(['Injected Current amplitude : ' FrameAll(1).Amplitude])
%amplitude = str2double(regexp(FrameAll(1).Amplitude,'\d*','Match'));
amplitude = 0.001; % A

for k=1:190
    VoltageRef(:,:,k)=reshape(VoltageRef_temp(:,k),16,16);
end

for k=1:302
    VoltageAnoMoving(:,:,k)=reshape(VoltageAnoMany_temp(:,k),16,16);
end


%% Convert Sciospec data to EIT data
NChannel=16;
NSkip=0;


for k=1:190
V=VoltageRef(:,:,k);
Veit=func_ConvertSciospecToEIT(V',NChannel,NSkip,false);
VeitRef(:,k)=Veit;
end


for k=1:302
V=VoltageAnoMoving(:,:,k);
Veit=func_ConvertSciospecToEIT(V',NChannel,NSkip,false);
VeitAnoMoving(:,k)=Veit;
end


v_all = real(VeitAnoMoving);
v_ref = real(VeitRef(:,2)); % reference frame for time-difference imaging
%% Make foward and invese models using eidors
imdl2D = mk_common_model('b2c',16); % use this model for dummy

fmdl2D = imdl2D.fwd_model;
figure(1);show_fem(imdl2D.fwd_model);
title('FEM model')
[fmdl2D.stimulation,fmdl2D.meas_select] = mk_stim_patterns(16,1,'{ad}','{ad}',{'no_meas_current'},amplitude);
%% Image reconstruction
img2D= inv_solve_diff_GN_one_step(imdl2D, v_ref, v_all);

img2D.calc_colours.ref_level=0;
img2D.type='image';
img2D.show_slices.img_cols=0;
img2D.calc_colours.ref_level=0;

figure(2);show_slices(img2D);
title('Difference Imaging')
%% Image reconstruction (Movie)
img2D= inv_solve_diff_GN_one_step(imdl2D, v_ref, v_all);
elem_data_all=-img2D.elem_data;

figure(3);
for i = 30:size(v_all,2)    
    img2D.elem_data=elem_data_all(:,i);
    figure(3);
    show_slices(img2D);        
    title({['Difference Imaging'] ; [num2str(i) 'th frame']})    
    pause(0.01)
end