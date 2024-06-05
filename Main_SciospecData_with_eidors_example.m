clear all; close all
%% Start-up eidors 
run '../../../eidors-v3.10/eidors/startup.m'; % start-up eidors 
%% Load Sciospec EIT data

% Reference
fpath='Reference/20200708 14.32.08';


for i=1:30
    fname=['Frame_' num2str(i) '.eit'];
    FrameAll(i)=fnc_read_SciospecData(fullfile(fpath,fname));
    
    VoltageRef_temp(:,i)=FrameAll(i).Voltages.voltage(:);
end


% Anomaly
fpath='Anomaly/20200708 14.32.28';

for i=1:180
    fname=['Frame_' num2str(i) '.eit'];
    FrameAll(i)=fnc_read_SciospecData(fullfile(fpath,fname));
    
    VoltageAnoMany_temp(:,i)=FrameAll(i).Voltages.voltage(:);
end

disp(['Injected Current amplitude : ' FrameAll(1).Amplitude])
amplitude = 0.01;


for k=1:30
    VoltageRef(:,:,k)=reshape(VoltageRef_temp(:,k),16,16);
end

for k=1:180
    VoltageAnoMoving(:,:,k)=reshape(VoltageAnoMany_temp(:,k),16,16);
end


%% Convert Sciospec data to EIT data
NChannel=16;
NSkip=0;


for k=1:30
V=VoltageRef(:,:,k);
Veit=func_ConvertSciospecToEIT(V',NChannel,NSkip,false);
VeitRef(:,k)=Veit;
end


for k=1:180
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
[fmdl2D.stimulation,fmdl2D.meas_select] = mk_stim_patterns(16,1,[1,0],[0,1],{'no_meas_current'},amplitude);
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