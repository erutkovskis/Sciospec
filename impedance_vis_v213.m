
clc;
close all;
clear all;

%% Get info from the setup file

fileID = fopen('./setup/setup.setUp');

% Skip useless lines
for iRead=1:6
    tline = fgetl(fileID);
end
tline = fgetl(fileID); % get software version
tline = strsplit(tline,':'); % split at colon
settings.softVersion = tline{2};
% Skip useless lines
for iRead = 8:22
    tline = fgetl(fileID);
end
tline = fgetl(fileID); % get channel order
tline = strsplit(tline,':');
chanOrder = strsplit(tline{2},',');
chanOrder = cell2mat(chanOrder);
chanOrder = str2double(regexp(chanOrder,'\d*','match'));
settings.chanOrder = chanOrder;
tline = fgetl(fileID); % Skip useless lines 
tline = fgetl(fileID);
tline = fgetl(fileID);
tline = fgetl(fileID);

% get current excitation pattern
for iRead = 1:length(chanOrder)
    tline = fgetl(fileID); 
    tline = strsplit(tline,',');
    tline = cell2mat(tline);
    CurrentExcitationPattern(iRead,:) = str2double(regexp(tline,'\d*','match'));
end
%CurrentExcitationPattern(:,2) = [];


fclose(fileID);
%% Read all voltage data from every injection file
files=dir('./setup/setup_*.eit');
%files={files.name};

% Measurement file structure: rows - file number, columns - injection pair
% number

injPair = cell(length(files),length(CurrentExcitationPattern));
measVolt = cell(length(files),length(CurrentExcitationPattern));
% scan files
for iFiles = 1:length(files)
    
    iFiles % show processed file number
    
    filepath = strjoin({files(iFiles).folder files(iFiles).name},'\');
    fileID = fopen(filepath);
    % skip useless lines
    % for iRead = 1:11
    %     tline = fgetl(fileID);
    % end
    headerRows = fgetl(fileID);
    %for iRead = 1:headerRows-1
    [EIT{iFiles}.fileVersion] = str2double(fgetl(fileID));
    [EIT{iFiles}.datasetName] = fgetl(fileID);
    [EIT{iFiles}.timeStamp] = fgetl(fileID);
    [EIT{iFiles}.freqMin] = str2double(fgetl(fileID));
    [EIT{iFiles}.freqMax] = str2double(fgetl(fileID));
    [EIT{iFiles}.freqScale] = str2double(fgetl(fileID)); % 0 lin, 1 log
    [EIT{iFiles}.freqCount] = str2double(fgetl(fileID));
    [EIT{iFiles}.amplitude] = str2double(fgetl(fileID));
    [EIT{iFiles}.framerate] = str2double(fgetl(fileID));
    [EIT{iFiles}.phaseCorrection] = str2double(fgetl(fileID));
    [EIT{iFiles}.gain] = str2double(fgetl(fileID));
    [EIT{iFiles}.ADCrange] = str2double(fgetl(fileID));
    [EIT{iFiles}.measureMode] = str2double(fgetl(fileID));
    [EIT{iFiles}.boundary] = str2double(fgetl(fileID));
    [EIT{iFiles}.switchType] = str2double(fgetl(fileID));
    [EIT{iFiles}.measChannels] = str2double(fgetl(fileID));
    [EIT{iFiles}.measChannels_Independent] = str2double(fgetl(fileID));

    iInj = 1;
    while ~feof(fileID) % scan the file until the end
        % injPair{iFiles,iInj} = fgetl(fileID); % get the injection pair
        % injPair{iFiles,iInj} = str2double(regexp(injPair{iFiles,iInj},'\d*','match')); 
        % 
        % measVolt{iFiles,iInj} = fgetl(fileID); % Get measured voltages per injected pair
        % measVolt{iFiles,iInj} = str2double(strsplit(measVolt{iFiles,iInj},'\t'));
        
        EIT{iFiles}.injPair(iInj,:) = str2double(regexp(fgetl(fileID),'\d*','match'));
        EIT{iFiles}.measVolt(iInj,:) = str2double(strsplit(fgetl(fileID),'\t'));
        EIT{iFiles}.real(iInj,:) = EIT{iFiles}.measVolt(iInj,[1:2:64]);
        EIT{iFiles}.imag(iInj,:) = EIT{iFiles}.measVolt(iInj,[2:2:64]);
        EIT{iFiles}.abs(iInj,:) = sqrt(EIT{iFiles}.real(iInj,:).^2 + EIT{iFiles}.imag(iInj,:).^2);
        iInj = iInj +1;
    end

    fclose(fileID);
end

%% Visualise voltage

figure(1);
bar(1:16,EIT{1,1}.abs(:,1:16));
xlabel("Injection"); ylabel("Voltage, V")

%% Calculate noise