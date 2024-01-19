
clc;
close all;
clear all;

%% Get info from setup file

fileID = fopen('setup.setUp');

% Skip useless lines
for iRead=1:7
    tline = fgetl(fileID);
end

% Measurement time in milliseconds - each injection
tline = fgetl(fileID);
tline=strsplit(tline,':');
settings.measTime=str2double(tline{2});

tline = fgetl(fileID); % Skip useless lines

% Using EIT
tline = fgetl(fileID);
tline=strsplit(tline,':');
settings.useEIT=strcmp(strtrim(tline{2}),'true');

% Using AUX data
tline = fgetl(fileID);
tline=strsplit(tline,':');
settings.useAUX=strcmp(strtrim(tline{2}),'true');

tline = fgetl(fileID); % Skip useless lines

% Using Time data
tline = fgetl(fileID);
tline=strsplit(tline,':');
settings.useTD=strcmp(strtrim(tline{2}),'true');

tline = fgetl(fileID); % Skip useless lines
tline = fgetl(fileID); % Skip useless lines

% Using Second Stimulator (i.e. bipolar port)
tline = fgetl(fileID);
tline=strsplit(tline,':');
settings.useStim2=strcmp(strtrim(tline{2}),'true');

if settings.useStim2
    
    % Second Stimulator Parameters
    
    tline = fgetl(fileID);
    tline=strsplit(tline,':');
    tline=strsplit(tline{2},',');
    tline=str2double(tline);
    
    settings.stim2.pulsewidth=tline(1)*1e6;   % Pulse width in microseconds
    settings.stim2.amplitude=tline(2);   % Pulse amplitude in A
    settings.stim2.frequency=tline(3);   % Repetition frequency in Hz
    
end

% Skip useless lines
for iRead=1:11
    tline = fgetl(fileID);
end

% Read injection pairs until you get to EIT settings line
isSetting=false;
prt=[];
prtLine=0;

while (isSetting==false)
    
    tline = fgetl(fileID);
    temp=strsplit(tline,':');
    
    if strcmp(temp{1},'Settings')
        settingsLine1=temp{2};
        settingsLine2=fgetl(fileID);
        isSetting=true;
    else
        tline = strsplit(tline,',');        
        prtLine=prtLine+1;        
        prt(prtLine,:)=[str2double(tline{1}) str2double(tline{2})] ;        
    end
    
end

settings.prt=prt;
settings.EITsettings=strsplit([settingsLine1 ',' settingsLine2],',');

fclose(fileID);

save expsetup prt settings


%% Read data from every injection file

files=dir('setup_Inj(*)_*.neitb');
files={files.name};

filesOrder=files;

numChAUX=8;
numChTD=32;

AUX=[];
TD=[];
EIT=[];

for iFiles=1:length(files)
    
    iFiles

    filename=files{iFiles};
    tempstr=strsplit(filename,'_');
    tempstr=tempstr{end};
    tempstr=strsplit(tempstr,'.neitb');
    tempstr=tempstr{1};
    iInj=str2double(tempstr)+1;

    AUX{iInj}.data=[];
    TD{iInj}.data=[];
    EIT{iInj}.real=[];
    EIT{iInj}.imag=[];
    
    fileID = fopen(filename,'r','b');   % Big endian, ID=3 or higher, open succeeded
    
    for i=1:1e10
        
        [header, byteCount]=fread(fileID, 1, '*uint8');
        
        if byteCount>0
            
            switch header
                
                case 1 % AUX Header
                    fread(fileID, 1, '*uint8');
                    AUX{iInj}.freq = fread(fileID, 1, 'float32');
                    fread(fileID, 1, '*uint8');
                    
                case 3  % TD Header
                    fread(fileID, 1, '*uint8');
                    TD{iInj}.freq = fread(fileID, 1, 'float32');
                    fread(fileID, 1, '*uint8');
                    
                case 2  % AUX Data
                    tempHeader=fread(fileID, 1, '*uint8');
                    numFramesAUX=floor(tempHeader/32);  % Divide by 8 chs, 4 bytes each
                    frameAUX = fread(fileID, numFramesAUX*numChAUX, 'float32'); % Do something with this
                    AUX{iInj}.data=[AUX{iInj}.data;frameAUX'];
                    fread(fileID, 1, '*uint8');  % Just read the final "02" value of auxiliary data
                    
                case 4  % TD Data
                    tempHeader=fread(fileID, 1, '*uint8');
                    numFramesTD=floor(tempHeader/4/numChTD);
                    frameTD = fread(fileID, numFramesTD*numChTD, 'float32');
                    TD{iInj}.data=[TD{iInj}.data;frameTD'];
                    fread(fileID, 1, '*uint8');
                    
                case 5   % EIT Header
                    fread(fileID, 1, '*uint8');
                    [EIT{iInj}.fmin] = fread(fileID, 1, 'float32');
                    [EIT{iInj}.fmax] = fread(fileID, 1, 'float32');
                    [EIT{iInj}.fscale] = fread(fileID, 1, 'float32');
                    [EIT{iInj}.fcount] = fread(fileID, 1, 'float32');
                    [EIT{iInj}.amplitude] = fread(fileID, 1, 'float32');
                    [EIT{iInj}.framerate] = fread(fileID, 1, 'float32');
                    [EIT{iInj}.phaseCorrection] = fread(fileID, 1, 'float32');    
                    [EIT{iInj}.gain] = fread(fileID, 1, 'float32');  
                    [EIT{iInj}.adcRange] = fread(fileID, 1, 'float32');  
                    [EIT{iInj}.measureMode] = fread(fileID, 1, 'float32');       
                    [EIT{iInj}.boundary] = fread(fileID, 1, 'float32');         
                    [EIT{iInj}.switchType] = fread(fileID, 1, 'float32');       
                    [EIT{iInj}.injectPlus] =  fread(fileID, 1, 'float32');       
                    [EIT{iInj}.injectMinus] = fread(fileID, 1, 'float32');       
                    fread(fileID, 1, '*uint8');
                  
                case 6 % EIT Data
                    
                    % Chs 1-16
                    
                    fread(fileID, 1, '*uint8');
                    fread(fileID, 5, '*uint8'); % Skip channel group and time stamp
                    
                    tempData=zeros(32,2);
                    for iCh=1:16
                        tempData(iCh,1)=fread(fileID, 1, 'float32');
                        tempData(iCh,2)=fread(fileID, 1, 'float32');
                    end
                    fread(fileID, 1, '*uint8');
                    
                    % Chs 17-32
                    
                    fread(fileID, 2, '*uint8');
                    fread(fileID, 5, '*uint8'); % Skip channel group and time stamp
                    
                    for iCh=17:32
                        tempData(iCh,1)=fread(fileID, 1, 'float32');
                        tempData(iCh,2)=fread(fileID, 1, 'float32');
                    end
                    fread(fileID, 1, '*uint8');
                    
                    % Add to sample list
                    
                    EIT{iInj}.real=[EIT{iInj}.real;tempData(:,1)'];
                        % Try to subtract noise from channel 1
                    EIT{iInj}.imag=[EIT{iInj}.imag;tempData(:,2)'];
                    
            end
            
            
        else
            break
            
        end
        
    end
    
    %% Close file
    
    fclose(fileID);     % Returns 0 if close succeeded
    
    %% Other stuff


  
    EIT{iInj}.abs=sqrt(EIT{iInj}.real.^2+EIT{iInj}.imag.^2);
    
    % Filename for this injection
    EIT{iInj}.filename=filename;
       
    
end


%% Get figures of merit

tWindow=1000:5000; 

BV=[];
noise=[];
Prt=[];

chs=[14	8 3	13 26 21 12 7 6	11 20 25 30	2 1];

for iInj=1:length(EIT)
   
    inj_BV=mean(EIT{iInj}.abs(tWindow,chs),1);        
    
    % Get noise in uV in the 100-1000 Hz band (makes sense for nerve)    
    fc=100;
    fs = EIT{iInj}.framerate;
    [b,a] = butter(5,fc/(EIT{iInj}.framerate/2),'high');   
    inj_noise=std(filtfilt(b,a,EIT{iInj}.abs(tWindow,chs)))*1e6;   
    
    injPos=[ find(EIT{iInj}.injectPlus==chs,1) find(EIT{iInj}.injectMinus==chs,1) ];    
    inj_Prt=[ repmat(injPos, 14, 1)  [1:14]' repmat(15,14,1) ];
    
    BV=[BV; inj_BV'];
    noise=[noise; inj_noise'];
    Prt=[Prt;inj_Prt];
    
end

save EIT_data EIT

save BV_and_noise BV noise Prt

%% Data visualisation for Martin - plot first injection
figure(1);
tWindow=1000:4000;
chs=[14	8 3	13 26 21 12 7 6	11 20 25 30	2 1];
plot(detrend(EIT{1,1}.abs(tWindow,chs)));

figure(1); title("Injection 14 - 22")
xlabel("Sample number");
ylabel("Voltage (V)");
%% Visualise standing voltages
figure(2); 
%time = 1:length(EIT{1,1}.abs);
tWindow=1000:4000;
chs=[14	8 3	13 26 21 12 7 6	11 20 25 30	2 1];
for i = 1:length(chs)
    title('Standing voltages');
    plot(tWindow,detrend(EIT{1,i}.abs(tWindow,chs)));

end
%% Visualise noise
figure(3);
tWindow=1000:4000;
chs=[14	8 3	13 26 21 12 7 6	11 20 25 30	2 1];
% for i = 1:length(chs)
%     title('Noise level');
%     bar(noise_12102023);
% 
% end
bar(noise);

figure(4);
for i = 1:length(chs)
    title('Standing voltages');
    bar(BV);

end

%% FFT 
%% Plot the amplitude spectrum


figure(5); 
hold on
%time = 1:length(EIT{1,1}.abs);
tWindow=1000:4000;
chs=[14	8 3	13 26 21 12 7 6	11 20 25 30	2 1];

Fs = fs; % Hz
T = 1/Fs; % s
L = length(EIT{1,i}.abs(tWindow,chs)); % samples
t = (0:L-1)*T; % time vector, s

[b,a] = butter(5,1/(Fs/2),'high');

for i = 1:length(chs)
    
    
    Y = fft(filtfilt(b,a,EIT{1,i}.abs(tWindow,chs)));
    P2 = abs(Y/L);
    P1 = 2*P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    f = Fs*(0:(L/2))/L;
    semilogy(f,P1) 


end

title("Single-Sided Amplitude Spectrum of X(t)")
xlabel("f (Hz)"); xlim([0 500])
ylabel("|P1(f)|")

hold off


%% PWelch

figure(6)
hold on



for i = 1:length(chs)
    
    pwelch(filtfilt(b,a,EIT{1,i}.abs(tWindow,chs))) 

end
title("PWelch")

hold off