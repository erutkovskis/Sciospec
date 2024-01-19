function SciospecData=fnc_read_SciospecData(fname)

% Updated on 2024.01.18
% parse fname as a filename of the *.eit file
% Edvards Rutkovskis: data format for EIT v2.1.3 - 18 header lines 

% Updated on 2020.06.06

formatIn = 'yyyy.mm.dd. HH:MM:SS.FFF';

% Date and time are saved as 'datenum'. To get exppress 'datenum' as
% seconds, multiply it by '24*3600'

%i=sqrt(-1);
% iter=0;
iter_Ia=0;
iter_Ib=0;

fid = fopen(fname, 'r');
tline = fgetl(fid);
numSkip=str2double(tline);

switch numSkip
    case 9
        numStart = 2;
        SciospecData.Version = 1;
    case 11
        numStart = 3;
        tline = fgetl(fid);
        SciospecData.Version = str2double(tline);
    case 18
        tline = fgetl(fid); % get data format version - most likely "2"
        SciospecData.Version = str2double(tline);
        numStart = 3; % line number to start???
end


iter = numStart;
tline = fgetl(fid); % get setup name

SciospecData.Name=tline; 
iter=iter+1; tline = fgetl(fid); % get recording date

SciospecData.Date=tline;
try
    SciospecData.Datenum = datenum(tline,formatIn);
catch ME
    SciospecData.Datenum = datenum(tline,'yyyy.mm.dd. HH:MM:SS');
end
iter=iter+1; tline = fgetl(fid); % get first injection frequency

Freq_min=str2double(tline);
SciospecData.Min_Freq=[tline ' Hz'];
iter=iter+1; tline = fgetl(fid); % get last injection frequency

Freq_max=str2double(tline);
SciospecData.Max_Freq=[tline ' Hz'];
iter=iter+1; tline = fgetl(fid); % frequency scale log/lin
        
Freq_isLog=str2double(tline);
iter=iter+1; tline = fgetl(fid); % frequency count

NofFreq=str2double(tline);
NoS=1+NofFreq;
iter=iter+1; tline = fgetl(fid); % injection amplitude, Amps

SciospecData.Amplitude=[tline ' A'];
iter=iter+1; tline = fgetl(fid); % framerate

SciospecData.FrameRate=[tline ' Frames/s'];
iter=iter+1; tline = fgetl(fid); % phase correction parameter

if(numSkip==11)
    SciospecData.PhaseCorrection=str2double(tline);
    iter=iter+1; tline = fgetl(fid);    
end

if (numSkip==18)
    % additional fields for software v2.1.3
    iter=iter+1; tline = fgetl(fid); % get gain setting, assume just a single number for now
    SciospecData.GainSetting = str2double(tline);

    iter=iter+1; tline = fgetl(fid); % get ADC range
    SciospecData.ADCRange = str2double(tline);

    iter=iter+1; tline = fgetl(fid); % get measure mode
    SciospecData.MeasureMode = str2double(tline);

    iter=iter+1; tline = fgetl(fid); % get boundary (?)
    SciospecData.Boundary = str2double(tline);

    iter=iter+1; tline = fgetl(fid); % get switch type
    SciospecData.SwitchType = str2double(tline);

    iter=iter+1; tline = fgetl(fid); % get measurements channels
    SciospecData.MeasChannels = str2double(regexp(tline,'\d*','match'));

    iter=iter+1; tline = fgetl(fid); % get MeasurementChannelsIndependentfromInjectionPattern (?) 
    SciospecData.MeasChanIndepFromInjPat = str2double(regexp(tline,'\d*','match'));
    
    iter=iter+1; tline = fgetl(fid); % proceed to the main part of the .eit file
end

num = iter-1;
while ischar(tline)
   
    NofV=iter-num;
    
    if mod(NofV,NoS)==1
        iter_Ia=iter_Ia+1;
        temp_in=str2num(tline);
        Injection_setting(iter_Ia,1:2)=temp_in;
    else
        iter_Ib=iter_Ib+1;
        temp_V=str2num(tline);
        VoltageData(iter_Ib,:)=temp_V(1:2:end)+1i*temp_V(2:2:end);
    end
    iter=iter+1; tline = fgetl(fid);    
    %disp(iter)
end
fclose(fid);


if ~Freq_isLog
    SciospecData.Frequencies=linspace(Freq_min,Freq_max,NofFreq);
else
    % disp('Add the codes for log!!')
%     a=NofFreq^(1/(Freq_max-Freq_min)); b=Freq_min;
%     SciospecData.Frequencies=log(1:NofFreq)/log(a)+b;
%     a=NofFreq^(1/(Freq_max-Freq_min)); b=Freq_min;
    SciospecData.Frequencies=(Freq_max/Freq_min).^(linspace(0,1,NofFreq))*Freq_min;
end

Nof_ij=size(Injection_setting,1);
for kk=1:NofFreq
    idxV=kk:NofFreq:NofFreq*Nof_ij;
    SciospecData.Voltages(kk).voltage=VoltageData(idxV,:);
    SciospecData.Voltages(kk).frequency=SciospecData.Frequencies(kk);
end

SciospecData.Injection_setting=Injection_setting;


% SciospecData.Voltages(kk).voltage(i,j) : voltage measured between j-th electrode and the ground electrdoe subject to the i-th current injection.
% SciospecData.Voltages(kk).voltage(i,:) : i-th row of the matrix is the measured voltages with respect to i-th injection current


