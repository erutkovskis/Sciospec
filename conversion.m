clear all

files=dir('setup_*.eit');
files={files.name};
EIT_rearranged = zeros(224,length(files));

for iFiles=1:length(files)
    iFiles
    filename=files{iFiles};
    SciospecData = fnc_read_SciospecData(filename);
    %EIT(:,iFiles) = func_ConvertSciospecToEIT(SciospecData.Voltages.voltage,16,4,1);
    numChannels = size(SciospecData.Voltages.voltage,1);
    EIT_noinj_temp = zeros(numChannels);
    EIT_noinj = [];
    % Delete injection voltages
    %EIT_noinj = zeros(224,1);
    for iMeasure = 1:numChannels
        EIT_noinj_temp = SciospecData.Voltages.voltage(:,iMeasure);
        EIT_noinj_temp([SciospecData.Injection_setting(iMeasure,:)]) = [];
        EIT_noinj = vertcat(EIT_noinj,EIT_noinj_temp);
    end

    % Append EIT voltages without injection channels
    EIT_rearranged(:,iFiles) = EIT_noinj;
end
% Convert to an array of real values for Kai
EIT_rearranged_real = real(EIT_rearranged);