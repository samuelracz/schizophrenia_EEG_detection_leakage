
addpath(genpath('functions'))

fnames = dir('Dataset1_processed/*.gdf');
ns = length(fnames);

subjIDs = {fnames.name};
n_hc = sum(contains(subjIDs,'hc'));
n_sz = sum(contains(subjIDs,'sz'));

subj_list = cell(ns,1);
for n = 1:ns
    subj_list{n} = fnames(n).name(1:6);
end
subj_hc = subj_list(contains(subj_list,'hc')==1);
subj_sz = subj_list(contains(subj_list,'sz')==1);

%% Specify EEG data preparation

nch = 16;
fs = 128;
chlab = {'F7','F3','F4','F8','T7','C3','Cz','C4','T8','P7','P3','Pz','P4','P8','O1','O2'};

fs_resample = 200;
TE = 5;
th_SD = 15;

T = 60;
Nt = floor(T/TE);

%% prepare and export data (per subject)

for subj = 1:ns
    [data,h] = sload([fnames(subj).folder '/' fnames(subj).name],'OVERFLOWDETECTION','off');

    % select equal amount of data for all subjects
    data = data(1:T*fs,:);
    
    if contains(fnames(subj).name,'hc')
        y_true = 0;
    elseif contains(fnames(subj).name,'sz')
        y_true = 1;
    else
        error('wrong data loaded')
    end

    % 1. upsampling to 200 Hz
    eeg_res = resample(data,fs_resample,fs);

    % 2. z-scoreing
    eeg_z = zscore(eeg_res);

    % 3. epoching
    eeg_mat = zeros(Nt,nch,TE*fs_resample);
    flag_AF = zeros(Nt,1);
    for t = 1:Nt
        t_start = (t-1)*TE*fs_resample + 1;
        t_stop = t_start + TE*fs_resample - 1;

        eeg_seg = eeg_z(t_start:t_stop,:);
        if any(abs(eeg_seg(:))>th_SD)
            flag_AF(t) = 1;
        else
            eeg_mat(t,:,:) = eeg_seg';
        end
    end

    % 4. removing outliers
    if any(flag_AF)
        eeg_mat(flag_AF==1,:,:) = [];
    end

    % structuring data
    eeg_reve_struct = struct(...
        'data', eeg_mat,...
        'electrode_names', {chlab'},...
        'subject_id', subj_list{subj},...
        'label', y_true,...
        'n_epochs', size(eeg_mat,1),...
        'fs', fs_resample,...
        'epoch_duration', TE);

    % saving data
    fname_out = [fnames(subj).name(1:end-14) '_reve_eeg.mat'];
    if ~exist('Dataset1_REVE_eeg','dir')
        mkdir('Dataset1_REVE_eeg')
    end
    save(['Dataset1_REVE_eeg/' fname_out],'eeg_reve_struct','-v7.3')
    disp([fnames(subj).name(1:end-14) ' done.'])


end
