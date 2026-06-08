
addpath(genpath('functions'))

fnames = dir('Dataset2_processed/*.gdf');
ns = length(fnames);

subjIDs = {fnames.name};
n_hc = sum(contains(subjIDs,'h0') | contains(subjIDs,'h1'));
n_sz = sum(contains(subjIDs,'s0') | contains(subjIDs,'s1'));

subj_list = cell(ns,1);
for n = 1:ns
    subj_list{n} = fnames(n).name(1:6);
end
subj_hc = subj_list(contains(subj_list,'h0') | contains(subj_list,'h1'));
subj_sz = subj_list(contains(subj_list,'h0') | contains(subj_list,'h1'));

%% Specify EEG data preparation

nch = 19;
fs = 250;
chlab = {'Fp2','F8','T8','P8','O2','Fp1','F7','T7','P7','O1','F4','C4','P4','F3','C3','P3','Fz','Cz','Pz'};

fs_resample = 200;
TE = 5;
th_SD = 15;

T = 12*60;
Nt = floor(T/TE);


%% prepare and export data (per subject)

for subj = 1:ns
    [data,h] = sload([fnames(subj).folder '/' fnames(subj).name],'OVERFLOWDETECTION','off');

    % select equal amount of data for all subjects
    data = data(1:T*fs,:);
    
    if contains(fnames(subj).name,'h0') || contains(fnames(subj).name,'h1')
        y_true = 0;
    elseif contains(fnames(subj).name,'s0') || contains(fnames(subj).name,'s1')
        y_true = 1;
    else
        error('wrong data loaded')
    end

    % 1. downsampling to 200 Hz
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
    if ~exist('Dataset2_REVE_eeg','dir')
        mkdir('Dataset2_REVE_eeg')
    end
    save(['Dataset2_REVE_eeg/' fname_out],'eeg_reve_struct','-v7.3')
    disp([fnames(subj).name(1:end-14) ' done.'])

end
