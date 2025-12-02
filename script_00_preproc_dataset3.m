% NOTE: in some versions, the processMARA() might crash with options set to
% [0,0,0,0,1] (no visualization, automatic rejection)
% 
% Author: F. Samuel Racz, 2025, The University of Texas at Austin
% email: fsr324@austin.utexas.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define path to the SU-SZ dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_raw_eeg = []; % replace with path
fnames = dir([path_to_raw_eeg '/*.dap']);
ns = length(fnames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define output path to save pre-processed EEG and MARA statistics %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_pre_processed_eeg = []; % replace with path
path_to_mara_statistics = []; % replace with path

% EEGLAB history file generated on the 21-Oct-2025
% ------------------------------------------------
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for subj = 1:ns
    % import data from matlab workspace
    fname_in = fnames(subj).name;
    EEG = loadcurry([path_to_raw_eeg, '/', fname_in], 'KeepTriggerChannel', 'True', 'CurryLocations', 'False');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');
    EEG = eeg_checkset( EEG );
    EEG.setname = 'Raw EEG';

    % select 10-20 locations only
    EEG = pop_select( EEG, 'channel',{'FP1','FP2','F7','F3','FZ','F4','F8','T7','C3','CZ','C4','T8','P7','P3','PZ','P4','P8','O1','O2'});
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','10-20 EEG','gui','off');
    EEG = eeg_checkset( EEG );

    % band-pass filter data
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',45);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','Filtered EEG','gui','off');
    EEG = eeg_checkset( EEG );

    % apply 'cleanline' to notch-filter data
    EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:19] ,'computepower',0,'linefreqs',[50 100] ,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',0,'winsize',4,'winstep',1);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname','Notch-filtered EEG','gui','off');
    EEG = eeg_checkset( EEG );

    % resample data to 250 hz
    EEG = pop_resample( EEG, 250);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'gui','off');
    EEG = eeg_checkset( EEG );

    % run ICA
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'rndreset','yes','interrupt','on');
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG.setname = 'ICA-EEG';

    % run MARA
    options = [0,0,0,0,1];
    processMARA( ALLEEG,EEG,CURRENTSET,options);
    EEG = eeg_checkset( EEG );
    EEG.setname = 'MARA-EEG';

    % re-reference
    EEG = pop_reref( EEG, []);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'setname','Re-referenced EEG','gui','off');
    EEG = eeg_checkset( EEG );

    % create output folders if necessary
    if ~exist(path_to_pre_processed_eeg,'dir')
        mkdir(path_to_pre_processed_eeg)
    end

    if ~exist(path_to_mara_statistics,'dir')
        mkdir(path_to_mara_statistics)
    end

    % saving data
    fname_out = [fname_in(1:end-4), '_processed.gdf'];
    if contains(fname_out,'sch_1') || contains(fname_out,'sch_3')
        fname_out = ['hc' fname_out(4:end)];
    end
    pop_writeeeg(EEG, [path_to_pre_processed_eeg, '\', fname_out], 'TYPE','GDF');

    % saving MARA statistics
    fname_out = [fname_in(1:end-4), '_MARA.set'];
    if contains(fname_out,'sch_1') || contains(fname_out,'sch_3')
        fname_out = ['hc' fname_out(4:end)];
    end
    EEG = pop_saveset( EEG, 'filename' ,fname_out, 'filepath', path_to_mara_statistics);

end