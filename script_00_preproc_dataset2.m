% NOTE: in some versions, the processMARA() might crash with options set to
% [0,0,0,0,1] (no visualization, automatic rejection)
% 
% Author: F. Samuel Racz, 2025, The University of Texas at Austin
% email: fsr324@austin.utexas.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define path to the RepOD dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_raw_eeg = []; % replace with path
fnames = dir([path_to_raw_eeg '/*.edf']);
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
    EEG = pop_biosig([fnames(subj).folder '/' fnames(subj).name]);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');
    EEG = eeg_checkset( EEG );
    EEG.setname = 'Unlabeled EEG';

    % add channel labels and locations - COMPLETE WITH CORRESPONDING PATH
    EEG=pop_chanedit(EEG, {'lookup',[]},'changefield',{3,'labels','T8'},'changefield',{4,'labels','P8'},'changefield',{8,'labels','T7'},'changefield',{9,'labels','P7'},'lookup',[]);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG.setname = 'Raw EEG';

    % band-pass filter data
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',45);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','Filtered EEG','gui','off');
    EEG = eeg_checkset( EEG );

    % apply 'cleanline' to notch-filter data
    EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:19] ,'computepower',0,'linefreqs',[50 100] ,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',0,'winsize',4,'winstep',1);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname','Notch-filtered EEG','gui','off');
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


    if ~exist(path_to_pre_processed_eeg,'dir')
        mkdir(path_to_pre_processed_eeg)
    end

    if ~exist(path_to_mara_statistics,'dir')
        mkdir(path_to_mara_statistics)
    end

    % saving data
    fname_out = [fname_in(1:end-4), '_processed.gdf'];
    pop_writeeeg(EEG, [path_to_pre_processed_eeg, '\', fname_out], 'TYPE','GDF');

    % saving MARA statistics
    fname_out = [fname_in(1:end-4), '_MARA.set'];
    EEG = pop_saveset( EEG, 'filename' ,fname_out, 'filepath', path_to_mara_statistics);


end