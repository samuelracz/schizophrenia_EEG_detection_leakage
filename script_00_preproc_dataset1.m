% NOTE: in some versions, the processMARA() might crash with options set to
% [0,0,0,0,1] (no visualization, automatic rejection)
% 
% Author: F. Samuel Racz, 2025, The University of Texas at Austin
% email: fsr324@austin.utexas.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define path to the MSU dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_raw_eeg = []; % replace with path
fnames = dir([path_to_raw_eeg '/*.mat']);
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
    EEG = pop_importdata('dataformat','matlab','nbchan',16,'data',[fnames(subj).folder '/' fnames(subj).name],'setname','Unlabeled EEG','srate',128,'pnts',0,'xmin',0);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');
    EEG = eeg_checkset( EEG );
    EEG.setname = 'Unlabeled EEG';

    % add channel labels and locations - COMPLETE WITH CORRESPONDING PATH
    EEG=pop_chanedit(EEG, {'lookup',[]},'changefield',{1,'labels','F7'},'append',1,'changefield',{2,'labels','F3'},'changefield',{2,'datachan',1},'append',2,'changefield',{3,'labels','F4'},'changefield',{3,'datachan',1},'append',3,'changefield',{4,'labels','F8'},'changefield',{4,'datachan',1},'append',4,'changefield',{5,'labels','T7'},'changefield',{5,'datachan',1},'append',5,'changefield',{6,'labels','C3'},'changefield',{6,'datachan',1},'append',6,'changefield',{7,'labels','Cz'},'changefield',{7,'datachan',1},'append',7,'changefield',{8,'labels','C4'},'changefield',{8,'datachan',1},'append',8,'changefield',{9,'labels','T8'},'changefield',{9,'datachan',1},'append',9,'changefield',{10,'labels','P7'},'changefield',{10,'datachan',1},'append',10,'changefield',{11,'labels','P3'},'changefield',{11,'datachan',1},'append',11,'changefield',{12,'labels','Pz'},'changefield',{12,'datachan',1},'append',12,'changefield',{13,'labels','P4'},'changefield',{13,'datachan',1},'append',13,'changefield',{14,'labels','P8'},'changefield',{14,'datachan',1},'append',14,'changefield',{15,'labels','O1'},'changefield',{15,'datachan',1},'append',15,'changefield',{16,'labels','O2'},'changefield',{16,'datachan',1},'lookup',[]);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG.setname = 'Raw EEG';

    % band-pass filter data
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',45);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','Filtered EEG','gui','off');
    EEG = eeg_checkset( EEG );

    % apply 'cleanline' to notch-filter data
    EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:16] ,'computepower',0,'linefreqs',[50 100] ,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',0,'winsize',4,'winstep',1);
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