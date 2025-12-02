% Author: F. Samuel Racz, 2025, The University of Texas at Austin
% email: fsr324@austin.utexas.edu

fnames = dir('Dataset2_processed/*.gdf'); % replace with appropriate path to processed RepOD data
ns = length(fnames);

subjIDs = {fnames.name};
n_hc = sum(contains(subjIDs,'h0') | contains(subjIDs,'h1'));
n_sz = sum(contains(subjIDs,'s0') | contains(subjIDs,'s1'));

subj_list = cell(ns,1);
for n = 1:ns
    subj_list{n} = fnames(n).name(1:6);
end
subj_hc = subj_list(contains(subj_list,'hc')==1);
subj_sz = subj_list(contains(subj_list,'sz')==1);


%% Specify feature extraction parameters

img_size = [224, 224]; % ResNET-18/50

load('ws_16chanlocs.mat')
nch = 19;
fs = 250;

TE = 5; % data segmentation epoch length
sTE = TE; % step size for epoching - non-overlapping windows
W = 1; % spectrogram window length
ws = 0.9; % spectrogram window overlap
freq = (1:0.25:45)'; % frequency resolution
Nf = length(freq);
Nxt = floor((TE-W)/(1-ws));

T = 2*60;
Nt = floor(T/TE);

feat_struct = struct(...
    'subjID', [],...
    'y', [],...
    'X', []);

%% obtain features
for subj = 1:ns
    [data,h] = sload([fnames(subj).folder '/' fnames(subj).name],'OVERFLOWDETECTION','off');

    if contains(fnames(subj).name,'h0') || contains(fnames(subj).name,'h1')
        y_true = 0;
    elseif contains(fnames(subj).name,'s0') || contains(fnames(subj).name,'s1')
        y_true = 1;
    else
        error('wrong data loaded')
    end

    for t = 1:Nt
        % get data segment
        tstart = (t-1)*floor(sTE*fs) + 1;
        tstop = tstart + TE*fs - 1;
        eeg = data(tstart:tstop,:);

        % z-score data
        eeg = zscore(eeg);

        % obtain spectrograms
        S_ch = cell(nch,1);
        for ch = 1:nch
            [S,~,xt] = spectrogram(eeg(:,ch),W*fs,floor(ws*fs),freq,fs);
            S_ch{ch} = 20*log10(abs(S));
        end

        % obtain scalograms
        CW_ch = cell(nch,1);
        for ch = 1:nch
            [CW,cf] = cwt(eeg(:,ch),'morse',fs);
            CW_ch{ch} = flipud(log(abs(CW(cf>=1 & cf<=45,:))));
        end
        cf = cf(cf>=1 & cf<=45);
        xt_cwt = linspace(0,TE,size(CW,2));

        % obtain Wigner-Ville distribution
        DWV_ch = cell(nch,1);
        for ch = 1:nch
            [DWV,bf,xt_dwv] = wvd(eeg(:,ch),fs,"smoothedPseudo");
            DWV_ch{ch} = log(abs(DWV(bf>=1 & bf<=45,:)));
        end
        bf = bf(bf>=1 & bf<=45);

        % create unscaled 'image'
        I1 = func_organize_spectrograms_data2(S_ch,freq,xt,[0, 45]);
        I2 = func_organize_spectrograms_data2(CW_ch,cf,xt_cwt,[1-eps, 45]);
        I3 = func_organize_spectrograms_data2(DWV_ch,bf,xt_dwv,[1-eps, 45]);

        % scale images to preferred size
        I1_res = imresize(I1,img_size);
        I2_res = imresize(I2,img_size);
        I3_res = imresize(I3,img_size);

        % collect images into 3D image
        I = cat(3,I1_res,I2_res);
        I = cat(3,I,I3_res);

        % store data
        feat_struct(t).subjID = subj;
        feat_struct(t).y = y_true;
        feat_struct(t).X = I;
        
    end

    % save data
    fname_out = [fnames(subj).name(1:end-14) '_features.mat'];
    if ~exist('Dataset2_CNN_features','dir')
        mkdir('Dataset2_CNN_features')
    end
    save(['Dataset2_CNN_features/' fname_out],'feat_struct')
    disp([fnames(subj).name(1:end-14) ' done.'])

end