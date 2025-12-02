% Author: F. Samuel Racz, 2025, The University of Texas at Austin
% email: fsr324@austin.utexas.edu

addpath(genpath('functions'))

fnames = dir('Dataset3_processed/*.gdf'); % replace with appropriate path to processed SU-SZ data
ns = length(fnames);

subjIDs = {fnames.name};
n_hc = sum(contains(subjIDs,'hc_'));
n_sz = sum(contains(subjIDs,'sch_'));

subj_list = cell(ns,1);
for n = 1:ns
    subj_list{n} = fnames(n).name(1:end-17);
end
subj_hc = subj_list(contains(subj_list,'hc_')==1);
subj_sz = subj_list(contains(subj_list,'sch_')==1);

%% Specify feature extraction parameters

nch = 19;
fs = 250;
chlab = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','C3','Cz','C4','T8','P7','P3','Pz','P4','P8','O1','O2'};
bands = {'delta','theta','alpha','beta','gamma'};
freqs = [2, 4; 4, 7; 7, 12; 12, 30; 30, 45];
Nb = length(bands);


% spectrogram parameters
TE = 10; % data segmentation epoch length
sTE = TE; % step size for epoching - non-overlapping windows
W = 1; % spectrogram window length
ws = 0.75; % spectrogram window overlap
freq = (1:0.25:45)'; % frequency resolution
Nf = length(freq);
Nxt = floor((TE-W)/(1-ws));

vnames_BLP_mean = cell(nch*Nb,1);
vnames_BLP_var = cell(nch*Nb,1);
for b = 1:Nb
    for ch = 1:nch
        vnames_BLP_mean{(b-1)*nch+ch} = [chlab{ch} '_' bands{b} '_mean'];
        vnames_BLP_var{(b-1)*nch+ch} = [chlab{ch} '_' bands{b} '_var'];
    end
end

% time domain parameters
Hjorth_params = {'activity', 'mobility', 'complexity'};
s = [8, 16, 32, 64, 128, 256];
M = TE*fs;

vnames_Hjorth = cell(nch*length(Hjorth_params),1);
for h = 1:length(Hjorth_params)
    for ch = 1:nch
        vnames_Hjorth{(h-1)*nch+ch} = [chlab{ch} '_' Hjorth_params{h}];
    end
end

vnames_DFA = cell(nch,1);
for ch = 1:nch
    vnames_DFA{ch} = [chlab{ch} '_DFA'];
end

% non-linear parameters
m = zeros(Nb,1);
L = zeros(Nb,1);
for b = 1:Nb
    [mb, Lb] = func_PE_params(  fs, freqs(b,1), freqs(b,2));
    m(b) = mb;
    L(b) = Lb;
end

vnames_PE = cell(nch*Nb,1);
for b = 1:Nb
    for ch = 1:nch
        vnames_PE{(b-1)*nch+ch} = [chlab{ch} '_' bands{b} '_PE'];
    end
end

% connectivity parameters
th_vec = [0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5];
Nth = length(th_vec);

vnames_PLI = cell(nch*Nb,1);
for b = 1:Nb
    for ch = 1:nch
        vnames_PLI{(b-1)*nch+ch} = [chlab{ch} '_' bands{b} '_PLI'];
    end
end

T = 100;
Nt = floor(T/TE);

vnames = [vnames_BLP_mean; vnames_BLP_var; vnames_Hjorth; vnames_DFA; vnames_PE; vnames_PLI];
Nv = length(vnames);

%% obtain features
for subj = 1:ns
    [data,h] = sload([fnames(subj).folder '/' fnames(subj).name],'OVERFLOWDETECTION','off');

    % select equal amount of data for all subjects
    data = data(1:T*fs,:);

    if contains(fnames(subj).name,'hc_')
        y_true = 0;
    elseif contains(fnames(subj).name,'sch_')
        y_true = 1;
    else
        error('wrong data loaded')
    end

    % store band-pass filtered data for connectivity analysis
    eeg_band = cell(Nb,1);
    for b = 1:Nb
        eeg_band{b} = bp_filter(data,fs,freqs(b,1),freqs(b,2),4);
    end

    % create data structure
    table_feat = table();
    table_feat.subjID = repmat(subj_list(subj),Nt,1);
    table_feat.y = repmat(y_true,Nt,1);
    for n = 1:Nv
        table_feat.(vnames{n}) = zeros(Nt,1);
    end

    for t = 1:Nt
        % get data segment
        tstart = (t-1)*floor(sTE*fs) + 1;
        tstop = tstart + TE*fs - 1;
        eeg = data(tstart:tstop,:);

        % store unscaled variance
        eeg_var = var(eeg,[],1);

        % z-score data
        eeg = zscore(eeg);

        % obtain spectrograms
        S_ch = cell(nch,1);
        for ch = 1:nch
            [S,~,xt] = spectrogram(eeg(:,ch),W*fs,floor(ws*fs),freq,fs);
            S_ch{ch} = 20*log10(abs(S));
        end

        % extract spectrogram features
        for b = 1:Nb
            for ch = 1:nch
                spec = S_ch{ch};
                spec_b = spec(freq>=freqs(b,1) & freq<=freqs(b,2),:);
                spec_BLP = mean(spec_b,1);
                spec_BLP_mean = mean(spec_BLP);
                spec_BLP_var = std(spec_BLP);

                vname_mean = [chlab{ch} '_' bands{b} '_mean'];
                vname_var = [chlab{ch} '_' bands{b} '_var'];

                table_feat.(vname_mean)(t) = spec_BLP_mean;
                table_feat.(vname_var)(t) = spec_BLP_var;
            end
        end

        % extract time domain features
        [h_comp, h_mob, ~] = func_Hjorth(eeg);
        Fc = DCCC_online_matrix(eeg,s,M);

        for ch = 1:nch
            p_DFA = polyfit(log(s),log(squeeze(Fc.DCCC(ch,ch,:))),1);

            vname_act = [chlab{ch} '_' Hjorth_params{1}];
            vname_mob = [chlab{ch} '_' Hjorth_params{2}];
            vname_comp = [chlab{ch} '_' Hjorth_params{3}];
            vname_DFA = [chlab{ch} '_DFA'];

            table_feat.(vname_act)(t) = eeg_var(ch); % use unscaled variance, as activity=1 for z-scored dat
            table_feat.(vname_mob)(t) = h_mob(ch);
            table_feat.(vname_comp)(t) = h_comp(ch);
            table_feat.(vname_DFA)(t) = p_DFA(1);
        end

        % extract nonlinear features
        for b = 1:Nb
            for ch = 1:nch
                eeg_b = eeg_band{b}(tstart:tstop,:);
                PE = PermEn(eeg_b(:,ch),m(b),L(b));
                vname_PE = [chlab{ch} '_' bands{b} '_PE'];

                table_feat.(vname_PE)(t) = PE.Hd_norm;
            end
        end

        % extract connectivity features
        for b = 1:Nb
            eeg_b = eeg_band{b}(tstart:tstop,:);
            [~,pli,~] = PLI_mat(eeg_b);

            PLI_th = zeros(nch,Nth);
            for th = 1:Nth
                mat_th = func_sparsity_th_matrix(pli,th_vec(th));
                mat_th(eye(nch)==1) = 0;
                PLI_th(:,th) = sum(mat_th,2);
            end

            PLI_vec = mean(PLI_th,2);

            for ch = 1:nch
                vname_PLI = [chlab{ch} '_' bands{b} '_PLI'];
                table_feat.(vname_PLI)(t) = PLI_vec(ch);
            end
        end



    end

    % save data
    fname_out = [fnames(subj).name(1:end-14) '_features.mat'];
    if ~exist('Dataset3_ML_features','dir')
        mkdir('Dataset3_ML_features')
    end
    save(['Dataset3_ML_features/' fname_out],'table_feat')
    disp([fnames(subj).name(1:end-14) ' done.'])

end