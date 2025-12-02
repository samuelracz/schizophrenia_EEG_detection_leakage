function [ I ] = func_organize_spectrograms_data2(S_cell, freq, xt, frange)

% channel map to output
Nf = length(freq(freq>frange(1) & freq<=frange(2)));
Nt = length(xt);

N1 = 5*Nf;
N2 = 5*Nt;


% channel mapping
ch_map = [1, 4;... Fp2
    2, 5;... F8
    3, 5;... T8 
    4, 5;... P8
    5, 4;... O2
    1, 2;... Fp1
    2, 1;... F7
    3, 1;... T7
    4, 1;... P7
    5, 2;... O1
    2, 4;... F4
    3, 4;... C4
    4, 4;... P4
    2, 2;... F3
    3, 2;... C3
    4, 2;... P3
    2, 3;... Fz
    3, 3;... Cz
    4, 3]; % Pz

I = zeros(N1,N2);

for ch = 1:length(ch_map)
    idx1 = ch_map(ch,1);
    idx2 = ch_map(ch,2);

    I((idx1-1)*Nf+1:idx1*Nf,(idx2-1)*Nt+1:idx2*Nt) = S_cell{ch}(freq>frange(1) & freq<=frange(2),:);
end


% filling up 'missing spots' with averages
avg_af7_idx = [6, 7, 14];
af7_pos = [1, 1];
S_af7 = zeros(Nf,Nt,length(avg_af7_idx));
for i = 1:length(avg_af7_idx)
    S_af7(:,:,i) = S_cell{avg_af7_idx(i)}(freq>frange(1) & freq<=frange(2),:);
end
I((af7_pos(1)-1)*Nf+1:af7_pos(1)*Nf,(af7_pos(2)-1)*Nt+1:af7_pos(2)*Nt) = mean(S_af7,3);

avg_fpz_idx = [1, 6, 17];
fpz_pos = [1, 3];
S_fpz = zeros(Nf,Nt,length(avg_fpz_idx));
for i = 1:length(avg_fpz_idx)
    S_fpz(:,:,i) = S_cell{avg_fpz_idx(i)}(freq>frange(1) & freq<=frange(2),:);
end
I((fpz_pos(1)-1)*Nf+1:fpz_pos(1)*Nf,(fpz_pos(2)-1)*Nt+1:fpz_pos(2)*Nt) = mean(S_fpz,3);

avg_af8_idx = [1, 2, 11];
af8_pos = [1, 5];
S_af8 = zeros(Nf,Nt,length(avg_af8_idx));
for i = 1:length(avg_af8_idx)
    S_af8(:,:,i) = S_cell{avg_af8_idx(i)}(freq>frange(1) & freq<=frange(2),:);
end
I((af8_pos(1)-1)*Nf+1:af8_pos(1)*Nf,(af8_pos(2)-1)*Nt+1:af8_pos(2)*Nt) = mean(S_af8,3);

avg_po7_idx = [9, 10, 16];
po7_pos = [5, 1];
S_po7 = zeros(Nf,Nt,length(avg_po7_idx));
for i = 1:length(avg_po7_idx)
    S_po7(:,:,i) = S_cell{avg_po7_idx(i)}(freq>frange(1) & freq<=frange(2),:);
end
I((po7_pos(1)-1)*Nf+1:po7_pos(1)*Nf,(po7_pos(2)-1)*Nt+1:po7_pos(2)*Nt) = mean(S_po7,3);

avg_oz_idx = [5, 10, 19];
oz_pos = [5, 3];
S_oz = zeros(Nf,Nt,length(avg_oz_idx));
for i = 1:length(avg_oz_idx)
    S_oz(:,:,i) = S_cell{avg_oz_idx(i)}(freq>frange(1) & freq<=frange(2),:);
end
I((oz_pos(1)-1)*Nf+1:oz_pos(1)*Nf,(oz_pos(2)-1)*Nt+1:oz_pos(2)*Nt) = mean(S_oz,3);

avg_po8_idx = [4, 5, 13];
po8_pos = [5, 5];
S_po8 = zeros(Nf,Nt,length(avg_po8_idx));
for i = 1:length(avg_po8_idx)
    S_po8(:,:,i) = S_cell{avg_po8_idx(i)}(freq>frange(1) & freq<=frange(2),:);
end
I((po8_pos(1)-1)*Nf+1:po8_pos(1)*Nf,(po8_pos(2)-1)*Nt+1:po8_pos(2)*Nt) = mean(S_po8,3);

% normalization
glob_min = min(I(:));
glob_max = max(I(:));

% scale between [0, 1]
I = (I - glob_min)./(glob_max-glob_min);

% upscale to [0, 255]
I = I*255;

end