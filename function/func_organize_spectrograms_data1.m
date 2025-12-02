function [ I ] = func_organize_spectrograms_data1(S_cell, freq, xt, frange)

% channel map to output
Nf = length(freq(freq>frange(1) & freq<=frange(2)));
Nt = length(xt);

N1 = 4*Nf;
N2 = 5*Nt;


% channel mapping
ch_map = [1, 1;... F7
    1, 2;... F3
    1, 4;... F4
    1, 5;... F8
    2, 1;... T7
    2, 2;... C3
    2, 3;... Cz
    2, 4;... C4
    2, 5;... T8
    3, 1;... P7
    3, 2;... P3
    3, 3;... Pz
    3, 4;... P4
    3, 5;... P8
    4, 2;... O1
    4, 4]; % O2

I = zeros(N1,N2);

for ch = 1:length(ch_map)
    idx1 = ch_map(ch,1);
    idx2 = ch_map(ch,2);

    I((idx1-1)*Nf+1:idx1*Nf,(idx2-1)*Nt+1:idx2*Nt) = S_cell{ch}(freq>frange(1) & freq<=frange(2),:);
end


% filling up 'missing spots' with averages
avg_fz_idx = [3, 4, 7];
fz_pos = [1, 3];
S_fz = zeros(Nf,Nt,length(avg_fz_idx));
for i = 1:length(avg_fz_idx)
    S_fz(:,:,i) = S_cell{avg_fz_idx(i)}(freq>frange(1) & freq<=frange(2),:);
end
I((fz_pos(1)-1)*Nf+1:fz_pos(1)*Nf,(fz_pos(2)-1)*Nt+1:fz_pos(2)*Nt) = mean(S_fz,3);

avg_po7_idx = [10, 11, 15];
po7_pos = [4, 1];
S_po7 = zeros(Nf,Nt,length(avg_po7_idx));
for i = 1:length(avg_po7_idx)
    S_po7(:,:,i) = S_cell{avg_po7_idx(i)}(freq>frange(1) & freq<=frange(2),:);
end
I((po7_pos(1)-1)*Nf+1:po7_pos(1)*Nf,(po7_pos(2)-1)*Nt+1:po7_pos(2)*Nt) = mean(S_po7,3);

avg_oz_idx = [12, 15, 16];
oz_pos = [4, 3];
S_oz = zeros(Nf,Nt,length(avg_oz_idx));
for i = 1:length(avg_oz_idx)
    S_oz(:,:,i) = S_cell{avg_oz_idx(i)}(freq>frange(1) & freq<=frange(2),:);
end
I((oz_pos(1)-1)*Nf+1:oz_pos(1)*Nf,(oz_pos(2)-1)*Nt+1:oz_pos(2)*Nt) = mean(S_oz,3);

avg_po8_idx = [13, 14, 16];
po8_pos = [4, 5];
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