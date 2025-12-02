function x_filtered = bp_filter(ts, sf, lf, hf, ord)
%parameters:
%ts: time series (last column should be time)
%output:
%xx: the filtered signal
%written by Frigyes Samuel Racz

if size(ts,1) < size(ts,2)
    ts = ts';
end

HalfSf = sf/2;
Wn = [lf hf]./HalfSf;

[bh,ah] = butter(ord,Wn,'bandpass');

x_filtered = filtfilt(bh,ah,ts);
