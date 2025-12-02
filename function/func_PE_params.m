function [ m, L ] = func_PE_params(fs,f_lo,f_hi)

L = ceil(fs/(3*f_hi));
m = ceil((3*f_hi)/f_lo) + 1;

end