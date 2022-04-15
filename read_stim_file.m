function [stim_amp] = read_stim_file (stim_dir)

stim_file = readmatrix(stim_dir);

stim_amp = stim_file(2:end-1,1);                                                                              % corp out useful part
stim_NaN = find(isnan(stim_amp));
stim_amp(stim_NaN) = [];                                                                                      % remove NaN
%stim_zeros = find(~stim_amp);
%stim_amp(stim_zeros) = []; 

end