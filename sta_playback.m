% freq_Hz = 25
% pulse_length_ms = 1
% repetition = 5
% 
% all_sta = [70 60 63]
% dir_sta = [10 1 5]
% indir_sta = [90 60 80]
% mean = 63

function []=sta_playback(freq_Hz,pulse_length_ms,repetition,all_sta,dir_sta,indir_sta, stim_mean)

cd D:\5-Mingsong

%% Frequency
peroid_dur_ms = (1/freq_Hz)*1000;

inter_pulse_ms = peroid_dur_ms - pulse_length_ms*2;

%% Clips
baseline = repmat(stim_mean,length(indir_sta), repetition);

gap = zeros(freq_Hz,1);


all_sta = repmat(all_sta, repetition,1)';                                    

dir_sta = repmat(dir_sta, repetition,1)';  

indir_sta = repmat(indir_sta, repetition,1)';


rev_all_sta = stim_mean-(all_sta-stim_mean);                                   

rev_dir_sta = stim_mean-(dir_sta-stim_mean);

rev_indir_sta = stim_mean -(indir_sta-stim_mean);


%% Format Data
amp_array_uA = [baseline(:);gap; all_sta(:);gap; dir_sta(:);gap; indir_sta(:);gap; rev_all_sta(:);gap; rev_dir_sta(:);gap; rev_indir_sta(:)];

Header={'Multi Channel Systems MC_Stimulus II','ASCII import Version 1.10',...
'channels: 8','output mode: current','format: ','3','channel: 1',...
'value time value time value time repeat'};

str = strcat(pwd,'\STA_playback_Freq=',num2str(freq_Hz),'_STA=',mat2str(round(indir_sta(end-5:end),1)),'.txt');
fid=fopen(str,'w');
fprintf(fid,'%s\r\n%s\r\n%s\r\n%s\r\n%s\t%s\r\n\r\n%s\r\n%s\t%s\r\n',Header{:},'');
fclose(fid);

stim=[0 0 0 0 0 0 1];

for i1=1:numel(amp_array_uA)
stim = [stim;[amp_array_uA(i1),pulse_length_ms*1000,-amp_array_uA(i1),pulse_length_ms*1000,0,inter_pulse_ms*1000,1] ];    
end

stim(end+1,:)=[0 0 0 0 0 0 1];

dlmwrite(str,stim,'roffset',0,'coffset',1,'-append','delimiter',' ','newline','pc','precision','%.3f');

end









