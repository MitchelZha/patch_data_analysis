% clc
% %clear all
% close all
% 
% a = 54.4
% b = 27.55
% 
% stm_duration_min = 4
% pulse_length_ms = 1
% freq_Hz = 25
% max_amp_uA = 200
% min_amp_uA= 0


% function [] = wn_generator(frequency, phase_width, stim_duration, max_amp_uA, min_amp_uA, a, b, c)

cd D:\5-Mingsong
    
%% get pulse number
    peroid_dur_ms = (1/freq_Hz)*1000
    
    inter_pulse_ms = peroid_dur_ms - pulse_length_ms*2
    
    pulse_num = (stm_duration_min*60*1000) / peroid_dur_ms
    
%% get mean and sd based on sigmoid
    syms x
    sig75 = double(vpa(solve(1/(1+exp(-(x-a)/b)) == 0.75, x)))
    sig50 = double(vpa(solve(1/(1+exp(-(x-a)/b)) == 0.5, x)))
    sig25 = double(vpa(solve(1/(1+exp(-(x-a)/b)) == 0.25, x)))

    amp_sd_uA = sig75-sig50

    amp_mean_uA = [sig25; sig75; sig75+2*amp_sd_uA]
   
%% restrain mean in a range & Generate raw stim amp array
    for i = 3
    
        if amp_mean_uA(i) > max_amp_uA - 1.645*amp_sd_uA
        
            amp_mean_uA(i) = max_amp_uA - 1.645*amp_sd_uA;
        
        elseif amp_mean_uA(i) < min_amp_uA + 1.645*amp_sd_uA
        
            amp_mean_uA(i) = min_amp_uA + 1.645*amp_sd_uA;
        
        end

        amp_array_uA = amp_mean_uA(i) + amp_sd_uA*randn(pulse_num,1);

        amp_array_uA = (amp_array_uA)* -1  

%% put stim array into a txt file
        Header={'Multi Channel Systems MC_Stimulus II','ASCII import Version 1.10',...
        'channels: 8','output mode: current','format: ','3','channel: 1',...
        'value time value time value time repeat'};
        
    
        str = strcat('wn_Pw=',num2str(pulse_length_ms),'Freq=',num2str(freq_Hz),'_Mean=',num2str(amp_mean_uA(i)),'_SD=',num2str(amp_sd_uA),'_Duration=',num2str(stm_duration_min),'.txt')
        
        while exist(str,"file") > 0                                                                                                                                 % in case new file replaced the file with the same name

            str = [(2), str]
            
        end

        str = strcat(pwd,'\',str)
        fid=fopen(str,'w');
        fprintf(fid,'%s\r\n%s\r\n%s\r\n%s\r\n%s\t%s\r\n\r\n%s\r\n%s\t%s\r\n',Header{:},'');
        fclose(fid);
        
        stim=[0 0 0 0 0 0 1];
        
        for i1=1:numel(amp_array_uA)                                                                                                                                % restain each pulse in the ideal range

            if amp_array_uA(i1) < - max_amp_uA
            
                amp_array_uA(i1) = max_amp_uA;
            
            elseif amp_array_uA(i1) > min_amp_uA
            
                amp_array_uA(i1) = 0.1;
            
            end

            stim = [stim;[amp_array_uA(i1),pulse_length_ms*1000,-amp_array_uA(i1),pulse_length_ms*1000,0,inter_pulse_ms*1000,1] ];

        end
        
        stim(end+1,:)=[0 0 0 0 0 0 1];
        
        dlmwrite(str,stim,'roffset',0,'coffset',1,'-append','delimiter',' ','newline','pc','precision','%.3f');
    
    end

% end