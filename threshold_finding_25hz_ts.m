%% Test
% clc
% clear all
% 


%%
function [max_fitresult_a, max_fitresult_b, max_fitresult_c]=threshold_finding_25hz_ts(pulsesintrain,repetition, phase_width_ms, peak_threshold_mV, peak_distance_sr, recording_dir, stim_dir, name, bin_left_sr, bin_right_sr)
%% fomular

    close all

    freq_Hz = 25;
    sample_rate = 50000;

    stim_amp = read_stim_file([stim_dir '.txt']);
    [trace] = abfload(recording_dir);
    
    peroid_dur_s = 1/freq_Hz;
    peroid_dur_ms = peroid_dur_s*1000;
    peroid_dur_sr = peroid_dur_s*sample_rate;
    phase_width_sr = phase_width_ms*0.001*sample_rate;
    pluse_width_sr = phase_width_sr*2;
    
    if isempty(bin_left_sr)
    bin_left_sr = sample_rate * 0.008 + pluse_width_sr;                                       %  ML are within 8ms following a pulse Jeson
    end
    
    if isempty(bin_right_sr)
        bin_right_sr = 10000000000000;                                       
    end


%% processing recording

    ttls = find(trace(:,2)>2);                                                 % pulse trigger has the height of 3, this line find all the samples that in the pulses
    
    trgs_on = ttls(find(diff(ttls)>pluse_width_sr*2)+1);                       % find the end of each pulse, +1 make this line find the start of pulses                                   
    
    trgs_on = [ttls(1); trgs_on];                                              % adding the onset sample of first pulse, it was not included
    
    amp_trigs_mat = [abs(stim_amp) trgs_on];
    
    end_trig = amp_trigs_mat(end,:);
    
    zeros_trgs = find(amp_trigs_mat(:,1) == 0.1);                              % deleting fake 0 triggers
    
    amp_trigs_mat(zeros_trgs, :) = [];
    
    amp_trigs_mat = [amp_trigs_mat; end_trig];                                 % add aditional ending trigger to make for loop easier

%% raster

    train_on = amp_trigs_mat(find(diff(amp_trigs_mat(:,2)) > pulsesintrain * peroid_dur_sr) + 1,:);       % find the end of each train, +1 make this line find the start of pulses  
    
    train_on = [amp_trigs_mat(1,:); train_on];                                  % adding the onset sample of first pulse, it was not included
    
    
    all_spks_timing = [];
    
    for i = 1:length(train_on)-1
    
        trace_clip =  trace(train_on(i, 2):train_on(i+1, 2),:);
    
        [spks_amp,spks_timing] =  findpeaks(trace_clip(:,1),1,'MINPEAKHEIGHT',peak_threshold_mV, 'MinPeakDistance',peak_distance_sr);
    
        spks_timing(:, 2) = train_on(i, 1);
    
        spks_timing = [spks_timing spks_amp];
    
        all_spks_timing = [all_spks_timing;spks_timing];
        
    end

    fig(3) = figure;
    scatter3(all_spks_timing(:,1),all_spks_timing(:,2),all_spks_timing(:,3),70,'|');
    view(2)
    
    t_peroid_sr = train_on(2,2)-train_on(1,2);
    t_peroid_ms = t_peroid_sr/sample_rate;
    xlim([0 t_peroid_sr])
    ylim([0 -min(stim_amp)])

    xline([pluse_width_sr, pluse_width_sr + peroid_dur_sr, pluse_width_sr + peroid_dur_sr*2, pluse_width_sr + peroid_dur_sr*3, pluse_width_sr + peroid_dur_sr*4,pluse_width_sr + peroid_dur_sr*5],'r')
    xline([bin_left_sr, bin_left_sr + peroid_dur_sr, bin_left_sr + peroid_dur_sr*2, bin_left_sr + peroid_dur_sr*3, bin_left_sr + peroid_dur_sr*4,bin_left_sr + peroid_dur_sr*5])
    xline(bin_right_sr)                                                    %%%%%% need adjust for each protocol

    xticks([0:peroid_dur_sr:t_peroid_sr])
    xticklabels([0:peroid_dur_ms:t_peroid_ms])

    xlabel('Time (ms)')
    ylabel('Amplitude (uA)')
    title(name)

%% BC sigmoid
    dir_spks_count = [];
    
    for i = 1:length(amp_trigs_mat) -1

        if amp_trigs_mat(i,2) + bin_left_sr > amp_trigs_mat(i+1,2)
            
            bin_left_i =  amp_trigs_mat(i+1,2) - amp_trigs_mat(i,2);

        else
            bin_left_i = bin_left_sr;

        end
    
        trace_clip =  trace(amp_trigs_mat(i,2) + pluse_width_sr:amp_trigs_mat(i,2) + bin_left_i,:);
    
        [spks_amp,spks_timing] =  findpeaks(trace_clip(:,1),1,'MINPEAKHEIGHT',peak_threshold_mV, 'MinPeakDistance',peak_distance_sr);
    
        spks_count = [length(spks_timing) amp_trigs_mat(i, 1)];
    
        dir_spks_count = [dir_spks_count; spks_count];
        
    
    end
    
    [i,j,k]=unique(dir_spks_count(:,2));
    dir_sum_spks_num=[i accumarray(k,dir_spks_count(:,1),[],@sum)];
    dir_mean_spks_num = [dir_sum_spks_num(:,1) dir_sum_spks_num(:,2)/repetition];
    
    fig(2) = figure;
    [dir_fitresult, dir_gof] = FitMitch(dir_mean_spks_num(:,1),dir_mean_spks_num(:,2));
    xlabel('Amplitude (uA)')
    ylabel('Spikes number')
    xlim([0 160])
    title(['BC ',name],['r^2=',num2str(dir_gof.rsquare), '  a=',num2str(dir_fitresult.a),'  b=', num2str(dir_fitresult.b), '  c=',num2str(dir_fitresult.c)])
    legend('OFF')

%% Indirect sigmoid
    all_spks_count = [];
    
    for i = 1:length(amp_trigs_mat) -1

        if amp_trigs_mat(i,2)+bin_right_sr > amp_trigs_mat(i+1,2)
            
            bin_right_i =  amp_trigs_mat(i+1,2) - amp_trigs_mat(i,2);

        else
            bin_right_i = bin_right_sr;

        end
    
        trace_clip =  trace(amp_trigs_mat(i,2)+bin_left_sr:amp_trigs_mat(i,2)+bin_right_i,:);
    
        [spks_amp,spks_timing] =  findpeaks(trace_clip(:,1),1,'MINPEAKHEIGHT',peak_threshold_mV, 'MinPeakDistance',peak_distance_sr);
    
        spks_count = [length(spks_timing) amp_trigs_mat(i, 1)];
    
        all_spks_count = [all_spks_count; spks_count];
        
    
    end
    
    [i,j,k]=unique(all_spks_count(:,2));

% probablidity Sigmoid
    % mean_spks_num=[i accumarray(k,all_spks_count(:,1),[],@mean)];
    % 
    % spks_prob = [(all_spks_count(:,1) - mean_spks_num(1,2)) > 0 (all_spks_count(:,2))];
    % 
    % [ii,jj,kk]=unique(spks_prob(:,2));
    % 
    % mean_spks_prob=[ii accumarray(kk,spks_prob(:,1),[],@mean)];
    % 
    % fig(3) = figure;
    % [prob_fitresult, prob_gof] = FitBoltzmann(mean_spks_prob(:,1),mean_spks_prob(:,2));
    % xlabel('Amplitude (uA)')
    % ylabel('Firing or not probability')
    % xlim([0 160])
    % title(['Probablity ',name],['r^2=',num2str(prob_gof.rsquare), '  a=',num2str(prob_fitresult.a),'  b=', num2str(prob_fitresult.b)])
    % legend('OFF')

% Maximun sigmoid
    sum_spks_num=[i accumarray(k,all_spks_count(:,1),[],@sum)];
    mean_spks_num = [sum_spks_num(:,1) sum_spks_num(:,2)/repetition];
    
    fig(1) = figure;
    [max_fitresult, max_gof] = FitMitch(mean_spks_num(:,1),mean_spks_num(:,2));
    xlabel('Amplitude (uA)')
    ylabel('Spikes number')
    xlim([0 160])
    title(['Indirect ',name],['r^2=',num2str(max_gof.rsquare), '  a=',num2str(max_fitresult.a),'  b=', num2str(max_fitresult.b), '  c=',num2str(max_fitresult.c)])
    legend('OFF')

%% plot for adjusting findpeak threshold
    figure;
    [spks_amp,spks_timing] =  findpeaks(trace(:,1),1,'MINPEAKHEIGHT',peak_threshold_mV, 'MinPeakDistance',peak_distance_sr);
    plot(trace)
    hold on
    plot(spks_timing,spks_amp,'o')
    xticks([0:sample_rate:length(trace)])
    xticklabels([0:1:(length(trace)/sample_rate)])
    xlim([amp_trigs_mat(end-pulsesintrain,2) amp_trigs_mat(end,2)])
    title(name)
    xlabel('Time (s)')
    ylabel('Memberine potential (mV)')
    hold off

%%
    saveas(fig,name)
    
    max_fitresult_a = max_fitresult.a;
    max_fitresult_b = max_fitresult.b;
    max_fitresult_c = max_fitresult.c; 

end

