%% Test
% clc
% clear all
% 
% freq_Hz = 2;
% phase_width_ms = 1;
% peak_threshold_mV = -15;
% peak_distance_sr = 30;
% addpath 'D:\5-Mingsong'; addpath 'D:\5-Mingsong\Code';
% 
% cd 'D:\5-Mingsong\220406';
% recording_dir = 'Clampex\2022_04_06_0020.abf';
% stim_dir = 'RRampup_Pulseshuffled_Pw=1_Freq=2_Amp=0.1-120_Req=10';
% name = ['220406 OFFT ' stim_dir '.fig'];
% bin_left_sr = 400; bin_right_sr = 7700;
% 
% % [a, b, c] = threshold_finding(freq_Hz, phase_width_ms, peak_threshold_mV, peak_distance_sr, recording_dir, stim_dir, name, bin_left_sr,bin_right_sr );
% % wn_generator(25, 1, 5, 200, 0, a, b, c)


%%
function [max_fitresult_a, max_fitresult_b, max_fitresult_c]=threshold_finding(freq_Hz, phase_width_ms, peak_threshold_mV, peak_distance_sr, recording_dir, stim_dir, name, bin_left_sr, bin_right_sr)

%% Fomular
    close all

    sample_rate = 50000;

    stim_amp = read_stim_file([stim_dir '.txt']);
    [trace] = abfload(recording_dir);
    
    peroid_dur_s = 1/freq_Hz;
    peroid_dur_ms = peroid_dur_s*1000;
    peroid_dur_sr = peroid_dur_ms*0.001*sample_rate;
    phase_width_sr = phase_width_ms*0.001*sample_rate;
    pluse_width_sr = phase_width_sr*2;
    
    if isempty(bin_left_sr)
        bin_left_sr = sample_rate * 0.008 + pluse_width_sr;                                       %  ML are within 8ms Jeson
    end
    
    if isempty(bin_right_sr)
        bin_right_sr = peroid_dur_sr;                                       
    end


%% Processing recording
    ttls = find(trace(:,2)>2);                                                 % pulse trigger has the height of 3, this line find all the samples that in the pulses
    
    trgs_on = ttls(find(diff(ttls)>pluse_width_sr*2)+1);                       % find the end of each pulse, +1 make this line find the start of pulses                                   
    
    trgs_on = [ttls(1); trgs_on];                                              % adding the onset sample of first pulse, it was not included
    
    missing_trgs = find(diff(trgs_on) > peroid_dur_sr+5);
    
    for i = 1:length(missing_trgs)
        
        missing_trgs = find(diff(trgs_on) > peroid_dur_sr+5);
    
        trgs_on = [trgs_on(1:missing_trgs(1)); trgs_on(missing_trgs(1)) + peroid_dur_sr+5; trgs_on(missing_trgs(1)+1:end)];
    
    end
    
    amp_trigs_mat = [abs(stim_amp) trgs_on];

%% Raster
    all_spks_timing = [];
    
    for i = 1:length(amp_trigs_mat)
    
        trace_clip =  trace(amp_trigs_mat(i, 2):amp_trigs_mat(i, 2)+peroid_dur_sr-1,:);
    
        [spks_amp,spks_timing] =  findpeaks(trace_clip(:,1),1,'MINPEAKHEIGHT',peak_threshold_mV, 'MinPeakDistance',peak_distance_sr);
    
        spks_timing(:, 2) = amp_trigs_mat(i, 1);
        spks_timing = [spks_timing spks_amp];
    
        all_spks_timing = [all_spks_timing;spks_timing];
        
    end
    
    fig = figure
    subplot(2,3,[1,2,4,5]);
    patch([bin_left_sr,bin_right_sr,bin_right_sr,bin_left_sr],[-1000,-1000,1000,1000],'black','FaceAlpha',.05)
    hold on
    scatter3(all_spks_timing(:,1),all_spks_timing(:,2),all_spks_timing(:,3),5,'filled');
    hold off
    view(2)
    xlim([0 peroid_dur_sr])
    ylim([0 -min(stim_amp)])

    xline(pluse_width_sr,'--r',{'pulse offset'})

    xticks([0:peroid_dur_sr/40:peroid_dur_sr])
    xticklabels([0:peroid_dur_ms/40:peroid_dur_ms])
    xlabel('Time (ms)')
    ylabel('Amplitude (uA)')
    title(name)

%% BC Sigmoid
    dir_spks_count = [];
    
    for i = 1:length(amp_trigs_mat)
    
        trace_clip =  trace(amp_trigs_mat(i, 2) + pluse_width_sr :amp_trigs_mat(i, 2)+bin_left_sr,:);
    
        [spks_amp,spks_timing] =  findpeaks(trace_clip(:,1),1,'MINPEAKHEIGHT',peak_threshold_mV, 'MinPeakDistance',peak_distance_sr);
    
        spks_count = [length(spks_timing) amp_trigs_mat(i, 1)];
    
        dir_spks_count = [dir_spks_count; spks_count];
        
    end
    
    [i,j,k]=unique(dir_spks_count(:,2));
    
    dir_spks_num=[i accumarray(k,dir_spks_count(:,1),[],@mean)];

    subplot(2,3,3);
    [dir_fitresult, dir_gof] = FitMitch(dir_spks_num(:,1),dir_spks_num(:,2));
    xlabel('Amplitude (uA)')
    ylabel('Spikes number')
    xlim([0 160])
    title('BC',['r^2=',num2str(dir_gof.rsquare), '  a=',num2str(dir_fitresult.a),'  b=', num2str(dir_fitresult.b), '  c=',num2str(dir_fitresult.c)])
    legend('OFF')


%% PR Sigmoid
    all_spks_count = [];
    
    for i = 1:length(amp_trigs_mat)
    
        trace_clip =  trace(amp_trigs_mat(i, 2)+bin_left_sr:amp_trigs_mat(i, 2)+bin_right_sr,:);
    
        [spks_amp,spks_timing] =  findpeaks(trace_clip(:,1),1,'MINPEAKHEIGHT',peak_threshold_mV, 'MinPeakDistance',peak_distance_sr);
    
        spks_count = [length(spks_timing) amp_trigs_mat(i, 1)];
    
        all_spks_count = [all_spks_count; spks_count];
        
    end
    
    [i,j,k]=unique(all_spks_count(:,2));
    
    mean_spks_num=[i accumarray(k,all_spks_count(:,1),[],@mean)];

% Probablidity sigmoid
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
    % title(['Probablity ',name],['r^2=',num2str(prob_gof.rsquare), '  a=',num2str(prob_fitresult.a),'  b=', num2str(prob_fitresult.b)]))
    % legend('OFF')


% Maximun sigmoid
   
    subplot(2,3,6);
    [max_fitresult, max_gof] = FitMitch(mean_spks_num(:,1),mean_spks_num(:,2));
    xlabel('Amplitude (uA)')
    ylabel('Spikes number')
    xlim([0 160])
    title('PR',['r^2=',num2str(max_gof.rsquare), '  a=',num2str(max_fitresult.a),'  b=', num2str(max_fitresult.b), '  c=',num2str(max_fitresult.c)])
    legend('OFF')

%% Plot for adjusting findpeak threshold
    figure;
    [spks_amp,spks_timing] =  findpeaks(trace(:,1),1,'MINPEAKHEIGHT',peak_threshold_mV, 'MinPeakDistance',peak_distance_sr);
    plot(trace)
    hold on
    plot(spks_timing,spks_amp,'o')
    xticks([0:sample_rate:length(trace)])
    xticklabels([0:1:(length(trace)/sample_rate)])
    xlim([trgs_on(end-5) trgs_on(end)])
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

