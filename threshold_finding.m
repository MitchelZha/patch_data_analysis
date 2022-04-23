% %% Test
% clc
% clear all
%
% addpath 'D:\OneDrive - UNSW\ephys'; addpath 'D:\OneDrive - UNSW\ephys\patch_data_analysis';
% pulsesintrain = 5;
% repetition = 2;
% phase_width_ms = 1;
% 
% cd 'D:\OneDrive - UNSW\ephys\220404';
% recording_dir = 'Clampex\2022_04_04_0029.abf';
% stim_dir = 'RRampup_Pulseshuffled_Pw=1_Freq=2_Amp=0.1-180_Req=10';
% name = ['220404 OFFT' stim_dir '.mat'];
% 
% peak_threshold_mV = -15;
% peak_distance_sr = 30;
% bin_left_sr = []; 
% bin_right_sr = [];
% % [a, b, c] = threshold_finding(freq_Hz, phase_width_ms, peak_threshold_mV, peak_distance_sr, recording_dir, stim_dir, name, bin_left_sr, bin_right_sr);
% % wn_generator(25, 1, 5, 200, 0, a, b, c)

%%
function [PR_fitresult_a, PR_fitresult_b, PR_fitresult_c]=threshold_finding(freq_Hz, phase_width_ms, peak_threshold_mV, peak_distance_sr, recording_dir, stim_dir, name, bin_left_sr, bin_right_sr)
%% Fomular
    close all

    sample_rate = 50000;

    stim_amp = read_stim_file([stim_dir '.txt']);
    [trace] = abfload(recording_dir);
    
    peroid_dur_s = 1/freq_Hz;
    peroid_dur_ms = peroid_dur_s*1000;
    peroid_dur_sr = peroid_dur_ms*0.001*sample_rate;
    phase_width_sr = phase_width_ms*0.001*sample_rate;
    pulse_width_sr = phase_width_sr*2;
    
    if isempty(bin_left_sr)
        bin_left_sr = sample_rate * 0.008 + pulse_width_sr;                                       %  ML are within 8ms Jeson
    end
    
    if isempty(bin_right_sr)
        bin_right_sr = peroid_dur_sr;                                       
    end


%% Processing recording
    ttls = find(trace(:,2)>2);                                                 % pulse trigger has the height of 3, this line find all the samples that in the pulses
    
    trgs_on = ttls(find(diff(ttls)>pulse_width_sr*3)+1);                       % find the end of each pulse, +1 make this line find the start of pulses                                   
    
    trgs_on = [ttls(1); trgs_on];                                              % adding the onset sample of first pulse, it was not included
    
%     missing_trgs = find(diff(trgs_on) > peroid_dur_sr+5);                    % Turn on when there is missing trigger caused by 0
%     for i = 1:length(missing_trgs) 
%         missing_trgs = find(diff(trgs_on) > peroid_dur_sr+5);
%         trgs_on = [trgs_on(1:missing_trgs(1)); trgs_on(missing_trgs(1)) + peroid_dur_sr+5; trgs_on(missing_trgs(1)+1:end)];
%     end

    stim_amp = abs(stim_amp);

    [spks_amp, spks_timing] = findpeaks(trace(1:end,1),'MINPEAKHEIGHT',peak_threshold_mV, 'MinPeakDistance',peak_distance_sr);

%% BC Sigmoid
    BC_spks_count = zeros(length(trgs_on),2);

    for i = 1:length(trgs_on) 
    
        spks = find(spks_timing > trgs_on(i) + pulse_width_sr + (0.0005*sample_rate) & spks_timing < trgs_on(i) + bin_left_sr);     %  adding 0.5ms to avoid delayed direct activaiton
    
        BC_spks_count(i,1) = length(spks);

        BC_spks_count(i,2) = stim_amp(i);

    end
    
    [i,j,k]=unique(BC_spks_count(:,2));
    
    BC_spks_sum=[i accumarray(k,BC_spks_count(:,1),[],@mean)];

    fig = figure;
    subplot(2,3,3);
    [dir_fitresult, dir_gof] = FitMitch(BC_spks_sum(:,1),BC_spks_sum(:,2));
    xlabel('Amplitude (uA)')
    ylabel('Spikes number')
    xlim([0 130])
    title('BC',['r^2=',num2str(dir_gof.rsquare), '  a=',num2str(dir_fitresult.a),'  b=', num2str(dir_fitresult.b), '  c=',num2str(dir_fitresult.c)])
    legend('OFF')

%% PR Sigmoid
    PR_spks_count = zeros(length(trgs_on),2);

    for i = 1:length(trgs_on) 
    
        spks = find(spks_timing > trgs_on(i) + bin_left_sr & spks_timing < trgs_on(i) +bin_right_sr);     % +10 to avoid delayed direct activaiton
    
        PR_spks_count(i,1) = length(spks);

        PR_spks_count(i,2) = stim_amp(i);

    end

    [i,j,k]=unique(PR_spks_count(:,2));
    
    PR_spks_sum=[i accumarray(k,PR_spks_count(:,1),[],@mean)];

% Probablidity sigmoid
    % spks_prob = [(PR_spks_count(:,1) - mean_spks_num(1,2)) > 0 (PR_spks_count(:,2))];
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
    [PR_fitresult, PR_gof] = FitMitch(PR_spks_sum(:,1),PR_spks_sum(:,2));
    xlabel('Amplitude (uA)')
    ylabel('Spikes number')
    xlim([0 130])
    title('PR',['r^2=',num2str(PR_gof.rsquare), '  a=',num2str(PR_fitresult.a),'  b=', num2str(PR_fitresult.b), '  c=',num2str(PR_fitresult.c)])
    legend('OFF')

%% Raster
    raster = [];

    for i = 1:length(trgs_on)

         spks = spks_timing(spks_timing > trgs_on(i) & spks_timing < trgs_on(i) + peroid_dur_sr - 1 );

         spks = spks - trgs_on(i);

         spks(:, 2) = stim_amp(i);

         raster = [raster; spks];

    end
    
    subplot(2,3,[1,2,4,5]);
    patch([bin_left_sr,bin_right_sr,bin_right_sr,bin_left_sr],[-1000,-1000,1000,1000],'black','FaceAlpha',.05)
    hold on
    scatter(raster(:,1),raster(:,2),5,'filled');
    hold off
%     view(2)
    xlim([0 peroid_dur_sr])
    ylim([0 max(stim_amp)])

    xline(pulse_width_sr,'--r',{'pulse offset'})

    xticks([0:peroid_dur_sr/40:peroid_dur_sr])
    xticklabels([0:peroid_dur_ms/40:peroid_dur_ms])
    xlabel('Time (ms)')
    ylabel('Amplitude (uA)')
    title(name)

%% Plot for adjusting findpeak threshold
    figure;
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

%%  Saving
    save(name)
%    saveas(fig,name)

    PR_fitresult_a = PR_fitresult.a;
    PR_fitresult_b = PR_fitresult.b;
    PR_fitresult_c = PR_fitresult.c;                                      

end

%% Discarded
%     amp_trigs_mat = [abs(stim_amp) trgs_on];

%     BC_spks_count = [];
%     for i = 1:length(amp_trigs_mat)
%         trace_clip =  trace(amp_trigs_mat(i, 2) + pulse_width_sr :amp_trigs_mat(i, 2)+bin_left_sr,:);
%         [spks_amp,spks_timing] =  findpeaks(trace_clip(:,1),1,'MINPEAKHEIGHT',peak_threshold_mV, 'MinPeakDistance',peak_distance_sr);
%         spks_count = [length(spks_timing) amp_trigs_mat(i, 1)];
%         BC_spks_count = [BC_spks_count; spks_count];
%     end

%     PR_spks_count = [];
%     for i = 1:length(amp_trigs_mat)
%         trace_clip =  trace(amp_trigs_mat(i, 2)+bin_left_sr:amp_trigs_mat(i, 2)+bin_right_sr,:);
%         [spks_amp,spks_timing] =  findpeaks(trace_clip(:,1),1,'MINPEAKHEIGHT',peak_threshold_mV, 'MinPeakDistance',peak_distance_sr);
%         spks_count = [length(spks_timing) amp_trigs_mat(i, 1)];
%         PR_spks_count = [PR_spks_count; spks_count];
%     end

%     all_spks_timing = [];
%     for i = 1:length(amp_trigs_mat)
%         trace_clip =  trace(amp_trigs_mat(i, 2):amp_trigs_mat(i, 2)+peroid_dur_sr-1,:);
%         [spks_amp,spks_timing] =  findpeaks(trace_clip(:,1),1,'MINPEAKHEIGHT',peak_threshold_mV, 'MinPeakDistance',peak_distance_sr);
%         spks_timing(:, 2) = amp_trigs_mat(i, 1);
%         spks_timing = [spks_timing spks_amp];
%         all_spks_timing = [all_spks_timing;spks_timing];
%     end
    