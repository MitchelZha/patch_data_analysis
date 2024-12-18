%             %% Test
%             clc
%             clear all
%  
%             %%
%             addpath 'D:\OneDrive - UNSW\ephys'; addpath 'D:\OneDrive - UNSW\ephys\patch_data_analysis';
%             freq_Hz = 25;
%             phase_width_ms = 1;
%             nkt = 50;
%            
%             cd 'D:\OneDrive - UNSW\ephys\220406';
%             recording_dir = 'Clampex\2022_04_06_0007.abf';
%             stim_dir = 'Mitch_Fixedwn_Freq=25_Mean=52_contrast=32';
%             name = ['220406 ONS 1 last min ' stim_dir '.mat'];
%            
%             crop_ratio = [0,1];
%             peak_threshold_mV = -15;     
%             peak_distance_sr = 35;
%             bin_left_sr = []; 
%             bin_right_sr = [];
%             % [nt_sta, BC_sta, PR_sta, stim_mean] = fix_fre_sta(freq_Hz, phase_width_ms, nkt, peak_threshold_mV, peak_distance_sr, recording_dir, stim_dir, name, crop_ratio, bin_left_sr, bin_right_sr);
%             % sta_playback(25, 1, 10, nt_sta, BC_sta, PR_sta, stim_mean);



%             %%
%             freq_Hz = 25;
%             phase_width_ms = 0.5;
%             nkt = 50;
%             peak_threshold_mV = -40;     
%             peak_distance_sr = 35;
%             addpath 'D:\OneDrive - UNSW\ephys'; addpath 'D:\OneDrive - UNSW\ephys\patch_data_analysis';
%             cd 'D:\OneDrive - UNSW\ephys\220304';
%             recording_dir = 'Clampex\2022_03_04_0035.abf';
%             stim_dir = 'Mitch_Fixedwn_Freq=25_Mean=50_contrast=17.05';
%             bin_left_sr = []%100;
%             bin_right_sr = []%500;
%             name = ['220304  ONS ' stim_dir '.mat'];
%             crop_ratio = [0,1];             
%             %[nt_sta, BC_sta, PR_sta, stim_mean] = fix_fre_sta(freq_Hz, phase_width_ms, nkt, peak_threshold_mV, peak_distance_sr, recording_dir, stim_dir, name, crop_ratio, bin_left_sr, bin_right_sr);             
%             %sta_playback(25, 1, 10, nt_sta, BC_sta, PR_sta, stim_mean);
% 
%             stim_amp = stim_amp(1:2811);
%             trgs_on = trgs_on(1:2811);
    



%             %%
%             addpath 'D:\OneDrive - UNSW\ephys'; addpath 'D:\OneDrive - UNSW\ephys\patch_data_analysis';
%             freq_Hz = 25;
%             phase_width_ms = 1;
%             nkt = 50;
%            
%             cd 'D:\OneDrive - UNSW\ephys\220404';
%             recording_dir = 'Clampex\2022_04_04_0020.abf';
%             stim_dir = 'Mitch_Fixedwn_Freq=25_Mean=85_contrast=32';
%             name = ['220404 OFF 1 ' stim_dir '.mat'];
%            
%             crop_ratio = [0,1];
%             bin_left_sr = 330; 
%             bin_right_sr = [];
%             peak_threshold_mV = -35;     
%             peak_distance_sr = 35;
%             % [nt_sta, BC_sta, PR_sta, stim_mean] = fix_fre_sta(freq_Hz, phase_width_ms, nkt, peak_threshold_mV, peak_distance_sr, recording_dir, stim_dir, name, crop_ratio, bin_left_sr, bin_right_sr);
%             % sta_playback(25, 1, 10, nt_sta, BC_sta, PR_sta, stim_mean);
% 
%             stim_amp = stim_amp(1:4000);
%             trgs_on = trgs_on(1:4000);




%             %%
%             freq_Hz = 25;
%             phase_width_ms = 1;
%             nkt = 50;
%             crop_ratio = [0,1];
%             peak_threshold_mV = -20;     
%             peak_distance_sr = 35;
%             
%             addpath 'D:\OneDrive - UNSW\ephys'; addpath 'D:\OneDrive - UNSW\ephys\patch_data_analysis'
%             cd 'D:\OneDrive - UNSW\ephys\220406';
%             recording_dir = 'Clampex\2022_04_06_0014.abf';
%             stim_dir = 'Mitch_Fixedwn_Freq=25_Mean=52_contrast=20';
%             name = ['220406 ONS ' stim_dir '.mat'];
%             bin_left_sr = []; bin_right_sr =[];
%             
%             % [nt_sta, BC_sta, PR_sta, stim_mean] = fix_fre_sta(freq_Hz, phase_width_ms, nkt, peak_threshold_mV, peak_distance_sr, recording_dir, stim_dir, name, crop_ratio);
%             % sta_playback(25, 1, 5, nt_sta, BC_sta, PR_sta, stim_mean);
% 
%             stim_amp = [stim_amp(1:1375); stim_amp(1700:end)]
%             trgs_on = [trgs_on(1:1375); trgs_on(1700:end)]





%             %%
%             freq_Hz = 25;
%             phase_width_ms = 1;
%             nkt = 50;
%             crop_ratio = 0.67;
%             peak_threshold_mV = -20;     
%             peak_distance_sr = 35;
%             
%             addpath 'D:\OneDrive - UNSW\ephys'; addpath 'D:\OneDrive - UNSW\ephys\patch_data_analysis'
%             cd 'D:\OneDrive - UNSW\ephys\220406';
%             recording_dir = 'Clampex\2022_04_06_0029.abf';
%             stim_dir = 'Mitch_Fixedwn_Freq=25_Mean=52_contrast=32';
%             name = ['220406 OFFT AD ' stim_dir '.mat'];
%             bin_left_sr = []; bin_right_sr =[];
%             
%             % [nt_sta, BC_sta, PR_sta, stim_mean] = fix_fre_sta(freq_Hz, phase_width_ms, nkt, peak_threshold_mV, peak_distance_sr, recording_dir, stim_dir, name, crop_ratio);
%             % sta_playback(25, 1, 5, nt_sta, BC_sta, PR_sta, stim_mean);
% 
% 
%             stim_amp = stim_amp(1:2563);
%             trgs_on = trgs_on(1:length(trgs_on)*crop_ratio);



%%
function [nt_sta, BC_sta, PR_sta, stim_mean]=fix_fre_sta(freq_Hz, phase_width_ms, nkt, peak_threshold_mV, peak_distance_sr, recording_dir, stim_dir, name, crop_ratio, bin_left_sr, bin_right_sr)
%% Fomular
    close all

    sample_rate = 50000;

    stim_amp = read_stim_file([stim_dir '.txt']);
    [trace] = abfload(recording_dir);
    
    peroid_dur_ms = (1/freq_Hz)*1000;
    peroid_dur_sr = peroid_dur_ms*0.001*sample_rate;
    phase_width_sr = phase_width_ms*0.001*sample_rate;
    pulse_width_sr = phase_width_sr *2;

    if isempty(bin_left_sr)
        bin_left_sr = sample_rate * 0.008 + pulse_width_sr;                                       %  ML are within 8ms Jeson and rizzo 
    end
    
    if isempty(bin_right_sr)
        bin_right_sr = peroid_dur_sr;                                       
    end



%% Processing recording
    ttls = find(trace(:,2)>2);
    
    trgs_on = ttls(find(diff(ttls)>pulse_width_sr*3)+1); 
    
    trgs_on = [ttls(1); trgs_on];
    
%     missing_trgs = find(diff(trgs_on) > peroid_dur_sr+5);                                  %  Turn on when there is missing trigger caused by 0 amplitude
%     for i = 1:length(missing_trgs)
%         missing_trgs = find(diff(trgs_on) > peroid_dur_sr+5);
%         trgs_on = [trgs_on(1:missing_trgs(1)); trgs_on(missing_trgs(1)) + peroid_dur_sr+5; trgs_on(missing_trgs(1)+1:end)];
%     end

    stim_amp = abs(stim_amp);


%% Crop the whitenoise
    stim_amp = stim_amp(1+length(stim_amp)*crop_ratio(1) : length(stim_amp)*crop_ratio(2));
    trgs_on = trgs_on(1+length(trgs_on)*crop_ratio(1) : length(trgs_on)*crop_ratio(2));


    stim_mean = mean(stim_amp);

    [spks_amp, spks_timing] = findpeaks(trace(1:end,1),'MINPEAKHEIGHT',peak_threshold_mV, 'MinPeakDistance',peak_distance_sr);
    
    
%%  Network(BC+PR) STA
    nt_spks_count = zeros(1,length(trgs_on));

    for i = 1:length(trgs_on)
    
        spks = find(spks_timing > trgs_on(i) + pulse_width_sr + (0.0005*sample_rate) & spks_timing < trgs_on(i) + peroid_dur_sr );   % adding 0.5ms to avoiding delayed direct activaiton
    
        nt_spks_count(i) = length(spks);
    
    end
    
    [nt_sta, nt_stc, mu, cov] = simpleSTC_hamed(stim_amp, nt_spks_count', nkt);
    
    tvec = (-nkt/2+1:nkt/2)'*1/freq_Hz-.5/freq_Hz;                                       % vector of time indices (in units of stim frames)
    
    fig = figure;
    subplot(2,2,2);
    patch([-2 0 0 -2],[1200 1200 -15  -15],'black','FaceAlpha',.05)
    hold on
    plot(tvec, nt_sta,'LineWidth',2)
    line([tvec(1),tvec(end)],[stim_mean ,stim_mean],'Color','k','LineStyle','--')
    xlabel('time before spike (sec)'); 
    ylabel('E-STA (A.U)');
    ylim([stim_mean-20, stim_mean+20])
    xlim([-1,.5])
    title('Network (BC+PR)', ['spike count=', num2str(sum(nt_spks_count))])
    hold off


%% BC STA
    BC_spks_count = zeros(1,length(trgs_on));

    for i = 1:length(trgs_on) 
    
        spks = find(spks_timing > trgs_on(i) + pulse_width_sr + 10 & spks_timing < trgs_on(i) + bin_left_sr);    % +10 to avoid delayed direct activaiton
    
        BC_spks_count(i) = length(spks);
    
    end
    
    [BC_sta, BC_stc, mu, cov] = simpleSTC_hamed(stim_amp, BC_spks_count', nkt);
    
    tvec = (-nkt/2+1:nkt/2)'*1/freq_Hz-.5/freq_Hz;                                       % vector of time indices (in units of stim frames)
    
    subplot(2,2,3);
    patch([-2 0 0 -2],[1200 1200 -15 -15],'black','FaceAlpha',.05)
    hold on
    plot(tvec, BC_sta,'LineWidth',2)
    line([tvec(1),tvec(end)],[stim_mean ,stim_mean],'Color','k','LineStyle','--')
    xlabel('time before spike (sec)'); 
    ylabel('E-STA (A.U)');
    ylim([stim_mean-20, stim_mean+20])
    xlim([-1,.5])
    title('BC', ['spike count=', num2str(sum(BC_spks_count))])
    hold off

%% PR STA
    PR_spks_count = zeros(1,length(trgs_on));

    for i = 1:length(trgs_on) 
    
        spks = find(spks_timing > trgs_on(i) + bin_left_sr & spks_timing < trgs_on(i) + bin_right_sr);
    
        PR_spks_count(i) = length(spks);
    
    end
    
    [PR_sta, PR_stc, mu, cov] = simpleSTC_hamed(stim_amp, PR_spks_count', nkt);
    
    tvec = (-nkt/2+1:nkt/2)'*1/freq_Hz-.5/freq_Hz;                                       % vector of time indices (in units of stim frames)
    
    subplot(2,2,4);
    patch([-2 0 0 -2],[1200 1200 -15  -15],'black','FaceAlpha',.05)
    hold on
    plot(tvec, PR_sta,'LineWidth',2)
    line([tvec(1),tvec(end)],[stim_mean ,stim_mean],'Color','k','LineStyle','--')
    xlabel('time before spike (sec)'); 
    ylabel('E-STA (A.U)');
    ylim([stim_mean-20, stim_mean+20])
    xlim([-1,.5])
    title('PR', ['spike count=', num2str(sum(PR_spks_count))])
    hold off


%% Raster
    %raster = zeros(length(spks_timing),2);
    raster = [];

    for i = 1:length(trgs_on)

         spks = spks_timing(spks_timing > trgs_on(i) & spks_timing < trgs_on(i) + peroid_dur_sr - 1 );

         spks = spks - trgs_on(i);

         spks(:, 2) = stim_amp(i);

         raster = [raster; spks];
    
%         trace_clip =  trace(trgs_on(i):trgs_on(i) + peroid_dur_sr - 1,:);
%     
%         [raster_spks_amp,raster_spks_timing] =  findpeaks(trace_clip(:,1),1,'MINPEAKHEIGHT',peak_threshold_mV, 'MinPeakDistance',peak_distance_sr);
%     
%         raster_spks_timing(:, 2) = stim_amp(i);
%     
%         raster = [raster; raster_spks_timing];
    
    end
    
    subplot(2,2,1);
    patch([bin_left_sr,bin_right_sr,bin_right_sr,bin_left_sr],[-1,-1,1000,1000],'black','FaceAlpha',0.05)
    hold on
    scatter(raster(:,1),raster(:,2),5,'filled')
    hold off
    xlim([0 peroid_dur_sr])
    ylim([0 max(stim_amp)])
    xline(pulse_width_sr,'--r',{'pulse offset'})
    xticks([0:peroid_dur_sr/40:peroid_dur_sr])
    xticklabels([0:peroid_dur_ms/40:peroid_dur_ms])
    xlabel('Time (ms)')
    ylabel('Amplitude(uA)')
    title(name)

%% multipulse Raster
    %       
    % mul_stim_amp =  movsum(stim_amp,2);
    % 
    % mul_raster_spks_count = [];
    % for i = 2:length(mul_stim_amp)
    % 
    %     trace_clip =  trace(trgs_on(i):trgs_on(i) + peroid_dur_sr - 1,:);
    % 
    %     [spks_amp,spks_timing] =  findpeaks(trace_clip(:,1),1,'MINPEAKHEIGHT',peak_threshold_mV, 'MinPeakDistance',peak_distance_sr);
    % 
    %     spks_timing(:, 2) = mul_stim_amp(i);
    % 
    %     mul_raster_spks_count = [mul_raster_spks_count; spks_timing];
    % 
    % end
    % 
    % figure;
    % 
    % %[C, ia] = unique(mul_raster_spks_count(:,2),'rows');
    % %scatter(mul_raster_spks_count(ia,1),mul_raster_spks_count(ia,2),5,'filled')
    % 
    % scatter(mul_raster_spks_count(:,1),mul_raster_spks_count(:,2),5,'filled')
    % hold on
    % patch([bin_left_sr,bin_right_sr,bin_right_sr,bin_left_sr],[-1,-1,1000,1000],'black','FaceAlpha',0.05)
    % hold off
    % xlim([0 peroid_dur_sr])
    % ylim([0 max(mul_stim_amp)])
    % xline(pulse_width_sr,'--r')
    % xticks([0:0.001*sample_rate:peroid_dur_sr])
    % xticklabels([0:1:peroid_dur_ms])
    % xlabel('Time (ms)')
    % ylabel('Amplitude(uA)')
    % title(['multipulse Raster',name])


%% For adjusting threshold
    figure;
    plot(trace)
    hold on
    plot(spks_timing,spks_amp,'o')
    xticks([0:0.01*sample_rate:length(trace)])
    xticklabels([0:0.01:(length(trace)/sample_rate)])
    xlim([trgs_on(end-5) trgs_on(end)])
    title(name)
    xlabel('Time (sec)')
    ylabel('Memberine potential (mV)')
    hold off

%% Saving
    save(name)
%   saveas(fig, name)
    
    nt_sta= nt_sta(1:nkt/2+2);
    BC_sta= BC_sta(1:nkt/2+2); 
    PR_sta= PR_sta(1:nkt/2+2);

end
