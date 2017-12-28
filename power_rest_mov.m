close all

clear all
root=fullfile(mdf,'stn_rotameter\STN_Rotameter\');
cd(root);

addpath(fullfile(mdf,'matlab_scripts\wjn_toolbox'));
addpath(fullfile(mdf,'stn_rotameter\rotameter_scripts'));

pat=patienten;
n=0;

for a=1:length(pat);

    if ~isfield(pat{a},'bad') && pat{a}.med;
        
        n=n+1;
        
        filename=pat{a}.id;

        D=spm_eeg_load(['tf_efspmeeg_' filename '.mat']);
        D=wjn_tf_smooth(D.fullfile,1,.3);
        
        hand=stringsplit(filename,'_'); 
        hand = hand{2};

        if hand == 'L' ;
            chand = 'R';
        else
            chand= 'L';
        end
    
        channel= ci (['STN' chand], D.chanlabels);
        
        rest_time=[wjn_sc(D.time,-2.0) wjn_sc(D.time,-1.5)];
        mov_time=[wjn_sc(D.time,0) wjn_sc(D.time,.5)];
        svel_conds=ci('onset-klein',D.conditions);
        lvel_conds=ci('onset-gross',D.conditions);
        
        rest_pow=squeeze(nanmean(nanmean(nanmean(D(channel,:,rest_time(1):rest_time(2),svel_conds),4),3),1));
        smov_pow=squeeze(nanmean(nanmean(nanmean(D(channel,:,mov_time(1):mov_time(2),svel_conds),4),3),1));
        lmov_pow=squeeze(nanmean(nanmean(nanmean(D(channel,:,mov_time(1):mov_time(2),lvel_conds),4),3),1));
        
%         rest_pow(46:54)=nanmean(rest_pow([45 55]));
%         mov_pow(46:54)=nanmean(mov_pow([45 55]));

        rest_npow=rest_pow./std(rest_pow).*100;
        smov_npow=smov_pow./std(rest_pow).*100;
        lmov_npow=lmov_pow./std(rest_pow).*100;
        
%         
%         rest_npow_sum=rest_pow./nansum(rest_pow'.*100);
%         smov_npow_sum=smov_pow./nansum(rest_pow'.*100);
%         lmov_npow_sum=lmov_pow./nansum(rest_pow'.*100);
%         
        log_all_rest_npow(n,:)=log(rest_npow);
        log_all_smov_npow(n,:)=log(smov_npow);
        log_all_lmov_npow(n,:)=log(lmov_npow);
        
        [freqmax,log_freqpeak]=max(log(lmov_npow(50:90)));
        log_freqpeak = log_freqpeak+50;
        
        [freqmax,freqpeak]=max(lmov_npow(50:90));
        freqpeak = freqpeak+50;
        
        all_freqpeak(n)= freqpeak;
        all_log_freqpeak(n)= log_freqpeak;
%         all_rest_pow(n,:)=rest_pow;
%         all_mov_pow(n,:)=mov_pow;
        
        all_rest_npow(n,:)=rest_npow;
        all_smov_npow(n,:)=smov_npow;
        all_lmov_npow(n,:)=lmov_npow;
    end
end

%% exemplary patient
filename=pat{6}.id;
D=spm_eeg_load(['tf_efspmeeg_' filename '.mat']);
channel= ci ('STNL', D.chanlabels);
rest_time=[wjn_sc(D.time,-2.0) wjn_sc(D.time,-1.5)];
mov_time=[wjn_sc(D.time,0) wjn_sc(D.time,.5)];
svel_conds=ci('onset-klein',D.conditions);
lvel_conds=ci('onset-gross',D.conditions);

pat_rest_pow=squeeze(nanmean(nanmean(D(channel,:,rest_time(1):rest_time(2),svel_conds),3),1));
pat_mov_pow=squeeze(nanmean(nanmean(D(channel,:,mov_time(1):mov_time(2),svel_conds),3),1));

% pat_rest_npow=log(pat_rest_pow./std(pat_rest_pow,0,2).*100);
% pat_mov_npow=log(pat_mov_pow./std(pat_rest_pow,0,2).*100);

%%
delta_lpow=log_all_lmov_npow-log_all_rest_npow;
delta_spow=log_all_smov_npow-log_all_rest_npow;

powr=log_all_rest_npow;
powm=log_all_lmov_npow;

cc=colorlover(6);
color=cc(4,:);

figure,
subplot(2,2,1)
        rest=mypower(1:400,powr(6,:)),
    hold on
        mov=mypower(1:400,powm(6,:),color)
        xlim([5 150])
        ylim([0 5])
        xlabel('Frequency [Hz]')
        ylabel('Spectral Power [a.u.]')
        legend([rest mov],'rest','mov','location','northeast')
        title('Exemplary patient')
subplot(2,2,3)
        mypower(delta_spow(6,:))
        xlim([33 100])
        ylim([0 0.45])
        xlabel('Frequency [Hz]')
        ylabel('? Mov-Rest Power [a.u.]')
subplot(2,2,2)
        rest=mypower(log_all_rest_npow),
    hold on
        mov=mypower(1:400,log_all_lmov_npow,color)
        xlim([5 150])
        ylim([1.5 4.5])
        xlabel('Frequency [Hz]')
        ylabel('Spectral Power [a.u.]')
        legend([rest mov],'rest','mov','location','northeast')
        title('Across patients')
subplot(2,2,4)
        mypower(delta_spow)
        xlim([33 100])
        ylim([0 0.23])
        xlabel('Frequency [Hz]')
        ylabel('? Mov-Rest Power [a.u.]')
figone(14,17)
myprint('allvsex_powerspectra_restmov_onset_0-5')

% all_rest_npow(:,50)=nan;
% all_mov_npow(:,50)=nan;
% 
% log_all_rest_pow=log(all_rest_npow(1,:));
% [Y,I] = max(log_all_mov_npow(2, 40:90))
% 
% figure,
% plot(log_all_rest_pow)
% 
% figure,
% mypower(log_all_rest_npow),
% hold on
% mypower(log_all_mov_npow),
% xlim([40 150])

%% BL corrected

n=0;

for a=1:length(pat);

    if ~isfield(pat{a},'bad') && ~pat{a}.med;
        
        n=n+1;
        
        filename=pat{a}.id;
        
        D=spm_eeg_load(['bgmtf_efspmeeg_' filename '.mat']);
        
        mov_time=[wjn_sc(D.time,0) wjn_sc(D.time,.5)];
        vel_conds=ci('velocity-gross',D.conditions);
        
        bl_mov_pow=squeeze(nanmean(nanmean(nanmean(D(:,:,mov_time(1):mov_time(2),vel_conds),4),3),1));
        all_bl_mov_pow_off(n,:)=bl_mov_pow;
        
    end
end

cc=colorlover(6);

log_power_rest=log(all_rest_npow);
log_power_mov=log(all_mov_npow);

figure,
mypower(log_power_rest), 
hold on, 
mypower(log_power_mov), 

figure,
mypower(all_bl_mov_pow),
hold on
mypower(all_bl_mov_pow_off),
xlim([30 120])
ylim([0 25])

figure,
mypower(f,all_rest_npow,cc(4,:),'sem'), 
hold on, 
mypower(f,all_mov_npow,cc(2,:),'sem'),
xlim([40 120])

% figure,
% plot(rest_pow)
% hold on
% figure, plot(rest_npow)