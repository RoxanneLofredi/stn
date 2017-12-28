
close all

clear all
root=fullfile(mdf,'stn_rotameter\STN_Rotameter\');
cd(root);

addpath(fullfile(mdf,'matlab_scripts\wjn_toolbox'));
addpath(fullfile(mdf,'stn_rotameter\rotameter_scripts'));

pat=patienten;
n=0;
m=0;

% a=17

for a=1:length(pat)
    
   if ~pat{a}.med && ~isfield(pat{a},'bad') && isfield(pat{a},'offseven') && a~=15;
    filename_off=pat{a}.id;
    Doff=spm_eeg_load(['rstf_efspmeeg_' filename_off '.mat']);
%     Doff=wjn_tf_smooth(Doff.fullfile,1,.3);
    
%     try
    filename_on=[pat{a}.id(1:end-3) 'ON'];
    Don=spm_eeg_load(['rstf_efspmeeg_' filename_on '.mat']);
%     Don=wjn_tf_smooth(Don.fullfile,1,.3);
    
    hand=stringsplit(filename_on,'_'); hand = hand{2};
    if hand == 'L' ;
       chand = 'R';
    else
       chand= 'L';
    end
    
    channels_on= ci (['STN' chand], Don.chanlabels);
    channels_off= ci (['STN' chand], Doff.chanlabels);
    
    conds_on=ci('onset', Don.conditions);
    conds_off=ci('onset', Doff.conditions);
    
    rest_time=[wjn_sc(Don.time,-2.0) wjn_sc(Don.time,-1.5)];
    mov_time=[wjn_sc(Don.time,0) wjn_sc(Don.time,.5)];
    
    for b=1:length(channels_on)
    
    std_rawsignal_on(n,b)=std(squeeze(nanmean(nanmean(nanmean(Don(channels_on(b),:,:,conds_on),1),4),3)));
    std_rawsignal_off(n,b)=std(squeeze(nanmean(nanmean(nanmean(Doff(channels_off(b),:,:,conds_off),1),4),3)));
    
    std_rawsignal_rest_on(n,b)=std(squeeze(nanmean(nanmean(nanmean(Don(channels_on(b),:,rest_time(1):rest_time(2),conds_on),1),4),3)));
    std_rawsignal_rest_off(n,b)=std(squeeze(nanmean(nanmean(nanmean(Doff(channels_off(b),:,rest_time(1):rest_time(2),conds_off),1),4),3)));
    
    std_rawsignal_mov_on(n,b)=std(squeeze(nanmean(nanmean(nanmean(Don(channels_on(b),:,mov_time(1):mov_time(2),conds_on),1),4),3)));
    std_rawsignal_mov_off(n,b)=std(squeeze(nanmean(nanmean(nanmean(Doff(channels_off(b),:,mov_time,conds_off),1),4),3)));

    nrawsignal_on(n,b,:)=(squeeze(nanmean(nanmean(nanmean(Don(channels_on(b),:,:,conds_on),1),4),3)))./std_rawsignal_on(n,b);
    nrawsignal_off(n,b,:)=(squeeze(nanmean(nanmean(nanmean(Doff(channels_off(b),:,:,conds_off),1),4),3)))./std_rawsignal_off(n,b);
    
    nrawsignal_rest_on(n,b,:)=(squeeze(nanmean(nanmean(nanmean(Don(channels_on(b),:,rest_time(1):rest_time(2),conds_on),1),4),3)))./std_rawsignal_rest_on(n,b);
    nrawsignal_rest_off(n,b,:)=(squeeze(nanmean(nanmean(nanmean(Don(channels_off(b),:,rest_time(1):rest_time(2),conds_off),1),4),3)))./std_rawsignal_rest_off(n,b);
    
    nrawsignal_mov_on(n,b,:)=(squeeze(nanmean(nanmean(nanmean(Don(channels_on(b),:,mov_time(1):mov_time(2),conds_on),1),4),3)))./std_rawsignal_mov_on(n,b);
    nrawsignal_mov_off(n,b,:)=(squeeze(nanmean(nanmean(nanmean(Don(channels_off(b),:,mov_time(1):mov_time(2),conds_off),1),4),3)))./std_rawsignal_mov_off(n,b);
    
    end
    
    beta.on(n,:)=squeeze(nanmean(nanmean(nanmean(Don(channel_on,13:30,:,conds_on),1),2),4));
    beta.off(n,:)=squeeze(nanmean(nanmean(nanmean(Doff(channel_off,13:30,:,conds_off),1),2),4));
    
    lowbeta.on(n,:)=squeeze(nanmean(nanmean(nanmean(Don(channel_on,13:20,:,conds_on),1),2),4));
    lowbeta.off(n,:)=squeeze(nanmean(nanmean(nanmean(Doff(channel_off,13:20,:,conds_off),1),2),4));
    
    highbeta.on(n,:)=squeeze(nanmean(nanmean(nanmean(Don(channel_on,20:30,:,conds_on),1),2),4));
    highbeta.off(n,:)=squeeze(nanmean(nanmean(nanmean(Doff(channel_off,20:30,:,conds_off),1),2),4));
    
    nbeta.on(n,:)=beta.on(n,:)./std_rawsignal_on(n).*100;
    nbeta.off(n,:)=beta.off(n,:)./std_rawsignal_off(n).*100;
    
    nbeta.low.on(n,:)=lowbeta.on(n,:)./std_rawsignal_on(n).*100;
    nbeta.low.off(n,:)=lowbeta.off(n,:)./std_rawsignal_off(n).*100;
    
    nbeta.high.on(n,:)=highbeta.on(n,:)./std_rawsignal_on(n).*100;
    nbeta.high.off(n,:)=highbeta.off(n,:)./std_rawsignal_off(n).*100;
    
%     catch
%         m=m+1;
%         notfound(m)=a;
    end
 end
% end

for c=1:size(nbeta.on,1)
    
    log_nrawsignal_on(c,:)=log(nrawsignal_on(c,:));
    log_nrawsignal_off(c,:)=log(nrawsignal_off(c,:));
    
    log_nrawsignal_rest_on(c,:)=log(nrawsignal_rest_on(c,:));   
    log_nrawsignal_rest_off(c,:)=log(nrawsignal_rest_off(c,:));  
    
    log_nrawsignal_mov_on(c,:)=log(nrawsignal_mov_on(c,:));  
    log_nrawsignal_mov_off(c,:)=log(nrawsignal_mov_off(c,:));  
end
    
nbeta.time=Doff.time;
beta.time=Doff.time;

cc=colorlover(6);
color=cc(4,:);

figure,
ON=mypower(nbeta.time,nbeta.on)
hold on,
OFF=mypower(nbeta.time,nbeta.off, color)
xlim([-3 4])
xlabel('Time [s]')
ylabel('Power [a.u.]')
legend([ON OFF],'ON','OFF','location','norteast')
title('Beta Dynamics ON OFF')
figone(7)
myprint('beta_dynamics_on_off_std')

figure,
ON=mypower(nbeta.time,nbeta.low.on)
hold on,
OFF=mypower(nbeta.time,nbeta.low.off, color)
xlim([-3 4])
xlabel('Time [s]')
ylabel('Power [a.u.]')
legend([ON OFF],'ON','OFF','location','northwest')
title('Low Beta (13-20 Hz) ON OFF')
figone(7)
myprint('lowbeta_dynamics_on_off_std')


figure,
ON=mypower(nbeta.time,nbeta.high.on)
hold on,
OFF=mypower(nbeta.time,nbeta.high.off, color)
xlim([-3 4])
xlabel('Time [s]')
ylabel('Power [a.u.]')
legend([ON OFF],'ON','OFF','location','northeast')
title('High Beta (20-30 Hz) ON OFF')
figone(7)
myprint('highbeta_dynamics_on_off_std')

figure,
ON=mypower(1:400,log_nrawsignal_on)
hold on,
OFF=mypower(1:400,log_nrawsignal_off, color)
xlim([3 30])
xlabel('Frequency [Hz]')
ylabel('Power [a.u.]')
legend([ON OFF],'ON','OFF','location','northeast')
title('Power ON OFF over trl')
figone(7)
myprint('log_power_onoff_trl')

figure,
ON=mypower(1:400,log_nrawsignal_rest_on)
hold on,
OFF=mypower(1:400,log_nrawsignal_rest_off, color)
xlim([3 30])
xlabel('Frequency [Hz]')
ylabel('Power [a.u.]')
legend([ON OFF],'ON','OFF','location','northeast')
title('Power ON OFF rest')
figone(7)
myprint('log_power_onoff_rest')

figure,
ON=mypower(1:400,nrawsignal_mov_on)
hold on,
OFF=mypower(1:400,nrawsignal_mov_off, color)
xlim([3 30])
xlabel('Frequency [Hz]')
ylabel('Power [a.u.]')
legend([ON OFF],'ON','OFF','location','northeast')
title('Power ON OFF mov')
figone(7)
myprint('power_onoff_mov')