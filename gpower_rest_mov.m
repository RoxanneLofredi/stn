
close all

clear all
root=fullfile(mdf,'stn_rotameter\STN_Rotameter\');
cd(root);

addpath(fullfile(mdf,'matlab_scripts\wjn_toolbox'));
addpath(fullfile(mdf,'stn_rotameter\rotameter_scripts'));

pat=patienten;
n=0;

%a=4;

for a=1:length(pat)
   if ~isfield(pat{a},'bad') && pat{a}.med;
    
   n=n+1;
   filename=pat{a}.id;
   
   try
   D=spm_eeg_load(['rstf_efspmeeg_' filename '.mat']);
   catch   
   D=spm_eeg_load(['tf_efspmeeg_' filename '.mat']);
   D=wjn_tf_smooth(D.fullfile,1,.3);
   end
        
    hand=stringsplit(filename,'_'); hand = hand{2};
    
    if hand == 'L' ;
       chand = 'R';
    else
       chand= 'L';
    end
    
    channels= ci (['STN' chand], D.chanlabels);
    lconds=ci('onset-gross', D.conditions);
    sconds=ci('onset-klein', D.conditions);
    
    rest_time=[wjn_sc(D.time,-2.0) wjn_sc(D.time,-1.5)];
    mov_time=[wjn_sc(D.time,0) wjn_sc(D.time,.5)];
    
    clear std_all std_rest std_mov pow_all pow_rest pow_mov pow_mov_alltrls
    
    for b=1:length(channels)
        
        std_all=std(squeeze(nanmean(nanmean(D(channels(b),[5:45 55:95],:,lconds),3),4)));
        std_rest=std(squeeze(nanmean(nanmean(D(channels(b),[5:45 55:95],rest_time(1):rest_time(2),lconds),3),4)));
        std_mov=std(squeeze(nanmean(nanmean(D(channels(b),[5:45 55:95],mov_time(1):mov_time(2),lconds),3),4)));
       
        pow_all(b,:)=(squeeze(nanmean(nanmean(D(channels(b),:,:,lconds),3),4))./std_all).*100;
        pow_rest(b,:)=(squeeze(nanmean(nanmean(D(channels(b),:,rest_time(1):rest_time(2),lconds),3),4))./std_rest).*100;
        pow_mov(b,:)=(squeeze(nanmean(nanmean(D(channels(b),:,mov_time(1):mov_time(2),lconds),3),4))./std_mov).*100;
        
        pow_mov_alltrls(b,:,:)=(squeeze(nanmean(D(channels(b),:,mov_time(1):mov_time(2),lconds),3))./std_mov).*100;

        std_sall=std(squeeze(nanmean(nanmean(D(channels(b),[5:45 55:95],:,sconds),3),4)));
        std_srest=std(squeeze(nanmean(nanmean(D(channels(b),[5:45 55:95],rest_time(1):rest_time(2),sconds),3),4)));
        std_smov=std(squeeze(nanmean(nanmean(D(channels(b),[5:45 55:95],mov_time(1):mov_time(2),sconds),3),4)));
        
        pow_sall(b,:)=(squeeze(nanmean(nanmean(D(channels(b),:,:,sconds),3),4))./std_sall).*100;
        pow_srest(b,:)=(squeeze(nanmean(nanmean(D(channels(b),:,rest_time(1):rest_time(2),sconds),3),4))./std_srest).*100;
        pow_smov(b,:)=(squeeze(nanmean(nanmean(D(channels(b),:,mov_time(1):mov_time(2),sconds),3),4))./std_smov).*100;
    end
    
    gpower.large.all(n,:)=log(nanmean(pow_all));
    gpower.large.rest(n,:)=log(nanmean(pow_rest));
    gpower.large.mov(n,:)=log(nanmean(pow_mov));
    
    gpower.large.movalltrl(n,:,:)=log(squeeze(nanmean(pow_mov_alltrls,1)));
    
    gpower.small.all(n,:)=log(nanmean(pow_sall));
    gpower.small.rest(n,:)=log(nanmean(pow_srest));
    gpower.small.mov(n,:)=log(nanmean(pow_smov));
   
   end
end

keep gpower
save gpower

gpower.large.delta=gpower.large.mov-gpower.large.rest;
gpower.small.delta=gpower.small.mov-gpower.small.rest;

for c=1:size(gpower.large.mov,2)
    gpower.large.pvalues(c)=wjn_ppt(gpower.large.mov(:,c),gpower.large.rest(:,c));
end


%% Gamma Figures

cc=colorlover(6);
color=cc(4,:);

signstat=gpower.large.pvalues<0.05;
signstat=signstat.*10;

figure,
rest=mypower(gpower.large.rest)
hold on,
mov=mypower(1:400,gpower.large.mov, color)
sigbar(gpower.frequency,gpower.large.pvalues<0.05)
xlim([2 10])
ylim([5 8])
xlabel('Frequency [Hz]')
ylabel('Power [a.u.]')
legend([rest mov],'rest','mov','location','northeast')
title('Gamma Power rest mov p<.05')
figone(7)
myprint('gamma_dynamics_rest_mov_large')

figure,
rest=mypower(gpower.large.rest(6,:))
hold on,
mov=mypower(1:400,gpower.large.mov(6,:), color)
% sigbar(gpower.frequency,gpower.large.pvalues<0.05)
xlim([10 100])
ylim([0 5])
xlabel('Frequency [Hz]')
ylabel('Power [a.u.]')
legend([rest mov],'rest','mov','location','northeast')
title('Gamma Power rest mov pat6')
figone(7)
myprint('gamma_dynamics_rest_mov_large_pat6')

figure,
rest=mypower(gpower.small.rest)
hold on,
mov=mypower(1:400,gpower.small.mov, color)
xlim([10 150])
xlabel('Frequency [Hz]')
ylabel('Power [a.u.]')
legend([rest mov],'rest','mov','location','northeast')
title('Gamma Power rest mov small')
figone(7)
myprint('gamma_dynamics_rest_mov_small')

figure,
rest=mypower(gpower.large.delta(6,:))
xlim([30 100])
ylim([0 1.7])
xlabel('Frequency [Hz]')
ylabel(' ? Power Mov-Rest [a.u.]')
title(' ? Gamma Power rest mov pat6')
figone(7)
myprint('gamma_delta_large_pat6')
