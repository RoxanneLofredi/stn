
close all

clear all
root=fullfile(mdf,'stn_rotameter\STN_Rotameter\');
cd(root);

addpath(fullfile(mdf,'matlab_scripts\wjn_toolbox'));
addpath(fullfile(mdf,'stn_rotameter\rotameter_scripts'));

pat=patienten;
n=0;

for a=1:length(pat)
   if ~pat{a}.med && ~isfield(pat{a},'bad') && isfield(pat{a},'offseven') && a~=15;
    
    n=n+1;
    
    filename_off=pat{a}.id;
    Doff=spm_eeg_load(['stf_efspmeeg_' filename_off '.mat']);

    if a~=55
    filename_on=[pat{a}.id(1:end-3) 'ON'];
    Don=spm_eeg_load(['stf_efspmeeg_' filename_on '.mat']);
    else
    filename_on=[pat{a}.id(1:end-9) 'ON'];
    Don=spm_eeg_load(['stf_efspmeeg_' filename_on '_merge.mat']);
    end
    
    hand=stringsplit(filename_on,'_'); hand = hand{2};
    if hand == 'L' ;
       chand = 'R';
    else
       chand= 'L';
    end
    
    channels_on= ci (['STN' chand], Don.chanlabels);
    channels_off= ci (['STN' chand], Doff.chanlabels);
    
    conds_on=ci('velocity', Don.conditions);
    conds_off=ci('velocity', Doff.conditions);
    
    rest_time=[wjn_sc(Don.time,-2.0) wjn_sc(Don.time,-1)];
    mov_time=[wjn_sc(Don.time,-.25) wjn_sc(Don.time,.25)];
    
    clear std_all_on std_rest_on std_mov.on pow_all_on pow_rest_on pow_mov_on beta.on beta.high.on beta.low.on ...
            std_all_off std_rest_off std_mov.off pow_all_off pow_rest_off pow_mov_off beta.off beta.high.off beta.low.off
    
    for b=1:length(channels_on)
        
        std_all_on=std(squeeze(nanmean(nanmean(Don(channels_on(b),[5:45 55:95],:,conds_on),3),4)));
        std_rest_on=std(squeeze(nanmean(nanmean(Don(channels_on(b),[5:45 55:95],rest_time(1):rest_time(2),conds_on),3),4)));
        std_mov_on=std(squeeze(nanmean(nanmean(Don(channels_on(b),[5:45 55:95],mov_time(1):mov_time(2),conds_on),3),4)));
        
        pow_all_on(b,:)=(squeeze(nanmean(nanmean(Don(channels_on(b),:,:,conds_on),3),4))./std_all_on).*100;
        pow_rest_on(b,:)=(squeeze(nanmean(nanmean(Don(channels_on(b),:,rest_time(1):rest_time(2),conds_on),3),4))./std_rest_on).*100;
        pow_mov_on(b,:)=(squeeze(nanmean(nanmean(Don(channels_on(b),:,mov_time(1):mov_time(2),conds_on),3),4))./std_mov_on).*100;
        
        beta.on(b,:)=(squeeze(nanmean(nanmean(Don(channels_on(b),13:30,:,conds_on),2),4))./std_all_on).*100;
        beta.high.on(b,:)=(squeeze(nanmean(nanmean(Don(channels_on(b),20:30,:,conds_on),2),4))./std_all_on).*100;
        beta.low.on(b,:)=(squeeze(nanmean(nanmean(Don(channels_on(b),13:20,:,conds_on),2),4))./std_all_on).*100;
    
    end
    
    for b=1:length(channels_off)
        
        std_all_off=std(squeeze(nanmean(nanmean(Doff(channels_off(b),[5:45 55:95],:,conds_off),3),4)));
        std_rest_off=std(squeeze(nanmean(nanmean(Doff(channels_off(b),[5:45 55:95],rest_time(1):rest_time(2),conds_off),3),4)));
        std_mov_off=std(squeeze(nanmean(nanmean(Doff(channels_off(b),[5:45 55:95],mov_time(1):mov_time(2),conds_off),3),4)));
        
        pow_all_off(b,:)=(squeeze(nanmean(nanmean(Doff(channels_off(b),:,:,conds_off),3),4))./std_all_off).*100;
        pow_rest_off(b,:)=(squeeze(nanmean(nanmean(Doff(channels_off(b),:,rest_time(1):rest_time(2),conds_off),3),4))./std_rest_off).*100;
        pow_mov_off(b,:)=(squeeze(nanmean(nanmean(Doff(channels_off(b),:,mov_time(1):mov_time(2),conds_off),3),4))./std_mov_off).*100;
        
        beta.off(b,:)=(squeeze(nanmean(nanmean(Doff(channels_off(b),13:30,:,conds_off),2),4))./std_all_off).*100;
        beta.high.off(b,:)=(squeeze(nanmean(nanmean(Doff(channels_off(b),20:30,:,conds_off),2),4))./std_all_off).*100;
        beta.low.off(b,:)=(squeeze(nanmean(nanmean(Doff(channels_off(b),13:20,:,conds_off),2),4))./std_all_off).*100;
    end
    
    rnpower.on(n,:)=log(nanmean(pow_all_on));
    rnpower.off(n,:)=log(nanmean(pow_all_off));
    
    rnpower.rest.on(n,:)=log(nanmean(pow_rest_on));
    rnpower.rest.off(n,:)=log(nanmean(pow_rest_off));    
    
    rnpower.mov.on(n,:)=log(nanmean(pow_mov_on));
    rnpower.mov.off(n,:)=log(nanmean(pow_mov_off));
    
    rnpower.beta.on(n,:)=nanmean(beta.on);
    rnpower.beta.low.on(n,:)=nanmean(beta.low.on);
    rnpower.beta.high.on(n,:)=nanmean(beta.high.on);
    
    rnpower.beta.off(n,:)=nanmean(beta.off);
    rnpower.beta.low.off(n,:)=nanmean(beta.low.off);
    rnpower.beta.high.off(n,:)=nanmean(beta.high.off); 
    
   end
end

rnpower.time=Doff.time;

rnpower.mbeta.rest.on=nanmean(rnpower.beta.on(:,21:41),2);
rnpower.mbeta.rest.off=nanmean(rnpower.beta.off(:,21:41),2);

rnpower.mbeta.mov.on=nanmean(rnpower.beta.on(:,51:71),2);
rnpower.mbeta.mov.off=nanmean(rnpower.beta.off(:,51:71),2);

rnpower.mbeta.low.rest.on=nanmean(rnpower.beta.low.on(:,21:41),2);
rnpower.mbeta.low.rest.off=nanmean(rnpower.beta.low.off(:,21:41),2);

rnpower.mbeta.low.mov.on=nanmean(rnpower.beta.low.on(:,51:71),2);
rnpower.mbeta.low.mov.off=nanmean(rnpower.beta.low.off(:,51:71),2);

rnpower.mbeta.high.rest.on=nanmean(rnpower.beta.high.on(:,21:41),2);
rnpower.mbeta.high.rest.off=nanmean(rnpower.beta.high.off(:,21:41),2);

rnpower.mbeta.high.mov.on=nanmean(rnpower.beta.high.on(:,51:71),2);
rnpower.mbeta.high.mov.off=nanmean(rnpower.beta.high.off(:,51:71),2);

rnpower.rbeta.on=((rnpower.mbeta.mov.on-rnpower.mbeta.rest.on)./rnpower.mbeta.rest.on)*100;
rnpower.rbeta.off=((rnpower.mbeta.mov.off-rnpower.mbeta.rest.off)./rnpower.mbeta.rest.off)*100;

rnpower.rbeta.low.on=((rnpower.mbeta.low.mov.on-rnpower.mbeta.low.rest.on)./rnpower.mbeta.low.rest.on)*100;
rnpower.rbeta.low.off=((rnpower.mbeta.low.mov.off-rnpower.mbeta.low.rest.off)./rnpower.mbeta.low.rest.off)*100;

rnpower.rbeta.high.on=((rnpower.mbeta.high.mov.on-rnpower.mbeta.high.rest.on)./rnpower.mbeta.high.rest.on)*100;
rnpower.rbeta.high.off=((rnpower.mbeta.high.mov.off-rnpower.mbeta.high.rest.off)./rnpower.mbeta.high.rest.off)*100;

keep rnpower 
save rnpower

prest_onoff=wjn_ppt(rnpower.mbeta.rest.on,rnpower.mbeta.rest.off);
pmov_onoff=wjn_ppt(rnpower.mbeta.mov.on,rnpower.mbeta.mov.off);
prel_onoff=wjn_ppt(rnpower.rbeta.on,rnpower.rbeta.off);

prest=num2str(prest_onoff);
pmov=num2str(pmov_onoff);
prel=num2str(prel_onoff);

plowrest_onoff=wjn_ppt(rnpower.mbeta.low.rest.on,rnpower.mbeta.low.rest.off);
plowmov_onoff=wjn_ppt(rnpower.mbeta.low.mov.on,rnpower.mbeta.low.mov.off);
plowrel_onoff=wjn_ppt(rnpower.rbeta.low.on,rnpower.rbeta.low.off);

plowrest=num2str(plowrest_onoff);
plowmov=num2str(plowmov_onoff);
plowrel=num2str(plowrel_onoff);

phighrest_onoff=wjn_ppt(rnpower.mbeta.high.rest.on,rnpower.mbeta.high.rest.off);
phighmov_onoff=wjn_ppt(rnpower.mbeta.high.mov.on,rnpower.mbeta.high.mov.off);
phighrel_onoff=wjn_ppt(rnpower.rbeta.high.on,rnpower.rbeta.high.off);

phighrest=num2str(phighrest_onoff);
phighmov=num2str(phighmov_onoff);
phighrel=num2str(phighrel_onoff);

%% permutation

for c=1:size(rnpower.mov.on,2)
    rnpower.mov.pvalues(c)=wjn_ppt(rnpower.mov.on(:,c),rnpower.mov.off(:,c));
end

%% FIGURES


cc=colorlover(6);
color=cc(4,:);

    figure,
for c=1:size(rnpower.beta.off,1) 
    plot(rnpower.beta.on(c,:))
    hold on
end

figure,
ON=mypower(rnpower.time,rnpower.beta.on)
hold on,
OFF=mypower(rnpower.time,rnpower.beta.off, color)
xlim([-3 4])
xlabel('Time [s]')
ylabel('Power [a.u.]')
legend([ON OFF],'ON','OFF','location','northeast')
title('Beta Dynamics ON OFF')
xlim([-2 2.5])
ylim([80 220])
figone(7)
myprint('rbeta_on_off_std')

figure,
ON=mypower(rnpower.time,rnpower.beta.low.on)
hold on,
OFF=mypower(rnpower.time,rnpower.beta.off, color)
xlim([-2 3])
xlabel('Time [s]')
ylabel('Power [a.u.]')
legend([ON OFF],'ON','OFF','location','southeast')
title('Low Beta (13-20 Hz) ON OFF')
figone(7)
myprint('lowrbeta_on_off_std')

figure,
ON=mypower(rnpower.time,rnpower.beta.high.on)
hold on,
OFF=mypower(rnpower.time,rnpower.beta.high.off, color)
xlim([-2 3])
xlabel('Time [s]')
ylabel('Power [a.u.]')
legend([ON OFF],'ON','OFF','location','southeast')
title('High Beta (20-30 Hz) ON OFF')
figone(7)
myprint('highrbeta_on_off_std')

figure,
ON=mypower(1:400,rnpower.on)
hold on,
OFF=mypower(1:400,rnpower.off, color)
xlim([3 30])
xlabel('Frequency [Hz]')
ylabel('Power [a.u.]')
legend([ON OFF],'ON','OFF','location','northeast')
title('Power ON OFF over trl')
figone(7)
myprint('rpower_onoff_trl')

grey=[0.8 0.8 0.8];
rose=[170/255 17/255 17/255];

figure,
sigbar(1:400,rnpower.mov.pvalues<0.01)
ONmov=mypower(1:400,rnpower.mov.on)
hold on
OFFmov=mypower(1:400,rnpower.mov.off, color)
sigbar(1:400,rnpower.mov.pvalues<0.05)
xlim([8 30])
ylim([3.5 6.5])
xlabel('Frequency [Hz]')
ylabel('Power [a.u.]')
legend([ONmov OFFmov],'ON mov','OFF mov','location','northeast')
title('Power ON OFF mov p<.05')
figone(7)
myprint('rpower_onoff_mov_pval')

figure,
ON=mypower(1:400,rnpower.mov.on)
hold on,
OFF=mypower(1:400,rnpower.mov.off, color)
xlim([3 30])
xlabel('Frequency [Hz]')
ylabel('Power [a.u.]')
legend([ON OFF],'ON','OFF','location','northeast')
title('Power ON OFF mov')
figone(7)
myprint('rpower_onoff_mov')


cc = [0 0 0; 0.8078 0.0941 0.2118; 0 0 0; 0 0 0; 0.8078 0.0941 0.2118; 0 0 0; 0 0 0; 0.8078 0.0941 0.2118];

figure,
mybar([rnpower.mbeta.rest.on rnpower.mbeta.rest.off nan(size(rnpower.mbeta.mov.on,1),1) rnpower.mbeta.mov.on rnpower.mbeta.mov.off nan(size(rnpower.mbeta.mov.on,1),1) rnpower.rbeta.on.*100 rnpower.rbeta.off.*100],cc)
sigbracket (['P=' prest(1:4)],[1 2],200,9);
sigbracket(['P=' pmov(1:4)],[4 5],200,9)
sigbracket(['P=' prel(1:4)],[7 8],200,9);
set(gca,'XTick',[1 2 4 5 7 8],'XTickLabel',{'ON','OFF','ON','OFF','ON','OFF'},'XTickLabelRotation',45);
xlabel('Rest               Mov              Rel')
ylabel('Power [a.u.]')
title('Beta Power dynamics')
myprint('all_beta_bars')

figure,
mybar([rnpower.mbeta.low.rest.on rnpower.mbeta.low.rest.off nan(size(rnpower.mbeta.mov.on,1),1) rnpower.mbeta.low.mov.on rnpower.mbeta.low.mov.off nan(size(rnpower.mbeta.mov.on,1),1) rnpower.rbeta.low.on.*100 rnpower.rbeta.low.off.*100],cc)
sigbracket (['P=' plowrest(1:4)],[1 2],250,9);
sigbracket(['P=' plowmov(1:4)],[4 5],250,9)
sigbracket(['P=' plowrel(1:4)],[7 8],250,9);
set(gca,'XTick',[1 2 4 5 7 8],'XTickLabel',{'ON','OFF','ON','OFF','ON','OFF'},'XTickLabelRotation',45);
xlabel('Rest               Mov              Rel')
ylabel('Power [a.u.]')
ylim([0 310])
title('Low Beta (13-20) -/+.5 MO -2/-1')
myprint('low_beta_aroundMO_bars')

figure,
mybar([rnpower.mbeta.high.rest.on rnpower.mbeta.high.rest.off nan(size(rnpower.mbeta.mov.on,1),1) rnpower.mbeta.high.mov.on rnpower.mbeta.high.mov.off nan(size(rnpower.mbeta.mov.on,1),1) rnpower.rbeta.high.on.*100 rnpower.rbeta.high.off.*100],cc)
sigbracket (['P=' phighrest(1:4)],[1 2],250,9);
sigbracket(['P=' phighmov(1:4)],[4 5],250,9)
sigbracket(['P=' phighrel(1:4)],[7 8],250,9);
set(gca,'XTick',[1 2 4 5 7 8],'XTickLabel',{'ON','OFF','ON','OFF','ON','OFF'},'XTickLabelRotation',45);
xlabel('Rest               Mov              Rel')
ylabel('Power [a.u.]')
ylim([0 300])
title('high Beta (13-20) dynamics')
myprint('high_beta_bars')

