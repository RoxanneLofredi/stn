clear all
close all

root=fullfile(mdf,'stn_rotameter\STN_Rotameter\');
cd(root);

addpath(fullfile(mdf,'matlab_scripts\wjn_toolbox'));
addpath(fullfile(mdf,'stn_rotameter\rotameter_scripts'));

pat=patienten;

fhz=[40 90];
conds='onset';
type='brown';
outname='gamma_brown';

gfiles=ffind('gamma_brown_burst_r*ON*.mat');

%% Calculate velocity and burst properties

for a=1:length(gfiles)
    
    g=load(gfiles{a});
    filename=gfiles{a}(21:end);
    
    %load velocity values for each cond
    D=spm_eeg_load(['re' filename]);
    ti=ci('onset',D.conditions);
    
    ig = ci('onset-gross',D.conditions(ti));
    im = ci('onset-mittel',D.conditions(ti));
    ik = ci('onset-klein',D.conditions(ti));
    
    D=spm_eeg_load(['re' filename(2:end)]);
    vel = D.btrialvel(ti);
    
    %write velocity values (sep for conditions) in gamma burst file
    g.trialvel=vel;
    g.meanvel = nanmean(vel);

    g.gvel = nanmean(vel(ig));
    g.mvel = nanmean(vel(im));
    g.kvel = nanmean(vel(ik));
    
    %calculate bl-corr mean rate/amp/dur increases during mov across trls
    g.meanblratemove = nanmean(g.bratemov,2);
    g.meanblampmove = nanmean(g.bampmov,2);
    g.meanbldurmove = nanmean(g.bdurmov,2);
    
    %sep bl-corr mean rate/amp/dur increases for conditions across trials
    g.gblratemov = nanmean(g.bratemov(:,ig),2);
    g.mblratemov = nanmean(g.bratemov(:,im),2);
    g.sblratemov = nanmean(g.bratemov(:,ik),2);
    
    g.gblampmov = nanmean(g.bampmov(:,ig),2);
    g.mblampmov = nanmean(g.bampmov(:,im),2);
    g.sblampmov = nanmean(g.bampmov(:,ik),2);
    
    g.gbldurmov = nanmean(g.bdurmov(:,ig),2);
    g.mbldurmov = nanmean(g.bdurmov(:,im),2);
    g.sbldurmov = nanmean(g.bdurmov(:,ik),2);
        
    % calculate mean rate/amp/dur without bl corr during mov
    t_move = wjn_sc(g.t,0):wjn_sc(g.t,.5);
    
    for c = 1:size(g.burstraster,1)
        g.ratemov(c,:)=nansum(g.burstraster(c,:,t_move),3);
        g.ampmov(c,:)=nanmean(g.burstamplitude(c,:,t_move),3);
        g.durmov(c,:)=nanmean(g.burstduration(c,:,t_move),3);
    end
    
    % non-bl corr rate/amp/dur values sep for conds
    g.gratemov = g.ratemov(:,ig);
    g.mratemov = g.ratemov(:,im);
    g.sratemov = g.ratemov(:,ik);
    
    g.gampmov = g.ampmov(:,ig);
    g.mampmov = g.ampmov(:,im);
    g.sampmov = g.ampmov(:,ik);
    
    g.gdurmov = g.durmov(:,ig);
    g.mdurmov = g.durmov(:,im);
    g.sdurmov = g.durmov(:,ik);

    % mean non-bl corr rate/amp/dur sep for conds
    
    g.meanratemov=nanmean(g.ratemov,2);
    g.meanampmov=nanmean(g.ampmov,2);
    g.meandurmov=nanmean(g.durmov,2);
    
    g.gmeanratemov=nanmean(g.gratemov,2);
    g.gmeanampmov=nanmean(g.gampmov,2);
    g.gmeandurmov=nanmean(g.gdurmov,2);
    
    g.mmeanratemov=nanmean(g.mratemov,2);
    g.mmeanampmov=nanmean(g.mampmov,2);
    g.mmeandurmov=nanmean(g.mdurmov,2);
    
    g.smeanratemov=nanmean(g.sratemov,2);
    g.smeanampmov=nanmean(g.sampmov,2);
    g.smeandurmov=nanmean(g.sdurmov,2);
    
    % save new gamma burst file
    save([outname '_burst_re' filename],'-struct','g');
end

%% Calulate correlations

clear all
gfiles=ffind('gamma_brown_burst_r*ON*.mat');
outname='gamma_brown';

%% corr power and burst properties

for a=1:length(gfiles)
    
    filename=gfiles{a}(21:end);
    g=load(gfiles{a});
    chs=g.channels;
    
    for b=1:length(chs)
      g.corrpowerbratemov(b)=wjn_pc(g.rawampmov(b,:)',g.bratemov(b,:)');
      g.corrpowerbampmov(b)=wjn_pc(g.rawampmov(b,:)',g.bampmov(b,:)');
      g.corrpowerdurmov(b)=wjn_pc(g.rawampmov(b,:)',g.bdurmov(b,:)');
    end 
    
     save([outname '_burst_re' filename],'-struct','g');
end

% corr across patients

for a=1:length(gfiles)
    
    g=load(gfiles{a},'corrpowerbratemov','corrpowerbampmov','corrpowerdurmov','npt');
    fnames = fieldnames(g);
    fnames=fnames(1:3);
   
    np(a)=g.npt;

    for c = 1:length(fnames);
            dg.(fnames{c})(a) = nanmean(g.(fnames{c}),2);
    end
end

for a = 1:length(fnames)
    mg.(fnames{a}) = wjn_match_cases(dg.(fnames{a}),np');
end

%stats across patients

pamprate=wjn_ppt(mg.corrpowerbratemov,mg.ccorrpowerbampmov);
pdurrate=wjn_ppt(mg.corrpowerdurmov,mg.corrpowerbratemov);
pduramp=wjn_ppt(mg.corrpowerdurmov,mg.corrpowerbampmov);

samprate=signrank(mg.corrpowerbratemov,mg.ccorrpowerbampmov);
sdurrate=signrank(mg.corrpowerdurmov,mg.corrpowerbratemov);
sduramp=signrank(mg.corrpowerdurmov,mg.corrpowerbampmov);

%figure 

cc = colorlover(1);

figure
mybar([mg.corrbratemov,mg.corrbampmov,mg.corrbdurmov],cc([1 3 5],:),[1 2 3])
sigbracket('***',[1 2],0.75)
sigbracket('***',[1 3],0.85)
sigbracket('***',[2 3],0.75)
ylim([0 0.95])
set(gca,'XTick',[1 2 3],'XTickLabel',{'Rate','Amp','Dur'}) %,'XTickLabelRotation',45)
ylabel('Spearmans \rho')

myprint('burst_fig6_corr_power_burst')

%% corr trlvel and rate/amp/dur during mov (non bl corr)

for a=2:length(gfiles)
    
    filename=gfiles{a}(21:end);
    g=load(gfiles{a});
    chs=g.channels;
    
    for b=1:length(chs)
      g.corrvelratemov(b)=wjn_pc(g.trialvel',g.ratemov(b,:)');
      g.corrvelampmov(b)=wjn_pc(g.trialvel',g.ampmov(b,:)');
      g.corrveldurmov(b)=wjn_pc(g.trialvel',g.durmov(b,:)');
    end 
    
     save([outname '_burst_re' filename],'-struct','g');
end


for a=1:length(gfiles)
    
    g=load(gfiles{a},'corrvelratemov','corrvelampmov','corrveldurmov','npt');
    fnames = fieldnames(g);
    fnames=fnames(1:3);
   
    np(a)=g.npt;

    for c = 1:length(fnames);
            dg.(fnames{c})(a) = nanmean(g.(fnames{c}),2);
    end
end

for a = 1:length(fnames)
    mg.(fnames{a}) = wjn_match_cases(dg.(fnames{a}),np');
end

pvelrateamp=wjn_ppt(mg.corrvelratemov,mg.corrvelampmov);
pvelratedur=wjn_ppt(mg.corrvelratemov,mg.corrveldurmov);
pvelduramp=wjn_ppt(mg.corrveldurmov,mg.corrvelampmov);

svelrateamp=signrank(mg.corrvelratemov,mg.corrvelampmov);
svelratedur=signrank(mg.corrvelratemov,mg.corrveldurmov);
svelduramp=signrank(mg.corrveldurmov,mg.corrvelampmov);

figure,
notBoxPlot([mg.corrvelratemov mg.corrvelampmov mg.corrveldurmov])
title('Corr coeff Trial vel ~ burst prop')
set(gca,'XTick',[1 2 3],'XTickLabel',{'Rate','Amp','Dur'}) 
myprint('velcorrfig_trlvel_noblcorr_burstprop');

%% corr trlvel and bl corr rate/amp/dur during mov

for a=1:length(gfiles)
    
    filename=gfiles{a}(21:end);
    g=load(gfiles{a});
    chs=g.channels;
    
    for b=1:length(chs)
      g.corrvelblratemov(b)=wjn_pc(g.trialvel',g.bratemov(b,:)');
      g.corrvelblampmov(b)=wjn_pc(g.trialvel',g.bampmov(b,:)');
      g.corrvelbldurmov(b)=wjn_pc(g.trialvel',g.bdurmov(b,:)');
    end 
    
     save([outname '_burst_re' filename],'-struct','g');
end


for a=1:length(gfiles)
    
    g=load(gfiles{a},'corrvelblratemov','corrvelblampmov','corrvelbldurmov','npt');
    fnames = fieldnames(g);
    fnames=fnames(1:3);
   
    np(a)=g.npt;

    for c = 1:length(fnames);
            dg.(fnames{c})(a) = nanmean(g.(fnames{c}),2);
    end
end

for a = 1:length(fnames)
    mg.(fnames{a}) = wjn_match_cases(dg.(fnames{a}),np');
end

pvelblrateamp=wjn_ppt(mg.corrvelblratemov,mg.corrvelblampmov);
pvelblratedur=wjn_ppt(mg.corrvelblratemov,mg.corrvelbldurmov);
pvelblduramp=wjn_ppt(mg.corrvelbldurmov,mg.corrvelblampmov);

svelrateamp=signrank(mg.corrvelblratemov,mg.corrvelblampmov);
svelratedur=signrank(mg.corrvelblratemov,mg.corrvelbldurmov);
svelduramp=signrank(mg.corrvelbldurmov,mg.corrvelblampmov);

figure,
notBoxPlot([mg.corrvelblratemov mg.corrvelblampmov mg.corrvelbldurmov])
title('Corr coeff Trial vel ~ bl corr burst prop')
set(gca,'XTick',[1 2 3],'XTickLabel',{'Rate','Amp','Dur'}) 
myprint('velcorrfig_trlvel_blcorr_burstprop');

%% mean bl-corr burst prop values per cond corr with mean vel per cond

for a=1:length(gfiles)
    
    filename=gfiles{a}(21:end);
    g=load(gfiles{a});
    chs=g.channels;
    
    for b=1:length(chs)
      g.corrvelblratemovcond(b)=wjn_corr_scatter([g.gvel g.mvel g.kvel]',[g.gblratemov(b,:) g.mblratemov(b,:) g.sblratemov(b,:)]');
      g.corrvelblampmovcond(b)=wjn_pc([g.gvel g.mvel g.kvel]',[g.gblampmov(b,:) g.mblampmov(b,:) g.sblampmov(b,:)]');
      g.corrvelbldurmovcond(b)=wjn_pc([g.gvel g.mvel g.kvel]',[g.gbldurmov(b,:) g.mbldurmov(b,:) g.sbldurmov(b,:)]');
    end 
    
     save([outname '_burst_re' filename],'-struct','g');
end


for a=1:length(gfiles)
    
    g=load(gfiles{a},'corrvelblratemovcond','corrvelblampmovcond','corrvelbldurmovcond','npt');
    fnames = fieldnames(g);
    fnames=fnames(1:3);
   
    np(a)=g.npt;

    for c = 1:length(fnames);
            dg.(fnames{c})(a) = nanmean(g.(fnames{c}),2);
    end
end

for a = 1:length(fnames)
    mg.(fnames{a}) = wjn_match_cases(dg.(fnames{a}),np');
end

pvelblrateampcond=wjn_ppt(mg.corrvelblratemovcond,mg.corrvelblampmovcond);
pvelblratedurcond=wjn_ppt(mg.corrvelblratemovcond,mg.corrvelbldurmovcond);
pvelbldurampcond=wjn_ppt(mg.corrvelbldurmovcond,mg.corrvelblampmovcond);

figure,
notBoxPlot([mg.corrvelblratemovcond mg.corrvelblampmovcond mg.corrvelbldurmovcond])
title('Corr coeff Trial vel conds ~ bl corr burst prop conds')
set(gca,'XTick',[1 2 3],'XTickLabel',{'Rate','Amp','Dur'}) 
myprint('velcorrfig_trlvel_blcorr_burstprop_conds');

%% mean burst prop values per cond corr with mean vel per cond (not bl corr)

for a=1:length(gfiles)
    
    filename=gfiles{a}(21:end);
    g=load(gfiles{a});
    chs=g.channels;
    
    for b=1:length(chs)
      g.corrvelratemovcond(b)=wjn_corr_scatter([g.gvel g.mvel g.kvel]',[g.gmeanratemov(b,:) g.mmeanratemov(b,:) g.smeanratemov(b,:)]');
      g.corrvelampmovcond(b)=wjn_pc([g.gvel g.mvel g.kvel]',[g.gmeanampmov(b,:) g.mmeanampmov(b,:) g.smeanampmov(b,:)]');
      g.corrveldurmovcond(b)=wjn_pc([g.gvel g.mvel g.kvel]',[g.gmeandurmov(b,:) g.mmeandurmov(b,:) g.smeandurmov(b,:)]');
    end 
    
     save([outname '_burst_re' filename],'-struct','g');
end


for a=1:length(gfiles)
    
    g=load(gfiles{a},'corrvelratemovcond','corrvelampmovcond','corrveldurmovcond','npt');
    fnames = fieldnames(g);
    fnames=fnames(1:3);
   
    np(a)=g.npt;

    for c = 1:length(fnames);
            dg.(fnames{c})(a) = nanmean(g.(fnames{c}),2);
    end
end

for a = 1:length(fnames)
    mg.(fnames{a}) = wjn_match_cases(dg.(fnames{a}),np');
end

pvelrateampcond=wjn_ppt(mg.corrvelratemovcond,mg.corrvelampmovcond);
pvelratedurcond=wjn_ppt(mg.corrvelratemovcond,mg.corrveldurmovcond);
pveldurampcond=wjn_ppt(mg.corrveldurmovcond,mg.corrvelampmovcond);

figure,
notBoxPlot([mg.corrvelratemovcond mg.corrvelampmovcond mg.corrveldurmovcond])
title('Corr coeff Trial vel conds ~ not bl corr burst prop conds')
set(gca,'XTick',[1 2 3],'XTickLabel',{'Rate','Amp','Dur'}) 
myprint('velcorrfig_trlvel_notblcorr_burstprop_conds');

%% ohne mittlere cond

for a=1:length(gfiles)
    
    filename=gfiles{a}(21:end);
    g=load(gfiles{a});
    chs=g.channels;
    
    for b=1:length(chs)
      g.corrvelratemovcond2(b)=wjn_corr_scatter([g.gvel g.kvel]',[g.gmeanratemov(b,:) g.smeanratemov(b,:)]');
      g.corrvelampmovcond2(b)=wjn_pc([g.gvel g.kvel]',[g.gmeanampmov(b,:) g.smeanampmov(b,:)]');
      g.corrveldurmovcond2(b)=wjn_pc([g.gvel g.kvel]',[g.gmeandurmov(b,:) g.smeandurmov(b,:)]');
    end 
    
     save([outname '_burst_re' filename],'-struct','g');
end


for a=1:length(gfiles)
    
    g=load(gfiles{a},'corrvelratemovcond2','corrvelampmovcond2','corrveldurmovcond2','npt');
    fnames = fieldnames(g);
    fnames=fnames(1:3);
   
    np(a)=g.npt;

    for c = 1:length(fnames);
            dg.(fnames{c})(a) = nanmean(g.(fnames{c}),2);
    end
end

for a = 1:length(fnames)
    mg.(fnames{a}) = wjn_match_cases(dg.(fnames{a}),np');
end

pvelrateampcond2=wjn_ppt(mg.corrvelratemovcond2,mg.corrvelampmovcond2);
pvelratedurcond2=wjn_ppt(mg.corrvelratemovcond2,mg.corrveldurmovcond2);
pveldurampcond2=wjn_ppt(mg.corrveldurmovcond2,mg.corrvelampmovcond2);

figure,
notBoxPlot([mg.corrvelratemovcond2 mg.corrvelampmovcond2 mg.corrveldurmovcond2])
title('Corr coeff Trial vel conds ~ not bl corr burst prop 2 conds')
set(gca,'XTick',[1 2 3],'XTickLabel',{'Rate','Amp','Dur'}) 
myprint('velcorrfig_trlvel_notblcorr_burstprop_2conds');

%% corr vel burst prop across patients

clear all
gfiles=ffind('gamma_brown_burst_r*ON*.mat');
outname='gamma_brown';

for a=1:length(gfiles)
    
    g=load(gfiles{a},'gvel','gblratemov','kvel','sblratemov','npt');
    fnames = fieldnames(g);
    fnames=fnames(1:4);
   
    np(a)=g.npt;

    for c = 1:length(fnames);
            dg.(fnames{c})(a) = nanmean(g.(fnames{c}));
    end
end

sub_increase_vel=(mg.gvel./mg.kvel);
sub_increase_rate=(mg.gblratemov - mg.sblratemov);


figure, wjn_corr_scatter(sub_increase_vel, sub_increase_rate)

for a = 1:length(fnames)
    mg.(fnames{a}) = wjn_match_cases(dg.(fnames{a}),np');
end

figure, wjn_corr_scatter(mg.gvel, mg.gblratemove)

pvelrateampcond2=wjn_ppt(mg.corrvelratemovcond2,mg.corrvelampmovcond2);
pvelratedurcond2=wjn_ppt(mg.corrvelratemovcond2,mg.corrveldurmovcond2);
pveldurampcond2=wjn_ppt(mg.corrveldurmovcond2,mg.corrvelampmovcond2);

figure,
notBoxPlot([mg.corrvelratemovcond2 mg.corrvelampmovcond2 mg.corrveldurmovcond2])
title('Corr coeff Trial vel conds ~ not bl corr burst prop 2 conds')
set(gca,'XTick',[1 2 3],'XTickLabel',{'Rate','Amp','Dur'}) 
myprint('velcorrfig_trlvel_notblcorr_burstprop_2conds');

%% best contact pairs correlation

clear all
gfiles=ffind('gamma_brown_burst_r*ON*.mat');
outname='gamma_brown';

for a=1:length(gfiles)
    
    g=load(gfiles{a},'gvel','kvel','gblratemov','sblratemov','npt');
    fnames = fieldnames(g);
    fnames=fnames(1:4);
   
    np(a)=g.npt;
  
    [maxv, maxi]= max (g.gblratemov);
    
    for c = 1:length(fnames);
        if c>2
            dg.(fnames{c})(a) = g.(fnames{c})(maxi,1);
        else
            dg.(fnames{c})(a) = nanmean(g.(fnames{c}));
        end
    end
end

for c=1:length(fnames);
dg.(fnames{c})(23)= [];
end

sub_increase_vel=dg.gvel./dg.kvel;
sub_increase_rate=dg.gblratemov./dg.sblratemov;

figure, wjn_corr_scatter(sub_increase_vel', sub_increase_rate')

for a = 1:length(fnames)
    mg.(fnames{a}) = wjn_match_cases(dg.(fnames{a}),np');
end

