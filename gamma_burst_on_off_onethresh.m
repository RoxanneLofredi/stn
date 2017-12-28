clear all
close all

root=fullfile(mdf,'stn_rotameter\STN_Rotameter\');
cd(root);

addpath(fullfile(mdf,'matlab_scripts\wjn_toolbox'));
addpath(fullfile(mdf,'stn_rotameter\rotameter_scripts'));

pat=patienten;

%%
fhz=[40 90];
conds='onset';
type='brown';
outname='newoneth_gamma_brown';
onethresh=1;

n=0;

for p=17:length(pat)
    
    if ~isfield(pat{p},'bad') && pat{p}.med;
       
        filename=(['refspmeeg_' pat{p}.id '.mat']);
        D=wjn_sl(filename); 
        
        if strcmp(pat{p}.side,'R');
            ch=channel_finder('STNL',D.chanlabels);
        else
            ch=channel_finder('STNR',D.chanlabels);
        end
        
        if isfield(pat{p},'notstnc')
            chs = ch(~ismember(ch,pat{p}.notstnc));   
        else
            chs=ch;
        end
        
        d=wjn_eeg_burst_analysis(filename,fhz,[],chs,conds,type,onethresh,[],outname);
        
        if d.noburst
            n=n+1;
            noburst(n,:)={filename, chs(d.noburst)};
        end
        
        try
            updrs = pat{p}.updrs(abs(pat{p}.med-2));
        catch
            updrs = nan;
        end
        
        d.npt = pat{p}.npatient;
        d.med = pat{p}.med;
        d.updrs = updrs;
        
        t_move = wjn_sc(d.t,0):wjn_sc(d.t,.5);
        t_base= wjn_sc(d.t,-2):wjn_sc(d.t,-1.5);
        
        for a=1:length(chs)
            for b=1:d.ntrials
                
        d.burstsummov(a,b)=nansum(squeeze(d.burstraster(a,b,t_move)));  
        d.mampmov(a,b)=nanmean(d.burstamplitude(a,b,t_move));
        d.mdurmov(a,b)=nanmean(d.burstduration(a,b,t_move));
                
        d.rawampmov(a,b)=(nanmean(d.rawamp(a,b,t_move),3)-nanmean(d.rawamp(a,b,t_base),3))./nanmean(d.rawamp(a,b,t_base),3);
        d.bratemov(a,b)=(nansum(d.burstraster(a,b,t_move),3)-nansum(d.burstraster(a,b,t_base),3))./nansum(d.burstraster(a,b,t_base),3);
        d.bampmov(a,b)=(nanmean(d.burstamplitude(a,b,t_move),3)-nanmean(d.burstamplitude(a,b,t_base),3))./nanmean(d.burstamplitude(a,b,t_base),3);
        d.bdurmov(a,b)=(nanmean(d.burstduration(a,b,t_move),3)-nanmean(d.burstduration(a,b,t_base),3))./nanmean(d.burstduration(a,b,t_base),3);
            end
        infrate=find(isinf(d.bratemov(a,:)));
        infamp=find(isinf(d.bampmov(a,:)));
        infdur=find(isinf(d.bdurmov(a,:)));
        d.bratemov(a,infrate)=nan;
        d.bampmov(a,infamp)=nan;
        d.bdurmov(a,infdur)=nan;
        
%         d.corrbratemov(a)=wjn_permcorr(d.rawampmov(a,:),d.bratemov(a,:));
%         d.corrbampmov(a)=wjn_permcorr(d.rawampmov(a,:),d.bampmov(a,:));
%         d.corrbdurmov(a)=wjn_permcorr(d.rawampmov(a,:),d.bdurmov(a,:));
        end
       save([outname '_burst_' filename],'-struct','d')
    end
end
%%
fhz=[40 90];
conds='onset';
outname='newoneth_gamma_brown';
type='brown';
onethresh=1;
n=0;

for p=1:length(pat)
    
    if ~isfield(pat{p},'bad') && ~pat{p}.med && isfield(pat{p},'offseven')
           
        filename=(['refspmeeg_' pat{p}.id '.mat']);
        D=wjn_sl(filename); 
        
        onfile=filename(1:20);
        
        try
        g=load(['onethresh_gamma_brown_burst_' onfile '_ON.mat'],'thresh');
        catch
        g=load(['onethresh_gamma_brown_burst_' onfile '_ON_merge.mat'],'thresh');
        end
        
        if strcmp(pat{p}.side,'R');
            ch=channel_finder('STNL',D.chanlabels);
        else
            ch=channel_finder('STNR',D.chanlabels);
        end
        
        if isfield(pat{p},'notstnc')
            chs = ch(~ismember(ch,pat{p}.notstnc));   
        else
            chs=ch;
        end
        
        type={g.thresh(:,1),'brown'};
        onethresh=1;
        
        d=wjn_eeg_burst_analysis(filename,fhz,[],chs,conds,type,onethresh,outname);
        
         if d.noburst
            n=n+1;
            noburst(n,:)={filename, chs(d.noburst)};
         end
        
        try
            updrs = pat{p}.updrs(abs(pat{p}.med-2));
        catch
            updrs = nan;
        end
        
        d.npt = pat{p}.npatient;
        d.med = pat{p}.med;
        d.updrs = updrs;
        
        t_move = wjn_sc(d.t,0):wjn_sc(d.t,.5);
        t_base= wjn_sc(d.t,-2):wjn_sc(d.t,-1.5);
        
        for a=1:length(chs)
            for b=1:d.ntrials
        d.burstsummov(a,b)=nansum(squeeze(d.burstraster(a,b,t_move)));  
        d.mampmov(a,b)=nanmean(d.burstamplitude(a,b,t_move));
        d.mdurmov(a,b)=nanmean(d.burstduration(a,b,t_move));
       
        d.rawampmov(a,b)=(nanmean(d.rawamp(a,b,t_move),3)-nanmean(d.rawamp(a,b,t_base),3))./nanmean(d.rawamp(a,b,t_base),3);
        d.bratemov(a,b)=(nansum(d.burstraster(a,b,t_move),3)-nansum(d.burstraster(a,b,t_base),3))./nansum(d.burstraster(a,b,t_base),3);
        d.bampmov(a,b)=(nanmean(d.burstamplitude(a,b,t_move),3)-nanmean(d.burstamplitude(a,b,t_base),3))./nanmean(d.burstamplitude(a,b,t_base),3);
        d.bdurmov(a,b)=(nanmean(d.burstduration(a,b,t_move),3)-nanmean(d.burstduration(a,b,t_base),3))./nanmean(d.burstduration(a,b,t_base),3);
            end
        infrate=find(isinf(d.bratemov(a,:)));
        infamp=find(isinf(d.bampmov(a,:)));
        infdur=find(isinf(d.bdurmov(a,:)));
        d.bratemov(a,infrate)=nan;
        d.bampmov(a,infamp)=nan;
        d.bdurmov(a,infdur)=nan;
        
%         d.corrbratemov(a)=wjn_corr_scatter(d.rawampmov(a,:),d.bratemov(a,:));
%         d.corrbampmov(a)=wjn_corr_scatter(d.rawampmov(a,:),d.bampmov(a,:));
%         d.corrbdurmov(a)=wjn_corr_scatter(d.rawampmov(a,:),d.bdurmov(a,:));
        end
        
        save([outname '_burst_' filename],'-struct','d')
   end
        

end