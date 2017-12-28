clear all
close all
% seven = {'r032DA55'    'r035EE53'    'r091NF42'    'r094RE60'    'r011AI40'    'r041EC34'  'r093AL56'};
% iids = ci(seven,ids);


root='';
cd(root);

pat=patienten;

prefix = 'sbgmtf_efspmeeg_';
behfix = 'rtf_efspmeeg_';

conds={'onset-klein','onset-mittel','onset-gross'};

n=0

for a=1:length(pat);
    
    id=strsplit(pat{a}.id, '_');
    
    if ~isfield(pat{a},'bad') && pat{a}.med && isfield(pat{a},'updrs'); %&& ~strcmp(pat{a}.id,'r004IR48_L_ON')
    
        n=n+1;
        
        D=spm_eeg_load([prefix pat{a}.id '.mat']);
        E=spm_eeg_load([behfix pat{a}.id '.mat']);
        rvs(n) = nanmean(E.btrialrvel(ci(conds{1},E.conditions)));
        rvl(n) = nanmean(E.btrialrvel(ci(conds{2},E.conditions)));
        rvm(n) = nanmean(E.btrialrvel(ci(conds{3},E.conditions)));
        
        rts(n) = nanmean(E.btrialrt(ci(conds{1},E.conditions)));
        rtl(n) = nanmean(E.btrialrt(ci(conds{2},E.conditions)));
        rtm(n) = nanmean(E.btrialrt(ci(conds{3},E.conditions)));
        
        rtall(n)=nanmean(E.btrialrt);
        
        if id{2}=='L'  
            cl=ci('STNR', D.chanlabels);
        else
            cl=ci('STNL', D.chanlabels);
        end
        
        med(n) = pat{a}.med;
     
        nid(n) = pat{a}.npatient;
        sv=ci(conds{1},D.conditions);
        gv=ci(conds{3},D.conditions);
        mv=ci(conds{2},D.conditions);
        
        
        if pat{a}.med
            updrs(n) = pat{a}.updrs(1);
        else
            updrs(n) = pat{a}.updrs(2);
        end

        
        gs(n,:) = squeeze(nanmean(nanmean(nanmean(nanmean(D(cl,40:90,sc(D.time,0):sc(D.time,0.5),sv),1),2),3),4));
        gm(n,:) = squeeze(nanmean(nanmean(nanmean(nanmean(D(cl,40:90,sc(D.time,0):sc(D.time,0.5),mv),1),2),3),4));
        gl(n,:) = squeeze(nanmean(nanmean(nanmean(nanmean(D(cl,40:90,sc(D.time,0):sc(D.time,0.5),gv),1),2),3),4));
   
        bs(n,:) = squeeze(nanmean(nanmean(nanmean(nanmean(D(cl,13:30,sc(D.time,0):sc(D.time,0.5),sv),1),2),3),4));
        bm(n,:) = squeeze(nanmean(nanmean(nanmean(nanmean(D(cl,13:30,sc(D.time,0):sc(D.time,0.5),mv),1),2),3),4));
        bl(n,:) = squeeze(nanmean(nanmean(nanmean(nanmean(D(cl,13:30,sc(D.time,0):sc(D.time,0.5),gv),1),2),3),4));
   
        ts(n,:) = squeeze(nanmean(nanmean(nanmean(nanmean(D(cl,2:8,sc(D.time,0):sc(D.time,0.5),sv),1),2),3),4));
        tm(n,:) = squeeze(nanmean(nanmean(nanmean(nanmean(D(cl,2:8,sc(D.time,0):sc(D.time,0.5),mv),1),2),3),4));
        tl(n,:) = squeeze(nanmean(nanmean(nanmean(nanmean(D(cl,2:8,sc(D.time,0):sc(D.time,0.5),gv),1),2),3),4));
         
    end
    
end

ngs=wjn_match_cases(gs,nid);
ngm=wjn_match_cases(gm,nid);
ngl=wjn_match_cases(gl,nid);
nbs=wjn_match_cases(bs,nid);
nbm=wjn_match_cases(bm,nid);
nbl=wjn_match_cases(bl,nid);
nts=wjn_match_cases(ts,nid);
ntm=wjn_match_cases(tm,nid);
ntl=wjn_match_cases(tl,nid);
nvs=wjn_match_cases(rvs,nid);
nvm=wjn_match_cases(rvm,nid);
nvl=wjn_match_cases(rvl,nid);
nupdrs = wjn_match_cases(updrs,nid);
nrts=wjn_match_cases(rts,nid);
nrtm=wjn_match_cases(rtm,nid);
nrtl=wjn_match_cases(rtl,nid);
nrtall=wjn_match_cases(rtall,nid);


%%
 pg = ((ngl-ngs)./ngl)*100 ;
% pg = ngl./ngs*100;
pv = ((nvl-nvs)./nvl)*100;
% pg(16) = [];
% pv(16)=[];
% ngl(16)=[];
% nvl(16)= [];
% nupdrs(16)=[];
% ng = ngl-ngs;

nmeang=nanmean([ngs; ngm; ngl]);
nmeant=nanmean([nts; ntm; ntl]);
nmeanb=nanmean([nbs; nbm; nbl]);
nmeanv=nanmean([nvs;nvm;nvl]);
nmeanrt=nanmean([nrts;nrtm;nrtl]);

figure,
wjn_corr_scatter(nmeang,nmeant);
xlabel('Relative power change [%] - Gamma')
ylabel('Relative power change [%] - Theta');

figure,
subplot(1,3,1)
wjn_corr_scatter(nupdrs,nmeant);
xlabel('UPDRS-ON')
ylabel('Relative power change [%] - Theta')
subplot(1,3,2)
wjn_corr_scatter(nupdrs,nmeanb);
xlabel('UPDRS-ON')
ylabel('Relative power change [%] - Beta')
subplot(1,3,3)
wjn_corr_scatter(nupdrs,nmeang);
xlabel('UPDRS-ON')
ylabel('Relative power change [%] - Gamma')
figone(9,35)
myprint('figures/correlations_tbg_ON_spearman_025maxV')


figure,
subplot(1,3,1)
wjn_corr_scatter(nmeanrt,nmeant);
xlabel('Reaction time [s]')
ylabel('Relative power change [%] - Theta')
subplot(1,3,2)
wjn_corr_scatter(nmeanrt,nmeanb);
xlabel('Reaction time [s]')
ylabel('Relative power change [%] - Beta')
subplot(1,3,3)
wjn_corr_scatter(nmeanrt,nmeang);
xlabel('Reaction time [s]')
ylabel('Relative power change [%] - Gamma')
figone(9,35)
myprint('figures/correlations_tbg_rt_ON_spearman_025maxV')

%%

figure,
subplot(1,3,1)
wjn_corr_scatter(nmeanv,nmeant);
subplot(1,3,2)
wjn_corr_scatter(nmeanv,nmeanb);
subplot(1,3,3)
wjn_corr_scatter(nmeanv,nmeang);
figone(9,35)

close all

figure
subplot(1,4,1)
wjn_corr_scatter(nupdrs',pg')
xlabel('UPDRS')
ylabel('% Gamma small vs. large')
subplot(1,4,2)
wjn_corr_scatter(ngl',pg')
xlabel('% Gamma large')
ylabel('% Gamma small vs. large')
subplot(1,4,3)
wjn_corr_scatter(nupdrs',nvl')
xlabel('UPDRS')
ylabel('% velocity large')
subplot(1,4,4)
wjn_corr_scatter(pv',pg')
xlabel('% velocity large')
ylabel('% Gamma small vs. large')

figure
subplot(1,3,1)
wjn_corr_scatter(nupdrs',ngs')
xlabel('UPDRS')
ylabel('% Gamma small')
subplot(1,3,2)
wjn_corr_scatter(nupdrs',ngm')
xlabel('UPDRS')
ylabel('% Gamma medium')
subplot(1,3,3)
wjn_corr_scatter(nupdrs',ngl')
xlabel('UPDRS')
ylabel('% Gamma large')

figure
subplot(1,3,1)
wjn_corr_scatter(nvs',ngs')
xlabel('velocity small')
ylabel('% Gamma small')
subplot(1,3,2)
wjn_corr_scatter(nvm',ngm')
xlabel('velocity med')
ylabel('% Gamma medium')
subplot(1,3,3)
wjn_corr_scatter(nvl',ngl')
xlabel('velocity large')
ylabel('% Gamma large')

figure
subplot(1,3,1)
wjn_corr_scatter(nvs',nts')
xlabel('velocity small')
ylabel('% Theta small')
subplot(1,3,2)
wjn_corr_scatter(nvm',ntm')
xlabel('velocity med')
ylabel('% Theta medium')
subplot(1,3,3)
wjn_corr_scatter(nvl',ntl')
xlabel('velocity large')
ylabel('% Theta large')


figure
subplot(1,3,1)
wjn_corr_scatter(nupdrs',ngl'-ngs')
xlabel('UPDRS')
ylabel('% Gamma small')
subplot(1,3,2)
wjn_corr_scatter(nvm',ngm')
xlabel('UPDRS')
ylabel('% Gamma medium')
subplot(1,3,3)
wjn_corr_scatter(nvl',ngl')
xlabel('UPDRS')
ylabel('% Gamma large')


figure
subplot(1,3,1)
wjn_corr_scatter(nupdrs',nbs')
xlabel('UPDRS')
ylabel('% Beta small')
subplot(1,3,2)
wjn_corr_scatter(nupdrs',nbm')
xlabel('UPDRS')
ylabel('% Beta medium')
subplot(1,3,3)
wjn_corr_scatter(nupdrs',nbl')
xlabel('UPDRS')
ylabel('% Beta large')


figure
subplot(1,3,1)
wjn_corr_scatter(nupdrs',nts')
xlabel('UPDRS')
ylabel('% theta small')
subplot(1,3,2)
wjn_corr_scatter(nupdrs',ntm')
xlabel('UPDRS')
ylabel('% theta medium')
subplot(1,3,3)
wjn_corr_scatter(nupdrs',ntl')
xlabel('UPDRS')
ylabel('% theta large')

figure
subplot(1,3,1)
wjn_corr_scatter(ngs',nts')
xlabel('UPDRS')
ylabel('% theta small')
subplot(1,3,2)
wjn_corr_scatter(ngm',ntm')
xlabel('UPDRS')
ylabel('% theta medium')
subplot(1,3,3)
wjn_corr_scatter(ngl',ntl')
xlabel('UPDRS')
ylabel('% theta large')