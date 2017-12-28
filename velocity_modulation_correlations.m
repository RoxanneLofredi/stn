clear all
close all

root='';
cd(root);

pat=patienten;

prefix = 'sbgmtf_efspmeeg_';

conds={'velocity-klein','velocity-gross'};

n=0

for a=1:length(pat);
    
    id=strsplit(pat{a}.id, '_');
    
    if ~isfield(pat{a},'bad') && pat{a}.med && ~strcmp(pat{a}.id,'r004IR48_L_ON') && isfield(pat{a},'updrs');
    
        n=n+1;
        
        D=spm_eeg_load([prefix pat{a}.id '.mat']);
        
        if id{2}=='L'  
            cl=ci('STNR', D.chanlabels);
        else
            cl=ci('STNL', D.chanlabels);
        end
        
        med(n) = pat{a}.med;
     
        nid(n) = pat{a}.npatient;
        sv=ci('velocity-klein',D.conditions);
        gv=ci('velocity-gross',D.conditions);
        
        if pat{a}.med
            updrs(n) = pat{a}.updrs(1);
        else
            updrs(n) = pat{a}.updrs(2);
        end

        
        gs(n,:) = squeeze(nanmean(nanmean(nanmean(nanmean(D(cl,40:90,sc(D.time,-0.2):sc(D.time,0.2),sv),1),2),3),4));
  gl(n,:) = squeeze(nanmean(nanmean(nanmean(nanmean(D(cl,40:90,sc(D.time,-0.2):sc(D.time,0.2),gv),1),2),3),4));
    end
    
end

ngs=wjn_match_cases(gs,nid)
ngl=wjn_match_cases(gl,nid)
nupdrs = wjn_match_cases(updrs,nid);

plot([ngs,ngl])
pg = (ngs./ngl)*100;
pg(16) = [];
nupdrs(16)=[];
% ng = ngl-ngs;
close all
wjn_corr_scatter(nupdrs',pg')