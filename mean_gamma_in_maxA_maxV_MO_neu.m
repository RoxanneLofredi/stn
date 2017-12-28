clear all
close all


pat=patienten;

prefix = 'sbgmtf_efspmeeg_';

conds={'onset', 'velocity', 'amplitude'};

n=0

for a=1:length(pat);
    
    id=strsplit(pat{a}.id, '_');
    
    if ~isfield(pat{a},'bad') && pat{a}.med && ~strcmp(pat{a}.id,'r004IR48_L_ON');
   
        n=n+1;
        
        D=spm_eeg_load([prefix pat{a}.id '.mat']);
        
        if id{2}=='L'  
            cl=ci('STNR', D.chanlabels);
        else
            cl=ci('STNL', D.chanlabels);
        end
        
        onset=ci('onset',D.conditions);
        amp=ci('amplitude',D.conditions);
        velo=ci('velocity',D.conditions);
       
        damp(n,:) = squeeze(nanmean(nanmean(nanmean(nanmean(D(cl,40:90,sc(D.time,-0.2):sc(D.time,0.2),amp),1),2),3),4));
        dvelo(n,:) = squeeze(nanmean(nanmean(nanmean(nanmean(D(cl,40:90,sc(D.time,-0.2):sc(D.time,0.2),velo),1),2),3),4));
        donset(n,:)= squeeze(nanmean(nanmean(nanmean(nanmean(D(cl,40:90,sc(D.time,-0.2):sc(D.time,0.2),onset),1),2),3),4));
         
    end
    
end

dall=[damp, dvelo, donset];
veloamp=dall(:,2)-dall(:,1);
velonset=dall(:,2)-dall(:,3);
amponset=dall(:,3)-dall(:,1);

pveloamp=wjn_pt(veloamp);
pvelonset=wjn_pt(velonset);
pamponset=wjn_pt(amponset);

cmap = colorlover(6,0);
cmap = cmap([2 5 4 2 5 4],:);

figure,
mybar(dall,cmap)
set(gca,'XTick',[1 2 3],'XTickLabel',{'maxA','maxV','MO'})        
ylabel('Relative power change [%,SEM]');
title('Mean gamma aligned maxA, maxV, MO');
sigbracket('P<0.001',[1 2],20,9)
sigbracket(['P=' num2str(pvelonset,1)],[2 3],20,9)
sigbracket(['P=' num2str(pamponset,1)],[1 3],24,9);
ylim([0 28]);
figone(7,10)
