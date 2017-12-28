close all
clear all

root= ;
cd(root)
mkdir('TF_CORR_ipsilateral')
addpath E:\Roxanne\scripts;
pat = patienten;

for a=1:length(pat) 
    close all
    keep pat a
    if ~isfield(pat{a},'bad');

        
        D=spm_eeg_load(ffind(['grtf_ef*' pat{a}.id '*.mat'],0));
        

        s = pat{a}.side;
        
% D=wjn_remove_bad_trials(D.fullfile);
ti=ci('velocity',D.conditions);
vel = D.btrialvel(ti);
acc = D.btrialacc(ti);
try 
rt = D.trialonset(ci('onset',D.conditions))-D.trialonset(ci('cue',D.conditions));
% D.rt = rt;
% save(D)
% D=wjn_grandchamp_mean(D.fullfile,[-2000 -500],'gain','cue');
% rt = D.rt;
Dr=wjn_tf_corr_avg_chan(D.fullfile,rt,'rt-onset-ipsi','onset','spearman',['STN' s],[-1 2],[1 100]);
Da=wjn_tf_corr_avg_chan(D.fullfile,acc,'acc-amp-ipsi','amplitude','spearman',['STN' s],[-1 2],[1 100]);
Dc=wjn_tf_corr_avg_chan(D.fullfile,vel,'vel-vel-ipsi','velocity','spearman',['STN' s],[-1 2],[1 100]);

figure,
subplot(2,2,[1 3])
imagesc(Dc.time,Dc.frequencies,squeeze(Dc(1,:,:,1))),
axis xy,
hold on,
contour(Dc.time,Dc.frequencies,squeeze(Dc(1,:,:,2))<.05,1,'color','k')
caxis([-.5 .5])
xlabel('Time [s]')
ylabel('Frequency [Hz]')
title(['VELOCITY IPSI' strrep(pat{a}.id,'_',' ')])
subplot(2,2,2)
[Rg.(pat{a}.id),Pg.(pat{a}.id),x,y] =wjn_tf_corr_avg_tf_chan(D.fullfile,vel,'velocity','spearman',['STN' s],[-.1 .5],[40 85]);
[nx,ny]=mycorrline(x,y,min(x),max(x));
scatter(x,y,'ro','filled','markeredgecolor','w')
hold on
myline(nx,ny,'color',[.5 .5 .5])
title({'Gamma (40 - 85 Hz) -.1 - .5 sec',['\rho = ' num2str(Rg.(pat{a}.id)) ' P = ' num2str(Pg.(pat{a}.id))]})
subplot(2,2,4)
[Rb.(pat{a}.id),Pb.(pat{a}.id),x,y] =wjn_tf_corr_avg_tf_chan(D.fullfile,vel,'velocity','spearman',['STN' s],[-.5 1],[13 30]);
[nx,ny]=mycorrline(x,y,min(x),max(x));
scatter(x,y,'bo','filled','markeredgecolor','w')
hold on
myline(nx,ny,'color',[.5 .5 .5])
title({'Gamma (13 - 35 Hz) -.5 - 1 sec',['\rho = ' num2str(Rb.(pat{a}.id)) ' P = ' num2str(Pb.(pat{a}.id))]})
figone(13,20);
myprint(fullfile('TF_CORR_ipsilateral',['velvel_TF_CORR_' pat{a}.id])) 


figure,
subplot(2,2,[1 3])
imagesc(Da.time,Da.frequencies,squeeze(Da(1,:,:,1))),
axis xy,
hold on,
contour(Da.time,Da.frequencies,squeeze(Da(1,:,:,2))<.05,1,'color','k')
caxis([-.5 .5])
xlabel('Time [s]')
ylabel('Frequency [Hz]')
title(['ACCURACY IPSI' strrep(pat{a}.id,'_',' ')])
subplot(2,2,2)
[Rg.(pat{a}.id),Pg.(pat{a}.id),x,y] =wjn_tf_corr_avg_tf_chan(D.fullfile,acc,'amplitude','spearman',['STN' s],[-.1 .5],[40 85]);
[nx,ny]=mycorrline(x,y,min(x),max(x));
scatter(x,y,'ro','filled','markeredgecolor','w')
hold on
myline(nx,ny,'color',[.5 .5 .5])
title({'Gamma (40 - 85 Hz) -.1 - .5 sec',['\rho = ' num2str(Rg.(pat{a}.id)) ' P = ' num2str(Pg.(pat{a}.id))]})
subplot(2,2,4)
[Rb.(pat{a}.id),Pb.(pat{a}.id),x,y] =wjn_tf_corr_avg_tf_chan(D.fullfile,acc,'amplitude','spearman',['STN' s],[-.5 1],[13 30]);
[nx,ny]=mycorrline(x,y,min(x),max(x));
scatter(x,y,'bo','filled','markeredgecolor','w')
hold on
myline(nx,ny,'color',[.5 .5 .5])
title({'Gamma (13 - 35 Hz) -.5 - 1 sec',['\rho = ' num2str(Rb.(pat{a}.id)) ' P = ' num2str(Pb.(pat{a}.id))]})
figone(13,20);
myprint(fullfile('TF_CORR_ipsilateral',['accamp_TF_CORR_' pat{a}.id])) 


figure,
subplot(2,2,[1 3])
imagesc(Dr.time,Dr.frequencies,squeeze(Dr(1,:,:,1))),
axis xy,
hold on,
contour(Dr.time,Dr.frequencies,squeeze(Dr(1,:,:,2))<.05,1,'color','k')
caxis([-.5 .5])
xlabel('Time [s]')
ylabel('Frequency [Hz]')
title(['REACTION TIME IPSI' strrep(pat{a}.id,'_',' ')])
subplot(2,2,2)
[Rg.(pat{a}.id),Pg.(pat{a}.id),x,y] =wjn_tf_corr_avg_tf_chan(D.fullfile,rt,'onset','spearman',['STN' s],[-.1 .5],[40 85]);
[nx,ny]=mycorrline(x,y,min(x),max(x));
scatter(x,y,'ro','filled','markeredgecolor','w')
hold on
myline(nx,ny,'color',[.5 .5 .5])
title({'Gamma (40 - 85 Hz) -.1 - .5 sec',['\rho = ' num2str(Rg.(pat{a}.id)) ' P = ' num2str(Pg.(pat{a}.id))]})
subplot(2,2,4)
[Rb.(pat{a}.id),Pb.(pat{a}.id),x,y] =wjn_tf_corr_avg_tf_chan(D.fullfile,rt,'onset','spearman',['STN' s],[-.5 1],[13 30]);
[nx,ny]=mycorrline(x,y,min(x),max(x));
scatter(x,y,'bo','filled','markeredgecolor','w')
hold on
myline(nx,ny,'color',[.5 .5 .5])
title({'Gamma (13 - 35 Hz) -.5 - 1 sec',['\rho = ' num2str(Rb.(pat{a}.id)) ' P = ' num2str(Pb.(pat{a}.id))]})
figone(13,20);
myprint(fullfile('TF_CORR_ipsilateral',['RTonset_TF_CORR_' pat{a}.id])) 
catch
    disp([' Nummer ' num2str(a) ' hat nicht funktioniert.'])
end
    end
end


