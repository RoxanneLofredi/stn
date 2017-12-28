close all
clear all

root=;
cd(root)
mkdir('TF_CORR_smooth_BLtfplot')

pat = patienten;

side={'ipsi','contra'}; 

for d=1:length(side);
for a=1:length(pat) 
    close all
    
    if d==1 && a==3;
    D=spm_eeg_load(ffind(['grtf_ef*' pat{a}.id '*.mat'],0));
    D=wjn_tf_smooth(D.fullfile,4,100);
    s=pat{a}.side;
    end
    
    if ~isfield(pat{a},'bad');

        D=spm_eeg_load(ffind(['grtf_ef*' pat{a}.id '*.mat'],0));
        D=wjn_tf_smooth(D.fullfile,4,100);
        
if d==1;
    s=pat{a}.side;
else
    s = wjn_cl_side(pat{a}.side);
end
        
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
Dr=wjn_tf_corr_avg_chan(D.fullfile,rt,['rt-onset-' side{d}],'onset','spearman',['STN' s],[-1 2],[1 100]);
Da=wjn_tf_corr_avg_chan(D.fullfile,acc,['acc-amp-' side{d}],'amplitude','spearman',['STN' s],[-1 2],[1 100]);
Dc=wjn_tf_corr_avg_chan(D.fullfile,vel,['vel-vel-' side{d}],'velocity','spearman',['STN' s],[-1 2],[1 100]);

figure,
subplot(2,2,[1 3])
imagesc(Dc.time,Dc.frequencies,squeeze(Dc(1,:,:,1))),
axis xy,
hold on,
contour(Dc.time,Dc.frequencies,squeeze(Dc(1,:,:,2))<.05,1,'color','k')
caxis([-.5 .5])
xlabel('Time [s]')
ylabel('Frequency [Hz]')
title(['VELOCITY_' side{d} '_BLsmooth' strrep(pat{a}.id,'_',' ')])
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
myprint(fullfile('TF_CORR_smooth_BLtfplot',['velvel_sTF_CORR_' side{d} '_' pat{a}.id])) 


figure,
subplot(2,2,[1 3])
imagesc(Da.time,Da.frequencies,squeeze(Da(1,:,:,1))),
axis xy,
hold on,
contour(Da.time,Da.frequencies,squeeze(Da(1,:,:,2))<.05,1,'color','k')
caxis([-.5 .5])
xlabel('Time [s]')
ylabel('Frequency [Hz]')
title(['ACCURACY ' side{d} '_BLsmooth' strrep(pat{a}.id,'_',' ')])
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
myprint(fullfile('TF_CORR_smooth_BLtfplot',['accamp_sTF_CORR_' side{d} '_' pat{a}.id])) 


figure,
subplot(2,2,[1 3])
imagesc(Dr.time,Dr.frequencies,squeeze(Dr(1,:,:,1))),
axis xy,
hold on,
contour(Dr.time,Dr.frequencies,squeeze(Dr(1,:,:,2))<.05,1,'color','k')
caxis([-.5 .5])
xlabel('Time [s]')
ylabel('Frequency [Hz]')
title(['REACTION TIME ' side{d} '_BLsmooth' strrep(pat{a}.id,'_',' ')])
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
myprint(fullfile('TF_CORR_smooth_BLtfplot',['RTonset_sTF_CORR_' side{d} '_' pat{a}.id])) 
catch
    disp([' Nummer ' num2str(a) ' hat nicht funktioniert.'])
end
    end
end
end

cd('E:\Roxanne\STN_Rotameter\TF_CORR_smooth_BLtfplot');
save('TF_corr_velaccrt_figures_allpatON_ipsicontra_BLtfsmoothed');

%% TF group correlation

close all
clear all

root='E:\Roxanne\STN_Rotameter';
cd(root)

addpath E:\Roxanne\scripts;
pat = patienten;

n=0; 

cond = {'vel','rt','acc'};
align={'vel', 'onset', 'amp'};
side={'ipsi', 'contra'};

non = 0;
noff = 0;

for d=1:length(side);
for c = 1:length(cond);
    
    non = 0;
    noff = 0;
    
for a=1:length(pat);
    
    if ~isfield(pat{a},'bad') && pat{a}.med % && pat{a}.onoff;
        non=non+1;
        
        if c==1;
        D=spm_eeg_load([cond{c} '-' align{c} side{d} '_rp_sgrtf_efspmeeg_' pat{a}.id])
        else
        D=spm_eeg_load([cond{c} '-' align{c} '-' side{d} '_rp_sgrtf_efspmeeg_' pat{a}.id]);
        end
%         D=wjn_tf_smooth(D.fullfile,4,100);
        don.(side{d}).(cond{c})(non,:,:) = squeeze(D(1,:,:,1));
        pon.(side{d}).(cond{c})(non,:,:) = squeeze(D(1,:,:,2));
        
    elseif ~isfield(pat{a},'bad') && ~pat{a}.med;
        noff=noff+1;
        
        if c==1;
        D=spm_eeg_load([cond{c} '-' align{c} side{d} '_rp_sgrtf_efspmeeg_' pat{a}.id])
        else        
        D=spm_eeg_load([cond{c} '-' align{c} '-' side{d} '_rp_sgrtf_efspmeeg_' pat{a}.id]);
        end
%         D=wjn_tf_smooth(D.fullfile,4,100);
        doff.(side{d}).(cond{c})(noff,:,:) = squeeze(D(1,:,:,1));
        poff.(side{d}).(cond{c})(noff,:,:) = squeeze(D(1,:,:,2));
        
    end
end
end
end

f=D.frequencies;
t = D.time;

cond = {'vel','rt','acc'};
align={'vel', 'onset', 'amp'};
side={'ipsi', 'contra'};

for d=1:2;
    for c=1:3;
        for a = 1:length(f);
            for b = 1:length(t);
                
                gpon.(side{d}).(cond{c})(a,b) = permTest(500,don.(side{d}).(cond{c})(:,a,b));
%         gpoff(a,b) = permTest(500,doff.(cond{c})(:,a,b));
        
%         gpon(a,b) = signrank(don.(cond{c})(:,a,b));
%         gpoff(a,b) = signrank(doff.(cond{c})(:,a,b));
            end
        end

sgpon.(side{d}).(cond{c})=zeros(size(gpon.(side{d}).(cond{c})));
sgpon.(side{d}).(cond{c})(gpon.(side{d}).(cond{c})<=.05)=-.5;
sgpon.(side{d}).(cond{c})(gpon.(side{d}).(cond{c})<=fdr_bh(gpon.(side{d}).(cond{c}),0.05))=-1;
    end
end

cd('E:\Roxanne\STN_Rotameter\TF_CORR_smooth_BLtfplot');
save('TF_corr_velaccrt_figures_allpatON_ipsicontra_BLtfsmoothed_permtest');

%% Figure

cd('E:\Roxanne\STN_Rotameter\TF_CORR_smooth_BLtfplot');
load ('TF_corr_velaccrt_figures_allpatON_ipsicontra_BLtfsmoothed_permtest');

cond = {'vel','rt','acc'};
align={'vel', 'onset', 'amp'};
side={'ipsi', 'contra'};

% cmb= colormap([0 0 0; 1 1 1]);

for d=1:2;
    
for c=1:3;
 
sgpon.(side{d}).(cond{c})=zeros(size(gpon.(side{d}).(cond{c})));
sgpon.(side{d}).(cond{c})(gpon.(side{d}).(cond{c})<=fdr_bh(gpon.(side{d}).(cond{c}),0.05))=-1;

[sgponccluster,sgponciclusters,sgponccluster_sizes] = wjn_remove_clusters(sgpon.(side{d}).(cond{c}),150);

cma = colormap([1 1 1; 0 0 0]);

figure,

ax2=subplot(1,2,1)
imagesc(t,f,squeeze(nanmean(don.(side{d}).(cond{c}),1)))
axis xy,
% contour(t,f,gpon<=fdr_bh(gpon,0.05),1,'color','r')
% contour(t,f,gpon<=.05,1,'color','k')
caxis([-.04 .1])
%colorbar
xlim([-1 1])
ylim([1 100])
% title(['ON ' (side{d}) ' ' (cond{c})])

ax1=subplot(1,2,2)
imagesc(t,f,sgponccluster)
axis xy,
xlim([-1 1])
ylim([1 100])
% 
% if c==2 || c==3
%     colormap(ax1, cmb);
% end

colormap(ax1,cma);
colormap(ax2,'parula');
suptitle([side{d} '-' cond{c} '-' align{c} ', p-values<.05, FDR corrected, cluster 150'])

figone(9,20)

myprint(fullfile('E:\Roxanne\STN_Rotameter\figures', [side{d} '_' cond{c} '_allON(n=16)_BLsmoothed_cluster150']));

end
end



%% Figure

load ('TF_corr_velaccrt_figures_allpatON_ipsicontra_BLtfsmoothed');

cond = {'vel','rt','acc'};
align={'vel', 'onset', 'amp'};
side={'ipsi', 'contra'};

% cmb= colormap([0 0 0; 1 1 1]);

    
for c=1:3;
 
sgpon.ipsi.(cond{c})=zeros(size(gpon.ipsi.(cond{c})));
sgpon.ipsi.(cond{c})(gpon.ipsi.(cond{c})<=fdr_bh(gpon.ipsi.(cond{c}),0.05))=-1;

sgpon.contra.(cond{c})=zeros(size(gpon.contra.(cond{c})));
sgpon.contra.(cond{c})(gpon.contra.(cond{c})<=fdr_bh(gpon.contra.(cond{c}),0.05))=-1;


[sgponccluster,sgponciclusters,sgponccluster_sizes] = wjn_remove_clusters(sgpon.contra.(cond{c}),150);
[sgponicluster,sgponiiclusters,sgponicluster_sizes] = wjn_remove_clusters(sgpon.ipsi.(cond{c}),150);


cma = colormap([1 1 1; 0 0 0]);



figure,

ax2=subplot(2,2,1)
colormap(ax2,'parula');
imagesc(t,f,squeeze(nanmean(don.contra.(cond{c}),1)))
axis xy,
caxis([-.04 .1])
xlim([-1 1])
ylim([1 100])
xlabel('Time [s]');
ylabel ('Frequency [Hz]');
k = get(gca,'ylabel');
k.String = {'contralateral',k.String};
set(gca,'ylabel',k);


ax1=subplot(2,2,2)
colormap(ax1,cma);
imagesc(t,f,sgponccluster)
xlabel('Time [s]');
ylabel ('Frequency [Hz]');
axis xy,
xlim([-1 1])
ylim([1 100])


ax2=subplot(2,2,3)
colormap(ax2,'parula');
imagesc(t,f,squeeze(nanmean(don.ipsi.(cond{c}),1)))
axis xy,
caxis([-.04 .1])
xlim([-1 1])
ylim([1 100])
ylabel ('Frequency [Hz]');
k = get(gca,'ylabel');
k.String = {'ipsilateral',k.String};
set(gca,'ylabel',k);
xlabel('Time [s]');

ax1=subplot(2,2,4)
colormap(ax1,cma);
imagesc(t,f,sgponicluster)
axis xy,
xlim([-1 1])
ylim([1 100])
xlabel('Time [s]');
ylabel ('Frequency [Hz]');

suptitle([cond{c} '-' align{c} ', p-values<.05, FDR corrected, cluster 150'])

figone(13,15)

end
end
