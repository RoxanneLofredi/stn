%% Figure all ipsi vs contra
clear all
close all

load('TF_pvalues_figures_allpat_ONOFF_ipsicontra');

pat = patienten;

[smpcontraipsicluster,smpcontraipsiclusters,smpcontraipsicluster_sizes] = wjn_remove_clusters(smpcontraipsi,10);

cma = colormap([1 1 1; 0 0 0]);

figure,
imagesc(t,f,smpcontraipsicluster)
hold on
title('p<.05, FDR, ipsi-vs contra,clust10')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy; xlim([-2 3]);
xlim([-2 3])
figone (7,8);
colormap(cma);

%% Figure only ON recordings, between cond

clear all
close all

cd('E:\Roxanne\STN_Rotameter');
load('TF_pvalues_figures_allpat_ONOFF_ipsicontra');

addpath E:\Roxanne\scripts;
pat = patienten;

[skgponcluster,skgponiclusters,skgponcluster_sizes] = wjn_remove_clusters(skgpon,150);
[skmponcluster,skmponiclusters,skmponcluster_sizes] = wjn_remove_clusters(skmpon,150);
[smgponcluster,smgponiclusters,smgponcluster_sizes] = wjn_remove_clusters(smgpon,150);

[skgponicluster,skgponiiclusters,skgponicluster_sizes] = wjn_remove_clusters(skgponi,150);
[skmponicluster,skmponiiclusters,skmponicluster_sizes] = wjn_remove_clusters(skmponi,150);
[smgponicluster,smgponiiclusters,smgponicluster_sizes] = wjn_remove_clusters(smgponi,150);



cma = colormap([1 1 1; 0 0 0]);

figure,
subplot(2,3,1)
imagesc(t,f,skgponcluster)
hold on
title('small - large')
ylabel('Frequency [Hz]')
k = get(gca,'ylabel');
k.String = {'contralateral',k.String};
set(gca,'ylabel',k);
xlabel('Time [s]');axis xy; xlim([-2 3]);

subplot(2,3,2)
imagesc(t,f,skmponcluster)
title('small - medium')
axis xy; xlim([-2 3]);
ylabel('Frequency [Hz]')
xlabel('Time [s]');
k = get(gca,'title');
k.String = {'All ON (n=16),FDR p<.05, Cluster 150',k.String};
set(gca,'title',k);

subplot(2,3,3)
imagesc(t,f,smgponcluster)
axis xy; xlim([-2 3]);
title('medium - large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');
colormap(cma);

subplot(2,3,4)
imagesc(t,f,skgponicluster)
hold on
ylabel('Frequency [Hz]')
k = get(gca,'ylabel');
k.String = {'ipsilateral',k.String};
set(gca,'ylabel',k);
xlabel('Time [s]');axis xy; xlim([-2 3]);

subplot(2,3,5)
imagesc(t,f,skmponicluster)
axis xy; xlim([-2 3]);
ylabel('Frequency [Hz]')
xlabel('Time [s]');

subplot(2,3,6)
imagesc(t,f,smgponicluster)
axis xy; xlim([-2 3]);
ylabel('Frequency [Hz]')
xlabel('Time [s]');
colormap(cma);


figone(15,18)

save('TF_pvalues_figures_allpat_ONOFF_ipsicontra');


%% Figure ipsi & contralateral images
clear all
close all

pat = patienten;

[skponccluster,skponciclusters,skponccluster_sizes] = wjn_remove_clusters(skponc,150);
[smponccluster,smponciclusters,skmponcluster_sizes] = wjn_remove_clusters(smponc,150);
[sgponccluster,sgponciclusters,sgponccluster_sizes] = wjn_remove_clusters(sgponc,150);

[skponicluster,skponiiclusters,skponicluster_sizes] = wjn_remove_clusters(skponi,150);
[smponicluster,smponiiclusters,skmponicluster_sizes] = wjn_remove_clusters(smponi,150);
[sgponicluster,sgponiiclusters,sgponicluster_sizes] = wjn_remove_clusters(sgponi,150);

[nii,xaxis,yaxis]=wjn_disp_nii('average_contralateral_onset-gross_ON.nii',[-3 3],[1 100]);

gc=spm_vol('average_contralateral_onset-gross_ON.nii');
gcread=spm_read_vols(gc);
gcreadx = linspace(xaxis(1),xaxis(end),size(gcread,2));
gcready = linspace(yaxis(1),yaxis(end),size(gcread,1));

mc=spm_vol('average_contralateral_onset-mittel_ON.nii');
mcread=spm_read_vols(mc);
mcreadx = linspace(xaxis(1),xaxis(end),size(mcread,2));
mcready = linspace(yaxis(1),yaxis(end),size(mcread,1));

kc=spm_vol('average_contralateral_onset-klein_ON.nii');
kcread=spm_read_vols(kc);
kcreadx = linspace(xaxis(1),xaxis(end),size(kcread,2));
kcready = linspace(yaxis(1),yaxis(end),size(kcread,1));

cd('E:\Roxanne\STN_Rotameter\ipsilateral_images\averages');

gi=spm_vol('average_ipsilateral_onset-gross_ON.nii');
giread=spm_read_vols(gi);
gireadx = linspace(xaxis(1),xaxis(end),size(giread,2));
giready = linspace(yaxis(1),yaxis(end),size(giread,1));

mi=spm_vol('average_ipsilateral_onset-mittel_ON.nii');
miread=spm_read_vols(mi);
mireadx = linspace(xaxis(1),xaxis(end),size(miread,2));
miready = linspace(yaxis(1),yaxis(end),size(miread,1));

ki=spm_vol('average_ipsilateral_onset-klein_ON.nii');
kiread=spm_read_vols(ki);
kireadx = linspace(xaxis(1),xaxis(end),size(kiread,2));
kiready = linspace(yaxis(1),yaxis(end),size(kiread,1));

cma = colormap([1 1 1; 0 0 0]);


figure,

ax2=subplot(3,4,1)
imagesc(newt,kcready,kcread);
title('small')
k = get(gca,'title');
k.String = {'contralateral',k.String};
set(gca,'title',k);
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy; xlim([-2 3]);
caxis([-15 30]);

ax1=subplot(3,4,2)
imagesc(t,f,skponccluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy; xlim([-2 3]);

colormap(ax1,cma);
colormap(ax2,'parula');

ax2=subplot(3,4,3)
imagesc(newt,kiready,kiread);
title('small')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy; xlim([-2 3]);
k = get(gca,'title');
k.String = {'ipsilateral',k.String};
set(gca,'title',k);
caxis([-15 30]);

ax1=subplot(3,4,4)
imagesc(t,f,skponicluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy; xlim([-2 3]);

colormap(ax1,cma)
colormap(ax2,'parula');

ax2=subplot(3,4,5)
imagesc(newt,mcready,mcread);
title('medium')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy; xlim([-2 3]);
caxis([-15 30]);

ax1=subplot(3,4,6)
imagesc(t,f,smponccluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy; xlim([-2 3]);

colormap(ax1,cma)
colormap(ax2,'parula');

ax2=subplot(3,4,7)
imagesc(newt,miready,miread);
title('medium')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy; xlim([-2 3]);
caxis([-15 30]);

ax1=subplot(3,4,8)
imagesc(t,f,smponicluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy; xlim([-2 3]);

colormap(ax1,cma)
colormap(ax2,'parula');

ax2=subplot(3,4,9)
imagesc(newt,gcready,gcread);
title('large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy; xlim([-2 3]);
caxis([-15 30]);

ax1=subplot(3,4,10)
imagesc(t,f,sgponccluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy; xlim([-2 3]);

colormap(ax1,cma)
colormap(ax2,'parula');

ax2=subplot(3,4,11)
imagesc(newt,giready,giread);
title('large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy; xlim([-2 3]);
caxis([-15 30]);

ax1=subplot(3,4,12)
imagesc(t,f,sgponicluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy; xlim([-2 3]);

colormap(ax1,cma)
colormap(ax2,'parula');

suptitle('All ON (n=16), permTest, p<0.05, FDR corrected, Cluster 150');

figone(25,28)

%% Figures ipsi vs contralateral

clear all
close all
pat = patienten;

[smpcontraipsicluster,smpcontraipsiiclusters,smpcontraipsi_sizes] = wjn_remove_clusters(smpcontraipsi,150);

cma = colormap([1 1 1; 0 0 0]);

figure,
imagesc(t,f, smpcontraipsicluster);
title('Contra- vs ipsilateral (n=16), p<.05, FDR corrected, permTest, cluster 150')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy; xlim([-2 3]);
colormap(cma)

figure,
subplot(1,3,1)
imagesc(t,f,skponicl)
hold on
title('small')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy; xlim([-2 3]);

subplot(1,3,2)
imagesc(t,f,smponicl)
title('medium')
axis xy; xlim([-2 3]);
ylabel('Frequency [Hz]')
k = get(gca,'ylabel');
k.String = {'ON',k.String};
set(gca,'ylabel',k);
xlabel('Time [s]');
k = get(gca,'title');
k.String = {'Ipsi- vs contralateral, n=16',k.String};
set(gca,'title',k);

subplot(1,3,3)
imagesc(t,f,sgponicl)
axis xy; xlim([-2 3]);
title('large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');

colormap('gray')
figone(7,20)

%%
save('TFplot_pvalues_betweencond_allpat');