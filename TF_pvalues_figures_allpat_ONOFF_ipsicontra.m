%% Ipsilateral average images

close all
clear all

addpath E:\Roxanne\scripts;
pat = patienten;

cond= {'onset-klein', 'onset-mittel', 'onset-gross'};
imagefolder = 'E:\Roxanne\STN_Rotameter\ipsilateral_images\';

cd(imagefolder)

for a=1:length(pat);
    if ~isfield (pat{a}, 'bad');

       for b=1:length(cond);
       id= pat{a}.id(1:8);

       images = ffind(['*' cond{b} '*' id '*' 'ON' '*' '.nii']);
       
       if strcmp(pat{a}.id,'r004IR48') 
           images = ffind(['*' cond{b} '*' id '*' 'L_ON' '*' '.nii']);
       end
       
            if ~isempty(images);
                outimage = (fullfile(cd,'averages', ['average_ipsilateral_' id '_' cond{b} '_ON.nii']));
                if length(images) == 2; 
                    spm_imcalc(images,outimage,'(i1+i2)/2');
                else
                    spm_imcalc(images,outimage, 'i1/1');
                end
            end
    end
end
end

imagefolder = 'E:\Roxanne\STN_Rotameter\ipsilateral_images\averages';

cd(imagefolder)

konifiles = ffind('*average_ipsilateral_r*_onset*klein*ON*.nii');
monifiles = ffind('*average_ipsilateral_r*_onset*mittel*ON*.nii');
gonifiles = ffind('*average_ipsilateral_r*_onset*gross*ON*.nii');


for a = 1:length(konifiles)
    [koni(a,:,:),t,f] = wjn_disp_nii(konifiles{a},[-3 3],[1 100],0);
    [moni(a,:,:),t,f] = wjn_disp_nii(monifiles{a},[-3 3],[1 100],0);
    [goni(a,:,:),t,f] = wjn_disp_nii(gonifiles{a},[-3 3],[1 100],0);
end

%% Grand average, all conditions, ipsi

close all
clear all

addpath E:\Roxanne\scripts;
pat = patienten;

cond= {'onset-klein', 'onset-mittel', 'onset-gross'};
imagefolder = 'E:\Roxanne\STN_Rotameter\ipsilateral_images\';

cd(imagefolder)

for a=1:length(pat);
    if ~isfield (pat{a}, 'bad');
        
       id= pat{a}.id(1:8);
       
       if strcmp(id,'r004IR48') 
           images = ffind(['*' id '*' 'L_ON' '*' '.nii']);
       else
       images = ffind(['*' id '*' 'ON' '*' '.nii']);
       end
       
            if ~isempty(images);
                outimage = (fullfile(cd,'averages', ['average_ipsilateral_' id '_sml_ON.nii']));
                if length(images) == 6; 
                    spm_imcalc(images,outimage,'(i1+i2+i3+i4+i5+i6)/6');
                else
                    spm_imcalc(images,outimage, 'i1+i2+i3/3');
                end
            end
    end
end


imagefolder = '';

cd(imagefolder)

onifiles = ffind('*average_ipsilateral_r*sml*ON*.nii');

for a = 1:length(onifiles)
    [oni(a,:,:),t,f] = wjn_disp_nii(onifiles{a},[-3 3],[1 100],1);
end



%% Grand average, contralateral, all conditions

close all
clear all

addpath E:\Roxanne\scripts;
pat = patienten;

cond= {'onset-klein', 'onset-mittel', 'onset-gross'};
imagefolder = 'E:\Roxanne\STN_Rotameter\contralateral_images\';

cd(imagefolder)

for a=1:length(pat);
    if ~isfield (pat{a}, 'bad');
        
       id= pat{a}.id(1:8);
       
       if strcmp(id,'r004IR48') 
           images = ffind(['*' id '*' 'R_ON' '*' '.nii']);
       else
       images = ffind(['*' id '*' 'ON' '*' '.nii']);
       end
       
            if ~isempty(images);
                outimage = (fullfile(cd,'averages', ['average_contralateral_' id '_sml_ON.nii']));
                if length(images) == 6; 
                    spm_imcalc(images,outimage,'(i1+i2+i3+i4+i5+i6)/6');
                else
                    spm_imcalc(images,outimage, 'i1+i2+i3/3');
                end
            end
    end
end


imagefolder = 'E:\Roxanne\STN_Rotameter\contralateral_images\averages';

cd(imagefolder)

oncfiles = ffind('*average_contralateral_r*sml*ON*.nii');

for a = 1:length(oncfiles)
    [onc(a,:,:),t,f] = wjn_disp_nii(oncfiles{a},[-3 3],[1 100],1);
end


%% Grand average, contralateral, all conditions

close all
clear all

addpath E:\Roxanne\scripts;
pat = patienten;

cond= {'onset-klein', 'onset-mittel', 'onset-gross'};
imagefolder = 'E:\Roxanne\STN_Rotameter\contralateral_images\';

cd(imagefolder)

for a=1:length(pat);
    if strcmp(~isfield (pat{a}, 'bad');
        
       id= pat{a}.id(1:8);
       
       if strcmp(id,'r004IR48') 
           images = ffind(['*' id '*' 'R_ON' '*' '.nii']);
       else
       images = ffind(['*' id '*' 'ON' '*' '.nii']);
       end
       
            if ~isempty(images);
                outimage = (fullfile(cd,'averages', ['average_contralateral_' id '_sml_ON.nii']));
                if length(images) == 6; 
                    spm_imcalc(images,outimage,'(i1+i2+i3+i4+i5+i6)/6');
                else
                    spm_imcalc(images,outimage, 'i1+i2+i3/3');
                end
            end
    end
end


imagefolder = '';

cd(imagefolder)

oncfiles = ffind('*average_contralateral_r*sml*ON*.nii');

for a = 1:length(oncfiles)
    [onc(a,:,:),t,f] = wjn_disp_nii(oncfiles{a},[-3 3],[1 100],1);
end


%% frame for fdr correction

cd('');
load('TFplot_pvalues_betweencond_allpat');

i_fdr = zeros(size(gponi));
i_fdr(1:100,sc(t,-2):sc(t,3))=1;
i_fdr=find(i_fdr);

save('TFplot_pvalues_betweencond_allpat');

%% Grand average contra vs ipsilateral

pat = patienten;

cd('');
load('TF_pvalues_figures_allpat_ONOFF_ipsicontra');

for a = 1:length(t);
    for b = 1:length(f);

ponicl (b,a)= wjn_pt(onc(:,b,a)-oni(:,b,a));

    end
end


for a=1:length(t);
   for b=1:length(f);
    
mcontra(:,b,a) = mean([konc(:,b,a),monc(:,b,a),gonc(:,b,a)],2);
mipsi(:,b,a) = mean([koni(:,b,a),moni(:,b,a),goni(:,b,a)],2);

   end
end

for a=1:length(t);
   for b=1:length(f);
mpcontraipsi(b,a)=wjn_pt(mcontra(:,b,a)-mipsi(:,b,a));
   end
end


smpcontraipsi=zeros(size(mpcontraipsi));
smpcontraipsi(mpcontraipsi<=.05)=-.5;
smpcontraipsi(mpcontraipsi<=fdr_bh(mpcontraipsi(i_fdr),0.05))=-1;

skponicl=zeros(size(kponicl));
skponicl(kponicl<=.05)=-.5;
skponicl(kponicl<=fdr_bh(kponicl(i_fdr),0.05))=-1;

smponicl=zeros(size(mponicl));
smponicl(mponicl<=.05)=-.5;
smponicl(mponicl<=fdr_bh(mponicl(i_fdr),0.05))=-1;

sgponicl=zeros(size(gponicl));
sgponicl(gponicl<=.05)=-.5;
sgponicl(gponicl<=fdr_bh(gponicl(i_fdr),0.05))=-1;

cd('E:\Roxanne\STN_Rotameter');
save('TF_pvalues_figures_allpat_ONOFF_ipsicontra');

%% Ipsilateral ON permTest, all pat

clear all
close all

cd('E:\Roxanne\STN_Rotameter');
load('TF_pvalues_figures_allpat_ONOFF_ipsicontra');

addpath E:\Roxanne\scripts;
pat = patienten;
 
for a = 1:length(t);
    for b = 1:length(f);
        kmponi(b,a) = wjn_pt(koni(:,b,a)-moni(:,b,a));
        kgponi(b,a) = wjn_pt(koni(:,b,a)-goni(:,b,a));
        mgponi(b,a) = wjn_pt(moni(:,b,a)-goni(:,b,a));
    end
end


skgponi=zeros(size(kgponi));
skgponi(kgponi<=.05)=-.5;
skgponi(kgponi<=fdr_bh(kgponi(i_fdr),0.05))=-1;

skmponi=zeros(size(kmponi));
skmponi(kmponi<=.05)=-.5;
skmponi(kmponi<=fdr_bh(kmponi(i_fdr),0.05))=-1;

smgponi=zeros(size(mgponi));
smgponi(mgponi<=.05)=-.5;
smgponi(mgponi<=fdr_bh(mgponi(i_fdr),0.05))=-1;

save('TFplot_pvalues_betweencond_allpat');

 %% Contralateral average images

addpath E:\Roxanne\scripts;
pat = patienten;

cond= {'onset-klein', 'onset-mittel', 'onset-gross'};
imagefolder = 'E:\Roxanne\STN_Rotameter\contralateral_images\';

cd(imagefolder)

for a=1:length(pat);
    if ~isfield (pat{a}, 'bad');
        
       for b=1:length(cond);
       id= pat{a}.id(1:8);

       images = ffind(['*' cond{b} '*' 'sbgmtf_efspmeeg' '*' id '*' 'ON' '*' '.nii']);
       
       if strcmp(pat{a}.id,'r004IR48') %
           images = ffind(['*' cond{b} '*' 'sbgmtf_efspmeeg' '*' id '*' 'R_ON' '*' '.nii']);
       end
       
            if ~isempty(images);
                outimage = (fullfile(cd,'averages', ['average_contralateral_' id '_' cond{b} '_ON.nii']));
                if length(images) == 2; 
                    spm_imcalc(images,outimage,'(i1+i2)/2');
                else
                    spm_imcalc(images,outimage, 'i1/1');
                end
            end
    end
end
end

imagefolder = 'E:\Roxanne\STN_Rotameter\contralateral_images\averages';

cd(imagefolder)

koncfiles = ffind('*average_contralateral_r*_onset*klein*ON*.nii');
moncfiles = ffind('*average_contralateral_r*_onset*mittel*ON*.nii');
goncfiles = ffind('*average_contralateral_r*_onset*gross*ON*.nii');


for a = 1:length(koncfiles)
    [konc(a,:,:),t,f] = wjn_disp_nii(koncfiles{a},[-3 3],[1 100],0);
    [monc(a,:,:),t,f] = wjn_disp_nii(moncfiles{a},[-3 3],[1 100],0);
    [gonc(a,:,:),t,f] = wjn_disp_nii(goncfiles{a},[-3 3],[1 100],0);
end

cd('E:\Roxanne\STN_Rotameter');
save('TF_pvalues_figures_allpat_ONOFF_ipsicontra');
 
%% Contralateral ON permTest, all pat
clear all
close all

cd('E:\Roxanne\STN_Rotameter');
load('TF_pvalues_figures_allpat_ONOFF_ipsicontra');

addpath E:\Roxanne\scripts;
pat = patienten;
 
for a = 1:length(t);
    for b = 1:length(f);
        kmpon(b,a) = wjn_pt(konc(:,b,a)-monc(:,b,a));
        kgpon(b,a) = wjn_pt(konc(:,b,a)-gonc(:,b,a));
        mgpon(b,a) = wjn_pt(monc(:,b,a)-gonc(:,b,a));
        
        kponc(b,a)=wjn_pt(konc(:,b,a));
        mponc(b,a)=wjn_pt(monc(:,b,a));
        gponc(b,a)=wjn_pt(gonc(:,b,a));
        
    end
end

skponc=zeros(size(kponc));
skponc(kponc<=.05)=-.5;
skponc(kponc<=fdr_bh(kponc(i_fdr),0.05))=-1;

smponc=zeros(size(mponc));
smponc(mponc<=.05)=-.5;
smponc(mponc<=fdr_bh(mponc(i_fdr),0.05))=-1;

sgponc=zeros(size(gponc));
sgponc(gponc<=.05)=-.5;
sgponc(gponc<=fdr_bh(gponc(i_fdr),0.05))=-1;

skgpon=zeros(size(kgpon));
skgpon(kgpon<=.05)=-.5;
skgpon(kgpon<=fdr_bh(kgpon(i_fdr),0.05))=-1;

skmpon=zeros(size(kmpon));
skmpon(kmpon<=.05)=-.5;
skmpon(kmpon<=fdr_bh(kmpon(i_fdr),0.05))=-1;

smgpon=zeros(size(mgpon));
smgpon(mgpon<=.05)=-.5;
smgpon(mgpon<=fdr_bh(mgpon(i_fdr),0.05))=-1;

save('TFplot_pvalues_betweencond_allpat');
%% Ipsilateral Variables vs BL

addpath E:\Roxanne\scripts;
pat = patienten;

for a = 1:length(t);
    for b = 1:length(f);
        
        kponi(b,a)=wjn_pt(koni(:,b,a));
        mponi(b,a)=wjn_pt(moni(:,b,a));
        gponi(b,a)=wjn_pt(goni(:,b,a));

    end
end

skponi=zeros(size(kponi));
skponi(kponi<=.05)=-.5;
skponi(kponi<=fdr_bh(kponi(i_fdr),0.05))=-1;

smponi=zeros(size(mponi));
smponi(mponi<=.05)=-.5;
smponi(mponi<=fdr_bh(mponi(i_fdr),0.05))=-1;

sgponi=zeros(size(gponi));
sgponi(gponi<=.05)=-.5;
sgponi(gponi<=fdr_bh(gponi(i_fdr),0.05))=-1;

cd('E:\Roxanne\STN_Rotameter');
save('TFplot_pvalues_betweencond_allpat');
%% variables ipsi vs cl pvalues

addpath E:\Roxanne\scripts;
pat = patienten;


for a = 1:length(t);
    for b = 1:length(f);

kponicl (b,a)= wjn_pt(konc(:,b,a)-koni(:,b,a));
mponicl (b,a)= wjn_pt(monc(:,b,a)-moni(:,b,a));
gponicl (b,a)= wjn_pt(gonc(:,b,a)-goni(:,b,a));

    end
end


for a=1:length(t);
   for b=1:length(f);
    
mcontra(:,b,a) = mean([konc(:,b,a),monc(:,b,a),gonc(:,b,a)],2);
mipsi(:,b,a) = mean([koni(:,b,a),moni(:,b,a),goni(:,b,a)],2);

   end
end

for a=1:length(t);
   for b=1:length(f);
mpcontraipsi(b,a)=wjn_pt(mcontra(:,b,a)-mipsi(:,b,a));
   end
end

smpcontraipsi=zeros(size(mpcontraipsi));
smpcontraipsi(mpcontraipsi<=.05)=-.5;
smpcontraipsi(mpcontraipsi<=fdr_bh(mpcontraipsi(i_fdr),0.05))=-1;

skponicl=zeros(size(kponicl));
skponicl(kponicl<=.05)=-.5;
skponicl(kponicl<=fdr_bh(kponicl(i_fdr),0.05))=-1;

smponicl=zeros(size(mponicl));
smponicl(mponicl<=.05)=-.5;
smponicl(mponicl<=fdr_bh(mponicl(i_fdr),0.05))=-1;

sgponicl=zeros(size(gponicl));
sgponicl(gponicl<=.05)=-.5;
sgponicl(gponicl<=fdr_bh(gponicl(i_fdr),0.05))=-1;

cd('E:\Roxanne\STN_Rotameter');
save('TF_pvalues_figures_allpat_ONOFF_ipsicontra');

%% Figure all ipsi vs contra

[smpcontraipsicluster,smpcontraipsiclusters,smpcontraipsicluster_sizes] = wjn_remove_clusters(smpcontraipsi,150);

cma = colormap([1 1 1; 0 0 0]);

figure,
imagesc(t,f,smpcontraipsicluster)
hold on
title('p<.05, not corrected, all ipsi- vs contralateral')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
xlim([-2 3])
figone (7,11);
colormap(cma);

%% Figure only ON recordings, between cond

clear all
close all

pat = patienten;

[skgponcluster,skgponiclusters,skgponcluster_sizes] = wjn_remove_clusters(skgpon,150);
[skmponcluster,skmponiclusters,skmponcluster_sizes] = wjn_remove_clusters(skmpon,150);
[smgponcluster,smgponiclusters,smgponcluster_sizes] = wjn_remove_clusters(smgpon,150);

cma = colormap([1 1 1; 0 0 0]);

figure,
subplot(1,3,1)
imagesc(t,f,skgponcluster)
hold on
title('small - large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy

subplot(1,3,2)
imagesc(t,f,skmponcluster)
title('small - medium')
axis xy
ylabel('Frequency [Hz]')
xlabel('Time [s]');
k = get(gca,'title');
k.String = {'All ON (n=16),FDR p<.05, Cluster 150',k.String};
set(gca,'title',k);

subplot(1,3,3)
imagesc(t,f,smgponcluster)
axis xy
title('medium - large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');
colormap(cma);

figone(7,17)


%% Ipsi perm Test all conditions

clear all
close all

cd('E:\Roxanne\STN_Rotameter');
load('TF_pvalues_figures_allpat_ONOFF_ipsicontra');

addpath E:\Roxanne\scripts;
pat = patienten;

[skgponicluster,skgponiiclusters,skgponicluster_sizes] = wjn_remove_clusters(skgponi,150);
[skmponicluster,skmponiiclusters,skmponicluster_sizes] = wjn_remove_clusters(skmponi,150);
[smgponicluster,smgponiiclusters,smgponicluster_sizes] = wjn_remove_clusters(smgponi,150);

cma = colormap([1 1 1; 0 0 0]);

figure,
subplot(1,3,1)
imagesc(t,f,skgponicluster)
hold on
title('small - large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy

subplot(1,3,2)
imagesc(t,f,skmponicluster)
title('small - medium')
axis xy
ylabel('Frequency [Hz]')
xlabel('Time [s]');
k = get(gca,'title');
k.String = {'All ON ipsi (n=16),FDR p<.05, Cluster 150',k.String};
set(gca,'title',k);

subplot(1,3,3)
imagesc(t,f,smgponicluster)
axis xy
title('medium - large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');
colormap(cma);

figone(7,17)

%% Figure ipsi & contralateral images
clear all
close all

cd('E:\Roxanne\STN_Rotameter');
load('TF_pvalues_figures_allpat_ONOFF_ipsicontra');

addpath E:\Roxanne\scripts;
pat = patienten;

[skponccluster,skponciclusters,skponccluster_sizes] = wjn_remove_clusters(skponc,150);
[smponccluster,smponciclusters,skmponcluster_sizes] = wjn_remove_clusters(smponc,150);
[sgponccluster,sgponciclusters,sgponccluster_sizes] = wjn_remove_clusters(sgponc,150);

[skponicluster,skponiiclusters,skponicluster_sizes] = wjn_remove_clusters(skponi,150);
[smponicluster,smponiiclusters,skmponicluster_sizes] = wjn_remove_clusters(smponi,150);
[sgponicluster,sgponiiclusters,sgponicluster_sizes] = wjn_remove_clusters(sgponi,150);

cd('E:\Roxanne\STN_Rotameter\contralateral_images\averages');

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
imagesc(tx,kcready,kcread);
title('small')
k = get(gca,'title');
k.String = {'contralateral',k.String};
set(gca,'title',k);
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
caxis([-15 30]);

ax1=subplot(3,4,2)
imagesc(tx,f,skponccluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy

colormap(ax1,cma);
colormap(ax2,'parula');

ax2=subplot(3,4,3)
imagesc(tx,kiready,kiread);
title('small')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
k = get(gca,'title');
k.String = {'ipsilateral',k.String};
set(gca,'title',k);
caxis([-15 30]);

ax1=subplot(3,4,4)
imagesc(kcreadx,f,skponicluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy

colormap(ax1,cma)
colormap(ax2,'parula');

ax2=subplot(3,4,5)
imagesc(kcreadx,mcready,mcread);
title('medium')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
caxis([-15 30]);

ax1=subplot(3,4,6)
imagesc(kcreadx,f,smponccluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy

colormap(ax1,cma)
colormap(ax2,'parula');

ax2=subplot(3,4,7)
imagesc(kcreadx,miready,miread);
title('medium')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
caxis([-15 30]);

ax1=subplot(3,4,8)
imagesc(kcreadx,f,smponicluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy

colormap(ax1,cma)
colormap(ax2,'parula');

ax2=subplot(3,4,9)
imagesc(kcreadx,gcready,gcread);
title('large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
caxis([-15 30]);

ax1=subplot(3,4,10)
imagesc(kcreadx,f,sgponccluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy

colormap(ax1,cma)
colormap(ax2,'parula');

ax2=subplot(3,4,11)
imagesc(kcreadx,giready,giread);
title('large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
caxis([-15 30]);

ax1=subplot(3,4,12)
imagesc(kcreadx,f,sgponicluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy

colormap(ax1,cma)
colormap(ax2,'parula');

suptitle('All ON (n=16), permTest, p<0.05, FDR corrected, Cluster 250');

figone(25,28)

myprint(fullfile('E:\Roxanne\STN_Rotameter\figures', 'TFplot_pvalues_gegBL_permTest_ON_ipsicontra'));

cd('E:\Roxanne\STN_Rotameter');
save('TF_pvalues_figures_allpat_ONOFF_ipsicontra');

%% Figures ipsi vs contralateral

figure,
imagesc(t,f, smpcontraipsi);
title('Contra- vs ipsilateral (n=16), p<.05, FDR corrected, permTest')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
colormap('gray')

myprint(fullfile('E:\Roxanne\STN_Rotameter\figures', 'TFplot_pvalues_permTest_meanipsivscontra_all'));

figure,
subplot(1,3,1)
imagesc(t,f,skponicl)
hold on
title('small')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy

subplot(1,3,2)
imagesc(t,f,smponicl)
title('medium')
axis xy
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
axis xy
title('large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');

colormap('gray')
figone(7,20)

myprint(fullfile('E:\Roxanne\STN_Rotameter\figures', 'TFplot_pvalues_SML_permTest_ipsivscontra_all'));

%% mean gamma in all patients ipsi vs contra

gammafreq=40:90
gammatime=sc(t,0):sc(t,0.5)

mean_onc= [mean(mean(konc(:,gammafreq,gammatime),2),3), mean(mean(monc(:,gammafreq,gammatime),2),3), mean(mean(gonc(:,gammafreq,gammatime),2),3)];
mean_oni= [mean(mean(koni(:,gammafreq,gammatime),2),3), mean(mean(moni(:,gammafreq,gammatime),2),3), mean(mean(goni(:,gammafreq,gammatime),2),3)];

pkmonc=wjn_pt(mean_onc(:,1)-mean_onc(:,2));
pmgonc=wjn_pt(mean_onc(:,2)-mean_onc(:,3));
pkgonc=wjn_pt(mean_onc(:,1)-mean_onc(:,3));

pkmoni=wjn_pt(mean_oni(:,1)-mean_oni(:,2));
pmgoni=wjn_pt(mean_oni(:,2)-mean_oni(:,3));
pkgoni=wjn_pt(mean_oni(:,1)-mean_oni(:,3));

pkon=wjn_pt(mean_onc(:,1)-mean_oni(:,1));
pmon=wjn_pt(mean_onc(:,2)-mean_oni(:,2));
pgon=wjn_pt(mean_onc(:,3)-mean_oni(:,3));

ponic=wjn_pt((mean(mean_oni(:,:),2))-(mean(mean_onc(:,:),2)));


pthresh = fdr_bh([pkmonc pmgonc pkgonc pkmoni pmgoni pkgoni pkon pmon pgon ponic],.05);

cmaps = colorlover(5,0);
cmaponi = [0,0,0;0,0,0;0,0,0]; %([0 0 0],:); %cmaps([1 1 1],:);
cmaponc = [0.8 0.8 0.8; 0.8 0.8 0.8; 0.8 0.8 0.8]; %cmaps([3 3 3],:);

figure
boni=mybar(mean_oni,cmaponi,.8:1:2.8)
hold on
bonc=mybar(mean_onc,cmaponc,1.2:1:3.2)
title('Mean Gamma [40-90,0-0.5]')
set(boni,'barwidth',.4)
set(bonc,'barwidth',.4)
set(gca,'XTick',[1:3],'XTickLabel',{'small','medium','large'})
legend([boni(1) bonc(1)],{'Ipsi','Contra'},'location','northeastoutside')
sigbracket(['P = ' num2str(pkon,1)],[.8 1.2],27,9)
sigbracket(['P = ' num2str(pmon,1)],[1.8 2.2],27,9)
sigbracket(['P = ' num2str(pgon,1)],[2.8 3.2],27,9)
ylabel('Relative power change [%, SEM]');
figone(7,14)

myprint('E:\Roxanne\STN_Rotameter\figures\mean_gamma_ipsicontra_grey_allON16');


cmap = colorlover(6);
cmap = [cmap([2,4,5],:);cmap([2,4,5],:)];

figure
    mybar([mean_oni mean_onc],cmap)
    title('Mean Gamma 40-90, 0-.5')
    set(gca,'XTickLabel',{'small','medium','large','small','medium','large'})
    xlabel('ipsilateral                          contralateral')
    ylabel('Relative power change [%, SEM]');

    sigbracket (['P=' num2str(pkmonc,2)],[4 5],29,9);
    sigbracket(['P=' num2str(pmgonc,2)],[5 6],29,9)
    sigbracket(['P=' num2str(pkgonc,2)],[4 6],34,9);
    
    sigbracket(['P=' num2str(pkgoni,2)],[1 2],29,9);
    sigbracket(['P=' num2str(pmgoni,2)],[2 3],29,9);
    sigbracket(['P=' num2str(pkgoni,2)],[1 3],34,9);
    
    sigbracket(['P=' num2str(ponic,2)],[2 5],39,9);
    
% sigbracket('P<0.01',[1 2],29,9)
% sigbracket('P=0.2',[2 3],29,9)
% sigbracket('P<0.01',[1 3],34,9);
% sigbracket('P<0.001',[4 5],29,9);
% sigbracket('P<0.05',[5 6],29,9);
% sigbracket('P<0.001',[4 6],34,9);
% sigbracket('P<0.05',[2 5],39,9);

figone(7,12)


%% mean beta in all patients ipsi vs contra

mean_onc_beta= [mean(mean(konc(:,13:30,sc(t,0.5):sc(t,1)),2),3), mean(mean(monc(:,13:30,sc(t,0.5):sc(t,1)),2),3), mean(mean(gonc(:,13:30,sc(t,0.5):sc(t,1)),2),3)];
mean_oni_beta= [mean(mean(koni(:,13:30,sc(t,0.5):sc(t,1)),2),3), mean(mean(moni(:,13:30,sc(t,0.5):sc(t,1)),2),3), mean(mean(goni(:,13:30,sc(t,0.5):sc(t,1)),2),3)];

pkmonc=wjn_pt(mean_onc_beta(:,1)-mean_onc_beta(:,2));
pmgonc=wjn_pt(mean_onc_beta(:,2)-mean_onc_beta(:,3));
pkgonc=wjn_pt(mean_onc_beta(:,1)-mean_onc_beta(:,3));

pkmoni=wjn_pt(mean_oni_beta(:,1)-mean_oni_beta(:,2));
pmgoni=wjn_pt(mean_oni_beta(:,2)-mean_oni_beta(:,3));
pkgoni=wjn_pt(mean_oni_beta(:,1)-mean_oni_beta(:,3));

pkon=wjn_pt(mean_onc_beta(:,1)-mean_oni_beta(:,1));
pmon=wjn_pt(mean_onc_beta(:,2)-mean_oni_beta(:,2));
pgon=wjn_pt(mean_onc_beta(:,3)-mean_oni_beta(:,3));

ponic=wjn_pt((mean(mean_oni_beta(:,:),2))-(mean(mean_onc_beta(:,:),2)));

% cmaps = colorlover(5,0);
% cmaponi = cmaps([1 1 1],:);
% cmaponc = cmaps([3 3 3],:);

cmaponi = [0,0,0;0,0,0;0,0,0]; %([0 0 0],:); %cmaps([1 1 1],:);
cmaponc = [0.8 0.8 0.8; 0.8 0.8 0.8; 0.8 0.8 0.8]; %cmaps([3 3 3],:);


figure
boni=mybar(mean_oni_beta,cmaponi,.8:1:2.8)
hold on
bonc=mybar(mean_onc_beta,cmaponc,1.2:1:3.2)
title('Mean Beta (13-30 Hz, .5-1)')
ylabel('Mean power change [%, SEM]');
set(boni,'barwidth',.4)
set(bonc,'barwidth',.4)
set(gca,'YDir','reverse')
set(gca,'XTick',[1:3],'XTickLabel',{'small','medium','large'})
legend([boni(1) bonc(1)],{'Ipsi','Contra'},'location','northeastoutside')
sigbracket(['P = ' num2str(pkon,1)],[.8 1.2],-19,9)
sigbracket(['P = ' num2str(pmon,1)],[1.8 2.2],-19,9)
sigbracket(['P = ' num2str(pgon,1)],[2.8 3.2],-19,9)
ylim([-23 0])
figone(7,14)

myprint('E:\Roxanne\STN_Rotameter\figures\mean_beta_ipsi_vs_contra_permTest_grey');

cmap=colorlover(6);
cmap = [cmap([2,4,5],:);cmap([2,4,5],:)];

figure
    mybar([mean_oni_beta mean_onc_beta],cmap)
    title('Mean Beta (13-30 Hz, .5-1)')
    set(gca,'XTickLabel',{'small','medium','large','small','medium','large'})
    xlabel('ipsilateral                          contralateral')
    ylabel('Mean power change [%,SEM]');
    set(gca,'YDir','reverse')
    
sigbracket(['P=' num2str(pkmoni,1)],[1 2],-18,9)
sigbracket(['P=' num2str(pmgoni,1)],[2 3],-18,9)
sigbracket(['P=' num2str(pkgoni,1)],[1 3],-21,9);
sigbracket(['P=' num2str(pkmonc,1)],[4 5],-18,9)
sigbracket(['P=' num2str(pmgonc,1)],[5 6],-18,9)
sigbracket(['P=' num2str(pkgonc,1)],[4 6],-21,9);
sigbracket(['P=' num2str(ponic,1)],[2 5],-24.5,9);

ylim([-30 0]);

figone(7,12)


%% mean theta in all patients ipsi vs contra

thetafreq=2:8
gammatime=sc(t,0):sc(t,0.5)

mean_onctheta= [mean(mean(konc(:,thetafreq,gammatime),2),3), mean(mean(monc(:,thetafreq,gammatime),2),3), mean(mean(gonc(:,thetafreq,gammatime),2),3)];
mean_onitheta= [mean(mean(koni(:,thetafreq,gammatime),2),3), mean(mean(moni(:,thetafreq,gammatime),2),3), mean(mean(goni(:,thetafreq,gammatime),2),3)];

pkmonc=wjn_pt(mean_onctheta(:,1)-mean_onctheta(:,2));
pmgonc=wjn_pt(mean_onctheta(:,2)-mean_onctheta(:,3));
pkgonc=wjn_pt(mean_onctheta(:,1)-mean_onctheta(:,3));

pkmoni=wjn_pt(mean_onitheta(:,1)-mean_onitheta(:,2));
pmgoni=wjn_pt(mean_onitheta(:,2)-mean_onitheta(:,3));
pkgoni=wjn_pt(mean_onitheta(:,1)-mean_onitheta(:,3));

pkon=wjn_pt(mean_onctheta(:,1)-mean_onitheta(:,1));
pmon=wjn_pt(mean_onctheta(:,2)-mean_onitheta(:,2));
pgon=wjn_pt(mean_onctheta(:,3)-mean_onitheta(:,3));

ponic=wjn_pt((mean(mean_onitheta(:,:),2))-(mean(mean_onctheta(:,:),2)));


pthresh = fdr_bh([pkmonc pmgonc pkgonc pkmoni pmgoni pkgoni pkon pmon pgon ponic],.05);

cmaps = colorlover(5,0);
cmaponi = [0,0,0;0,0,0;0,0,0]; %([0 0 0],:); %cmaps([1 1 1],:);
cmaponc = [0.8 0.8 0.8; 0.8 0.8 0.8; 0.8 0.8 0.8]; %cmaps([3 3 3],:);

figure
boni=mybar(mean_onitheta,cmaponi,.8:1:2.8)
hold on
bonc=mybar(mean_onctheta,cmaponc,1.2:1:3.2)
title('Mean Theta [2-8,0-0.5]')
set(boni,'barwidth',.4)
set(bonc,'barwidth',.4)
set(gca,'XTick',[1:3],'XTickLabel',{'small','medium','large'})
legend([boni(1) bonc(1)],{'Ipsi','Contra'},'location','northeastoutside')
sigbracket(['P = ' num2str(pkon,1)],[.8 1.2],52,9)
sigbracket(['P = ' num2str(pmon,1)],[1.8 2.2],52,9)
sigbracket(['P = ' num2str(pgon,1)],[2.8 3.2],52,9)
ylabel('Relative power change [%, SEM]');
figone(7,14)

myprint('E:\Roxanne\STN_Rotameter\figures\mean_theta_ipsicontra_grey_allON16');


cmap = colorlover(6);
cmap = [cmap([2,4,5],:);cmap([2,4,5],:)];

figure
    mybar([mean_onitheta mean_onctheta],cmap)
    title('Mean Theta 2-8, 0-.5')
    set(gca,'XTickLabel',{'small','medium','large','small','medium','large'})
    xlabel('ipsilateral                          contralateral')
    ylabel('Relative power change [%, SEM]');

    sigbracket (['P=' num2str(pkmonc,2)],[4 5],55,9);
    sigbracket(['P=' num2str(pmgonc,2)],[5 6],55,9)
    sigbracket(['P=' num2str(pkgonc,2)],[4 6],65,9);
    
    sigbracket(['P=' num2str(pkgoni,2)],[1 2],55,9);
    sigbracket(['P=' num2str(pmgoni,2)],[2 3],55,9);
    sigbracket(['P=' num2str(pkgoni,2)],[1 3],65,9);
    
    sigbracket(['P=' num2str(ponic,2)],[2 5],75,9);
    
figone(7,12)
