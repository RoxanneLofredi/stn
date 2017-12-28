%% Ipsilateral average images

close all
clear all

pat = patienten;

cond= {'onset-klein', 'onset-mittel', 'onset-gross'};
imagefolder = ;

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
                
%                 wjn_disp_nii(['average' id '_' cond{a} '.nii'],[-3 3],[1 100],1);
%                 caxis ([-15 30]); colorbar; title ([id ' allmed contralateral onset']);
%         
%                 myprint(fullfile(cd,'averaged_TF_plots',['averagesidecondmed_' id]));
            end
    end
end
end


cd(imagefolder)

 for b=1:length(cond);
%     if ~isfield (pat{a}, 'bad');
%        id= pat{a}.id(1:8); 
       images = ffind(['average_ipsilateral_r' '*' cond{b} '*' 'ON.nii']);
            if ~isempty(images);
                outimage = ['average_ipsilateral' cond{b} 'ON.nii'];
                
                    spm_imcalc(images,outimage,'(i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12+i13+i14+i15+i16)/16');

%                 wjn_disp_nii(['average' id '_' cond{a} '.nii'],[-3 3],[1 100],1);
%                 caxis ([-15 30]); colorbar; title ([id ' allmed contralateral onset']);
%         
%                 myprint(fullfile(cd,'averaged_TF_plots',['averagesidecondmed_' id]));
            end
%     end 
 end
 
 %% Contralateral average images
 
close all
clear all

cond= {'onset-klein', 'onset-mittel', 'onset-gross'};

cd(imagefolder)

for a=1:length(pat);
    if ~isfield (pat{a}, 'bad');
        
       for b=1:length(cond);
       id= pat{a}.id(1:8);

       images = ffind(['*' cond{b} '*' 'sbgmtf_efspmeeg' '*' id '*' 'ON' '*' '.nii']);
       
       if strcmp(pat{a}.id,'r004IR48')
           images = ffind(['*' cond{b} '*' 'sbgmtf_efspmeeg' '*' id '*' 'R_ON' '*' '.nii']);
       end
       
            if ~isempty(images);
                outimage = (fullfile(cd,'averages', ['average_contralateral_' id '_' cond{b} '_ON.nii']));
                if length(images) == 2; 
                    spm_imcalc(images,outimage,'(i1+i2)/2');
                else
                    spm_imcalc(images,outimage, 'i1/1');
                end
                
%                 wjn_disp_nii(['average' id '_' cond{a} '.nii'],[-3 3],[1 100],1);
%                 caxis ([-15 30]); colorbar; title ([id ' allmed contralateral onset']);
%         
%                 myprint(fullfile(cd,'averaged_TF_plots',['averagesidecondmed_' id]));
            end
    end
end
end


% Für OFF Files

for a=1:length(pat);
    if ~isfield (pat{a}, 'bad') && ~pat{a}.med;
        
       for b=1:length(cond);
       id= pat{a}.id(1:8);

       images = ffind(['*' cond{b} '*' 'sbgmtf_efspmeeg' '*' id '*' 'OFF' '*' '.nii']);
       
       if strcmp(pat{a}.id,'r004IR48')
           images = ffind(['*' cond{b} '*' 'sbgmtf_efspmeeg' '*' id '*' 'R_OFF' '*' '.nii']);
       end
       
            if ~isempty(images);
                outimage = (fullfile(cd,'averages', ['average_contralateral_' id '_' cond{b} '_OFF.nii']));
                if length(images) == 2; 
                    spm_imcalc(images,outimage,'(i1+i2)/2');
                else
                    spm_imcalc(images,outimage, 'i1/1');
                end
                
%                 wjn_disp_nii(['average' id '_' cond{a} '.nii'],[-3 3],[1 100],1);
%                 caxis ([-15 30]); colorbar; title ([id ' allmed contralateral onset']);
%         
%                 myprint(fullfile(cd,'averaged_TF_plots',['averagesidecondmed_' id]));
            end
    end
end
end


cd(imagefolder)

 for b=1:length(cond);
%     if ~isfield (pat{a}, 'bad');
%        id= pat{a}.id(1:8); 
       images = ffind(['average_contralateral_r' '*' cond{b} '*' 'ON.nii']);
            if ~isempty(images);
                outimage = ['average_contralateral' cond{b} 'ON.nii'];
                
                    spm_imcalc(images,outimage,'(i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12+i13+i14+i15+i16)/16');

%                 wjn_disp_nii(['average' id '_' cond{a} '.nii'],[-3 3],[1 100],1);
%                 caxis ([-15 30]); colorbar; title ([id ' allmed contralateral onset']);
%         
%                 myprint(fullfile(cd,'averaged_TF_plots',['averagesidecondmed_' id]));
            end
%     end 
 end
 
 % Für OFF Files
 
  for b=1:length(cond);
%     if ~isfield (pat{a}, 'bad');
%        id= pat{a}.id(1:8); 
       images = ffind(['average_contralateral_r' '*' cond{b} '*' 'OFF.nii']);
            if ~isempty(images);
                outimage = ['average_contralateral' cond{b} 'OFF.nii'];
                
                    spm_imcalc(images,outimage,'(i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12)/12');

%                 wjn_disp_nii(['average' id '_' cond{a} '.nii'],[-3 3],[1 100],1);
%                 caxis ([-15 30]); colorbar; title ([id ' allmed contralateral onset']);
%         
%                 myprint(fullfile(cd,'averaged_TF_plots',['averagesidecondmed_' id]));
            end
%     end 
  end


 
%% Contralateral ON and OFF permTest, all pat
clear all
close all

root = ;

cd(root);

konfiles = ffind('*average_contralateral_r*_onset*klein*ON*.nii');
monfiles = ffind('*average_contralateral_r*_onset*mittel*ON*.nii');
gonfiles = ffind('*average_contralateral_r*_onset*gross*ON*.nii');

% konfiles([1;3;4;6;9;10;11;15;16])=[];
% monfiles([1;3;4;6;9;10;11;15;16])=[];
% gonfiles([1;3;4;6;9;10;11;15;16])=[];

kofffiles = ffind('*average_contralateral_r*_onset*klein*OFF*.nii');
mofffiles = ffind('*average_contralateral_r*_onset*mittel*OFF*.nii');
gofffiles = ffind('*average_contralateral_r*_onset*gross*OFF*.nii');
% 
% kofffiles([1;4;7;11;12])=[];
% mofffiles([1;4;7;11;12])=[];
% gofffiles([1;4;7;11;12])=[];
seven = {'r032DA55'    'r035EE53'    'r091NF42'    'r094RE60'    'r011AI40'    'r041EC34'  'r093AL56'};

offi=ci(seven,kofffiles)
oni=ci(seven,konfiles)

konfiles = konfiles(oni);
monfiles = monfiles(oni);
gonfiles = gonfiles(oni);

kofffiles = kofffiles(offi);
mofffiles = mofffiles(offi);
gofffiles = gofffiles(offi);

for a = 1:length(konfiles)
    [kon(a,:,:),t,f] = wjn_disp_nii(konfiles{a},[-3 3],[1 100],0);
    [mon(a,:,:),t,f] = wjn_disp_nii(monfiles{a},[-3 3],[1 100],0);
    [gon(a,:,:),t,f] = wjn_disp_nii(gonfiles{a},[-3 3],[1 100],0);
end

for a = 1:length(t);
    for b = 1:length(f);
        kmpon(b,a) = wjn_pt(kon(:,b,a)-mon(:,b,a));
        kgpon(b,a) = wjn_pt(kon(:,b,a)-gon(:,b,a));
        mgpon(b,a) = wjn_pt(mon(:,b,a)-gon(:,b,a));
        
        kpon(b,a)=wjn_pt(kon(:,b,a));
        mpon(b,a)=wjn_pt(mon(:,b,a));
        gpon(b,a)=wjn_pt(gon(:,b,a));
        
%         kpon(b,a) = signrank(kon(:,b,a));
%         mpon(b,a) = signrank(mon(:,b,a));
%         gpon(b,a) = signrank(gon(:,b,a));
%         kmpon(b,a) = signrank(kon(:,b,a),mon(:,b,a));
%         kgpon(b,a) = signrank(kon(:,b,a),gon(:,b,a));
%         mgpon(b,a) = signrank(mon(:,b,a),gon(:,b,a));
    end
end

for a = 1:length(kofffiles)
    [koff(a,:,:),t,f] = wjn_disp_nii(kofffiles{a},[-3 3],[1 100],0);
    [moff(a,:,:),t,f] = wjn_disp_nii(mofffiles{a},[-3 3],[1 100],0);
    [goff(a,:,:),t,f] = wjn_disp_nii(gofffiles{a},[-3 3],[1 100],0);
end

for a = 1:length(t);
    for b = 1:length(f);
        
        kmpoff(b,a) = wjn_pt(koff(:,b,a)-moff(:,b,a));
        kgpoff(b,a) = wjn_pt(koff(:,b,a)-goff(:,b,a));
        mgpoff(b,a) = wjn_pt(moff(:,b,a)-goff(:,b,a));
        
        kpoff(b,a)=wjn_pt(koff(:,b,a));
        mpoff(b,a)=wjn_pt(moff(:,b,a));
        gpoff(b,a)=wjn_pt(goff(:,b,a));

%         kpoff(b,a) = signrank(koff(:,b,a));
%         mpoff(b,a) = signrank(moff(:,b,a));
%         gpoff(b,a) = signrank(goff(:,b,a));
%         kmpoff(b,a) = signrank(koff(:,b,a),moff(:,b,a));
%         kgpoff(b,a) = signrank(koff(:,b,a),goff(:,b,a));
%         mgpoff(b,a) = signrank(moff(:,b,a),goff(:,b,a));
    end
end

skpon=zeros(size(kpon));
skpon(kpon<=.01)=-.5;
% skpon(kpon<=fdr_bh(kpon,0.05))=-1;

smpon=zeros(size(mpon));
smpon(mpon<=.01)=-.5;
% smpon(mpon<=fdr_bh(mpon,0.05))=-1;

sgpon=zeros(size(gpon));
sgpon(gpon<=.01)=-.5;
% sgpon(gpon<=fdr_bh(gpon,0.05))=-1;

skgpon=zeros(size(kgpon));
skgpon(kgpon<=.01)=-.5;
% skgpon(kgpon<=fdr_bh(kgpon,0.05))=-1;

skmpon=zeros(size(kmpon));
skmpon(kmpon<=.01)=-.5;
% skmpon(kmpon<=fdr_bh(kmpon,0.05))=-1;

smgpon=zeros(size(mgpon));
smgpon(mgpon<=.01)=-.5;
% smgpon(mgpon<=fdr_bh(mgpon,0.05))=-1;

skpoff=zeros(size(kpoff));
skpoff(kpoff<=.01)=-.5;
% skpoff(kpoff<=fdr_bh(kpoff,0.05))=-1;

smpoff=zeros(size(mpoff));
smpoff(mpoff<=.01)=-.5;
% smpoff(mpoff<=fdr_bh(mpoff,0.05))=-1;

sgpoff=zeros(size(gpoff));
sgpoff(gpoff<=.01)=-.5;
% sgpoff(gpoff<=fdr_bh(gpoff,0.05))=-1;

skgpoff=zeros(size(kgpoff));
skgpoff(kgpoff<=.01)=-.5;
% skgpoff(kgpoff<=fdr_bh(kgpoff,0.05))=-1;

skmpoff=zeros(size(kmpoff));
skmpoff(kmpoff<=.01)=-.5;
% skmpoff(kmpoff<=fdr_bh(kmpoff,0.05))=-1;

smgpoff=zeros(size(mgpoff));
smgpoff(mgpoff<=.01)=-.5;
% smgpoff(mgpoff<=fdr_bh(mgpoff,0.05))=-1;



% skgpoff=zeros(size(kgpoff));
% skgpoff(kgpoff<=.001)=-.5;
% skgpoff(kgpoff<=fdr_bh(kgpoff,0.001))=-1;

%% Ipsilateral Variables

cd(rooti);

konifiles = ffind('*average_ipsilateral*_onset*klein*ON*.nii');
monifiles = ffind('*average_ipsilateral*_onset*mittel*ON*.nii');
gonifiles = ffind('*average_ipsilateral*_onset*gross*ON*.nii');

for a = 1:length(konifiles)
    [koni(a,:,:),t,f] = wjn_disp_nii(konifiles{a},[-3 3],[1 100],0);
    [moni(a,:,:),t,f] = wjn_disp_nii(monifiles{a},[-3 3],[1 100],0);
    [goni(a,:,:),t,f] = wjn_disp_nii(gonifiles{a},[-3 3],[1 100],0);
end

for a = 1:length(t);
    for b = 1:length(f);
        
        kponi(b,a)=wjn_pt(koni(:,b,a));
        mponi=wjn_pt(moni(:,b,a));
        gponi=wjn_pt(goni(:,b,a));

    end
end

skponi=zeros(size(kponi));
skponi(kponi<=.001)=-.5;
skponi(kponi<=fdr_bh(kponi,0.05))=-1;

smponi=zeros(size(mponi));
smponi(mponi<=.001)=-.5;
smponi(mponi<=fdr_bh(mponi,0.05))=-1;

sgponi=zeros(size(gponi));
sgponi(gponi<=.001)=-.5;
sgponi(gponi<=fdr_bh(gponi,0.05))=-1;


%% variables ipsi vs cl pvalues

t=[-3 3];
f=[1 100];

for a = 1:length(t);
    for b = 1:length(f);
        
% kponoff (b,a)= signrank(koff(:,b,a),kon(:,b,a));
% mponoff (b,a)= signrank(moff(:,b,a), mon(:,b,a));
% gponoff (b,a)= signrank(goff(:,b,a),gon(:,b,a));

kponicl (b,a)= wjn_pt(kon(:,b,a)-koni(:,b,a));
mponicl (b,a)= wjn_pt(mon(:,b,a)-moni(:,b,a));
gponicl (b,a)= wjn_pt(gon(:,b,a)-goni(:,b,a));

    end
end

skponicl=zeros(size(kponicl));
skponicl(kponicl<=.05)=-.5;
skponicl(kponicl<=fdr_bh(kponicl,0.05))=-1;

smponicl=zeros(size(mponicl));
smponicl(mponicl<=.05)=-.5;
smponicl(mponicl<=fdr_bh(mponicl,0.05))=-1;

sgponicl=zeros(size(gponicl));
sgponicl(gponicl<=.05)=-.5;
sgponicl(gponicl<=fdr_bh(gponicl,0.05))=-1;

%% Figures ON and OFF, between conditions, contralateral, all pat

figure,
subplot(2,3,1)
imagesc(t,f,skgpon)
hold on
title('klein - gross')
ylabel('ON')
axis xy
subplot(2,3,2)
imagesc(t,f,skmpon)
title('klein - mittel')
axis xy
subplot(2,3,3)
imagesc(t,f,smgpon)
axis xy
title('mittel - gross')
colormap('gray')
subplot(2,3,4)
imagesc(t,f,skgpoff)
ylabel('OFF')
hold on
axis xy
subplot(2,3,5)
imagesc(t,f,skmpoff)
axis xy
subplot(2,3,6)
imagesc(t,f,smgpoff)
axis xy
colormap('gray')

%% Figure only ON recordings, between cond
t=[-3 3];
f=[1 100];

figure,
subplot(1,3,1)
imagesc(t,f,skgpon)
hold on
title('small - large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
subplot(1,3,2)
imagesc(t,f,-1*(kmpon<=.001))
title('small - medium')
axis xy
ylabel('Frequency [Hz]')
xlabel('Time [s]');
k = get(gca,'title');
k.String = {'All(n=17),ON,p<.001',k.String};
set(gca,'title',k);
subplot(1,3,3)
imagesc(t,f,-1*(mgpon<=.001))
axis xy
title('medium - large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');
colormap('gray')
figone(12,22)
myprint(fullfile('E:\Roxanne\STN_Rotameter\figures', 'TFplot_pvalues_SML_permTest_ON'));

%% Figure ON vs OFF contralateral

% imagefolderi = 'E:\Roxanne\STN_Rotameter\ipsilateral_images\averages';
% imagefoldercl = 'E:\Roxanne\STN_Rotameter\contralateral_images\averages';

[skponcluster,skponiclusters,skponccluster_sizes] = wjn_remove_clusters(skpon,150);
[smponcluster,smponiclusters,skmponcluster_sizes] = wjn_remove_clusters(smpon,150);
[sgponcluster,sgponiclusters,sgponccluster_sizes] = wjn_remove_clusters(sgpon,150);

[skpoffcluster,skpofficlusters,skpoffcluster_sizes] = wjn_remove_clusters(skpoff,150);
[smpoffcluster,smpofficlusters,skmpoffcluster_sizes] = wjn_remove_clusters(smpoff,150);
[sgpoffcluster,sgpofficlusters,sgpoffcluster_sizes] = wjn_remove_clusters(sgpoff,150);

cma = colormap([1 1 1; 0 0 0]);

figure,

ax2=subplot(3,4,1)
imagesc(t,f, squeeze(nanmean(kon(:,:,:),1)));
title('small')
k = get(gca,'title');
k.String = {'ON',k.String};
set(gca,'title',k);
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
caxis([-15 30]);
xlim([-2 3])

ax1=subplot(3,4,2)
imagesc(t,f,skponcluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
xlim([-2 3])

colormap(ax1,cma);
colormap(ax2,'parula');

ax2=subplot(3,4,3)
imagesc(t,f, squeeze(nanmean(koff(:,:,:),1)));
title('small')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
k = get(gca,'title');
k.String = {'OFF',k.String};
set(gca,'title',k);
caxis([-15 30]);
xlim([-2 3])

ax1=subplot(3,4,4)
imagesc(t,f,skpoffcluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
xlim([-2 3])

colormap(ax1,cma)
colormap(ax2,'parula');

ax2=subplot(3,4,5)
imagesc(t,f, squeeze(nanmean(mon(:,:,:),1)));
title('medium')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
caxis([-15 30]);
xlim([-2 3])

ax1=subplot(3,4,6)
imagesc(t,f,smponcluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
xlim([-2 3])

colormap(ax1,cma)
colormap(ax2,'parula');

ax2=subplot(3,4,7)
imagesc(t,f, squeeze(nanmean(moff(:,:,:),1)));
title('medium')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
caxis([-15 30]);
xlim([-2 3])

ax1=subplot(3,4,8)
imagesc(t,f,smpoffcluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
xlim([-2 3])

colormap(ax1,cma)
colormap(ax2,'parula');

ax2=subplot(3,4,9)
imagesc(t,f, squeeze(nanmean(gon(:,:,:),1)));
title('large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
caxis([-15 30]);
xlim([-2 3])

ax1=subplot(3,4,10)
imagesc(t,f,sgponcluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
xlim([-2 3])

colormap(ax1,cma)
colormap(ax2,'parula');

ax2=subplot(3,4,11)
imagesc(t,f, squeeze(nanmean(goff(:,:,:),1)));
title('large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
caxis([-15 30]);
xlim([-2 3])

ax1=subplot(3,4,12)
imagesc(t,f,sgpoffcluster)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
xlim([-2 3])

colormap(ax1,cma)
colormap(ax2,'parula');

suptitle('ON and OFF (n=17), permTest, p<0.05, FDR corrected, Cluster 150');

figone(25,28)


%% Figure ipsi & contralateral images

imagefolderi = '';
imagefoldercl = '';

figure,

subplot(2,6,1)
wjn_disp_nii(fullfile(cd(imagefoldercl),'average_contralateralonset-kleinON.nii'),[-3 3],[1 100],1); caxis ([-15 30]); 
title('contralateral small')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
colorbar;
hold on

subplot(2,6,2)
imagesc(t,f,skpon)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
hold on

subplot(2,6,3)
wjn_disp_nii(fullfile(cd(imagefoldercl),'average_contralateralonset-mittelON.nii'),[-3 3],[1 100],1); caxis ([-15 30]); 
title('contralateral medium')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
hold on

subplot(2,6,4)
imagesc(t,f,smpon)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
hold on

subplot(2,6,5)
wjn_disp_nii(fullfile(cd(imagefoldercl),'average_contralateralonset-grossON.nii'),[-3 3],[1 100],1); caxis ([-15 30]); 
title('contralateral medium')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
hold on

subplot(2,6,6)
imagesc(t,f,sgpon)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
colorbar;
hold on

subplot(2,6,7)
wjn_disp_nii(fullfile(cd(imagefolderi),'average_ipsilateralonset-kleinON.nii'),[-3 3],[1 100],1); caxis ([-15 30]); 
title('ipsilateral small')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
colorbar;
hold on

subplot(2,6,8)
imagesc(t,f,skponi)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
hold on

subplot(2,6,9)
wjn_disp_nii(fullfile(cd(imagefolderi),'average_ipsilateralonset-mittelON.nii'),[-3 3],[1 100],1); caxis ([-15 30]); 
title('ipsilateral medium')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
hold on

subplot(2,6,10)
imagesc(t,f,smponi)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
hold on

subplot(2,6,11)
wjn_disp_nii(fullfile(cd(imagefolderi),'average_ipsilateralonset-grossON.nii'),[-3 3],[1 100],1); caxis ([-15 30]); 
title('ipsilateral large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
hold on

subplot(2,6,12)
imagesc(t,f,sgponi)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy

hold on

colormap('gray')
figone(12,22)


%% Figures ipsi vs contralateral

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
figone (10,20)

%%
save('TFplot_pvalues_betweencond_allpat');