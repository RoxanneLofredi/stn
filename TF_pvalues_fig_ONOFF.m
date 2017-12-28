clear all
close all

cd(root);

pat = patienten;

cond= {'onset-klein', 'onset-mittel', 'onset-gross'};
med={'ON','OFF'};


for b=1:length(cond);
    for c=1:length(med)

       images = ffind(['*' cond{b} '*' med{c} '*' '.nii']);
       outimage=['contralateral_' cond{b} '_' med{c} '.nii'];
       spm_imcalc(images,outimage,'(i1+i2+i3+i4+i5+i6+i7)/7');
       
       [data,x,y]=wjn_disp_nii(['contralateral_' cond{b} '_' med{c} '.nii'],[-3 3],[1 100],1);
       imagesc(x,y,data),axis xy
       caxis ([-15 30]); colorbar; title ([cond{b} '_ ' med{c}]);
        
        myprint([cond{b} '_' med{c}]);
        
     end
end

%%
[kON,x,y]=wjn_disp_nii(['contralateral_onset-klein_ON.nii'],[-3 3],[1 100],0);
[mON,x,y]=wjn_disp_nii(['contralateral_onset-mittel_ON.nii'],[-3 3],[1 100],0);      
[gON,x,y]=wjn_disp_nii(['contralateral_onset-gross_ON.nii'],[-3 3],[1 100],0);

[kOFF,x,y]=wjn_disp_nii(['contralateral_onset-klein_OFF.nii'],[-3 3],[1 100],0);
[mOFF,x,y]=wjn_disp_nii(['contralateral_onset-mittel_OFF.nii'],[-3 3],[1 100],0);      
[gOFF,x,y]=wjn_disp_nii(['contralateral_onset-gross_OFF.nii'],[-3 3],[1 100],0);


%% Figures
load ('TF-pvalues_betweencond_ONOFF.mat');

figure,
subplot(2,3,1)
imagesc(t,f,skgpon)
hold on
title('small - large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
k = get(gca,'ylabel');
k.String = {'ON',k.String};
set(gca,'ylabel',k);
subplot(2,3,2)
imagesc(t,f,skmpon)
title('small - medium')
axis xy
ylabel('Frequency [Hz]')
xlabel('Time [s]');
k = get(gca,'title');
k.String = {'ONOFF (n=7),FDR p<.05',k.String};
set(gca,'title',k);
subplot(2,3,3)
imagesc(t,f,smgpon)
axis xy
title('medium - large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');
subplot(2,3,4)
imagesc(t,f,skgpoff)
hold on
title('small - large')
ylabel('Frequency [Hz]')
k = get(gca,'ylabel');
k.String = {'OFF',k.String};
set(gca,'ylabel',k);
xlabel('Time [s]');axis xy
subplot(2,3,5)
imagesc(t,f,skmpoff)
title('small - medium')
axis xy
ylabel('Frequency [Hz]')
xlabel('Time [s]');
subplot(2,3,6)
imagesc(t,f,smgpoff)
axis xy
title('medium - large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');
colormap('gray')
figone (20,35)

%% Figures ON vs OFF
 
figure,
subplot(3,3,1)
imagesc(x,y,kON);
title('small')
ylabel('Frequency [Hz]')
k = get(gca,'ylabel');
k.String = {'ON',k.String};
set(gca,'ylabel',k);
xlabel('Time [s]');axis xy
caxis([-15 30]);
hold on

subplot(3,3,2)
imagesc(x,y,mON);
title('medium')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
caxis([-15 30]);
hold on

subplot(3,3,3)
imagesc(x,y,gON);
title('large')
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
caxis([-15 30]);
hold on

subplot(3,3,4)
imagesc(x,y,kOFF);
ylabel('Frequency [Hz]')
k = get(gca,'ylabel');
k.String = {'OFF',k.String};
set(gca,'ylabel',k);
xlabel('Time [s]');axis xy
caxis([-15 30]);
hold on

subplot(3,3,5)
imagesc(x,y,mOFF);
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
caxis([-15 30]);
hold on

subplot(3,3,6)
imagesc(x,y,gOFF);
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
caxis([-15 30]);
hold on

ax1=subplot(3,3,7)
imagesc(t,f,skponoff)
ylabel('Frequency [Hz]')
xlabel('Time [s]');axis xy
k = get(gca,'ylabel');
k.String = {'ON vs OFF',k.String};
set(gca,'ylabel',k);
colormap(ax1, 'gray')
hold on

ax1=subplot(3,3,8)
imagesc(t,f,smponoff)
axis xy
ylabel('Frequency [Hz]')
xlabel('Time [s]');
colormap(ax1, 'gray')
hold on

ax1=subplot(3,3,9)
imagesc(t,f,sgponoff)
axis xy
ylabel('Frequency [Hz]')
xlabel('Time [s]');

colormap(ax1, 'gray')
figone(25,25)



save('TF_pvalues_fig_ONOFF')