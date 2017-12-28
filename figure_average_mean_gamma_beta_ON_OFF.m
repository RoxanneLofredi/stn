clear all
close all

pat = patienten;

timerange = [-0.5 1];
freqrange = [13 35; 40 90]
franges = {'beta','gamma'};
med={'ON','OFF'};
conds = {'onset-klein'; 'onset-mittel'; 'onset-gross'};

gammafreq= 40:90
betafreq= 13:30

n=0

%% IMAGE TO TF PLOT

for a=1:length(pat);
    
    if isfield(pat{a}, 'dupdrs');
        n=n+1;
        id= stringsplit(pat{a}.id,'_');
        id=id{1};
        
        for b=1:length(med);
        for c=1:length(conds);

        filename=['average_contralateral_' id '_' conds{c} '_' med{b}];
        
       [tfdata,x] = wjn_disp_nii(fullfile(root,[filename '.nii']),[-3 3],[1 100],0);
    
        dgamma(n,b,c,:,:)=tfdata(gammafreq,:);
        dbeta(n,b,c,:,:)=tfdata(betafreq,:);
        end
        end
    end
end

%% MEAN GAMMA

gamma_mean_on=squeeze(squeeze(squeeze(squeeze(mean(mean(mean(dgamma(:,1,:,:,:),1),3),4)))));
gamma_mean_off=squeeze(squeeze(squeeze(squeeze(mean(mean(mean(dgamma(:,2,:,:,:),1),3),4)))));

gamma_small_on=squeeze(squeeze(squeeze(squeeze(mean(mean(dgamma(:,1,1,:,:),1),4)))));
gamma_med_on=squeeze(squeeze(squeeze(squeeze(mean(mean(dgamma(:,1,2,:,:),1),4)))));
gamma_large_on=squeeze(squeeze(squeeze(squeeze(mean(mean(dgamma(:,1,3,:,:),1),4)))));

gamma_small_off=squeeze(squeeze(squeeze(squeeze(mean(mean(dgamma(:,2,1,:,:),1),4)))));
gamma_med_off=squeeze(squeeze(squeeze(squeeze(mean(mean(dgamma(:,2,2,:,:),1),4)))));
gamma_large_off=squeeze(squeeze(squeeze(squeeze(mean(mean(dgamma(:,2,3,:,:),1),4)))));

%% MEAN BETA

beta_mean_on=squeeze(squeeze(squeeze(squeeze(mean(mean(mean(dbeta(:,1,:,:,:),1),3),4)))));
beta_mean_off=squeeze(squeeze(squeeze(squeeze(mean(mean(mean(dbeta(:,2,:,:,:),1),3),4)))));

beta_small_on=squeeze(squeeze(squeeze(squeeze(mean(mean(dbeta(:,1,1,:,:),1),4)))));
beta_med_on=squeeze(squeeze(squeeze(squeeze(mean(mean(dbeta(:,1,2,:,:),1),4)))));
beta_large_on=squeeze(squeeze(squeeze(squeeze(mean(mean(dbeta(:,1,3,:,:),1),4)))));

beta_small_off=squeeze(squeeze(squeeze(squeeze(mean(mean(dbeta(:,2,1,:,:),1),4)))));
beta_med_off=squeeze(squeeze(squeeze(squeeze(mean(mean(dbeta(:,2,2,:,:),1),4)))));
beta_large_off=squeeze(squeeze(squeeze(squeeze(mean(mean(dbeta(:,2,3,:,:),1),4)))));

%% FIGURE

dgamma_small=squeeze(nanmean(dgamma(:,1,1,:,:),4))-squeeze(nanmean(dgamma(:,2,1,:,:),4));
dgamma_med=squeeze(nanmean(dgamma(:,1,2,:,:),4))-squeeze(nanmean(dgamma(:,2,2,:,:),4));
dgamma_large=squeeze(nanmean(dgamma(:,1,3,:,:),4))-squeeze(nanmean(dgamma(:,2,3,:,:),4));

clear dgamma_sm dgamma_sl dgamma_ml 

for a=1:length(x)
    dgamma_sm(a,:)=wjn_pt(dgamma_small(:,a)-dgamma_med(:,a));
    dgamma_sl(a,:)=wjn_pt(dgamma_small(:,a)-dgamma_large(:,a));
    dgamma_ml(a,:)=wjn_pt(dgamma_med(:,a)-dgamma_large(:,a));
end

% dummy = nan(size(x));    
%              gamma_sm=dummy;
%              gamma_sm(find(dgamma_sm<=fdr_bh(dgamma_sm)))=1;
%              gamma_sl=dummy;
%              gamma_sl(find(dgamma_sl<=fdr_bh(dgamma_sl)))=1;
%              gamma_ml=dummy;
%              gamma_ml(find(dgamma_ml<=fdr_bh(dgamma_ml)))=1;

             
             dummy = nan(size(x));    
             gamma_sm=dummy;
             gamma_sm(find(dgamma_sm<0.05))=1;
             gamma_sl=dummy;
             gamma_sl(find(dgamma_sl<0.05))=1;
             gamma_ml=dummy;
             gamma_ml(find(dgamma_ml<0.05))=1;
             
             
c=colorlover(6)
c = [c(2,:);c(4,:);c(5,:)];

figure,
for a = 3:-1:1,
    plot(x,squeeze(nanmean(nanmean(nanmean(dgamma(:,1,a,:,:),1),3),4))-squeeze(nanmean(nanmean(nanmean(dgamma(:,2,a,:,:),1),3),4)),'color',c(a,:),'linewidth',3),
%     set(gca,'Ydir','reverse'),
    hold on,
end
% ylim([-25 13]);
xlabel('Time [s]');
ylabel('\Delta relative power change [%]');
title('\Delta between gamma power ON-OFF');
% legend([h1,h2,h3],'small', 'medium', 'large','location','northeast')

legend('large', 'medium', 'small')
ylim([-7 25])
figone(7,11)
myprint('Difference_gamma_ON_OFF','E:\Roxanne\STN_Rotameter\figures');
% 
% 
%             line(x,gamma_sm*7,'linewidth',0.8, 'LineStyle','-', 'Color', [0 0 0]);
%             line(x,gamma_sl*9,'linewidth',0.8, 'LineStyle','-', 'Color', [0 0 0]);
%             line(x,gamma_ml*11,'linewidth',0.8, 'LineStyle','-', 'Color', [0 0 0]);
%            
%             text(-2.5,7,'M vs S');
%             text(-2.5,9,'L vs M');
%             text(-2.5,11,'L vs S');
            
            

figure,
c=colorlover(6)
h1=mypower(x,squeeze(squeeze(squeeze(squeeze(mean(mean(mean(dgamma(:,1,:,:,:),1),3),4))))),c(2,:),'sem',1)
hold on
h2=mypower(x,squeeze(squeeze(squeeze(squeeze(mean(mean(mean(dgamma(:,2,:,:,:),1),3),4))))),'sem',1)
xlim([-2 3])
legend([h1,h2],'ON', 'OFF','location','northwest') 

title('Gamma')
%title(sub1, 'OFF state')
xlabel('time [s]')
ylabel('power change [%, SEM]') 

figone(25,30)

annotation('textbox',[.035 .77 .2 .2],'string','OFF','linestyle','none','fontsize',20)
annotation('textbox',[.035 .3 .2 .2],'string','ON','linestyle','none','fontsize',20)