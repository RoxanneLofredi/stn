close all
clear all

load LEAD_groupanalysis.mat;

cond= {'velocity-klein', 'velocity-mittel', 'velocity-gross'};

channels_r={'STNR01','STNR12','STNR23'};
channels_l={'STNL01','STNL12','STNL23'};
contact_r=[1 2 3];
contact_l=[4 5 6];

n=[];
n=0;

for a=1:length(pat);
    
    if ~isfield(pat{a},'bad')  && pat{a}.med==1;
    
    pat_id=strsplit(pat{a}.id,'_');
    id=pat_id{1}(2:end);
    pat_pos=strmatch(id,M.patient.list(:,3));
    
    filename=['sbgmtf_efspmeeg_' pat{a}.id '.mat'];
    D=spm_eeg_load(filename);
    
    gammafreq= 40:90; 
    
    if pat{a}.side=='L';
    
    try
    mni=M.elstruct(pat_pos).coords_mm{1}(1:3,:)+0.5*diff(M.elstruct(pat_pos).coords_mm{1});
    catch
    disp ([id ' not found...']);
    mni=[nan nan nan; nan nan nan; nan nan nan];
    end
    
    for b=1:length(channels_r);
    
    if ci(channels_r(b),D.chanlabels);
    
    n=n+1;
    
    gamma_tf(n,:,:,:) = squeeze(nanmean(D(ci(channels_r(b),D.chanlabels), :, :, ci('velocity',D.conditions)),4));
   
    gamma_all(n,:) = squeeze(nanmean(nanmean(nanmean(D(ci(channels_r(b),D.chanlabels), gammafreq, sc(D.time,-0.25):sc(D.time,0.25), ci('velocity',D.conditions)),2),3),4));
    gamma_l(n,:)= squeeze(nanmean(nanmean(nanmean(D(ci(channels_r(b),D.chanlabels), gammafreq, sc(D.time,-0.25):sc(D.time,0.25), ci(cond(3),D.conditions)),2),3),4));
    gamma_m(n,:) = squeeze(nanmean(nanmean(D(ci(channels_r(b),D.chanlabels), gammafreq, sc(D.time,-0.25):sc(D.time,0.25), ci(cond(2),D.conditions)),2),3));
    gamma_s(n,:) = squeeze(nanmean(nanmean(D(ci(channels_r(b),D.chanlabels), gammafreq, sc(D.time,-0.25):sc(D.time,0.25), ci(cond(1),D.conditions)),2),3));
    
    %gamma_max(n,:)=squeeze(nanmax(nanmax(nanmean(D(ci(channels_r(b),D.chanlabels), gammafreq, sc(D.time,-0.25):sc(D.time,0.25), ci('velocity',D.conditions)),[],2),[],3),4));
    
    npat(n,:)=pat{a}.npatient;
    side(n,:)=contact_r(b);
    
    mni_x(n,:)=mni(b,1);
    mni_y(n,:)=mni(b,2);
    mni_z(n,:)= mni(b,3);
    
    else
    
    n=n+1;
    
    gamma_all(n,:)=nan;
    gamma_l(n,:)=nan;
    gamma_m(n,:)=nan;
    gamma_s(n,:)=nan;
    
   % gamma_max(n,:)=nan;
    
    npat(n,:)=pat{a}.npatient;
    side(n,:)=contact_r(b);
    
    mni_x(n,:)=mni(b,1);
    mni_y(n,:)=mni(b,2);
    mni_z(n,:)= mni(b,3);
    
    end
    end
    
    else
    try    
    mni=M.elstruct(pat_pos).coords_mm{2}(1:3,:)+0.5*diff(M.elstruct(pat_pos).coords_mm{2});
    catch
    disp ([id ' not found...'])
    mni=[nan nan nan; nan nan nan; nan nan nan];
    end
    
    hemi=2;
    
    for b=1:length(channels_l);
    
    if ci(channels_l(b),D.chanlabels)
    
    n=n+1;
    gamma_tf(n,:,:,:) = squeeze(nanmean(D(ci(channels_l(b),D.chanlabels), :, :, ci('velocity',D.conditions)),4));
   
    gamma_all(n,:) = squeeze(nanmean(nanmean(nanmean(D(ci(channels_l(b),D.chanlabels), gammafreq, sc(D.time,-0.25):sc(D.time,0.25), ci('velocity',D.conditions)),2),3),4));
    gamma_l(n,:)= squeeze(nanmean(nanmean(nanmean(D(ci(channels_l(b),D.chanlabels), gammafreq, sc(D.time,-0.25):sc(D.time,0.25), ci(cond(3),D.conditions)),2),3),4));
    gamma_m(n,:) = squeeze(nanmean(nanmean(D(ci(channels_l(b),D.chanlabels), gammafreq, sc(D.time,-0.25):sc(D.time,0.25), ci(cond(2),D.conditions)),2),3));
    gamma_s(n,:) = squeeze(nanmean(nanmean(D(ci(channels_l(b),D.chanlabels), gammafreq, sc(D.time,-0.25):sc(D.time,0.25), ci(cond(1),D.conditions)),2),3));
    
    %gamma_max(n,:)=squeeze(nanmax(nanmax(nanmean(D(ci(channels_l(b),D.chanlabels), gammafreq, sc(D.time,-0.25):sc(D.time,0.25), ci('velocity',D.conditions)),2),3),4));
    
    npat(n,:)=pat{a}.npatient;
    side(n,:)=contact_l(b);
    mni_x(n,:)=mni(b,1);
    mni_y(n,:)=mni(b,2);
    mni_z(n,:)= mni(b,3);
    
    else
    
    n=n+1;
    
    gamma_all(n,:)=nan;
    gamma_l(n,:)=nan;
    gamma_m(n,:)=nan;
    gamma_s(n,:)=nan;
    %gamma_max(n,:)=nan;
    
    npat(n,:)=pat{a}.npatient;
    side(n,:)=contact_l(b);
    mni_x(n,:)=mni(b,1);
    mni_y(n,:)=mni(b,2);
    mni_z(n,:)= mni(b,3);
    
    end
    end
    end
    end
end

mni_gamma=[npat side mni_x mni_y mni_z gamma_all gamma_l gamma_m gamma_s ]; %gamma_max

scatter3(abs(mni_x),mni_y,mni_z,40,gamma_all,'filled');
ax = gca;
ax.XDir = 'reverse';
view(-31,14)
xlabel('lat-med')
ylabel('postero-ant')
zlabel('vent-dors')
caxis([-10,50]);
cb=colorbar;
cb.Label.String='Rel. Gamma Power change %';
title('Gamma Heat Map');
% myprint('gamma_heat_scatter');


