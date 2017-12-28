close all

clear all

pat=patienten;

load TF_master_results_badout_sbgmtf_efspmeeg;

condition = {'onset', 'velocity', 'amplitude'}; %'cue',
vars = {'gres'};

for x =1:length(vars);
for y=1:length(condition);

sourcevar = eval(vars{x}); 
med = {'OFF','ON'};
side = {'R','L'};
iside = {'L','R'};
hemisphere = {'contralateral','ipsilateral'};
conds = {[condition{y} '-gross'], [condition{y} '-mittel'], [condition{y} '-klein']};

data=[];

for a=1:length(med);
    fnames = fieldnames(sourcevar.(med{a})); 
    for b = 1:length(fnames)   
        for c = 1:length(side) 
            for d = 1:length(conds)
                try
                ic = ci(conds{d},sourcevar.(med{a}).(fnames{b}).(side{c}).conditions);
                ccc = ci(['STN' iside{c}],sourcevar.(med{a}).(fnames{b}).(side{c}).channels);
                cci = ci(['STN' side{c}],sourcevar.(med{a}).(fnames{b}).(side{c}).channels);
               
                data.(med{a}).contralateral(b,c,d,:,:) = squeeze(nanmean(nanmean(sourcevar.(med{a}).(fnames{b}).(side{c}).data(ccc,:,:,ic),1),4));
                data.(med{a}).ipsilateral(b,c,d,:,:) = squeeze(nanmean(nanmean(sourcevar.(med{a}).(fnames{b}).(side{c}).data(cci,:,:,ic),1),4));
               
                catch
                    disp([med{a} ' ' fnames{b} ' ' side{c} ]);
                end
            end
        end
    end
end
end
end


data.f = sourcevar.OFF.r004IR48.L.f; 
data.t = sourcevar.OFF.r004IR48.L.time;


%%
for a=1:length(med)
    for b=1:length(hemisphere)
        data.(med{a}).(hemisphere{b})(data.(med{a}).(hemisphere{b})==0)=nan; 
    end
end

condition_gamma(y,:) =squeeze(nanmean(nanmean(nanmean(nanmean(data.ON.contralateral(:,:,:,40:90,sc(data.t,0):sc(data.t,0.5)),2),3),4),5));
 m=squeeze(nanmean(nanmean(data.(med{a}).(hemisphere{b}),2),1));
 gamma_m=squeeze(nanmean(nanmean(m(:,40:90,:),2),1)); % 
 
 c=colorlover(6)
 c = [c(2,:);c(4,:);c(5,:)];

a=2
b=1

%     figure,
hold on,
    plot(data.t,gamma_m,'color',c(1,:),'linewidth',3),axis xy;
    xlim([-2 2])
    xlabel('Time [s]')
    ylabel ('Relative power change [%, SEM]')
    legend('onset', 'velocity','amplitude'); 
    title('Mean power change Gamma (40-90 Hz)');
    figone(9,14)
