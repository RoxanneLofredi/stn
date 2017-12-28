clear all
close all

pat=patienten;

cond = {'onset-klein','onset-mittel','onset-gross'};

noff=0;
non=0;

neuesrotaoff=[];
neuesrotaon=[];

cid = '';

for a=1:length(pat); 
    
   id=pat{a}.id(1:4);
   
      if ~isfield (pat{a}, 'bad') && ~strcmp(id,cid) && pat{a}.med;
          cid=pat{a}.id(1:4);
          files = ffind(['efspmeeg_' pat{a}.id(1:4) '*ON*.mat']);
        sneuesrotaon=[];
        sdneuesrotaon=[];
        srton=[];
          for b=1:length(files)
              cside = strsplit(files{b},'_')
              cside = cside{3};
               D=spm_eeg_load(files{b});
               ic = ci('rota',D.chanlabels);
               if isempty(ic)
                    ic=ci('response', D.chanlabels);
               end
      

            for c=1:length (cond);
                sneuesrotaon(b,c,:) = mean(D(ic,:,ci(cond{c},D.conditions)),3);
                if strcmp(cside,'R')
                    sneuesrotaon(b,c,:) = sneuesrotaon(b,c,:)*-1;
                end
                sneuesrotaon(b,c,:) = (sneuesrotaon(b,c,:)-min(sneuesrotaon(b,c,:)));
                sdneuesrotaon(b,c,:) = smooth(mydiff(sneuesrotaon(b,c,:)),100);
                srton(b,c) = mean(D.trialonset(ci(cond{c},D.conditions))-D.trialonset(ci(cond{c},D.conditions)-1));
            end
        
          end
         
             non=non+1;
            idon{non,1} = pat{a}.id;
            
                neuesrotaon(non,:,:) = mean(sneuesrotaon,1);
                dneuesrotaon(non,:,:) = mean(sdneuesrotaon,1);
                rton(non,:) = mean(srton,1);
    end
end

%%
% neuesrotaon=neuesrotaon ./1000;
dneuesrotaon=dneuesrotaon*1000;
t=D.time;

for a =1:length(t);
    vpls(a) = wjn_pt(dneuesrotaon(:,1,a)-dneuesrotaon(:,3,a),[],500);
end

for a =1:length(t);
vplm(a)= wjn_pt(dneuesrotaon(:,1,a)-dneuesrotaon(:,2,a),[],500);
vpms(a)=wjn_pt(dneuesrotaon(:,2,a)-dneuesrotaon(:,3,a),[],500);
end


for a = 1:3
    for b = 1:size(neuesrotaon,3);
        for c = 1:size(neuesrotaon,1);
        neuesrotaon(c,a,b)=neuesrotaon(c,a,b)-nanmean(neuesrotaon(c,a,1:2800),3);
        end
    end
end

cc=colorlover(6);
 figure
 subplot(1,2,1);
    klein=mypower(D.time,squeeze(neuesrotaon(:,1,:)),cc(2,:),'sem');
    hold on
        mittel=mypower(D.time,squeeze(neuesrotaon(:,2,:)),cc(4,:),'sem');
            gross=mypower(D.time,squeeze(neuesrotaon(:,3,:)),cc(5,:),'sem');
   
   title ('Angular position');
    ylabel('mean rel. angular position');
    xlabel ('time [s]');
    xlim ([-1 3]);
    ylim([-0.5 2.5]);

 subplot(1,2,2);
     klein=mypower(D.time,squeeze(dneuesrotaon(:,1,:)),cc(2,:),'sem');
    hold on
        mittel=mypower(D.time,squeeze(dneuesrotaon(:,2,:)),cc(4,:),'sem');
            gross=mypower(D.time,squeeze(dneuesrotaon(:,3,:)),cc(5,:),'sem');
  
 
            
            dummy = nan(size(t));    
             vms=dummy;
             vms(find(vpms<=fdr_bh(vpms)))=1;
             vlm=dummy;
             vlm(find(vplm<=fdr_bh(vplm)))=1;
             vls=dummy;
             vls(find(vpls<=fdr_bh(vpls)))=1;
  
             
% [vmscluster,vmsiclusters,vmscluster_sizes] = wjn_remove_clusters(vms,5);
% [vlmcluster,vlmiclusters,vlmcluster_sizes] = wjn_remove_clusters(vlm,5);
% [vlscluster,vlsiclusters,vlscluster_sizes] = wjn_remove_clusters(vls,5);

            line(t,vms*-0.8,'linewidth',0.8, 'LineStyle','-', 'Color', [0 0 0]);
            line(t,vlm*-0.5,'linewidth',0.8, 'LineStyle','-', 'Color', [0 0 0]);
            line(t,vls*-0.2,'linewidth',0.8, 'LineStyle','-', 'Color', [0 0 0]);
           
            text(-1.8,-0.8,'M vs S');
            text(-1.8,-0.5,'L vs M');
            text(-1.8,-0.2,'L vs S');
            
    ylim([-1 5.5])
    xlim ([-2 2]);
    
    title ('Mean angular speed');
    ylabel('mean rel. angular speed [a.u.]');
    xlabel ('time [s]');
    legend ([klein mittel gross],'small', 'medium', 'large','location','northeast');

 
suptitle ('Experimental set-up, FDR corr, n=16');
    
figone(12,20)
   
%%

von=[max(abs(squeeze(dneuesrotaon(:,1,:)))');max(abs(squeeze(dneuesrotaon(:,2,:)))');max(squeeze(abs(dneuesrotaon(:,3,:)))')];
von=von';

%%

% idcond={'idoff', 'idon'};
% rtcond={'rtoff', 'rton'};
vcond={'voff', 'von'};
med={'OFF', 'ON'};

pat=patienten;

%%
mvon= mean(von,2);

ponvelsm=permTest(1000,von(:,1)-von(:,2),zeros(size(von,1)));
ponvelml=permTest(1000,von(:,2)-von(:,3),zeros(size(von,1)));
ponvelsl=permTest(1000,von(:,1)-von(:,3),zeros(size(von,1)));

figure
mybar(von)
if ponvelsm<.05;
sigbracket('*',[1 2])
end

if ponvelml<.05;
sigbracket('*',[2 3])
end

if ponvelsl<.05;
sigbracket('*',[1 3]);
end

set(gca,'XTickLabel',{'small','medium','large'})
title('ON n=16 (permTest, p=0)')
ylabel('Velocity [au]')

%%
v = nan(11,6);
v(:,1:3) = von;
v(1:11,4:6)=voff;

cmap = colorlover(5);
cmap = [cmap([1,3,4],:);cmap([1,3,4],:)];

figure,
mybar(v,cmap)
set(gca,'XTickLabel',{'small','medium','large','small','medium','large'})
ylabel('Velocity [au]')
xlabel('ON                               OFF')
sigbracket('*',[1 3],0.007)
sigbracket('*',[1 2],0.006)
sigbracket('*',[2 3],0.006)
sigbracket('*',[4 6],0.007)
sigbracket('*',[4 5],0.006)
sigbracket('*',[5 6],0.006)
if pvall<.05
    sigbracket(['P=' num2str(pvall)], [2 5], 0.008);
end

title('OnOff Patients (n=11)');
ylim([0 0.009])
figone(7,12)
myprint(fullfile(cd,'figures','Velocity_ON_OFF_bar_plot_onoffpat'));


rt = nan(27,6);
rt(:,1:3) = rton;
rt(1:18,4:6)=rtoff;

mrton=mean(rton,2);
mrtoff=mean(rtoff,2);
% drtonoff=mrton-mrtoff;
% prtonoff=permTest(1000,drtonoff,zeros(size(drtonoff)));

ponrtsm=permTest(1000,rton(:,1)-rton(:,2),zeros(size(rton,1)));
ponrtml=permTest(1000,rton(:,2)-rton(:,3),zeros(size(rton,1)));
ponrtsl=permTest(1000,rton(:,1)-rton(:,3),zeros(size(rton,1)));

poffrtsm=permTest(1000,rtoff(:,1)-rtoff(:,2),zeros(size(rtoff,1)));
poffrtml=permTest(1000,rtoff(:,2)-rtoff(:,3),zeros(size(rtoff,1)));
poffrtsl=permTest(1000,rtoff(:,1)-rtoff(:,3),zeros(size(rtoff,1)));

cmap = colorlover(5);
cmap = [cmap([1,3,4],:);cmap([1,3,4],:)];

figure,
mybar(rt,cmap)
set(gca,'XTickLabel',{'small','medium','large','small','medium','large'})
ylabel('Reaction time [s]')
xlabel('ON                               OFF')
if ponrtsm<.06;
sigbracket(['* P=' num2str(ponrtsm)],[1 2],0.8)
end

if ponrtml<.06;
sigbracket(['* P=' num2str(ponrtml)],[2 3],0.8)
end

if ponrtsl<.06;
sigbracket(['* P=' num2str(ponrtsl)],[1 3],1);
end

if poffrtsm<.06;
sigbracket(['* P=' num2str(poffrtsm)],[4 5],0.8)
end

if poffrtml<.06;
sigbracket(['* P=' num2str(poffrtml)],[5 6],0.8)
end

if poffrtsl<.06;
sigbracket(['* P=' num2str(poffrtsl)],[4 6],1);
end

% sigbracket('P=0.03',[1 3],0.8)
% sigbracket('P=0.05',[4 6],0.8)
% sigbracket ('P=0.5', [2 5], 1.2)
title('All');
ylim([0 1.2])
figone(7,12)
myprint(fullfile(cd,'figures','Reactiontime_ON_OFF_bar_plot_all'))