
clear all
close all
% seven = {'r032DA55'    'r035EE53'    'r091NF42'    'r094RE60'    'r011AI40'    'r041EC34'  'r093AL56'};
% iids = ci(seven,ids);

pat=patienten;

cond = {'onset-klein','onset-mittel','onset-gross'};

noff=0;
non=0;

neuesrotaoff=[];  
param = {'btrialrt','btrialvel','btrialacc'};
conds = {'klein','mittel','gross'};
med = {'OFF','ON'};
n=0;
non = 0;
noff = 0;

for a=1:length(pat);  
   
      if ~isfield (pat{a}, 'bad') && isfield(pat{a}, 'onoff');
         if pat{a}.med
          non=non+1;
          idson{non} = pat{a}.id(1:8);
         else
              noff = noff+1;
              idsoff{noff} = pat{a}.id(1:8);
          end
 
          files = ffind(['rtf_efspmeeg_' pat{a}.id(1:4) '*' med{pat{a}.med+1} '*.mat']);
          if isempty(files)
              disp(num2str(a))
          end
          for b=1:length(files)
                
            D=spm_eeg_load(files{b});
%             D.rt = D.trialonset(ci('onset',D.conditions))-D.trialonset(ci('cue',D.conditions));
%             D.btrialrt = [];
%             for c = 1:length(D.rt);
%                 D.btrialrt = [D.btrialrt ones(1,4)*D.rt(c)];
%             end
%             
%             save(D);
            for d = 1:length(param);
            for c = 1:length(conds)
      
                aic = ci(['onset-'],D.conditions);
                ic = ci(['onset-' conds{c}],D.conditions);
                temp.(param{d})(b,c) = nanmean(abs(D.(param{d})(ic)));
               
            end                     
            end
          end     
          
          for d = 1:length(param)
              if pat{a}.med
                  on.(param{d})(non,:) = nanmean(temp.(param{d}));
              else
                  off.(param{d})(noff,:) = nanmean(temp.(param{d}));
              end
          end
      end

end
            
%%

seven = {'r032DA55'    'r035EE53'    'r091NF42'    'r094RE60'    'r011AI40'    'r041EC34'  'r093AL56'};
iids_off = ci(seven,idsoff);
iids_on = ci(seven,idson);
param = {'btrialrt','btrialvel','btrialacc'};

cmap = colorlover(6);
cmap = [cmap([2,4,5],:);cmap([2,4,5],:)];

for a = 2
    
    d = off.(param{a})(iids_off,:)-on.(param{a})(iids_on,:);
    moff=off.(param{a})(iids_off,:);
    mon=on.(param{a})(iids_on,:);
    dk = d(:,1);
    dm = d(:,2);
    dg = d(:,3);
    
    onk = mon(:,1);
    onm = mon(:,2);
    ong = mon(:,3);
    offk = moff(:,1);
    offm = moff(:,2);
    offg = moff(:,3);

    ponoff = permTest(10000,mean(moff,2)-mean(mon,2),zeros(7,1));
    ponoffk= permTest(10000,mean(offk,2)-mean(onk,2),zeros(7,1));
    ponoffm= permTest(10000,mean(offm,2)-mean(onm,2),zeros(7,1));
    ponoffg= permTest(10000,mean(offg,2)-mean(ong,2),zeros(7,1));
    
    ponkm = permTest(10000,mean(onk,2)-mean(onm,2),zeros(7,1));
    ponmg = permTest(10000,mean(onm,2)-mean(ong,2),zeros(7,1));
    ponkg = permTest(10000,mean(onk,2)-mean(ong,2),zeros(7,1));
    
    poffkm = permTest(10000,mean(offk,2)-mean(offm,2),zeros(7,1));
    poffmg = permTest(10000,mean(offm,2)-mean(offg,2),zeros(7,1));
    poffkg = permTest(10000,mean(offk,2)-mean(offg,2),zeros(7,1));
    
    
    figure
    mybar([on.(param{a})(iids_on,:) off.(param{a})(iids_off,:)],cmap)
    title([param{a} ' n=7'])
    set(gca,'XTickLabel',{'small','medium','large','small','medium','large'})
    xlabel('ON                              OFF')
    ylabel('Mean velocity [a.u.]')
    
if ponkm<.001;
sigbracket('P<.001',[1 2],0.3,9)   
% sigbracket(['P<.001' num2str(ponkm,1)],[1 2],0.8)
end

if ponmg<.001;
sigbracket('P<0.001',[2 3],0.3,9)
% sigbracket(['P=' num2str(ponmg,1)],[2 3],0.8)
end

if ponkg<.001;
    sigbracket('P<0.001',[1 3],0.4,9);
% sigbracket(['P=' num2str(ponkg,1)],[1 3],1);
end

if poffkm<.001;
sigbracket('P<.001',[4 5],0.3,9)
end

if poffmg<.001;
sigbracket('P<.001',[5 6],0.3,9)
% sigbracket(['P=' num2str(poffmg,1)],[5 6],0.8)
end

if poffkg<.001;
sigbracket('P<.001',[4 6],0.4,9);
end

sigbracket(['P=' num2str(ponoff,1)],[2 5],0.5,9);


ylim([0 0.6])
figone(7,12)

end
myprint(fullfile(cd,'figures', [param{a} '_ON_OFF_n=7']));


cmaps = colorlover(5,0);
cmapoff = ones(3,3)*0;
cmapon = ones(3,3).*0.8;

figure
boff=mybar(moff,cmapoff,.8:1:2.8)
hold on
bon=mybar(mon,cmapon,1.2:1:3.2)
title([param{a} ' n=7'])
set(boff,'barwidth',.4)
set(bon,'barwidth',.4)
set(gca,'XTick',[1:3],'XTickLabel',{'small','medium','large'})
ylabel('Mean velocity [a.u.]');
legend([boff(1) bon(1)],{'OFF','ON'},'location','northeastoutside')
sigbracket(['P = ' num2str(ponoffk,1)],[.8 1.2],0.3,9)
sigbracket(['P = ' num2str(ponoffm,1)],[1.8 2.2],0.3,9)
sigbracket(['P = ' num2str(ponoffg,1)],[2.8 3.2],0.3,9)
figone(7,14)


for a = 2
    figure
    mybar([on.(param{a}) off.(param{a})])
    title([param{a} ' n=7'])
    
    d = off.(param{a})(iids,:)-on.(param{a})(iids,:);
    moff=off.(param{a})(iids,:);
    mon=on.(param{a})(iids,:);
    dk = d(:,1);
    dm = d(:,2);
    dg = d(:,3);
    
    onk = mon(:,1);
    onm = mon(:,2);
    ong = mon(:,3);
    offk = moff(:,1);
    offm = moff(:,2);
    offg = moff(:,3);

    ponoff = permTest(10000,mean(moff,2)-mean(mon,2),zeros(7,1));
    ponoffk= permTest(10000,mean(offk,2)-mean(onk,2),zeros(7,1));
    ponoffm= permTest(10000,mean(offm,2)-mean(onm,2),zeros(7,1));
    ponoffg= permTest(10000,mean(offg,2)-mean(ong,2),zeros(7,1));
    
    ponkm = permTest(10000,mean(onk,2)-mean(onm,2),zeros(7,1));
    ponmg = permTest(10000,mean(onm,2)-mean(ong,2),zeros(7,1));
    ponkg = permTest(10000,mean(onk,2)-mean(ong,2),zeros(7,1));
    
    poffkm = permTest(10000,mean(offk,2)-mean(offm,2),zeros(7,1));
    poffmg = permTest(10000,mean(offm,2)-mean(offg,2),zeros(7,1));
    poffkg = permTest(10000,mean(offk,2)-mean(offg,2),zeros(7,1));
    
    set(gca,'XTickLabel',{'small','medium','large','small','medium','large'})
    xlabel('ON                               OFF')
    
if ponkm<.06;
sigbracket(['P=' num2str(ponkm)],[1 2],0.4)
end

if ponmg<.06;
sigbracket(['P=' num2str(ponmg)],[2 3],0.4)
end

if ponkg<.06;
sigbracket(['P=' num2str(ponkg)],[1 3],0.5);
end

if poffkm<.06;
sigbracket(['P=' num2str(poffkm)],[4 5],0.4)
end

if poffmg<.06;
sigbracket(['P=' num2str(poffmg)],[5 6],0.4)
end

if poffkg<.06;
sigbracket(['P=' num2str(poffkg)],[4 6],0.5);
end

sigbracket(['P=' num2str(ponoff)],[2 5],0.6);


ylim([0 0.7])
figone(7,12)
myprint(fullfile(cd,'figures', [param{a} '_ONOFF_n=7']));


cmaps = colorlover(5,0);
cmapoff = cmaps([1 1 1],:);
cmapon = cmaps([3 3 3],:);

figure
boff=mybar(moff,cmapoff,.8:1:2.8)
hold on
bon=mybar(mon,cmapon,1.2:1:3.2)
title([param{a} ' n=7'])
set(boff,'barwidth',.4)
set(bon,'barwidth',.4)
set(gca,'XTick',[1:3],'XTickLabel',{'small','medium','large'})
legend([boff(1) bon(1)],{'OFF','ON'},'location','northeastoutside')
sigbracket(['P = ' num2str(ponoffk)],[.8 1.2],0.4,9)
sigbracket(['P = ' num2str(ponoffm)],[1.8 2.2],0.4,9)
sigbracket(['P = ' num2str(ponoffg)],[2.8 3.2],0.4,9)
figone(7,14)

end


for a = 3
    figure
    mybar([on.(param{a}) off.(param{a})])
    title([param{a} ' n=7'])
    
    d = off.(param{a})(iids,:)-on.(param{a})(iids,:);
    moff=off.(param{a})(iids,:);
    mon=on.(param{a})(iids,:);
    dk = d(:,1);
    dm = d(:,2);
    dg = d(:,3);
    
    onk = mon(:,1);
    onm = mon(:,2);
    ong = mon(:,3);
    offk = moff(:,1);
    offm = moff(:,2);
    offg = moff(:,3);

    ponoff = permTest(10000,mean(moff,2)-mean(mon,2),zeros(7,1));
    ponoffk= permTest(10000,mean(offk,2)-mean(onk,2),zeros(7,1));
    ponoffm= permTest(10000,mean(offm,2)-mean(onm,2),zeros(7,1));
    ponoffg= permTest(10000,mean(offg,2)-mean(ong,2),zeros(7,1));
    
    ponkm = permTest(10000,mean(onk,2)-mean(onm,2),zeros(7,1));
    ponmg = permTest(10000,mean(onm,2)-mean(ong,2),zeros(7,1));
    ponkg = permTest(10000,mean(onk,2)-mean(ong,2),zeros(7,1));
    
    poffkm = permTest(10000,mean(offk,2)-mean(offm,2),zeros(7,1));
    poffmg = permTest(10000,mean(offm,2)-mean(offg,2),zeros(7,1));
    poffkg = permTest(10000,mean(offk,2)-mean(offg,2),zeros(7,1));
    
    set(gca,'XTickLabel',{'small','medium','large','small','medium','large'})
    xlabel('ON                               OFF')
    
if ponkm<.06;
sigbracket(['P=' num2str(ponkm)],[1 2],5)
end

if ponmg<.06;
sigbracket(['P=' num2str(ponmg)],[2 3],5)
end

if ponkg<.06;
sigbracket(['P=' num2str(ponkg)],[1 3],6);
end

if poffkm<.06;
sigbracket(['P=' num2str(poffkm)],[4 5],5)
end

if poffmg<.06;
sigbracket(['P=' num2str(poffmg)],[5 6],5)
end

if poffkg<.06;
sigbracket(['P=' num2str(poffkg)],[4 6],6);
end

sigbracket(['P=' num2str(ponoff)],[2 5],7);


ylim([0 8])
figone(7,12)
myprint(fullfile(cd,'figures', [param{a} '_ONOFF_n=7']));

cmaps = colorlover(5,0);
cmapoff = cmaps([1 1 1],:);
cmapon = cmaps([3 3 3],:);

figure
boff=mybar(moff,cmapoff,.8:1:2.8)
hold on
bon=mybar(mon,cmapon,1.2:1:3.2)
title([param{a} ' n=7'])
set(boff,'barwidth',.4)
set(bon,'barwidth',.4)
set(gca,'XTick',[1:3],'XTickLabel',{'small','medium','large'})
legend([boff(1) bon(1)],{'OFF','ON'},'location','northeastoutside')
sigbracket(['P = ' num2str(ponoffk)],[.8 1.2],5,9)
sigbracket(['P = ' num2str(ponoffm)],[1.8 2.2],5,9)
sigbracket(['P = ' num2str(ponoffg)],[2.8 3.2],5,9)
figone(7,14)

myprint(fullfile(cd,'figures', [param{a} '_ONOFFcondcomp_n=7']));
end
