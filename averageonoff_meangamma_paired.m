close all
clear all

root= 'E:\Roxanne\STN_Rotameter';
cd(root)

addpath E:\Roxanne\scripts;
pat = patienten;

med={'ON', 'OFF'};
cond= {'onset-klein', 'onset-mittel', 'onset-gross'};

% cluster=bwlabeln(~isnan(wjn_disp_nii('E:\Roxanne\STN_Rotameter\models\contralateral_fullfactorial\alles_baseline.nii',[-3 3],[1 100],0)));
% igamma = cluster==3;
% igamma(1:40,:) = 0;
% ibeta = cluster == 2;
% itheta = cluster ==1;
imagefolder = 'E:\Roxanne\STN_Rotameter\contralateral_images';

non=0;
noff=0;



for a=1:length(pat)
    if ~isfield(pat{a},'bad')  && isfield (pat{a}, 'onoff') && isfield(pat{a}, 'updrs')
        
            if ~pat{a}.med %&& numel(pat{a}.updrs)==2;
            noff=noff+1;
            idoff{noff} = pat{a}.id;
            %nmessungoff(noff)=a;
             uoff(noff) = pat{a}.updrs(2);

            for b=1:length(cond);

            filename=['contralateral_' cond{b} '_sbgmtf_efspmeeg_' pat{a}.id ];

            [tfdata,x] = wjn_disp_nii(fullfile(imagefolder,[filename '.nii']),[-3 3],[1 100],0);
            gammafreq= 40:90 %[pat{a}.gammafreq(1):pat{a}.gammafreq(2)];
            gammatime= sc(x,0):sc(x,0.5); %[sc(x,pat{a}.gammatime(1)):sc(x,pat{a}.gammatime(2))]
            
            betafreq=13:30;
            betatime=sc(x,0):sc(x,0.5);
            
            thetafreq=2:8;
     
            
            data=tfdata(gammafreq,gammatime);
            betadata=tfdata(betafreq, betatime);
            thetadata=tfdata(thetafreq, betatime);
            %data(data<0)=nan;
            offmgamma(noff,b) = squeeze(mean(mean(data,1),2));
            offmbeta (noff,b)=  squeeze(mean(mean(betadata,1),2));
            offmtheta (noff,b)=  squeeze(mean(mean(thetadata,1),2));
%             offmtheta (noff,b) =  mean(tfdata(itheta));
            
            end

            elseif pat{a}.med %&& numel(pat{a}.updrs)==2;
            non=non+1
            idon{non} = pat{a}.id;
            %nmessungon(non)=a;
            uon(non) = pat{a}.updrs(1);
            for b=1:length(cond);
            
            filename=['contralateral_' cond{b} '_sbgmtf_efspmeeg_' pat{a}.id ];

            [tfdata,x] = wjn_disp_nii(fullfile(imagefolder,[filename '.nii']),[-3 3],[1 100],0);
            gammafreq= 40:90 %[pat{a}.gammafreq(1):pat{a}.gammafreq(2)];
            gammatime= sc(x,-0.1):sc(x,1) 
%             gammafreq = 70:90;
%             gammatime = sc(x,.25):sc(x,1.5);
            betafreq=13:30;
            betatime=sc(x,0):sc(x,0.5);
            
            thetafreq=2:8;

            data=tfdata(gammafreq,gammatime);
            betadata=tfdata(betafreq, betatime);
            thetadata=tfdata(thetafreq, betatime);
            %data(data<0)=nan;
            onmgamma(non,b) = squeeze(mean(mean(data,1),2));
            onmbeta (non,b)=  squeeze(mean(mean(betadata,1),2));
            onmtheta(non,b)=  squeeze(mean(mean(thetadata,1),2));
%             onmbeta (non,b)=  mean(tfdata(ibeta));
%             onmtheta (non,b) =  mean(tfdata(itheta));
            
            end

            end
        end
end

uon = uon';
uoff = uoff';
dupdrs=(uoff-uon)./uoff*100;
f=[];
f=dupdrs<30;
dupdrs(f)=[];
onmgamma(f,:)=[];
offmgamma(f,:)=[];

onmbeta(f,:)=[];
offmbeta(f,:)=[];

onmtheta(f,:)=[];
offmtheta(f,:)=[];
%% Mean Beta

close all
clear all

root= 'E:\Roxanne\STN_Rotameter';
cd(root)

addpath E:\Roxanne\scripts;
pat = patienten;

med={'ON', 'OFF'};
cond= {'onset-klein', 'onset-mittel', 'onset-gross'};

% cluster=bwlabeln(~isnan(wjn_disp_nii('E:\Roxanne\STN_Rotameter\models\contralateral_fullfactorial\alles_baseline.nii',[-3 3],[1 100],0)));
% igamma = cluster==3;
% igamma(1:40,:) = 0;
% ibeta = cluster == 2;
% itheta = cluster ==1;
imagefolder = 'E:\Roxanne\STN_Rotameter\contralateral_images';

non=0;
noff=0;



for a=1:length(pat)
    if ~isfield(pat{a},'bad')  && isfield (pat{a}, 'onoff') && isfield(pat{a}, 'updrs')
        
            if ~pat{a}.med %&& numel(pat{a}.updrs)==2;
            noff=noff+1;
            idoff{noff} = pat{a}.id;
            %nmessungoff(noff)=a;
             uoff(noff) = pat{a}.updrs(2);

            for b=1:length(cond);

            filename=['contralateral_' cond{b} '_sbgmtf_efspmeeg_' pat{a}.id ];

            [tfdata,x] = wjn_disp_nii(fullfile(imagefolder,[filename '.nii']),[-3 3],[1 100],0);
            betafreq= 13:30 %[pat{a}.gammafreq(1):pat{a}.gammafreq(2)];
            gammatime= sc(x,-0.1):sc(x,1) %[sc(x,pat{a}.gammatime(1)):sc(x,pat{a}.gammatime(2))];
            
            data=tfdata(gammafreq,gammatime);
            %data(data<0)=nan;
            offmgamma(noff,b) = squeeze(mean(mean(data,1),2));
%             offmbeta (noff,b)=  mean(tfdata(ibeta));
%             offmtheta (noff,b) =  mean(tfdata(itheta));
            
            end

            elseif pat{a}.med %&& numel(pat{a}.updrs)==2;
            non=non+1
            idon{non} = pat{a}.id;
            %nmessungon(non)=a;
            uon(non) = pat{a}.updrs(1);
            for b=1:length(cond);
            
            filename=['contralateral_' cond{b} '_sbgmtf_efspmeeg_' pat{a}.id ];

            [tfdata,x] = wjn_disp_nii(fullfile(imagefolder,[filename '.nii']),[-3 3],[1 100],0);
            gammafreq= 40:90 %[pat{a}.gammafreq(1):pat{a}.gammafreq(2)];
            gammatime= sc(x,-0.1):sc(x,1) 
%             gammafreq = 70:90;
%             gammatime = sc(x,.25):sc(x,1.5);

            data=tfdata(gammafreq,gammatime);
            %data(data<0)=nan;
            onmgamma(non,b) = squeeze(mean(mean(data,1),2));
%             onmbeta (non,b)=  mean(tfdata(ibeta));
%             onmtheta (non,b) =  mean(tfdata(itheta));
            
            end

            end
        end
end

uon = uon';
uoff = uoff';
dupdrs=(uoff-uon)./uoff*100;
f=[];
f=dupdrs<30;
dupdrs(f)=[];
onmgamma(f,:)=[];
offmgamma(f,:)=[];


%% p values

monmgamma= mean(onmgamma,2);
moffmgamma = mean(offmgamma,2);
dgamma = monmgamma-moffmgamma;

monmbeta= mean(onmbeta,2);
moffmbeta = mean(offmbeta,2);
dbeta = monmbeta-moffmbeta;

monmtheta= mean(onmtheta,2);
moffmtheta = mean(offmtheta,2);
dtheta = monmtheta-moffmtheta;

pall = permTest(100000,dgamma, zeros(size(dgamma)));
pallbeta=permTest(100000,dbeta, zeros(size(dbeta)));
palltheta=permTest(100000,dtheta, zeros(size(dtheta)));

% [r,p]=mypermCorr(dgamma,dupdrs,'pearson');
% [l1,l2]=mycorrline(dupdrs, dgamma, 30, 90);
% 
% 
% [ron,pon]=mypermCorr(monmgamma,uon,'spearman');
% % [ron,pon]=corr(log(monmgamma),log(uon),'type','pearson');
% [lon1,lon2]=mycorrline(monmgamma,uon, 1, 5);
% 
% [roff,poff]=mypermCorr(moffmgamma,uoff,'spearman');
% [loff1,loff2]=mycorrline(moffmgamma,uoff, 1, 5);

%% Figure

figure,
subplot(2,1,1)
scatter (log(monmgamma), log(uon)), 
xlabel ('Mean Gamma Change [%]', 'FontSize', 10), ylabel('UPDRS score', 'FontSize', 10), 
title (['Rho=' num2str(ron) ' P=' num2str(pon)], 'FontSize', 12);
hold on
plot(lon1,lon2);

subplot(2,1,2)
scatter (log(moffmgamma), log(uoff)), 
xlabel ('Mean Gamma Change [%]', 'FontSize', 10), ylabel('UPDRS score', 'FontSize', 10), 
title (['Rho=' num2str(roff) ' P=' num2str(poff)], 'FontSize', 12);
hold on
plot(loff1,loff2);

figure, 
subplot(3,1,1)
figure,
scatter (dupdrs, dgamma), 
xlabel ('Delta UPDRS [%]', 'FontSize', 18), ylabel('Mean Gamma Change [%]', 'FontSize', 18), 
title (['Rho=' num2str(r) ' P=' num2str(p)], 'FontSize', 18);
hold on
plot(l1,l2);




%% Figure
close

cmapoff = [0,0,0;0,0,0;0,0,0]; %([0 0 0],:); %cmaps([1 1 1],:);
cmapon = [0.8 0.8 0.8; 0.8 0.8 0.8; 0.8 0.8 0.8]; %cmaps([3 3 3],:);
cmap= [0 0 0; 0.8 0.8 0.8];

mybar({moffmgamma,monmgamma}, cmap)
figone(7,8)

if pall <= .05 
   sigbracket('P<0.001',[1 2],23,9)  
%    sigbracket(['P = ' num2str(pall,1)],[1 2],34,9)
end

set(gca,'XTicklabel',{'OFF' 'ON'})
ylabel('Relative power change [%, SEM]') 
title('Mean Gamma(40-90 Hz,0-.5 s)')



mybar({moffmbeta,monmbeta}, cmap)
figone(7,8)

   sigbracket(['P = ' num2str(pallbeta,1)],[1 2],-37,9)

set(gca,'XTicklabel',{'OFF' 'ON'})
ylabel('Relative power change [%, SEM]') 
ylim([-43 0 ])
set(gca,'Ydir','reverse')
title('Mean Beta(13-30, 0 bis .5)')


mybar({moffmtheta,monmtheta}, cmap)
figone(7,8)

   sigbracket(['P = ' num2str(palltheta,1)],[1 2],90,9)

set(gca,'XTicklabel',{'OFF' 'ON'})
ylabel('Relative power change [%, SEM]') 
title('Mean Theta(2-8, 0 bis .5)')


%% PermTest

for a = 1:3,pcondonoff(a) = permTest(10000,onmgamma(:,a)-offmgamma(:,a), zeros(size(dgamma)));end

for a = 1:3,pcondonoffbeta(a) = permTest(10000,onmbeta(:,a)-offmbeta(:,a), zeros(size(dbeta)));end


% if pcondonoff <= .05 
%    sigbracket('*',[1 3],bylims(a,2)-2)
% end
% 
% if mplmon <= .05 
%    sigbracket('*',[1 2],bylims(a,2)-5)
% end
% 
% if mpmson <= .05 
%    sigbracket('*',[2 3],bylims(a,2)-8)
% end
% 
% set(gca,'XTicklabel',{'small' 'medium' 'large'})
% ylabel('power change [%, SEM]') 
% ylim(bylims(a,:))

% for a = 1:3,pcondonoff(a) = permTest(10000,onmgamma(:,a)-offmgamma(:,a), zeros(size(dgamma)));end

pcondonsm= permTest(10000,onmgamma(:,1)-onmgamma(:,2),zeros(size(dgamma)));
pcondonsl= permTest(10000,onmgamma(:,1)-onmgamma(:,3),zeros(size(dgamma)));
pcondonml= permTest(10000,onmgamma(:,2)-onmgamma(:,3),zeros(size(dgamma)));

pcondoffsm= permTest(10000,offmgamma(:,1)-offmgamma(:,2),zeros(size(dgamma)));
pcondoffsl= permTest(10000,offmgamma(:,1)-offmgamma(:,3),zeros(size(dgamma)));
pcondoffml= permTest(10000,offmgamma(:,2)-offmgamma(:,3),zeros(size(dgamma)));

pcondonsmbeta= permTest(10000,onmbeta(:,1)-onmbeta(:,2),zeros(size(dbeta)));
pcondonslbeta= permTest(10000,onmbeta(:,1)-onmgamma(:,3),zeros(size(dgamma)));
pcondonmlbeta= permTest(10000,onmbeta(:,2)-onmgamma(:,3),zeros(size(dgamma)));

pcondoffsmbeta= permTest(10000,offmbeta(:,1)-offmbeta(:,2),zeros(size(dbeta)));
pcondoffslbeta= permTest(10000,offmbeta(:,1)-offmbeta(:,3),zeros(size(dbeta)));
pcondoffmlbeta= permTest(10000,offmbeta(:,2)-offmbeta(:,3),zeros(size(dbeta)));


%% Figure 
cmap = colorlover(6,0);
cmap = cmap([2 4 5 2 4 5],:);

close all
figure
mybar([offmgamma onmgamma],cmap)
xlabel('OFF                                           ON')
set(gca,'XTick',[1 2 3 4 5 6 ],'XTickLabel',{'small','medium','large','small','medium','large'}) %,'XTickLabelRotation',45)
ylabel('Relative power change [%,SEM]')
ylim([0 53])
title ('Mean Gamma (40 90 0 0.5)')
sigbracket(['P = ' num2str(pcondonsm,2)],[4 5],32,9) %sigbracket('P<0.01',[4 5],48,9) 
sigbracket(['P =' num2str(pcondonml,2)],[5 6],32,9) %sigbracket('P <0.05 ',[5 6],48,9) %
sigbracket(['P = ' num2str(pcondonsl,2)],[4 6],38,9) %sigbracket(['P<0.0001'],[4 6],55,9) %
sigbracket(['P = ' num2str(pall,2)],[2 5],44,9) %sigbracket('P<0.05',[2 5],64,9) % 

sigbracket(['P = ' num2str(pcondoffsm,2)],[1 2],32,9) %sigbracket('P<0.05',[1 2],48,9) 
sigbracket(['P= ' num2str(pcondoffml,2)],[2 3],32,9) %sigbracket('P =0.5 ',[2 3],48,9)
sigbracket(['P = ' num2str(pcondoffsl,2)],[1 3],38,9) %sigbracket('P<0.05',[1 3],55,9)

ylim([0 73])

figone(9,14)


%% Figure
% cmaps = colorlover(5,0);
% cmapoff = cmaps([1 1 1],:);
% cmapon = cmaps([3 3 3],:);
close all

cmapoff = [0,0,0;0,0,0;0,0,0]; %([0 0 0],:); %cmaps([1 1 1],:);
cmapon = [0.8 0.8 0.8; 0.8 0.8 0.8; 0.8 0.8 0.8]; %cmaps([3 3 3],:);

figure

boff=mybar(offmgamma,cmapoff,.8:1:2.8)
hold on
bon=mybar(onmgamma,cmapon,1.2:1:3.2)
set(boff,'barwidth',.4)
set(bon,'barwidth',.4)
ylabel('Relative power change [%,SEM]')
set(gca,'XTick',[1:3],'XTickLabel',{'small','medium','large'})
title('Mean Gamma ONOFF 40-90 0-0.5')
legend([boff(1) bon(1)],{'OFF','ON'},'location','northeastoutside')

sigbracket(['P=' num2str(pcondonoff(1),2)],[.8 1.2],33,9)
sigbracket(['P=' num2str(pcondonoff(2),2)],[1.8 2.2],33,9)
sigbracket(['P=' num2str(pcondonoff(3),2)],[2.8 3.2],33,9)

ylim ([0 40]);

% sigbracket(['P = ' num2str(pcondonoff(1),2)],[.8 1.2],40,9)
% sigbracket(['P = ' num2str(pcondonoff(2),2)],[1.8 2.2],40,9)
% sigbracket(['P = ' num2str(pcondonoff(3),2)],[2.8 3.2],40,9)

% sigbracket(['P = ' num2str(pcondonsm,2)],[1.2 2.1],50,9)
% sigbracket(['P = ' num2str(pcondonml,2)],[2.3 3.2],50,9)
% sigbracket(['P = ' num2str(pcondonsl,2)],[1.2 3.2],58,9)

figone(7,14)

%nmessung(monmgamma<3)
%% write excel files

offfreq={'offmgamma'} %, 'offmbeta', 'offmtheta'};
onfreq={'onmgamma'} %, 'onmbeta', 'onmtheta'};
freq={'gamma'} %, 'beta', 'theta'};
 
for a=1:length(offfreq);

    off=eval(offfreq{a});
    on=eval(onfreq{a});
    
    xlswrite('mean_freqrange_value_allcond_meangamma_paired', {['ON_' freq{a} '_' 'small'], ['ON_' freq{a} '_' 'medium'], ['ON_' freq{a} '_' 'large'], ['OFF_' freq{a} '_' 'small'], ['OFF_' freq{a} '_' 'medium'], ['OFF_' freq{a} '_' 'large']},freq{a});
    
    xlswrite ('mean_freqrange_value_allcond_meangamma_paired', off, freq{a}, 'A2');
    xlswrite ('mean_freqrange_value_allcond_meangamma_paired', on, freq{a}, 'D2'); 

    
end
