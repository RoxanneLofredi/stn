
clear
pat = patienten;

% close all


timerange = [-0.5 1];
freqrange = [13 35; 40 90] %; 100 250;];
franges = {'beta','gamma'}; %,'hfo'};
med={'ON','OFF'};
conds = {'onset-klein'; 'onset-mittel'; 'onset-gross'};


for b=1:size(freqrange,1);
        non=0;
        noff=0;
          
for a = 1:length(pat);
   if ~isfield (pat{a}, 'bad') && isfield(pat{a}, 'dupdrs');       
        if ~pat{a}.med
            noff=noff+1;
            for b=1:length(conds);

            filename=['contralateral_' conds{b} '_sbgmtf_efspmeeg_' pat{a}.id ];

            [tfdata,x] = wjn_disp_nii(fullfile(imagefolder,[filename '.nii']),[-3 3],[1 100],0);
            
            gammafreq= 40:90
            gammatime= sc(x,-0.5):sc(x,3)
            gammadata=tfdata(gammafreq,gammatime);
            offmgamma(noff,b) = squeeze(mean(mean(gammadata,1),2));
            
            betafreq=13:30
            betatime= sc(x,-0.5):sc(x,3)
            betadata=tfdata(betafreq,betatime);
            offmbeta(noff,b) = squeeze(mean(mean(betadata,1),2));
            
            end

            elseif pat{a}.med;
            non=non+1
            idon{non} = pat{a}.id;
            uon(non) = pat{a}.updrs(1);
            
            for b=1:length(cond);
            
            filename=['contralateral_' cond{b} '_sbgmtf_efspmeeg_' pat{a}.id ];

             [tfdata,x] = wjn_disp_nii(fullfile(imagefolder,[filename '.nii']),[-3 3],[1 100],0);
            
            gammafreq= 40:90
            gammatime= sc(x,-0.5):sc(x,3)
            gammadata=tfdata(gammafreq,gammatime);
            offmgamma(noff,b) = squeeze(mean(mean(gammadata,1),2));
            
            betafreq=13:30
            betatime= sc(x,-0.5):sc(x,3)
            betadata=tfdata(betafreq,betatime);
            offmbeta(noff,b) = squeeze(mean(mean(betadata,1),2));
            
            end

            end
        end
end

end

    if ~pat{a}.med
 
    noff=noff+1;
    
    for c=1:length(conds);
         ic = ci(conds{c},D.conditions);
         
            doff(noff,c,b,:) = squeeze(squeeze((nanmean(nanmean(D(cci,freqrange(b,1):freqrange(b,2),:,ic),1),2))));
       
    end
    else
        non=non+1;
    
    for c=1:length(conds);
         ic = ci(conds{c},D.conditions);
            don(non,c,b,:) = squeeze(squeeze((nanmean(nanmean(D(cci,freqrange(b,1):freqrange(b,2),:,ic),1),2))));
       
    end
    end
   end
end
end
 


%% Figure

figure
c=colorlover(6)
freqbands ={'Beta-Band [13-35 Hz]','Gamma-Band [60-90 Hz]','HFO [100-250 Hz]'}

for a = 1:3
    
sub1=subplot(2,3,a);
h1=mypower(D.time, squeeze(doff(:,1,a,:)),c(2,:),'sem',1) 
hold on
h2=mypower(D.time, squeeze(doff(:,2,a,:)),c(3,:),'sem',1) 
h3=mypower(D.time, squeeze(doff(:,3,a,:)),c(5,:),'sem',1) 
xlim([-2 3])

if a==1
legend([h1,h2,h3],'small', 'medium', 'large','location','northwest') 
end

title([freqbands{a}])
%title(sub1, 'OFF state')
xlabel('time [s]')
ylabel('power change [%, SEM]') 

sub2=subplot(2,3,a+3);
h1=mypower(D.time, squeeze(don(:,1,a,:)),c(2,:),'sem',1) 
hold on
h2=mypower(D.time, squeeze(don(:,2,a,:)),c(3,:),'sem',1) 
h3=mypower(D.time, squeeze(don(:,3,a,:)),c(5,:),'sem',1) 
xlim([-2 3])

if a==1
legend([h1,h2,h3],'small', 'medium', 'large','location','northwest')
end

title(freqbands{a})
%title(sub2, 'ON state')
xlabel('time [s]')
ylabel('power change [%, SEM]') 

end
figone(25,30)

annotation('textbox',[.035 .77 .2 .2],'string','OFF','linestyle','none','fontsize',20)
annotation('textbox',[.035 .3 .2 .2],'string','ON','linestyle','none','fontsize',20)

myprint('mean power change all freq')

myprint(fullfile(cd,'powerchange_freqrange_filter_smooth_contralat_onset_5-1s_onoff','mean power change onoff'));


%%
close all

tranges = [0 1;0 1;0 1];

cm = [c(2,:);c(3,:);c(5,:)];
ylims = [-50 25; -10 36; -5 15];
bylims = [-30 0; 0 22; 0 12];

for a = 1:length(franges)
    trange = [sc(D.time,tranges(a,1)):sc(D.time,tranges(a,2))];
    for b=1:length(conds)
    
    if a==1
    
    mdon1(:,b) = mean(don(:,b,1,trange),4);
    mdoff1(:,b) = mean(doff(:,b,1,trange),4);
   
    elseif a==2
    
    mdon2(:,b) = mean(don(:,b,a,trange),4);
    mdoff2(:,b) = mean(doff(:,b,a,trange),4);
    
    elseif a==3
    mdon3(:,b) = mean(don(:,b,a,trange),4);
    mdoff3(:,b) = mean(doff(:,b,a,trange),4);
    end

    end
end

%% write excel files

offfreq={'mdoff1', 'mdoff2', 'mdoff3'};
onfreq={'mdon1', 'mdon2', 'mdon3'};
freq={'beta', 'gamma', 'hfo'};


    
for a=1:length(offfreq);

    off=eval(offfreq{a});
    on=eval(onfreq{a});
    
    xlswrite('powerchange_freqrange_filter_smooth_contralat_onset_5-1s_onoff', {['OFF_' freq{a} '_' 'small'], ['OFF_' freq{a} '_' 'medium'], ['OFF_' freq{a} '_' 'large'], ['ON_' freq{a} '_' 'small'], ['ON_' freq{a} '_' 'medium'], ['ON_' freq{a} '_' 'large']},freq{a});
    
    xlswrite ('powerchange_freqrange_filter_smooth_contralat_onset_5-1s_onoff', off, freq{a}, 'A2');
    xlswrite ('powerchange_freqrange_filter_smooth_contralat_onset_5-1s_onoff', on, freq{a}, 'D2'); 

    
end

%%

for a=1: length(freqbands)
 
 mdon= eval(['mdon' num2str(a)]);
 mdoff = eval(['mdoff' num2str(a)]);
    
figure

subplot (2,3,[1 2]);

h1=mypower(D.time, squeeze(don(:,1,a,:)),c(2,:),'sem',1) 
hold on
h2=mypower(D.time, squeeze(don(:,2,a,:)),c(3,:),'sem',1) 
h3=mypower(D.time, squeeze(don(:,3,a,:)),c(5,:),'sem',1) 

xlim([-2 3])
legend([h1,h2,h3],'small', 'medium', 'large','location','northeast')
title([freqbands{a} ' ' 'ON'])
xlabel('time [s]')
ylabel('power change [%, SEM]') 

if  a==1
    legend([h1,h2,h3],'small', 'medium', 'large','location','southeast')
end

ylim(ylims(a,:));

subplot (2,3,3);

mybar(mdon,cm)
mplson = permTest(10000,mdon(:,1),mdon(:,3));
mplmon = permTest(10000,mdon(:,1),mdon(:,2));
mpmson = permTest(10000,mdon(:,2),mdon(:,3));


if mplson <= .05 
   sigbracket('*',[1 3],bylims(a,2)-2)
end

if mplmon <= .05 
   sigbracket('*',[1 2],bylims(a,2)-5)
end

if mpmson <= .05 
   sigbracket('*',[2 3],bylims(a,2)-8)
end

set(gca,'XTicklabel',{'small' 'medium' 'large'})
ylabel('power change [%, SEM]') 
ylim(bylims(a,:))

subplot (2,3,[4 5]);

h1=mypower(D.time, squeeze(doff(:,1,a,:)),c(2,:),'sem',1) 
hold on
h2=mypower(D.time, squeeze(doff(:,2,a,:)),c(3,:),'sem',1) 
h3=mypower(D.time, squeeze(doff(:,3,a,:)),c(5,:),'sem',1) 
legend([h1,h2,h3],'small', 'medium', 'large','location','southeast')
if a==2
    ylim([-25 25])
    legend([h1,h2,h3],'small', 'medium', 'large','location','northeast')
elseif a==3
    ylim([-25 25])
    legend([h1,h2,h3],'small', 'medium', 'large','location','northeast')
end
ylim(ylims(a,:))
xlim([-2 3])

title([freqbands{a} ' ' 'OFF'])
xlabel('time [s]')
ylabel('power change [%, SEM]') 


subplot (2,3,6);

mybar(mdoff,cm)
mplsoff = permTest(10000,mdoff(:,1),mdoff(:,3));
mplmoff = permTest(10000,mdoff(:,1),mdoff(:,2));
mpmsoff = permTest(10000,mdoff(:,2),mdoff(:,3));

if mplsoff <= .05 
   sigbracket('*',[1 3],bylims(a,2)-2)
end

if mplmoff <= .05 
   sigbracket('*',[1 2],bylims(a,2)-5)
end

if mpmsoff <= .05 
   sigbracket('*',[2 3],bylims(a,2)-8)
end

set(gca,'XTicklabel',{'small' 'medium' 'large'})
ylabel('Power change [%,SEM]')
ylim(bylims(a,:))
figone(16,23)

myprint(['mean power change onoff' ' ' franges{a}]);

myprint(fullfile(cd,'powerchange_freqrange_filter_smooth_contralat_onset_5-1s_onoff', ['mean power change onoff' ' ' franges{a}]));

end

save powerchange_freqrange_allchannels_filter_onoff
