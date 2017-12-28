

pat = patienten;

% copy2otherfile = 0;

for a = 1:length(pat)
    
%     if ~copy2otherfile
    
    keep pat a %copy2otherfile
    
    if isfield (pat{a},'eichamp') && ~isfield (pat{a},'bad');
        
    D=spm_eeg_load(['fspmeeg_' pat{a}.id]);
    cD=D; 
    
    x = pat{a}.eichamp(:,1);
    i = ci('rota',D.chanlabels);
    
    if isempty(i)
       i = ci('response',D.chanlabels);
    end
    
    D= chanlabels(D,i,'rota'); 
    save(D)
    
    for b=1:length(x);
        tx(b) = sc(D.time,x(b)-.1); 
    end
    
    sd = [];
    
    for b = 1:length(x) 
        if x(b)
  
            sd = ones(1,3000);
            imax(b) = 1;
            c=0;
            
            while isempty(find(sd(imax(b):end)<=0.00001,1,'first'))
                c=c+500;
                sd = abs(smooth(mydiff(D(i,tx(b):-1:tx(b)-(3000+c),1)),100));
                [~,imax(b)]=max(sd);
            end
            

       
        imin(b) = find(sd(imax(b):end)<=0.00001,1,'first')+imax(b);     
        
        else
            imax(b) = nan;
            imin(b) =  nan;
        end
    end     
    % figure, plot(D(i,tx(b)-5000:tx(b)+5000,1))
    % figure, plot(sd),hold on,plot(ones(size(sd))*0.00001)
    
    figure
    subplot(1,2,1)
    plot(D.time,squeeze(D(i,:,1)))
    hold on
    plot(D.time,squeeze(D(ci('event',D.chanlabels),:,1)))
    
    if sum(x)
    hold on
    scatter(D.time(tx(tx>1)),squeeze(D(i,tx(tx>1),1)))
    scatter(D.time(tx(tx>1)-imin(~isnan(imin))),squeeze(D(i,tx(tx>1)-imin(~isnan(imin)),1)))
    end  
    title(strrep(pat{a}.id,'_',' '))
    
    for b = 1:length(tx);
        if tx(b)>1
         maxamp(b) = squeeze(D(i,tx(b),1));
         minamp(b) = D(i,tx(b)-imin(b),1);
        else
            maxamp(b) = nan;
            minamp(b) = nan;
        end
    end
    
        if a==54;
        minamp= [-1.024 -1.022];
        end
 
        if a==55;
        minamp= [-0.9789 1.388];
        maxamp= [1.579 -0.9573];
        x=[52.9 912.4];
        end
        
        if a==56;
        minamp= [0];
        maxamp=[4.648];
        end
        
        if a==53;
        minamp= [-0.9286 1.294 1.256];
        maxamp= [1.524 -0.9889 -0.6932];
        end        
        
        if a==33;
        minamp= 1.635;
        maxamp= -1.058;
        x=12.8;
        end        
         
    D.maxamp = maxamp;
    D.minamp = minamp;
    

    save(D);
    
    D=spm_eeg_load(['e' D.fname]);
    D.eichamp = pat{a}.eichamp;
    
    if length(x) == 1;
        xtrials = [1 D.ntrials];
    else
    for b = 2:length(x);
        xtrials(b-1,2) = sc(D.trialonset,x(b))-1;
    end
    
    xtrials(b,2) = D.ntrials;
    
    for b = 1:length(xtrials)
        if b ==1
        xtrials(b,1) = 1;
        else
            xtrials(b,1) = xtrials(b-1,2)+1;
        end
    end
    end
    
    D.maxamp = maxamp; %epochiertes file
    D.minamp = minamp;
    D.xtrials = xtrials;
    save(D)
    
    i = ci('rota',D.chanlabels);
    
    if isempty(i)
    i=ci('response', D.chanlabels);
    end

    dr = squeeze(D(i,:,:));
    rdr=[];
    
   for b = 1:length(x);
       trials = xtrials(b,1):xtrials(b,2);
       
   if isnan(minamp(b)) || isnan(maxamp(b));
       
       ctrials = trials(find(ismember(trials,ci('cue-gross',D.conditions))));
       atrials = trials(find(ismember(trials,ci('amplitude-gross',D.conditions))));
       D.minamp(b) = nanmedian(dr(sc(D.time,0),ctrials));
       D.maxamp(b) = nanmedian(dr(sc(D.time,0),atrials));
      
   end
   end
   
    for b = 1:length(x);
 
            rdr(:,D.xtrials(b,1):D.xtrials(b,2))=100*((dr(:,D.xtrials(b,1):D.xtrials(b,2))-D.minamp(b))./(D.maxamp(b)-D.minamp(b)));

    end
    
    subplot(1,2,2)
    plot(D.time,rdr(:,ci('amplitude-gross',D.conditions)))
    hold on
    plot(D.time,rdr(:,ci('amplitude-mittel',D.conditions)));
    plot(D.time,rdr(:,ci('amplitude-klein',D.conditions)));
    figone(15,30)
    myprint(fullfile('eichpics',pat{a}.id))
    
    nrdr = [];
    nrdr(1,:,:) = rdr(:,:);
%     
    
    nD=clone(D,['r' D.fname]);
    nD(i,:,:) = nrdr;
    save(nD)
    D=nD;
    itrials = ci('amplitude',D.conditions);
    trialamp = [];
    amp = [];acc=[];trialacc =[];trialvel = [];
    vel = [];
    for b = 1:length(itrials);
        amp(b) = rdr(sc(D.time,0),itrials(b));
        trialamp = [trialamp ones(1,4)*amp(b)];
        if strcmp(D.conditions{itrials(b)}(end),'s');
            acc(b) = amp(b)-(3/3*100);
            bcond(b) = 3;
        elseif strcmp(D.conditions{itrials(b)}(end),'l');
            acc(b) = amp(b)-(2/3*100);
            bcond(b) = 2;
        elseif strcmp(D.conditions{itrials(b)}(end),'n');
            acc(b) = amp(b)-(1/3*100);
            bcond(b) = 1;
        end
        trialacc= [trialacc ones(1,4)*acc(b)];
        
   
        vel(b) =  max(abs(smooth(mydiff(D(i,sc(D.time,-.5):sc(D.time,.5),itrials(b)-1)),50)));
        trialvel= [trialvel ones(1,4)*vel(b)];
    end
    
    rvel = vel/max(vel)*100;
    trialrvel = trialvel/max(vel)*100;
    

    
    D.rvel = rvel;
    D.vel = vel;
    D.acc = acc;
    D.amp = amp;
    D.bcond = bcond;
    D.trialvel = trialvel;
    D.trialrvel = trialrvel;
    D.trialacc = trialacc;
    D.trialamp = trialamp;
    
    bt=badtrials(D);
    d=ones(size(D.trialvel));
    d(bt)=0;
    gt=find(d);
    D.btrialvel = D.trialvel(gt);
    D.btrialrvel = D.trialrvel(gt);
    D.btrialacc = D.trialacc(gt);
    D.btrialamp = D.trialamp(gt);
    
    save(D)
    end
%     else
%         
%       nD=spm_eeg_load(['espmeeg_' pat{a}.id ]);
%       D=spm_eeg_load(['efspmeeg_' pat{a}.id ])
%       
%       
%     D.eichamp = pat{a}.eichamp;
%     D.maxamp = nD.maxamp; %epochiertes file
%     D.minamp = nD.minamp;
%     D.xtrials = nD.xtrials;
%             D.rvel = nD.rvel;
%     D.vel = nD.vel;
%     D.acc = nD.acc;
%     D.amp = nD.amp;
%     D.bcond = nD.bcond;
%     D.trialvel = nD.trialvel;
%     D.trialrvel = nD.trialrvel;
%     D.trialacc = nD.trialacc;
%     D.trialamp = nD.trialamp;
%     
%     bt=badtrials(D);
%     d=ones(size(D.trialvel));
%     d(bt)=0;
%     gt=find(d);
%     D.btrialvel = D.trialvel(gt);
%     D.btrialrvel = D.trialrvel(gt);
%     D.btrialacc = D.trialacc(gt);
%     D.btrialamp = D.trialamp(gt);
%     
%     save(D)
%     end
end


    
    
    

    