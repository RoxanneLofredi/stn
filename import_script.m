clear
% close all
root = 'E:\Roxanne\STN_Rotameter';
cd(root)
dateiname = 'r015NI58_R';

hand=stringsplit(dateiname,'_'); 
hand = hand{2}; 
datei = dateiname;
data = load(datei);
cond = {'klein','mittel','gross'}; 

try
    event = data.event; 
catch
    event = data.Event;
end

try
    rota = data.rota;
catch
    try
    rota = data.Rota;
    catch
        rota=data.response;
    end
end

[punkte,wert]=hist(rota.values,3); 
[~,index_max_punkte]=max(punkte); 

% if index_max_punkte >1 && strcmp(hand,'R') % wenn index>1 bedeutet das, dass meisten Punkte der Rota values 
% %     && strcmp(datei,'r025EA64_L_OFF') == 0 && strcmp(datei,'r035EE53_L_ON') == 0 && strcmp(datei,'r046AÖ47_L_ON_merge') == 0
%     rota.values = rota.values*-1;
% end

if strcmp(datei,'r041EC34_R_OFF_merge')

rota.values(898601:1599958) = rota.values(898601:1599958)*-1-mean(rota.values(898601:1599958)*-1);
    
end

rota.values = rota.values-median(rota.values);
rota.values = rota.values*-1


t = event.times;
d = event.values;
r = rota.values;
% sr = smooth(rota.values,10);

[b,a] = butter(5,.04);
sr = filter(b,a,r); 
figure,
plot(t,r)
hold on 
plot(t,sr)

dsr = mydiff(sr);
% dsr = smooth(dsr,50);

x = mythresh(d,.5);

d = round(d);

try
    if dateiname == 'r093AL56_L_OFF'
        
        x(:,2)=round(r(x(:,1)+3000)); 
    else
        x(:,2)=d(x(:,1)+500); 
    end
catch
    x(end)=[]
    x(:,2)=d(x(:,1)+500); 
end

x(x(:,2)==0,:) = []; 

figure
plot(t,d);
hold on
plot(t(x(:,1)),d(x(:,1)),'r+')
plot(t,sr)
pause

nn=0;
trl = [];
for a = 1:length(x);
    if a < length(x) && x(a,1)+10000 < length(t);
        epoche = [x(a,1)+100:x(a,1)+10000]; 
    elseif a==length(x) || x(a,1)+10000 > length(t)  
        epoche = [x(a,1)+100:x(a,1)+200]; 
    end

try
%     edsr=smooth(dsr(epoche),10);
esr = sr(epoche);  
nesr = (esr-min(esr));
rnesr = nesr./max(nesr); 
    
edsr = filter(b,a,dsr(epoche));
rdsr = edsr./max(edsr);

    istart = mythresh(rdsr,.15);

istart = istart(1); 
[nvmax,ivmax] = max(rdsr(istart:istart+1000));
[namax,iamax] =  max(rnesr(istart:istart+ivmax+1000)); 
nn=nn+1; 

figure
if strcmp(dateiname,'r093AL56_L_OFF')
plot(epoche,round(r(epoche)))
else
    plot(epoche,d(epoche))
end
title([strrep(datei,'_',' ') ' Trial Nummer ' num2str(nn) ' C: ' num2str(x(a,2)) ]) % strrep (string replace) in der datei, das Zeichen '_' mit ' ' ersetzen
hold on %
plot(epoche,rnesr)%
plot(epoche,rdsr) %
plot(epoche(istart),rnesr(istart),'or') % 
plot(epoche(istart+ivmax),rdsr(istart+ivmax),'ob') %
plot(epoche(istart+iamax),rnesr(istart+iamax),'og') % 
legend('event','amplitude','speed','istart','ivmax','iamax')% 
dcmObj = datacursormode; set(dcmObj,'UpdateFcn',@updateFcn); % 
drawnow
pause % 
% keyboard
titlestring = {'N','Condition','Cue Onset','Movement Onset','Highest Velocity','Highest Amplitude'}; % eigene Tabelle die die Titel der Spalten der gespeicherten Variable enthält
trl(nn,:) = [nn x(a,2) x(a,1) x(a,1)+istart+100 x(a,1)+istart+ivmax+100 x(a,1)+istart+iamax+100]; % Variable die berechnete Werte enthält (an row nn -> reihe an der sich skript gerade befindet und in alle Spalten soll folgened eingetragen werden [ Trialnummer, Wert der Variable x an Stelle a und Zeile 2 (Gibt Kondition des events an: aufgerundeter Wert: 1 / 2 / 3 = klein/mittel/groß), Spalte 1 aus Variable x in reihe a (X-achsenwert der threshold Überschreitung); Thresholdüberschreitung + intervall bis istart + 100 (?) um absoluten x-achsenabschnittswert von istart zu erhalten, dasselbe für ivmax und iamax im anschluss]
conditionlabels{nn} = cond{x(a,2)}; %  


close 
catch ME
    disp(['Trial ' num2str(a) ' weg']) %
%     keyboard
    continue % 
end
end

save(['nntrl_' datei],'trl','conditionlabels') % 
