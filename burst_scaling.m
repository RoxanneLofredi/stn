clear all
close all

root=fullfile(mdf,'stn_rotameter\STN_Rotameter\');
cd(root);

addpath(fullfile(mdf,'matlab_scripts\wjn_toolbox'));
addpath(fullfile(mdf,'stn_rotameter\rotameter_scripts'));

pat=patienten;

gfiles=ffind('newoneth_gamma_brown_burst_r*.mat');

for c=1:length(gfiles);
    
    d=load(gfiles{c});
    filename=(gfiles{c});
    
    t_move = wjn_sc(d.t,0):wjn_sc(d.t,.5);
    t_base= wjn_sc(d.t,-2):wjn_sc(d.t,-1.5);
    
    ig=ci('gross',d.trials);
    im=ci('mittel',d.trials);
    is=ci('klein',d.trials);
    
        for a=1:length(d.gc)
            for b=1:length(ig)
        
        d.gburstsummov(a,b)=nansum(squeeze(d.burstraster(a,ig(b),t_move)));  
        d.gmampmov(a,b)=nanmean(d.burstamplitude(d.gc(a),ig(b),t_move));
        d.gmdurmov(a,b)=nanmean(d.burstduration(d.gc(a),ig(b),t_move));  
        
        d.grawampmov(a,b)=d.rawampmov(d.gc(a),ig(b));
        d.gbratemov(a,b)=d.bratemov(d.gc(a),ig(b));
        d.gbampmov(a,b)=d.bampmov(d.gc(a),ig(b));
        d.gbdurmov(a,b)=d.bdurmov(d.gc(a),ig(b));
            end
        end
    save(filename,'-struct','d') 
         
        for a=1:length(d.gc)
            for b=1:length(im)
        
        d.miburstsummov(a,b)=nansum(squeeze(d.burstraster(d.gc(a),im(b),t_move)));  
        d.mimampmov(a,b)=nanmean(d.burstamplitude(d.gc(a),im(b),t_move));
        d.mimdurmov(a,b)=nanmean(d.burstduration(d.gc(a),im(b),t_move));  
        
        d.mirawampmov(a,b)=d.rawampmov(d.gc(a),im(b));
        d.mibratemov(a,b)=d.bratemov(d.gc(a),im(b));
        d.mibampmov(a,b)=d.bampmov(d.gc(a),im(b));
        d.mibdurmov(a,b)=d.bdurmov(d.gc(a),im(b));
            end
        end
        
        save(filename,'-struct','d')
         
        for a=1:length(d.gc)
            for b=1:length(is)
        
        d.kburstsummov(a,b)=nansum(squeeze(d.burstraster(d.gc(a),is(b),t_move)));  
        d.kmampmov(a,b)=nanmean(d.burstamplitude(d.gc(a),is(b),t_move));
        d.kmdurmov(a,b)=nanmean(d.burstduration(d.gc(a),is(b),t_move));  
        
        d.krawampmov(a,b)=d.rawampmov(d.gc(a),is(b));
        d.kbratemov(a,b)=d.bratemov(d.gc(a),is(b));
        d.kbampmov(a,b)=d.bampmov(d.gc(a),is(b));
        d.kbdurmov(a,b)=d.bdurmov(d.gc(a),is(b));
            end
        end   
        save(filename,'-struct','d')
end
