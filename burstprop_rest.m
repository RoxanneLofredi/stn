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
    
    filename=gfiles{c};
    
    t_move = wjn_sc(d.t,0):wjn_sc(d.t,.5);
    t_base= wjn_sc(d.t,-2):wjn_sc(d.t,-1.5);
    
        for a=1:length(d.gc)
            for b=1:d.ntrials
                
        d.burstsumbase(a,b)=nansum(squeeze(d.burstraster(a,b,t_base)));  
        d.mampbase(a,b)=nanmean(d.burstamplitude(a,b,t_base));
        d.mdurbase(a,b)=nanmean(d.burstduration(a,b,t_base));
        
            end
        end
    save(filename,'-struct','d')    
end
