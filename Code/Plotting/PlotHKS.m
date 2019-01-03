function PlotHKS(folderpath)
    %PlotHKS plots HKS and saves as eps at different calculation times.
    %   PlotHKS(FOLDERPATH) Plots the HKS values for a part at 5 equally
    %   spaced time intervals (except start first interval to avoid t=0),
    %   and saves figures in FOLDERPATH as partname-HKS-PLOT_TIME#s.eps

    % Load global file separator fs
    global fs

    % Retreive filename from folderpath
    indsep = strfind(folderpath,fs);
    %tradename = folderpath(1:indsep(1)-1);
    partname = folderpath(indsep(1)+1:indsep(2)-1);

    % Load HKS and other variables for HKS plotting
    file=['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat'];
    load(file,'tri','coord','HKS','stepold','iterold');

    % Generate level 2 screen comment
    ScreenComment('',['Creating and saving ' partname ' HKS plots']);

    % Plot and save eps HKS figure at specified 
    plot_inds_frac = [0.025,0.25,0.5,0.75,1];
    plot_inds = floor(plot_inds_frac.*iterold);
    plot_times = plot_inds.*stepold;
    for i = 1:length(plot_inds)
        plot_ind = plot_inds(i);
        plot_name = [partname '-HKS-' num2str(plot_inds(i))];
        plot_dir = ['ShapeTerra' fs 'Output' fs folderpath fs plot_name];
        h=figure();
        trisurf(tri,coord(:,1),coord(:,2),coord(:,3),HKS(:,plot_ind))
        axis equal 
        title(plot_name)
        colorbar
        saveas(h, plot_dir, 'fig')
        close(h)
    end
end