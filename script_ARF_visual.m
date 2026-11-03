
%%% -- Modified by D.N.Rakhmatov, 27 September 2024 -- %%%

% IMPORTANT:
% Run main script to populate workspace first!


%% Phantom II

topSlice = 106;
botSlice = 262;

% representative z-x frame(s)
for it = [10 28 46 64 82]

    f1Amed = disp_data1A(:,:,it);
    f1AB = disp_data1A_B(:,:,it);
    f1A = disp_data1A_pre(:,:,it);
    f1B = disp_data1B_pre(:,:,it);


    figure; imagesc(x_axis*1e3,z_axis*1e3,f1Amed);
    shading flat; colormap gray; axis image; caxis([-1 1]*max(abs(f1Amed(:))));
    set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
    xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
    ylabel('z (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
    set(gca,'YDir','reverse'); set(gca,'YTick',[5 15 25 35 45]); set(gca,'XTick',[-16 -8 0 8 16]);
    saveas(gcf,sprintf('dataAmed_tIndex%d',it),'epsc');
    saveas(gcf,sprintf('dataAmed_tIndex%d.fig',it));

    figure; imagesc(x_axis*1e3,z_axis*1e3,f1AB);
    shading flat; colormap gray; axis image; caxis([-1 1]*max(abs(f1AB(:))));
    set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
    xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
    ylabel('z (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
    set(gca,'YDir','reverse'); set(gca,'YTick',[5 15 25 35 45]); set(gca,'XTick',[-16 -8 0 8 16]);
    saveas(gcf,sprintf('dataAB_tIndex%d',it),'epsc');
    saveas(gcf,sprintf('dataAB_tIndex%d.fig',it));

    figure; imagesc(x_axis*1e3,z_axis*1e3,f1A);
    shading flat; colormap gray; axis image; caxis([-1 1]*max(abs(f1A(:))));
    set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
    xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
    ylabel('z (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
    set(gca,'YDir','reverse'); set(gca,'YTick',[5 15 25 35 45]); set(gca,'XTick',[-16 -8 0 8 16]);
    saveas(gcf,sprintf('dataA_tIndex%d',it),'epsc');
    saveas(gcf,sprintf('dataA_tIndex%d.fig',it));

    figure; imagesc(x_axis*1e3,z_axis*1e3,f1B);
    shading flat; colormap gray; axis image; caxis([-1 1]*max(abs(f1B(:))));
    set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
    xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
    ylabel('z (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
    set(gca,'YDir','reverse'); set(gca,'YTick',[5 15 25 35 45]); set(gca,'XTick',[-16 -8 0 8 16]);
    saveas(gcf,sprintf('dataB_tIndex%d',it),'epsc');
    saveas(gcf,sprintf('dataB_tIndex%d.fig',it));

end


% representative x-t slice(s)
for iz = [topSlice botSlice]

    s1_before = squeeze(disp_data1A_pre(iz,:,:));
    s1_guide = squeeze(disp_data1B_pre(iz,:,:));
    s1_after = squeeze(disp_data1A_B(iz,:,:));


    figure; imagesc(x_axis*1e3,t_axis(first_frame:Nt)*1e3,s1_before(:,first_frame:Nt).');
    shading flat; colormap gray; axis image;
    set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
    xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
    ylabel('t (ms)','FontName','Courier','FontSize',14,'FontWeight','bold');
    saveas(gcf,sprintf('before_zIndex%d',iz),'epsc');
    saveas(gcf,sprintf('before_zIndex%d.fig',iz));

    figure; imagesc(x_axis*1e3,t_axis(first_frame:Nt)*1e3,s1_after(:,first_frame:Nt).');
    shading flat; colormap gray; axis image;
    set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
    xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
    ylabel('t (ms)','FontName','Courier','FontSize',14,'FontWeight','bold');
    saveas(gcf,sprintf('after_zIndex%d',iz),'epsc');
    saveas(gcf,sprintf('after_zIndex%d.fig',iz));

    figure; imagesc(x_axis*1e3,t_axis(first_frame:Nt)*1e3,s1_guide(:,first_frame:Nt).');
    shading flat; colormap gray; axis image;
    set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
    xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
    ylabel('t (ms)','FontName','Courier','FontSize',14,'FontWeight','bold');
    saveas(gcf,sprintf('guide_zIndex%d',iz),'epsc');
    saveas(gcf,sprintf('guide_zIndex%d.fig',iz));

end

