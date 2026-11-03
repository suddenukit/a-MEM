
%%% -- Modified by D.N.Rakhmatov, 27 September 2024 -- %%%

clear all; close all;
clc;

addpath(genpath('C:\Users\daler\Desktop\Reference_US\USTB'));

addpath(genpath('C:\Users\daler\Desktop\Reference_US\aMEM\ImpulseIQ'));


%% Read displacement data/parameters

%phantomID = 'ARF_phantom_I';
phantomID = 'ARF_phantom_II';
%phantomID = 'ARF_phantom_III';

caseID = [phantomID '_z2t1'];


load(['disp_auto_' caseID '.mat']);

% 2D spatial grid
x_axis = disp_auto.x_axis;
dx = disp_auto.dx;
Nx = disp_auto.Nx;

z_axis = disp_auto.z_axis;
dz = disp_auto.dz;
Nz = disp_auto.Nz;

% push/tracking info
x_push = disp_auto.x_push;
Nt = disp_auto.Nt;
dt = disp_auto.dt;
t_axis = dt*(0:Nt-1).';
tx_ratio = dt/dx;

ix_offset1 = floor((x_push(1)-x_axis(1))/dx) + 1;
ix_offset2 = floor((x_push(2)-x_axis(1))/dx) + 1;
ix_offset2 = Nx - ix_offset2;

% other parameters
izExt2 = disp_auto.izExt2;
itExt1 = disp_auto.itExt1;
first_frame = disp_auto.first_frame;

% RHS (1) and LHS (2) displacement data
disp_data1A = disp_auto.dataA{1};           % progressive referencing
disp_data2A = disp_auto.dataA{2};           % progressive referencing
disp_data1B = disp_auto.dataB{1};           % fixed referencing
disp_data2B = disp_auto.dataB{2};           % fixed referencing


%% Configure SWS mapping parameters

swsRange = [0.2 5.0];                   % SWS range: [min, max] m/s

swsRangeAdj = swsRange*tx_ratio;                       
angleMin = atand(swsRangeAdj(1));       % min SWS angle [deg]
angleMax = atand(swsRangeAdj(2));       % max SWS angle [deg]

x_window = 2.0*1e-3;                    % lateral window length [m]
ixExt2 = max(1,round(x_window/dx/2));   % lateral two-sided window extent

gfx_window = 1.0*1e-3;                  % guided filter: lateral span [m]
gft_window = 1.0*1e-3;                  % guided filter: temporal span [s]
ixt_filtG = [max(3,round(gfx_window/dx)) max(3,round(gft_window/dt))];

zx_filtM = [1.0 1.0]*1e-3;              % median filter size
izx_filtM = [max(3,round(zx_filtM(1)/dz)) max(3,round(zx_filtM(2)/dx))];

zx_filtW = [2.0 2.0]*1e-3;              % Wiener filter size
izx_filtW = [max(3,round(zx_filtW(1)/dz)) max(3,round(zx_filtW(2)/dx))];

sigma = 0.5;                            % Gaussian filter: standard dev.
ifkx_filtS = [3 3];                     % Gaussian filter: size

sws_res = 0.10;                         % SWS accuracy limit [m/s]
dtInterp = dx*sws_res;                  % interpolated t-axis step size [s]
t_axisInterp = (t_axis(first_frame):dtInterp:t_axis(end)).';
NtInterp = length(t_axisInterp);

t_offset = x_window/swsRange(2);
it_offset = round(t_offset/dtInterp);   % min lag at max SWS


%% Configure Fourier grid

NtFFT = pow2(nextpow2(Nt)+1);
NxFFT = pow2(nextpow2(Nx)+1);

% quadrant mask for [min max] velocity trapezoid 
kx_slope = 1 - angleMax/90;     % [0, 1] = [keep-all, reject-all]
f_slope = angleMin/90;          % [0, 1] = [keep-all, reject-all]
kxFFT = round((1:NxFFT/2)*kx_slope);
fFFT = round((1:NtFFT/2)*f_slope);
subMask = ones(NxFFT/2,NtFFT/2);
for kx = 1:NxFFT/2
    subMask(kx,1:kxFFT(kx)) = 0;
end
for f = 1:NtFFT/2
    subMask(1:fFFT(f),f) = 0;
end
% [0 1] fractional band --> [1 NtFFT/2] index range
f_bw = [0 1];	
f_min = max(1,round(f_bw(1)*NtFFT/2));
f_max = min(NtFFT/2,round(f_bw(2)*NtFFT/2));   
subMask(:,1:f_min-1) = 0;
subMask(:,f_max+1:NtFFT/2) = 0;

% full-spectrum mask
fftMask = ones(NxFFT,NtFFT);
% discard backward velocities
fftMask(1:NxFFT/2,1:NtFFT/2) = 0;
fftMask(NxFFT/2+1:NxFFT,NtFFT/2+1:NtFFT) = 0;
% keep [min max] velocities
fftMask(1:NxFFT/2,NtFFT/2+1:NtFFT) = ...
    fftMask(1:NxFFT/2,NtFFT/2+1:NtFFT).*flip(subMask,1);
fftMask(NxFFT/2+1:NxFFT,1:NtFFT/2) = ...
    fftMask(NxFFT/2+1:NxFFT,1:NtFFT/2).*flip(subMask,2);

%fftMask = imgaussfilt(fftMask,sigma,'FilterSize',ifkx_filtS);


%% Prepare displacement data

% filter (z,x) frames
for it = 1:Nt

    temp_frame = disp_data1A(:,:,it);
    temp_frame = medfilt2(temp_frame,izx_filtM);
    temp_frame = temp_frame/max(abs(temp_frame(:)));
    temp_frame(isnan(temp_frame)) = 0;
    disp_data1A(:,:,it) = temp_frame;

    temp_frame = disp_data2A(:,:,it);
    temp_frame = medfilt2(temp_frame,izx_filtM);
    temp_frame = temp_frame/max(abs(temp_frame(:)));
    temp_frame(isnan(temp_frame)) = 0;
    disp_data2A(:,:,it) = temp_frame;

    temp_frame = disp_data1B(:,:,it);
    temp_frame = medfilt2(temp_frame,izx_filtM);
    temp_frame = temp_frame/max(abs(temp_frame(:)));
    temp_frame(isnan(temp_frame)) = 0;
    disp_data1B(:,:,it) = temp_frame;

    temp_frame = disp_data2B(:,:,it);
    temp_frame = medfilt2(temp_frame,izx_filtM);
    temp_frame = temp_frame/max(abs(temp_frame(:)));
    temp_frame(isnan(temp_frame)) = 0;
    disp_data2B(:,:,it) = temp_frame;

end

% data versions before guided filtering (for comparison)
disp_data1A_pre = zeros(size(disp_data1A));     
disp_data2A_pre = zeros(size(disp_data2A));     
disp_data1B_pre = zeros(size(disp_data1B));
disp_data2B_pre = zeros(size(disp_data2B));

% data versions after guided filtering
disp_data1A_A = zeros(size(disp_data1A));       % filtered A guided by A
disp_data2A_A = zeros(size(disp_data2A));       % filtered A guided by A
disp_data1A_B = zeros(size(disp_data1A));       % filtered A guided by B
disp_data2A_B = zeros(size(disp_data2A));       % filtered A guided by B
disp_data1B_B = zeros(size(disp_data1B));       % filtered B guided by B
disp_data2B_B = zeros(size(disp_data2B));       % filtered B guided by B
disp_data1B_A = zeros(size(disp_data1B));       % filtered B guided by A
disp_data2B_A = zeros(size(disp_data2B));       % filtered B guided by A

% guided filter parameter: degree of smoothing parameter
degSmooth1A = max(abs(disp_data1A(:)));
degSmooth2A = max(abs(disp_data2A(:)));
degSmooth1B = max(abs(disp_data1B(:)));
degSmooth2B = max(abs(disp_data2B(:)));

% filter (x,t) frames
for iz = 1:Nz
    
    % RHS and LHS data slices
    s1A = squeeze(disp_data1A(iz,:,:));
    s1B = squeeze(disp_data1B(iz,:,:));
    s2A = flip(squeeze(disp_data2A(iz,:,:)),1);     % reverse x-axis
    s2B = flip(squeeze(disp_data2B(iz,:,:)),1);     % reverse x-axis

    % data A: reflection suppression + median filtering
    s1Afft = fftshift(fft2(s1A,NxFFT,NtFFT));
    s1Afft = s1Afft.*fftMask;
    s1Af = ifft2(ifftshift(s1Afft),'symmetric');
    s1Af = s1Af(1:Nx,1:Nt);
    s1Af = medfilt2(s1Af,ixt_filtG);

    s2Afft = fftshift(fft2(s2A,NxFFT,NtFFT));
    s2Afft = s2Afft.*fftMask;
    s2Af = ifft2(ifftshift(s2Afft),'symmetric');
    s2Af = s2Af(1:Nx,1:Nt);
    s2Af = medfilt2(s2Af,ixt_filtG);

    disp_data1A_pre(iz,:,:) = s1Af;
    disp_data2A_pre(iz,:,:) = s2Af;


    % data B: reflection suppression (?) + median filtering
    s1Bfft = fftshift(fft2(s1B,NxFFT,NtFFT));
    s1Bfft = s1Bfft.*fftMask;
    s1Bf = ifft2(ifftshift(s1Bfft),'symmetric');
    s1Bf = s1Bf(1:Nx,1:Nt);
    s1Bf = medfilt2(s1Bf,ixt_filtG);

    s2Bfft = fftshift(fft2(s2B,NxFFT,NtFFT));
    s2Bfft = s2Bfft.*fftMask;
    s2Bf = ifft2(ifftshift(s2Bfft),'symmetric');
    s2Bf = s2Bf(1:Nx,1:Nt);
    s2Bf = medfilt2(s2Bf,ixt_filtG);

    disp_data1B_pre(iz,:,:) = s1Bf;
    disp_data2B_pre(iz,:,:) = s2Bf;


%     % data A: guided filtering
%     disp_data1A_A(iz,:,:) = ...
%         imguidedfilter(s1Af,s1Af,'NeighborhoodSize',ixt_filtG,'DegreeOfSmoothing',degSmooth1A);
%     disp_data2A_A(iz,:,:) = ...
%         imguidedfilter(s2Af,s2Af,'NeighborhoodSize',ixt_filtG,'DegreeOfSmoothing',degSmooth2A);
%     disp_data1A_B(iz,:,:) = ...
%         imguidedfilter(s1Af,s1Bf,'NeighborhoodSize',ixt_filtG,'DegreeOfSmoothing',degSmooth1B);
%     disp_data2A_B(iz,:,:) = ...
%         imguidedfilter(s2Af,s2Bf,'NeighborhoodSize',ixt_filtG,'DegreeOfSmoothing',degSmooth2B);
% 
%     % data B: guided filtering
%     disp_data1B_B(iz,:,:) = ...
%         imguidedfilter(s1Bf,s1Bf,'NeighborhoodSize',ixt_filtG,'DegreeOfSmoothing',degSmooth1B);
%     disp_data2B_B(iz,:,:) = ...
%         imguidedfilter(s2Bf,s2Bf,'NeighborhoodSize',ixt_filtG,'DegreeOfSmoothing',degSmooth2B);
%     disp_data1B_A(iz,:,:) = ...
%         imguidedfilter(s1Bf,s1Af,'NeighborhoodSize',ixt_filtG,'DegreeOfSmoothing',degSmooth1A);
%     disp_data2B_A(iz,:,:) = ...
%         imguidedfilter(s2Bf,s2Af,'NeighborhoodSize',ixt_filtG,'DegreeOfSmoothing',degSmooth2A);
   
end
% IMPORTANT: 'data2' is now flipped (has x-axis reversed) to match 'data1'


% dataID = 'AB';
% %dataID = 'BA';
% %dataID = 'AA';
% %dataID = 'BB';
% 
% switch dataID
%     case 'AB'
%         disp_data1 = disp_data1A_B;
%         disp_data2 = disp_data2A_B;
%     case 'BA'
%         [~,~,disp_data1] = gradient(disp_data1B_A,dz,dx,dt);
%         [~,~,disp_data2] = gradient(disp_data2B_A,dz,dx,dt);
%     case 'AA'
%         disp_data1 = disp_data1A_A;
%         disp_data2 = disp_data2A_A;
%     case 'BB'
%         [~,~,disp_data1] = gradient(disp_data1B_B,dz,dx,dt);
%         [~,~,disp_data2] = gradient(disp_data2B_B,dz,dx,dt);
% end

dataID = 'noGF';
disp_data1 = disp_data1A_pre;
disp_data2 = disp_data2A_pre;


%% Visualize filtered data

filt_disp = uff.beamformed_data();
filt_disp.scan = uff.linear_scan();

filt_disp.scan.x_axis = x_axis;
filt_disp.scan.z_axis = z_axis;

filt_disp.data = reshape(disp_data1,[Nz*Nx,1,1,Nt]);
%filt_disp.data = reshape(disp_data1,[Nz*Nx,1,1,Nt]);

% display estimated displacements
fig_disp = figure;
filt_disp.plot(fig_disp,['Filtered Displacement'],[],'none');
colormap(gca(fig_disp),'gray');


%% Interpolate

s1 = zeros(Nz,Nx,NtInterp);
s2 = zeros(Nz,Nx,NtInterp);
for iz = 1:Nz
    s1f = squeeze(disp_data1(iz,:,:));
    s2f = squeeze(disp_data2(iz,:,:));
    itRange = (first_frame:Nt);
    for ix = 1:Nx
        s1(iz,ix,:) = interp1(t_axis(itRange),squeeze(s1f(ix,itRange)),t_axisInterp,'spline',0);
        s2(iz,ix,:) = interp1(t_axis(itRange),squeeze(s2f(ix,itRange)),t_axisInterp,'spline',0);
    end
end


%% Process displacement data

sws1 = zeros(Nz,Nx);
sws2 = zeros(Nz,Nx);

for iz = 1:Nz
    
    % RHS and LHS data slices
    s1I = squeeze(s1(iz,:,:));
    s2I = squeeze(s2(iz,:,:));

    for ix = 1:1:Nx

        ixRange = max(1,ix-ixExt2):min(ix+ixExt2,Nx);
        ixLen = length(ixRange);
        delta_x = (ixLen-1)*dx;

        itInterpRange = 1:NtInterp;
        itInterpLen = length(itInterpRange);

        % estimate SWS from RHS data
        template = squeeze(s1I(ixRange(1),itInterpRange));
        signal = squeeze(s1I(ixRange(end),itInterpRange));
        % cross-correlation lags: -(L-1),...,-1,0,+1,..,+(L-1)
        CC = normxcorr2(template,signal);  
        % find max starting from zero-lag (index L)
        [~,loc] = max(CC(itInterpLen+it_offset:end));   
        % compute SWS
        delta_t1 = (loc+it_offset-1)*dtInterp;
        sws1(iz,ix) = delta_x/delta_t1;

        % estimate SWS from LHS data
        template = squeeze(s2I(ixRange(1),itInterpRange));
        signal = squeeze(s2I(ixRange(end),itInterpRange));
        % cross-correlation lags: -(L-1),...,-1,0,+1,..,+(L-1)
        CC = normxcorr2(template,signal);  
        % find max starting from zero-lag (index L)
        [~,loc] = max(CC(itInterpLen+it_offset:end));   
        % compute SWS
        delta_t2 = (loc+it_offset-1)*dtInterp;
        sws2(iz,ix) = delta_x/delta_t2;

    end

end

sws2 = fliplr(sws2);


%% Combine and filter SWS maps

% weighted sum: (s1/(s1+s2))*s1 + (s2/(s1+s2))*s2
sws_map = (sws1.^2 + sws2.^2)./(sws1 + sws2);
sws_mapf = wiener2(sws_map,izx_filtW);


%% Display SWS maps

figure; imagesc(x_axis*1e3,z_axis*1e3,sws1);
axis image; colormap jet; colorbar; caxis([0,swsRange(2)]); 

figure; imagesc(x_axis*1e3,z_axis*1e3,sws2);
axis image; colormap jet; colorbar; caxis([0,swsRange(2)]); 

figure; imagesc(x_axis*1e3,z_axis*1e3,sws_mapf);
set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
axis image; shading flat; colormap jet; colorbar; caxis([0,swsRange(2)]); 
set(gca,'XTick',[-16 -8 0 8 16]); set(gca,'YTick',[5 15 25 35 45]);
xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
ylabel('z (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
saveas(gcf,sprintf('sws_%s_%s_[%d_%d].fig',caseID,dataID,round(swsRange(1)*10),round(swsRange(2)*10)));
%saveas(gcf,sprintf('sws_%s_%s_[%d_%d]',caseID,dataID,round(swsRange(1)*10),round(swsRange(2)*10)),'epsc');


%% Evaluate SWS maps

logFILE = fopen(sprintf('eval_%s_%s_[%d_%d].txt',...
    caseID,dataID,round(swsRange(1)*10),round(swsRange(2)*10)),'w');

% annotation and evaluation of Phantom I
if strcmp(phantomID,'ARF_phantom_I')
    mask_I = zeros(Nz,Nx);
    for ix = 1:Nx
        if (x_axis(ix) >= -15.42e-3) && (x_axis(ix) <= 15.42e-3)
            for iz = 1:Nz
                if (z_axis(iz) >= 4.0e-3) && (z_axis(iz) <= 46.51e-3)
                    mask_I(iz,ix) = 1;
                end
            end
        end
    end
    
    StDev = std(sws_mapf(mask_I>0))

    fprintf(logFILE,'\n\nStDev = %.4f\n\n',StDev);

    figure; imagesc(x_axis*1e3,z_axis*1e3,sws_mapf);
    set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
    axis image; shading flat; colormap jet; caxis([0,swsRange(2)]);
    set(gca,'XTick',[-16 -8 0 8 16]); set(gca,'YTick',[5 15 25 35 45]);
    xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
    ylabel('z (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
 
    rectangle('Position',[-15.42,4.0,30.84,42.51],...
        'EdgeColor',"k",'LineWidth',1,'LineStyle','--');

    saveas(gcf,sprintf('sws_%s_%s_[%d_%d]_marked',caseID,dataID,round(swsRange(1)*10),round(swsRange(2)*10)),'epsc');
    saveas(gcf,sprintf('sws_%s_%s_[%d_%d]_marked.fig',caseID,dataID,round(swsRange(1)*10),round(swsRange(2)*10)));

end


% annotation and evaluation of Phantom II
if strcmp(phantomID,'ARF_phantom_II')

    center_x = 0.327e-3;
    center_z = 16.478e-3;
    radius_i = 7.39e-3;
    radius_o1 = 9.39e-3;
    radius_o2 = 11.95e-3;

    mask_II_i = zeros(Nz,Nx);
    mask_II_o = zeros(Nz,Nx);
    for ix = 1:Nx
        dist = sqrt((x_axis(ix)-center_x).^2 + (z_axis-center_z).^2);
        for iz = 1:Nz
            if dist(iz) <= radius_i
                mask_II_i(iz,ix) = 1;
            elseif (dist(iz) >= radius_o1) && (dist(iz) <= radius_o2)
                mask_II_o(iz,ix) = 1;
            end            
        end
    end

    Mean_i = mean(sws_mapf(mask_II_i>0).^2);
    Mean_o = mean(sws_mapf(mask_II_o>0).^2);
    Var_i = var(sws_mapf(mask_II_i>0).^2);
    Var_o = var(sws_mapf(mask_II_o>0).^2);
    
    CNR = abs(Mean_i-Mean_o)/sqrt((Var_i+Var_o)/2)

    fprintf(logFILE,'\nMean / Var (inside)  = %.4f / %.4f\n',Mean_i,Var_i);
    fprintf(logFILE,'\nMean / Var (outside) = %.4f / %.4f\n',Mean_o,Var_o);
    fprintf(logFILE,'\nCNR = %.4f\n',CNR);

    figure; imagesc(x_axis*1e3,z_axis*1e3,sws_mapf);
    set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
    axis image; shading flat; colormap jet; caxis([0,swsRange(2)]);
    set(gca,'XTick',[-16 -8 0 8 16]); set(gca,'YTick',[5 15 25 35 45]);
    xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
    ylabel('z (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
 
    %viscircles([0.327,16.478],[7.39],'EnhanceVisibility',0,'Color',"w",'LineWidth',1,'LineStyle','--');
    %viscircles([0.327,16.478],[9.39],'EnhanceVisibility',0,'Color',"w",'LineWidth',1);
    %viscircles([0.327,16.478],[11.95],'EnhanceVisibility',0,'Color',"w",'LineWidth',1);

    annotation(gcf,'ellipse',[0.4321 0.5310 0.1875 0.2714],'LineWidth',1,'LineStyle','--'); 

    saveas(gcf,sprintf('sws_%s_%s_[%d_%d]_marked',caseID,dataID,round(swsRange(1)*10),round(swsRange(2)*10)),'epsc');
    saveas(gcf,sprintf('sws_%s_%s_[%d_%d]_marked.fig',caseID,dataID,round(swsRange(1)*10),round(swsRange(2)*10)));

end


% annotation and evaluation of Phantom III
if strcmp(phantomID,'ARF_phantom_III')

    center_x = 0.327e-3;
    center_z = 35.266e-3;
    radius_i = 8.36e-3;
    radius_o1 = 10.36e-3;
    radius_o2 = 13.31e-3;

    mask_III_i = zeros(Nz,Nx);
    mask_III_o = zeros(Nz,Nx);
    for ix = 1:Nx
        dist = sqrt((x_axis(ix)-center_x).^2 + (z_axis-center_z).^2);
        for iz = 1:Nz
            if dist(iz) <= radius_i
                mask_III_i(iz,ix) = 1;
            elseif (dist(iz) >= radius_o1) && (dist(iz) <= radius_o2)
                mask_III_o(iz,ix) = 1;
            end            
        end
    end

    Mean_i = mean(sws_mapf(mask_III_i>0).^2);
    Mean_o = mean(sws_mapf(mask_III_o>0).^2);
    Var_i = var(sws_mapf(mask_III_i>0).^2);
    Var_o = var(sws_mapf(mask_III_o>0).^2);
    
    CNR = abs(Mean_i-Mean_o)/sqrt((Var_i+Var_o)/2)

    fprintf(logFILE,'\nMean / Var (inside)  = %.4f / %.4f\n',Mean_i,Var_i);
    fprintf(logFILE,'\nMean / Var (outside) = %.4f / %.4f\n',Mean_o,Var_o);
    fprintf(logFILE,'\nCNR = %.4f\n',CNR);

        figure; imagesc(x_axis*1e3,z_axis*1e3,sws_mapf);
    set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
    axis image; shading flat; colormap jet; caxis([0,swsRange(2)]);
    set(gca,'XTick',[-16 -8 0 8 16]); set(gca,'YTick',[5 15 25 35 45]);
    xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
    ylabel('z (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
 
    %viscircles([0.327,35.266],[8.36],'EnhanceVisibility',0,'Color',"w",'LineWidth',1,'LineStyle','--');
    %viscircles([0.327,35.266],[10.36],'EnhanceVisibility',0,'Color',"w",'LineWidth',1);
    %viscircles([0.327,35.266],[13.31],'EnhanceVisibility',0,'Color',"w",'LineWidth',1);

    annotation(gcf,'ellipse',[0.3964 0.2262 0.2518 0.2667],'LineWidth',1,'LineStyle','--');

    saveas(gcf,sprintf('sws_%s_%s_[%d_%d]_marked',caseID,dataID,round(swsRange(1)*10),round(swsRange(2)*10)),'epsc');
    saveas(gcf,sprintf('sws_%s_%s_[%d_%d]_marked.fig',caseID,dataID,round(swsRange(1)*10),round(swsRange(2)*10)));

end

fclose(logFILE);
