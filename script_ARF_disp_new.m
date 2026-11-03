
%%% -- Modified by D.N.Rakhmatov, 27 September 2024 -- %%%

clear all; close all;
clc;

addpath(genpath('C:\Users\daler\Desktop\Reference_US\USTB'));

addpath(genpath('C:\Users\daler\Desktop\Reference_US\aMEM\ImpulseIQ'));


%% Read challenge data/parameters

%phantomID = 'ARF_phantom_I';
phantomID = 'ARF_phantom_II';
%phantomID = 'ARF_phantom_III';

load([phantomID '.mat']);

fsRF = 20*1e6;                  % RF sampling frequency (from specs)
first_frame = 6;                % first frame after push (from specs)
prepush_frame = 1;              % pre-push frame = index 1

Np = length(X_push);            % number of pushes
x_push = 1e-3*X_push;           % lateral push locations

c = SoS;
fs = 1e6*IQ_sample_rate;
fc = 1e6*IQ_demodulation_f;

f_ratio = IQ_demodulation_f/IQ_sample_rate;

lambda = c/fc;

dz = 1e-3*dZ;
dx = 1e-3*dX;
dt = 1/frame_rate;

Nz = length(Z_IQ);
Nx = length(X_IQ);
Nt = size(IQ_Data{1},3);        % number of frames

z_axis = 1e-3*Z_IQ.';
x_axis = 1e-3*X_IQ.';

disp_auto.dataA = cell(1,Np);
disp_auto.dataB = cell(1,Np);

disp_auto.Nz = Nz;
disp_auto.Nx = Nx;
disp_auto.Nt = Nt;
disp_auto.dz = dz;
disp_auto.dx = dx;
disp_auto.dt = dt;
disp_auto.z_axis = z_axis;
disp_auto.x_axis = x_axis;
disp_auto.x_push = x_push;
disp_auto.first_frame = first_frame;

disp_auto.c = c;
disp_auto.lambda = lambda;


%% Estimate displacement (Loupas autocorrelation method)

izExt2 = max(1,round(1.0*lambda/dz));   % axial two-sided extent
itExt1 = 1;                             % temporal one-sided extent

sws_ref = 5.0;                      % reference SWS [m/s]
x_tdist = sws_ref*dt*(itExt1+1);    % lateral travel distance at ref. SWS
fprintf('\nTemporal window --> Travel distance = %.3f mm at %.1f m/s\n\n',x_tdist*1000,sws_ref);

disp_auto.izExt2 = izExt2;
disp_auto.itExt1 = itExt1;

for ip = 1:Np

    IQ = IQ_Data{ip};

    dataA = zeros(Nz,Nx,Nt);
    dataB = zeros(Nz,Nx,Nt);

    for ix = 1:Nx

        IQframe = squeeze(IQ(:,ix,:));

        for it = first_frame:Nt

            itRange = it:min(it+itExt1,Nt);
            Nn = length(itRange);

            for iz = 1:Nz

                izRange = max(1,iz-izExt2):min(iz+izExt2,Nz);
                Nm = length(izRange);

                % data A: progressive referencing
                A = IQframe(izRange,itRange);
                Ac = conj(A);

                % Fixed x (z index = m, t index = n):
                % --> Im = I(m,n+1)Q(m,n) - I(m,n)Q(m,n+1)
                % --> Re = I(m,n)I(m,n+1) + Q(m,n)Q(m,n+1)
                numerA = 0;
                for m = 1:Nm
                    for n = 1:Nn-1
                        numerA = numerA + A(m,n).*Ac(m,n+1);
                    end
                end

                % Fixed x (z index = m, t index = n):
                % -- Im = I(m+1,n)Q(m,n) - I(m,n)Q(m+1,n)
                % -- Re = I(m,n)I(m+1,n) + Q(m,n)Q(m+1,n)
                denomA = 0;
                for n = 1:Nn
                    for m = 1:Nm-1
                        denomA = denomA + A(m,n).*Ac(m+1,n);
                    end
                end

                % data B: fixed (pre-push) referencing
                refVec = IQframe(izRange,prepush_frame);
                nowVec = A(:,1);  % current time --> column 1 of A
                B = sum(refVec).*conj(sum(nowVec));

                % estimated displacement data
                dataA(iz,ix,it) = (c/2)*(angle(numerA)/(2*pi*fc))./(1 + angle(denomA)/(2*pi*f_ratio));
                dataB(iz,ix,it) = (c/2)*(angle(B)/(2*pi*fc));

            end
        end
    end
    
    disp_auto.dataA{ip} = dataA;
    disp_auto.dataB{ip} = dataB;

end

save(sprintf('disp_auto_%s_z%dt%d.mat',phantomID,izExt2,itExt1), 'disp_auto');


%% Display estimated images for a given push (USTB)

est_disp = uff.beamformed_data();
est_disp.scan = uff.linear_scan();

est_disp.scan.x_axis = 1e-3*X_IQ.';
est_disp.scan.z_axis = 1e-3*Z_IQ.';

ip = 1;

%tmp = disp_auto.dataA{ip};
tmp = disp_auto.dataB{ip};
est_disp.data = reshape(tmp,[Nz*Nx,1,1,Nt]);


% display estimated displacements
fig_disp = figure;
est_disp.plot(fig_disp,['Estimated Displacement'],[],'none');
colormap(gca(fig_disp),'hot');
caxis([-1e-6 1e-6]);
%caxis([-lambda/2 lambda/2]);


% display B-mode image (first frame)
env = abs(IQ_Data{ip}(:,:,1));
bmode = 20*log10(env/max(env(:)));

figure; imagesc(x_axis*1000,z_axis*1000,bmode);
set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
ylabel('t (ms)','FontName','Courier','FontSize',14,'FontWeight','bold');
set(gca,'YTick',[5 15 25 35 45]); set(gca,'XTick',[-16 -8 0 8 16]);
set(gca,'YDir','reverse'); axis image; shading flat;
colormap gray; caxis([-60 0]);

saveas(gcf,sprintf('bmode_%s.fig',phantomID));
%saveas(gcf,sprintf('bmode_%s',phantomID),'epsc');
