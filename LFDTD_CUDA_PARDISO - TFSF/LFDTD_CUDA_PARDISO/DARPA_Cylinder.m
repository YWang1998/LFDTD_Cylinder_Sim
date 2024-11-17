%% 3-D Laguerre-FDTD Program
% Yifan Wang 
% School of Electrical and Computer Engineering
% Georgia Institute of Technology
% Atlanta, GA
%
% Description
% Simulation of DARPA defined cylinder structure with bottom apearture

clc;
clear;
fprintf('0. 3-D Laguerre-FDTD developed by Yifan Wang\n');

%% Define constants

eps0 = 1/36/pi*1e-9; % epsilon for vacuum
mu0 = 4*pi*1e-7; % mu for vacuum
ita0 = sqrt(mu0/eps0); % impedance for vacuum
v0 = 1/sqrt(eps0*mu0); % light speed in vacuum

%% Define simulation control parameters
dt = 10e-12; % time interval
tStep = 1000; % Number of time steps
s = 6e10; % Time scale factor (Needs to be carefully chosen)
tArray = 1:tStep;
tArray = tArray.*dt;
qStop = 150; % Number of basis functions used in the simulation
fac = 2; % Factor used to decrease the plot time step (1/fac*tStep)
recordNum = 4; % Number of Laguerre polynomials and basis functions used for plot
boundTypeIndex = 2; % Outmost boundary type: 1-PEC; 2-ABC; 3-PML;
TFSF_Update = 1; % TFSF Field Update
plotSec = 3; % 1-x cross section; 2-y cross section; 3-z cross section
plotLayer = 20; % Layer index

% Pulse related

pulseTypeIndex = 3; % 1-Gaussian; 2-Gaussian Derivative; 3-Modulated Gaussian

% Define plot Indicators
plotWave = 1; % Plot the generated source current waveform
plotMesh = 0; % Plot the mesh and locations of source and probes - require HUGE memory (use with caution)
plotMatrix = 0; % Plot the A matrix with spy(A)
plotLag = 0; % Plot the first n (recordNum) Laguerre polynomials and basis functions
plotRecJ = 0; % Plot the recoverd source current from q Laguerre basis function (require memory)
printQ = 1; % Print the Laguerre coefficient for the probe on the screen
plotProbe = 0; % Plot the time domain waveform of the probe voltage
plotFieldEx = 0; % Plot the Ex field
plotFieldEy = 0; % Plot the Ey field
plotFieldEz = 0; % Plot the Ez field
plotFieldHx = 0; % Plot the Hx field
plotFieldHy = 0; % Plot the Hy field
plotFieldHz = 0; % Plot the Hz field
PlotE_waveform = 0; % Plot 1D E field waveform
plot1DField = 0;

%% Define structurE

nx = 50; % Number of Cells in x direction
ny = 50; % Number of Cells in y direction
nz = 64; % Number of cells in z direction

nx_TFSF = 5;
ny_TFSF = 5;
nz_TFSF = 5;

dx = zeros(nx+1,1);
dy = zeros(ny+1,1);
dz = zeros(nz+1,1);

dx(:,:) = 0.009;
dy(:,:) = 0.009;
dz(:,:) = 0.01;

%% Define material - Cylinder Structure

eps = zeros(nx+1,ny+1,nz+1);
mu = zeros(nx+1,ny+1,nz+1);

eps(:,:,:) = eps0;
mu(:,:,:) = mu0;

nr = [15,14,2];
sigma = [3.56e+7,0,0];
cor = [1,1;1,1;1,1]; % offset for the bottom apearture location - [x_cor_start, x_cor_end; y_cor_start, y_cor_end; z_cor_start, z_cor_end]
[sigma_x,sigma_y,sigma_z] = ConformalGrid_3D_ML(nx, ny, nz, dx(1), dy(1), nr, sigma, cor);

%% Define Parameters (material related)

cex = zeros(nx+1,ny+1,nz+1);
cey = zeros(nx+1,ny+1,nz+1);
cez = zeros(nx+1,ny+1,nz+1);

chx = zeros(nx+1,ny+1,nz+1);
chy = zeros(nx+1,ny+1,nz+1);
chz = zeros(nx+1,ny+1,nz+1);

dxe = dx;
dye = dy;
dze = dz;

dxh = zeros(nx+1,1);
dyh = zeros(ny+1,1);
dzh = zeros(nz+1,1);


for i=1:nx+1
    
    if i==1
        dxh(i) = dx(i,1);
    else
        dxh(i) = (dx(i,1)+dx(i-1,1))/2;
    end
    
end

for i=1:ny+1
    
    if i==1
        dyh(i) = dy(i,1);
    else
        dyh(i) = (dy(i,1)+dy(i-1,1))/2;
    end
    
end

for i=1:nz+1
    
    if i==1
        dzh(i) = dz(i,1);
    else
        dzh(i) = (dz(i,1)+dz(i-1,1))/2;
    end
    
end

for i=1:nx+1
    for j=1:ny+1
        for k=1:nz+1
            
            cex(i,j,k) = 2/(s*eps(i,j,k)*dxh(i,1));
            cey(i,j,k) = 2/(s*eps(i,j,k)*dyh(j,1));
            cez(i,j,k) = 2/(s*eps(i,j,k)*dzh(k,1));
            chx(i,j,k) = 2/(s*mu(i,j,k)*dxe(i,1));
            chy(i,j,k) = 2/(s*mu(i,j,k)*dye(j,1));
            chz(i,j,k) = 2/(s*mu(i,j,k)*dze(k,1));
            
        end
    end
end

% Put a break point here and run the MatrixExport.m function if you want to
% convert the Cylinder Geometrical data into C++ code

%% Define source

% Define source location

% jCellIndex = [1+npml,nx-npml,1+npml,ny-npml,2+npml,2+npml]; % Indicate the Cells which contains the current source [xmin,xmax,ymin,ymax,zmin,zmax]
jCellIndex = [nx_TFSF+1,nx-nx_TFSF,ny_TFSF+1,ny+1-ny_TFSF,nz_TFSF+1,nz+1-nz_TFSF];
jCellIndex_Z = [nx_TFSF+1,nx+1-nx_TFSF,ny_TFSF+1,ny+1-ny_TFSF,nz_TFSF+1,nz-nz_TFSF];
jDirecIndex = [1,0,1]; % Indicate the direction of the current source [x direction,y direction,z direction]

mCellIndex = [nx_TFSF+1,nx-nx_TFSF,ny_TFSF+1,ny+1-ny_TFSF,nz_TFSF,nz-nz_TFSF+1];
mCellIndex_Z = [nx_TFSF+1,nx-nx_TFSF,ny_TFSF,ny+1-ny_TFSF,nz_TFSF+1,nz-nz_TFSF+1];
mDirecIndex = [0,1,1]; % Indicate the direction of the magnetic source [x direction,y direction,z direction]

% Current Source
Jx = zeros(nx+1,ny+1,nz+1);
Jy = zeros(nx+1,ny+1,nz+1);
Jz = zeros(nx+1,ny+1,nz+1);

% Magnetic Source
Mx = zeros(nx,ny,nz);
My = zeros(nx,ny,nz);
Mz = zeros(nx,ny,nz);


%% Define Pulse

fc = 4e9; % central frequency
mf = 1.5; % factor
td = 1/(mf*fc); % modulation determination
tc = 2.0e-9; % central time
omega = 2*pi*fc;

waveform = zeros(1,tStep);

switch pulseTypeIndex
    case 1
        for i=1:tStep
            waveform(i) = exp(-((i-tc)^2/td^2));
        end
    case 2
        for i=1:tStep
            waveform(i)=-2*(dt*i-tc)/td*exp(1)^(-((dt*i-tc)/td)^2);
        end
    case 3
        for i=1:tStep
            waveform(i)=sin(2*pi*fc*(dt*i-tc))*exp(1)^(-((dt*i-tc)/td)^2);
        end
    case 4
        for i=1:tStep
            waveform(i)=sin(2*pi*fc*(dt*(i-1)));
        end
    otherwise
        fprintf('Program terminated, please choose correct waveform\n');
        return;
end

% Frequency characteristic of the waveform

fs = 1/dt;
fx = fft(waveform);
df = fs/tStep;
nn = 0:tStep/2;
fn = nn*df/(1e+9);

if plotWave==1
    figure;
    subplot(2,1,1);
    plot(tArray/10^(-9),waveform);
    title('Time Domain');
    xlabel('Time (ns)');
    ylabel('Magnitude');
    grid;
    subplot(2,1,2);
    plot(fn,abs(fx(nn+1))*2/tStep); 
    xlim([0 5]);
    title('Frequency Domain');
    xlabel('Frequency (GHz)');
    ylabel('Magnitude');
    grid;
end

%% Define Probe

numProbe = 2;
probe = zeros(tStep,numProbe);
probeDirecIndex = [0,0,1]; % Indicate the direction of the probe E field [x direction,y direction,z direction]
probeCellIndex = zeros(6,numProbe);

% Probe 1
probeCellIndex(:,1) = [3,3,3,3,3,4]; % [xmin,xmax,ymin,ymax,zmin,zmax]

% Probe 2
probeCellIndex(:,2) = [6,6,6,6,3,4]; % [xmin,xmax,ymin,ymax,zmin,zmax]

%% Plot mesh, source and probe

if plotMesh==1
    
    LagPos.Plot_Mesh(nx,ny,nz,lx,ly,lz,dx,dy,dz,eps,sigma,eps0);

    return;
    
end

%% Numbering nodes

nodeNum = zeros(nx+1,ny+1,nz+1,3);
nodeTypeIndex = zeros(nx+1,ny+1,nz+1,3); % 0-Ordinary node; 1-Special node (inner boundary)
nnode = 1;

for i=1:nx+1
    for j=1:ny+1
        for k=1:nz+1

            if i==nx+1
                nodeNum(i,j,k,1) = 0;
            else
                nodeNum(i,j,k,1) = nnode;
                nnode = nnode+1;
            end

            if j==ny+1
                nodeNum(i,j,k,2) = 0;
            else
                nodeNum(i,j,k,2) = nnode;
                nnode = nnode+1;
            end

            if k==nz+1
                nodeNum(i,j,k,3) = 0;
            else
                nodeNum(i,j,k,3) = nnode;
                nnode = nnode+1;
            end

        end
    end
end

if TFSF_Update == 1
    nodeNum_1D = zeros(nz-nz_TFSF+2,1);
    nodeNum_1D(:) =(1:nz-nz_TFSF+2)';
    nnode_1D = nz-nz_TFSF+2;
end

nnode = nnode-1;

fprintf('1. Number of total nodes is %d\n',nnode);

%% Build A matrix

index = 0; % Final index should equal to number of non zero terms in A matrix
nmat = zeros(3,nnode*13); % Store the node x and y index and node value, size not exact nnode*13
imat = [1:nnode;1:nnode;ones(1,nnode)];

if TFSF_Update == 1
    
    index_1D = 0;
    nmat_1D = zeros(3,nnode_1D*13);
    
    % 1D Ex equation
    
    for i=1:1
        for j=1:1
            for k = 2:nz-nz_TFSF+1
                
                ex111 = nodeNum_1D(k);
                ex112 = nodeNum_1D(k+1);
                ex110 = nodeNum_1D(k-1);
                
                index_1D = index_1D+1; nmat_1D(:,index_1D) = [ex111,ex111,1+cez(i,j,k)*chz(i,j,k)+cez(i,j,k)*chz(i,j,k-1)];
                index_1D = index_1D+1; nmat_1D(:,index_1D) = [ex111,ex112,-cez(i,j,k)*chz(i,j,k)];
                index_1D = index_1D+1; nmat_1D(:,index_1D) = [ex111,ex110,-cez(i,j,k)*chz(i,j,k-1)];
                
            end
        end
    end
    
    % Outmost ABC boundary for Ex at z = Z
    
    for i=1:1
        for j = 1:1
            for k = nz-nz_TFSF+2:nz-nz_TFSF+2
            
                ex111 = nodeNum_1D(k);
                ex110 = nodeNum_1D(k-1);
                
                index_1D = index_1D+1; nmat_1D(:,index_1D) = [ex111,ex111,1/dze(k-1)+s/(4*v0)];
                index_1D = index_1D+1; nmat_1D(:,index_1D) = [ex111,ex110,-1/dze(k-1)+s/(4*v0)];
            
            end
        end
    end
    
    % Outmost TF/SF boundary for Ex
    
    for i=1:1
        for j=1:1
            for k = 1:1
                
                ex111 = nodeNum_1D(k);
                
                index_1D = index_1D+1; nmat_1D(:,index_1D) = [ex111,ex111,1];
                
            end
        end
    end
    
    nmat_1D = nmat_1D(:,1:index_1D);
    A_1D = sparse(nmat_1D(1,:),nmat_1D(2,:),nmat_1D(3,:),nnode_1D,nnode_1D);
    x_1D = zeros(nnode_1D,1);
    b_1D = zeros(nnode_1D,1);
    sumE_1D = zeros(nnode_1D,1);
    hy_1D = zeros(1,1,nz-nz_TFSF+1);
    sumHy_1D = zeros(1,1,nz-nz_TFSF+1);
end


% Field components except inner boundary and outmost boundary
    
% Ex equation
for i=1:nx
    for j=2:ny
        for k=2:nz

            ex111 = nodeNum(i,j,k,1);

            ex121 = nodeNum(i,j+1,k,1);
            ey211 = nodeNum(i+1,j,k,2);
            ey111 = nodeNum(i,j,k,2);

            ex101 = nodeNum(i,j-1,k,1);
            ey201 = nodeNum(i+1,j-1,k,2);
            ey101 = nodeNum(i,j-1,k,2);

            ez211 = nodeNum(i+1,j,k,3);
            ez111 = nodeNum(i,j,k,3);
            ex112 = nodeNum(i,j,k+1,1);

            ez210 = nodeNum(i+1,j,k-1,3);
            ez110 = nodeNum(i,j,k-1,3);
            ex110 = nodeNum(i,j,k-1,1);

            index = index+1; nmat(:,index) = [ex111,ex111,1+cey(i,j,k)*chy(i,j,k)+cey(i,j,k)*chy(i,j-1,k)+cez(i,j,k)*chz(i,j,k)+cez(i,j,k)*chz(i,j,k-1)+2*sigma_x(i,j,k)/(s*eps(i,j,k))];

            index = index+1; nmat(:,index) = [ex111,ex121,-cey(i,j,k)*chy(i,j,k)];
            index = index+1; nmat(:,index) = [ex111,ey211,cey(i,j,k)*chx(i,j,k)];
            index = index+1; nmat(:,index) = [ex111,ey111,-cey(i,j,k)*chx(i,j,k)];

            index = index+1; nmat(:,index) = [ex111,ex101,-cey(i,j,k)*chy(i,j-1,k)];
            index = index+1; nmat(:,index) = [ex111,ey201,-cey(i,j,k)*chx(i,j-1,k)];
            index = index+1; nmat(:,index) = [ex111,ey101,cey(i,j,k)*chx(i,j-1,k)];

            index = index+1; nmat(:,index) = [ex111,ez211,cez(i,j,k)*chx(i,j,k)];
            index = index+1; nmat(:,index) = [ex111,ez111,-cez(i,j,k)*chx(i,j,k)];
            index = index+1; nmat(:,index) = [ex111,ex112,-cez(i,j,k)*chz(i,j,k)];

            index = index+1; nmat(:,index) = [ex111,ez210,-cez(i,j,k)*chx(i,j,k-1)];
            index = index+1; nmat(:,index) = [ex111,ez110,cez(i,j,k)*chx(i,j,k-1)];
            index = index+1; nmat(:,index) = [ex111,ex110,-cez(i,j,k)*chz(i,j,k-1)];
        end
    end
end

% Ey equation
for i=2:nx
    for j=1:ny
        for k=2:nz

            ey111 = nodeNum(i,j,k,2);

            ey112 = nodeNum(i,j,k+1,2);
            ez121 = nodeNum(i,j+1,k,3);
            ez111 = nodeNum(i,j,k,3);

            ey110 = nodeNum(i,j,k-1,2);
            ez120 = nodeNum(i,j+1,k-1,3);
            ez110 = nodeNum(i,j,k-1,3);

            ex121 = nodeNum(i,j+1,k,1);
            ex111 = nodeNum(i,j,k,1);
            ey211 = nodeNum(i+1,j,k,2);

            ex021 = nodeNum(i-1,j+1,k,1);
            ex011 = nodeNum(i-1,j,k,1);
            ey011 = nodeNum(i-1,j,k,2);

            index = index+1; nmat(:,index) = [ey111,ey111,1+cez(i,j,k)*chz(i,j,k)+cez(i,j,k)*chz(i,j,k-1)+cex(i,j,k)*chx(i,j,k)+cex(i,j,k)*chx(i-1,j,k)+2*sigma_y(i,j,k)/(s*eps(i,j,k))];

            index = index+1; nmat(:,index) = [ey111,ey112,-cez(i,j,k)*chz(i,j,k)];
            index = index+1; nmat(:,index) = [ey111,ez121,cez(i,j,k)*chy(i,j,k)];
            index = index+1; nmat(:,index) = [ey111,ez111,-cez(i,j,k)*chy(i,j,k)];

            index = index+1; nmat(:,index) = [ey111,ey110,-cez(i,j,k)*chz(i,j,k-1)];
            index = index+1; nmat(:,index) = [ey111,ez120,-cez(i,j,k)*chy(i,j,k-1)];
            index = index+1; nmat(:,index) = [ey111,ez110,cez(i,j,k)*chy(i,j,k-1)];

            index = index+1; nmat(:,index) = [ey111,ex121,cex(i,j,k)*chy(i,j,k)];
            index = index+1; nmat(:,index) = [ey111,ex111,-cex(i,j,k)*chy(i,j,k)];
            index = index+1; nmat(:,index) = [ey111,ey211,-cex(i,j,k)*chx(i,j,k)];

            index = index+1; nmat(:,index) = [ey111,ex021,-cex(i,j,k)*chy(i-1,j,k)];
            index = index+1; nmat(:,index) = [ey111,ex011,cex(i,j,k)*chy(i-1,j,k)];
            index = index+1; nmat(:,index) = [ey111,ey011,-cex(i,j,k)*chx(i-1,j,k)];
        end
    end
end

% Ez equation
for i=2:nx
    for j=2:ny
        for k=1:nz

            ez111 = nodeNum(i,j,k,3);

            ez211 = nodeNum(i+1,j,k,3);
            ex112 = nodeNum(i,j,k+1,1);
            ex111 = nodeNum(i,j,k,1);

            ez011 = nodeNum(i-1,j,k,3);
            ex012 = nodeNum(i-1,j,k+1,1);
            ex011 = nodeNum(i-1,j,k,1);

            ey112 = nodeNum(i,j,k+1,2);
            ey111 = nodeNum(i,j,k,2);
            ez121 = nodeNum(i,j+1,k,3);

            ey102 = nodeNum(i,j-1,k+1,2);
            ey101 = nodeNum(i,j-1,k,2);
            ez101 = nodeNum(i,j-1,k,3);

            index = index+1; nmat(:,index) = [ez111,ez111,1+cex(i,j,k)*chx(i,j,k)+cex(i,j,k)*chx(i-1,j,k)+cey(i,j,k)*chy(i,j,k)+cey(i,j,k)*chy(i,j-1,k)+2*sigma_z(i,j,k)/(s*eps(i,j,k))];

            index = index+1; nmat(:,index) = [ez111,ez211,-cex(i,j,k)*chx(i,j,k)];
            index = index+1; nmat(:,index) = [ez111,ex112,cex(i,j,k)*chz(i,j,k)];
            index = index+1; nmat(:,index) = [ez111,ex111,-cex(i,j,k)*chz(i,j,k)];

            index = index+1; nmat(:,index) = [ez111,ez011,-cex(i,j,k)*chx(i-1,j,k)];
            index = index+1; nmat(:,index) = [ez111,ex012,-cex(i,j,k)*chz(i-1,j,k)];
            index = index+1; nmat(:,index) = [ez111,ex011,cex(i,j,k)*chz(i-1,j,k)];

            index = index+1; nmat(:,index) = [ez111,ey112,cey(i,j,k)*chz(i,j,k)];
            index = index+1; nmat(:,index) = [ez111,ey111,-cey(i,j,k)*chz(i,j,k)];
            index = index+1; nmat(:,index) = [ez111,ez121,-cey(i,j,k)*chy(i,j,k)];

            index = index+1; nmat(:,index) = [ez111,ey102,-cey(i,j,k)*chz(i,j-1,k)];
            index = index+1; nmat(:,index) = [ez111,ey101,cey(i,j,k)*chz(i,j-1,k)];
            index = index+1; nmat(:,index) = [ez111,ez101,-cey(i,j,k)*chy(i,j-1,k)];
        end
    end
end




% Outmost boundary

if boundTypeIndex==1
    
    % Outmost PEC boundary for Ex
    
    for i=1:nx
        
        % Edge
        
        for j=1:ny:ny+1
            for k=1:nz:nz+1 
                ex111 = nodeNum(i,j,k,1);
                index = index+1; nmat(:,index) = [ex111,ex111,1];    
            end
        end
        
        % Face
        
        for j=1:ny:ny+1
            for k=2:nz
                ex111 = nodeNum(i,j,k,1);
                index = index+1; nmat(:,index) = [ex111,ex111,1];
            end
        end
        
        for j=2:ny
            for k=1:nz:nz+1
                ex111 = nodeNum(i,j,k,1);
                index = index+1; nmat(:,index) = [ex111,ex111,1];
            end
        end
        
    end
    
    % Outmost PEC boundary for Ey
    
    for j=1:ny
        
        % Edge
        
        for i=1:nx:nx+1
            for k=1:nz:nz+1
                ey111 = nodeNum(i,j,k,2);
                index = index+1; nmat(:,index) = [ey111,ey111,1];
            end
        end
        
        % Face
        
        for i=1:nx:nx+1
            for k=2:nz
                ey111 = nodeNum(i,j,k,2);
                index = index+1; nmat(:,index) = [ey111,ey111,1];
            end
        end
        
        for i=2:nx
            for k=1:nz:nz+1
                ey111 = nodeNum(i,j,k,2);
                index = index+1; nmat(:,index) = [ey111,ey111,1];
            end
        end
        
    end
    
    % Outmost PEC boundary for Ez
    
    for k=1:nz
        
        % Edge
        
        for i=1:nx:nx+1
            for j=1:ny:ny+1             
                ez111 = nodeNum(i,j,k,3);
                index = index+1; nmat(:,index) = [ez111,ez111,1];                
            end
        end
        
        % Face
        
        for i=1:nx:nx+1
            for j=2:ny            
                ez111 = nodeNum(i,j,k,3);
                index = index+1; nmat(:,index) = [ez111,ez111,1];                
            end
        end
        
        for i=2:nx
            for j=1:ny:ny+1             
                ez111 = nodeNum(i,j,k,3);
                index = index+1; nmat(:,index) = [ez111,ez111,1];                
            end
        end
        
    end
    
elseif boundTypeIndex==2 % ABC boundary
    
    % Outmost ABC boundary for Ex
    
    for i=1:nx
        
        % Edge
        
        for j=1:ny:ny+1
            for k=1:nz:nz+1

                if j==1 && k==1 % Case (a-1)
                    ex111 = nodeNum(i,j,k,1);
                    ex122 = nodeNum(i,j+1,k+1,1);
                    index = index+1; nmat(:,index) = [ex111,ex111,1/(sqrt(dye(j)^2+dze(k)^2))+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ex111,ex122,-1/(sqrt(dye(j)^2+dze(k)^2))+s/(4*v0)];                                          
                elseif j==1 && k==nz+1 % Case (a-2)
                    ex111 = nodeNum(i,j,k,1);
                    ex120 = nodeNum(i,j+1,k-1,1);
                    index = index+1; nmat(:,index) = [ex111,ex111,1/(sqrt(dye(j)^2+dze(k-1)^2))+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ex111,ex120,-1/(sqrt(dye(j)^2+dze(k-1)^2))+s/(4*v0)];                       
                elseif j==ny+1 && k==1 % Case (a-3)
                    ex111 = nodeNum(i,j,k,1);
                    ex102 = nodeNum(i,j-1,k+1,1);
                    index = index+1; nmat(:,index) = [ex111,ex111,1/(sqrt(dye(j-1)^2+dze(k)^2))+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ex111,ex102,-1/(sqrt(dye(j-1)^2+dze(k)^2))+s/(4*v0)];
                elseif j==ny+1 && k==nz+1 % Case (a-4)    
                    ex111 = nodeNum(i,j,k,1);
                    ex100 = nodeNum(i,j-1,k-1,1);
                    index = index+1; nmat(:,index) = [ex111,ex111,1/(sqrt(dye(j-1)^2+dze(k-1)^2))+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ex111,ex100,-1/(sqrt(dye(j-1)^2+dze(k-1)^2))+s/(4*v0)];
                end

            end
        end
        
        % Face
            
        for j=1:ny:ny+1
            for k=2:nz

                if j==1 % Case (1-1)              
                    ex111 = nodeNum(i,j,k,1);
                    ex121 = nodeNum(i,j+1,k,1);
                    index = index+1; nmat(:,index) = [ex111,ex111,1/dye(j)+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ex111,ex121,-1/dye(j)+s/(4*v0)];
                elseif j==ny+1 % Case (1-2)
                    ex111 = nodeNum(i,j,k,1);
                    ex101 = nodeNum(i,j-1,k,1);
                    index = index+1; nmat(:,index) = [ex111,ex111,1/dye(j-1)+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ex111,ex101,-1/dye(j-1)+s/(4*v0)];
                end

            end
        end

        for j=2:ny
            for k=1:nz:nz+1

                if k==1 % Case (1-3)
                    ex111 = nodeNum(i,j,k,1);
                    ex112 = nodeNum(i,j,k+1,1);
                    index = index+1; nmat(:,index) = [ex111,ex111,1/dze(k)+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ex111,ex112,-1/dze(k)+s/(4*v0)];
                elseif k==nz+1 % Case (1-4)
                    ex111 = nodeNum(i,j,k,1);
                    ex110 = nodeNum(i,j,k-1,1);
                    index = index+1; nmat(:,index) = [ex111,ex111,1/dze(k-1)+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ex111,ex110,-1/dze(k-1)+s/(4*v0)];
                end

            end
        end
        
    end
    
    % Outmost ABC boundary for Ey
                    
    for j=1:ny

        % Edge

        for i=1:nx:nx+1
            for k=1:nz:nz+1

                if k==1 && i==1 % Case (b-1)
                    ey111 = nodeNum(i,j,k,2);
                    ey212 = nodeNum(i+1,j,k+1,2);
                    index = index+1; nmat(:,index) = [ey111,ey111,1/(sqrt(dze(k)^2+dxe(i)^2))+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ey111,ey212,-1/(sqrt(dze(k)^2+dxe(i)^2))+s/(4*v0)];
                elseif k==1 && i==nx+1 % Case (b-2)
                    ey111 = nodeNum(i,j,k,2);
                    ey012 = nodeNum(i-1,j,k+1,2);
                    index = index+1; nmat(:,index) = [ey111,ey111,1/(sqrt(dze(k)^2+dxe(i-1)^2))+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ey111,ey012,-1/(sqrt(dze(k)^2+dxe(i-1)^2))+s/(4*v0)];
                elseif k==nz+1 && i==1 % Case (b-3)
                    ey111 = nodeNum(i,j,k,2);
                    ey210 = nodeNum(i+1,j,k-1,2);
                    index = index+1; nmat(:,index) = [ey111,ey111,1/(sqrt(dze(k-1)^2+dxe(i)^2))+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ey111,ey210,-1/(sqrt(dze(k-1)^2+dxe(i)^2))+s/(4*v0)];
                elseif k==nz+1 && i==nx+1 % Case (b-4)
                    ey111 = nodeNum(i,j,k,2);
                    ey010 = nodeNum(i-1,j,k-1,2);
                    index = index+1; nmat(:,index) = [ey111,ey111,1/(sqrt(dze(k-1)^2+dxe(i-1)^2))+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ey111,ey010,-1/(sqrt(dze(k-1)^2+dxe(i-1)^2))+s/(4*v0)];
                end

            end
        end

        % Face

        for i=2:nx
            for k=1:nz:nz+1

                if k==1 % Case (2-1)
                    ey111 = nodeNum(i,j,k,2);
                    ey112 = nodeNum(i,j,k+1,2);
                    index = index+1; nmat(:,index) = [ey111,ey111,1/dze(k)+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ey111,ey112,-1/dze(k)+s/(4*v0)];
                elseif k==nz+1 % Case (2-2)
                    ey111 = nodeNum(i,j,k,2);
                    ey110 = nodeNum(i,j,k-1,2);
                    index = index+1; nmat(:,index) = [ey111,ey111,1/dze(k-1)+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ey111,ey110,-1/dze(k-1)+s/(4*v0)];
                end

            end
        end

        for i=1:nx:nx+1
            for k=2:nz

                if i==1 % Case (2-3)
                    ey111 = nodeNum(i,j,k,2);
                    ey211 = nodeNum(i+1,j,k,2);
                    index = index+1; nmat(:,index) = [ey111,ey111,1/dxe(i)+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ey111,ey211,-1/dxe(i)+s/(4*v0)];
                elseif i==nx+1 % Case (2-4)
                    ey111 = nodeNum(i,j,k,2);
                    ey011 = nodeNum(i-1,j,k,2);
                    index = index+1; nmat(:,index) = [ey111,ey111,1/dxe(i-1)+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ey111,ey011,-1/dxe(i-1)+s/(4*v0)];
                end

            end
        end

    end

    % Outmost ABC boundary for Ez

    for k=1:nz

        % Edge

        for i=1:nx:nx+1
            for j=1:ny:ny+1

                if i==1 && j==1 % Case (c-1)
                    ez111 = nodeNum(i,j,k,3);
                    ez221 = nodeNum(i+1,j+1,k,3);
                    index = index+1; nmat(:,index) = [ez111,ez111,1/(sqrt(dxe(i)^2+dye(j)^2))+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ez111,ez221,-1/(sqrt(dxe(i)^2+dye(j)^2))+s/(4*v0)];
                elseif i==1 && j==ny+1 % Case (c-2)
                    ez111 = nodeNum(i,j,k,3);
                    ez201 = nodeNum(i+1,j-1,k,3);
                    index = index+1; nmat(:,index) = [ez111,ez111,1/(sqrt(dxe(i)^2+dye(j-1)^2))+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ez111,ez201,-1/(sqrt(dxe(i)^2+dye(j-1)^2))+s/(4*v0)];
                elseif i==nx+1 && j==1 % Case (c-3)
                    ez111 = nodeNum(i,j,k,3);
                    ez021 = nodeNum(i-1,j+1,k,3);
                    index = index+1; nmat(:,index) = [ez111,ez111,1/(sqrt(dxe(i-1)^2+dye(j)^2))+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ez111,ez021,-1/(sqrt(dxe(i-1)^2+dye(j)^2))+s/(4*v0)];
                elseif i==nx+1 && j==ny+1 % Case (c-4)
                    ez111 = nodeNum(i,j,k,3);
                    ez001 = nodeNum(i-1,j-1,k,3);
                    index = index+1; nmat(:,index) = [ez111,ez111,1/(sqrt(dxe(i-1)^2+dye(j-1)^2))+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ez111,ez001,-1/(sqrt(dxe(i-1)^2+dye(j-1)^2))+s/(4*v0)];
                end

            end
        end

        % Face

        for i=1:nx:nx+1
            for j=2:ny

                if i==1 % Case (3-1)
                    ez111 = nodeNum(i,j,k,3);
                    ez211 = nodeNum(i+1,j,k,3);
                    index = index+1; nmat(:,index) = [ez111,ez111,1/dxe(i)+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ez111,ez211,-1/dxe(i)+s/(4*v0)];
                elseif i==nx+1 % Case (3-2)
                    ez111 = nodeNum(i,j,k,3);
                    ez011 = nodeNum(i-1,j,k,3);
                    index = index+1; nmat(:,index) = [ez111,ez111,1/dxe(i-1)+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ez111,ez011,-1/dxe(i-1)+s/(4*v0)];
                end

            end
        end

        for i=2:nx
            for j=1:ny:ny+1

                if j==1 % Case (3-3)
                    ez111 = nodeNum(i,j,k,3);
                    ez121 = nodeNum(i,j+1,k,3);
                    index = index+1; nmat(:,index) = [ez111,ez111,1/dye(j)+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ez111,ez121,-1/dye(j)+s/(4*v0)];
                elseif j==ny+1 % Case (3-4)
                    ez111 = nodeNum(i,j,k,3);
                    ez101 = nodeNum(i,j-1,k,3);
                    index = index+1; nmat(:,index) = [ez111,ez111,1/dye(j-1)+s/(4*v0)];
                    index = index+1; nmat(:,index) = [ez111,ez101,-1/dye(j-1)+s/(4*v0)];
                end

            end
        end

    end

else
    
    fprintf('Boundary setting error!\n');
    return;
    
end

% Rescale nmat

nmat = nmat(:,1:index);

% Transfer nmat into sparse A matrix

A = sparse(nmat(1,:),nmat(2,:),nmat(3,:),nnode,nnode);
fprintf('3. Number of non-zero terms in A matrix is %d\n',index);

% Spy A matrix

if plotMatrix==1
    figure;
    spy(A);
end

tol = 1e-6;
maxit = 40;

D = diag(diag(A));

%% Iteratively solve the problem

lagPoly = zeros(4,tStep);
q = 0;
scale = 230/(s*(tStep-1)*dt/2); % Similar as in Myunghyun's paper

recordLagPoly = zeros(recordNum+1,tStep);
recordLagBasisFunc = zeros(recordNum+1,tStep);

if plotRecJ==1
    recordJq = zeros(qStop+1,1);
    recordJqPhi = zeros(qStop+1,tStep);
    recordRecJ = zeros(tStep,1);
end

if plotFieldEx==1 || plotFieldEy == 1 || plotFieldEz == 1
    plotStep = floor(tStep/fac);
    recordE = zeros(nnode,plotStep);
    recordHx = zeros(nx,ny,nz,plotStep);
    recordHy = zeros(nx,ny,nz,plotStep);
    recordHz = zeros(nx,ny,nz,plotStep);
    recordHy_1D = zeros(1,1,nz-nz_TFSF+1,plotStep);
end

sumE = zeros(nnode,1);
hx = zeros(nx,ny,nz);
hy = zeros(nx,ny,nz);
hz = zeros(nx,ny,nz);
sumHx = zeros(nx,ny,nz);
sumHy = zeros(nx,ny,nz);
sumHz = zeros(nx,ny,nz);

b = zeros(nnode,1);
x = zeros(nnode,1);

fprintf('4. Begin iterative calculation...\n');

tic;
while q<=qStop
    % Calculate Laguerre polynomials
    
    if q==0
        
        for i=1:tStep
            lagPoly(3,i) = 1.0;
            lagPoly(4,i) = -s*i*dt/2;
        end
        
    elseif q==1
        
        for i=1:tStep
            lagPoly(2,i) = lagPoly(3,i);
            lagPoly(3,i) = 1-s*i*dt;
        end
        
    else
        
        for i=1:tStep
            lagPoly(1,i) = lagPoly(2,i);
            lagPoly(2,i) = lagPoly(3,i);
        end
        
        for i=1:tStep
            lagPoly(3,i) = (1/q)*((2*q-1-s*i*dt)*lagPoly(2,i)-(q-1)*lagPoly(1,i));
        end
        
        for i=1:tStep
            if lagPoly(3,i)>1e100 % Make sure that Laguerre polynomial does not go to infinity
                lagPoly(1,i) = lagPoly(1,i)*exp(-s*i*dt/2*scale);
                lagPoly(2,i) = lagPoly(2,i)*exp(-s*i*dt/2*scale);
                lagPoly(3,i) = lagPoly(3,i)*exp(-s*i*dt/2*scale);
                lagPoly(4,i) = lagPoly(4,i)+s*i*dt/2*scale;
            end
        end
    end
    
    % Record first n laguerre polynomials
    
    if plotLag==1

        if q<=recordNum
        
            recordLagPoly(q+1,:) = lagPoly(3,:);
            recordLagBasisFunc(q+1,:) = lagPoly(3,:).*exp(lagPoly(4,:));
            
            if (q==recordNum)
                figure;
                for i=1:recordNum+1
                    plot(tArray*s,recordLagPoly(i,:));
                    hold on;
                end
                title('Laguerre Polynomials');
        
                figure;
                for i=1:recordNum+1
                    plot(tArray*s,recordLagBasisFunc(i,:));
                    hold on;
                end
                title('Laguerre Basis Functions');
            end
        end
    end
    
    %% Compute Laguerre Coefficients for the source
    
    jq = 0;
    
    for i=1:tStep
        jq = jq+waveform(1,i)*(lagPoly(3,i).*exp(lagPoly(4,i)))*s*dt;
    end
    
    
    % Record Laguerre coefficients for the source and recover the source
    
    if plotRecJ==1
        
        recordJq(q+1,1) = jq;
        recordJqPhi(q+1,:) = recordJq(q+1,1)*(lagPoly(3,:).*exp(lagPoly(4,:)));
        
        if q==qStop
            for i=1:tStep
                recordRecJ(i,1) = sum(recordJqPhi(:,i));
            end
            
            % Time Domain
            figure;
            subplot(2,1,1);
            plot(tArray,waveform);
            title('Current source');
            subplot(2,1,2);
            plot(tArray,recordRecJ);
            title('Current source recovered from sum of q');
            
            % Frequency Domain
            figure;
            RecJ_Freq = fft(recordRecJ);
            plot(fn,abs(RecJ_Freq(nn+1))*2/tStep,'Linewidth',2);
            xlim([0 3]);
            title('Frequency Domain');
            xlabel('Frequency (GHz)');
            ylabel('Magnitude');
            grid;
            
        end
    end
    
    
    %% Build 1D b vector
    
    if TFSF_Update == 1
        
        b_1D(:,:) = 0;
        
        % 1D Ex equation
        for i=1:1
            for j=1:1
                for k = 2:nz-nz_TFSF+1
                    
                    ex111 = nodeNum_1D(k);

                    b_1D(ex111,1) = 2*cez(i,j,k)*(sumHy_1D(i,j,k)-sumHy_1D(i,j,k-1))-2*sumE_1D(ex111,1);
                    
                end
            end
        end
        
        % Outmost ABC boundary for Ex at z = Z
        
        for i=1:1
            for j = 1:1
                for k = nz-nz_TFSF+2:nz-nz_TFSF+2
                    
                    ex111 = nodeNum_1D(k);
                    ex110 = nodeNum_1D(k-1);
                    
                    b_1D(ex111,1) = -s/(2*v0)*(sumE_1D(ex111,1)+sumE_1D(ex110,1));
                    
                end
            end
        end
        
        % Outmost TF/SF boundary for
        
        for i=1:1
            for j=1:1
                for k = 1:1
                    
                    ex111 = nodeNum_1D(k);
                    
                    b_1D(ex111,1) = jq;
                    
                end
            end
        end
        
        % Solve 1D x=A\b
        
        x_1D = A_1D\b_1D;
        
        % Update 1D variable
        
        sumE_1D = sumE_1D + x_1D;
        
        for i=1:1
            for j=1:1
                for k = 1:nnode_1D-1
                    
                    ex112 = nodeNum_1D(k+1);
                    ex111 = nodeNum_1D(k);
                    
                    hy_1D(i,j,k) = -chz(i,j,k)*(x_1D(ex112,1)-x_1D(ex111,1))-2*sumHy_1D(i,j,k);
                    
                end
            end
        end
        
        sumHy_1D = sumHy_1D+hy_1D;
        
    end

    
    %% Current source
    
    % Jqx
    if jDirecIndex(1) == 1
        
        % Bottom current source
        Jx(jCellIndex(1):jCellIndex(2),jCellIndex(3):jCellIndex(4),jCellIndex(5)) = (-1/dzh(nz_TFSF,1))*hy_1D(1,1,1);
        
        % Top current source
        Jx(jCellIndex(1):jCellIndex(2),jCellIndex(3):jCellIndex(4),jCellIndex(6)) = (1/dzh(end-nz_TFSF+1,1))*hy_1D(1,1,end-nz_TFSF+1);
        
    end
    
    % Jqy
    if jDirecIndex(2) == 1
        for i = jCellIndex(1):jCellIndex(2)
            for j = jCellIndex(3):jCellIndex(4)
                for k = jCellIndex(5):jCellIndex(6)
                    Jy(i,j,k) = jq;
                end
            end
        end
    end
    
    % Jqz
    if jDirecIndex(3) == 1
        
        % Left Jz current source

        for i = jCellIndex_Z(1):jCellIndex_Z(1)
            for j = jCellIndex_Z(3):jCellIndex_Z(4)
                for k = jCellIndex_Z(5):jCellIndex_Z(6)
                    Jz(i,j,k) = (1/dxh(nz_TFSF,1))*hy_1D(1,1,k-nz_TFSF+1);
                end
            end
        end
        
        % Right Jz current source
        
        for i = jCellIndex_Z(2):jCellIndex_Z(2)
            for j = jCellIndex_Z(3):jCellIndex_Z(4)
                for k = jCellIndex_Z(5):jCellIndex_Z(6)
                    Jz(i,j,k) = (-1/dxh(nz_TFSF,1))*hy_1D(1,1,k-nz_TFSF+1);
                end
            end
        end
        
    end
    
    % Mqx
    if mDirecIndex(1) == 1
        for i = mCellIndex(1):mCellIndex(2)
            for j = mCellIndex(3):mCellIndex(4)
                for k = mCellIndex(5):mCellIndex(6)
                    Mx(i,j,k) = jq*ita0;
                end
            end
        end
    end
    
    % Mqy
    if mDirecIndex(2) == 1
        
        % Bottom magnetic source
        My(mCellIndex(1):mCellIndex(2),mCellIndex(3):mCellIndex(4),mCellIndex(5)) = (-1/dze(nz_TFSF,1))*x_1D(2,1);
        
        % Top magnetic source
        My(mCellIndex(1):mCellIndex(2),mCellIndex(3):mCellIndex(4),mCellIndex(6)) = (1/dze(end-nz_TFSF+1,1))*x_1D(end-nz_TFSF,1);

    end
    
    % Mqz
    if mDirecIndex(3) == 1
        
        % Front magnetic source
        for i = mCellIndex_Z(1):mCellIndex_Z(2)
            for j = mCellIndex_Z(3):mCellIndex_Z(3)
                for k = mCellIndex_Z(5):mCellIndex_Z(6)
                    Mz(i,j,k) = (1/dye(nz_TFSF,1))*x_1D(k-nz_TFSF+1,1);
                end
            end
        end
        
        % Back magnetic source
        for i = mCellIndex_Z(1):mCellIndex_Z(2)
            for j = mCellIndex_Z(4):mCellIndex_Z(4)
                for k = mCellIndex_Z(5):mCellIndex_Z(6)
                    Mz(i,j,k) = (-1/dye(nz_TFSF,1))*x_1D(k-nz_TFSF+1,1);
                end
            end
        end
        
    end
    %% Build b vector
    
    % Clear b vector
    
    b(:,:) = 0;
            
            
    % Ex equation except outmost PEC boundary
    % No re-assignment of b value for outmost PEC boundary
    for i=1:nx
        for j=2:ny
            for k=2:nz

                ex111 = nodeNum(i,j,k,1);
                b(ex111,1) = -2*cey(i,j,k)*(sumHz(i,j,k)-sumHz(i,j-1,k))+2*cez(i,j,k)*(sumHy(i,j,k)-sumHy(i,j,k-1))-2/(s*eps(i,j,k))*Jx(i,j,k)...
                    +cez(i,j,k)*(2/(s*mu(i,j,k)))*(My(i,j,k)-My(i,j,k-1))-cey(i,j,k)*(2/(s*mu(i,j,k)))*(Mz(i,j,k)-Mz(i,j-1,k))...
                    -2*sumE(ex111,1);

            end
        end
    end

    % Ey equation except outmost PEC boundary
    % No re-assignment of b value for outmost PEC boundary
    for i=2:nx
        for j=1:ny
            for k=2:nz

                ey111 = nodeNum(i,j,k,2);
                b(ey111,1) = -2*cez(i,j,k)*(sumHx(i,j,k)-sumHx(i,j,k-1))+2*cex(i,j,k)*(sumHz(i,j,k)-sumHz(i-1,j,k))-2/(s*eps(i,j,k))*Jy(i,j,k)...
                    -cez(i,j,k)*(2/(s*mu(i,j,k)))*(Mx(i,j,k)-Mx(i,j,k-1))+cex(i,j,k)*(2/(s*mu(i,j,k)))*(Mz(i,j,k)-Mz(i-1,j,k))...
                    -2*sumE(ey111,1);

            end
        end
    end

    % Ez equation except outmost PEC boundary
    % No re-assignment of b value for outmost PEC boundary
    for i=2:nx
        for j=2:ny
            for k=1:nz

                ez111 = nodeNum(i,j,k,3);
                b(ez111,1) = -2*cex(i,j,k)*(sumHy(i,j,k)-sumHy(i-1,j,k))+2*cey(i,j,k)*(sumHx(i,j,k)-sumHx(i,j-1,k))-2/(s*eps(i,j,k))*Jz(i,j,k)...
                    +cey(i,j,k)*(2/(s*mu(i,j,k)))*(Mx(i,j,k)-Mx(i,j-1,k))-cex(i,j,k)*(2/(s*mu(i,j,k)))*(My(i,j,k)-My(i-1,j,k))...
                    -2*sumE(ez111,1);

            end
        end
    end
        
    
    
    % Boundary for outmost PEC and ABC (PEC, b is 0, nothing to re-assign)
    
    if boundTypeIndex==2 % ABC boundary
        
       % Outmost ABC boundary for Ex
    
        for i=1:nx

            % Edge

            for j=1:ny:ny+1
                for k=1:nz:nz+1

                    if j==1 && k==1 % Case (a-1)
                        ex111 = nodeNum(i,j,k,1);
                        ex122 = nodeNum(i,j+1,k+1,1);
                        b(ex111,1) = -s/(2*v0)*(sumE(ex111,1)+sumE(ex122,1));       
                    elseif j==1 && k==nz+1 % Case (a-2)
                        ex111 = nodeNum(i,j,k,1);
                        ex120 = nodeNum(i,j+1,k-1,1);
                        b(ex111,1) = -s/(2*v0)*(sumE(ex111,1)+sumE(ex120,1));                      
                    elseif j==ny+1 && k==1 % Case (a-3)
                        ex111 = nodeNum(i,j,k,1);
                        ex102 = nodeNum(i,j-1,k+1,1);
                        b(ex111,1) = -s/(2*v0)*(sumE(ex111,1)+sumE(ex102,1));
                    elseif j==ny+1 && k==nz+1 % Case (a-4)    
                        ex111 = nodeNum(i,j,k,1);
                        ex100 = nodeNum(i,j-1,k-1,1);
                        b(ex111,1) = -s/(2*v0)*(sumE(ex111,1)+sumE(ex100,1));
                    end

                end
            end

            % Face

            for j=1:ny:ny+1
                for k=2:nz

                    if j==1 % Case (1-1)              
                        ex111 = nodeNum(i,j,k,1);
                        ex121 = nodeNum(i,j+1,k,1);
                        b(ex111,1) = -s/(2*v0)*(sumE(ex111,1)+sumE(ex121,1));
                    elseif j==ny+1 % Case (1-2)
                        ex111 = nodeNum(i,j,k,1);
                        ex101 = nodeNum(i,j-1,k,1);
                        b(ex111,1) = -s/(2*v0)*(sumE(ex111,1)+sumE(ex101,1));
                    end

                end
            end

            for j=2:ny
                for k=1:nz:nz+1

                    if k==1 % Case (1-3)
                        ex111 = nodeNum(i,j,k,1);
                        ex112 = nodeNum(i,j,k+1,1);
                        b(ex111,1) = -s/(2*v0)*(sumE(ex111,1)+sumE(ex112,1));
                    elseif k==nz+1 % Case (1-4)
                        ex111 = nodeNum(i,j,k,1);
                        ex110 = nodeNum(i,j,k-1,1);
                        b(ex111,1) = -s/(2*v0)*(sumE(ex111,1)+sumE(ex110,1));
                    end

                end
            end

        end

        % Outmost ABC boundary for Ey

        for j=1:ny

            % Edge

            for i=1:nx:nx+1
                for k=1:nz:nz+1

                    if k==1 && i==1 % Case (b-1)
                        ey111 = nodeNum(i,j,k,2);
                        ey212 = nodeNum(i+1,j,k+1,2);
                        b(ey111,1) = -s/(2*v0)*(sumE(ey111,1)+sumE(ey212,1));
                    elseif k==1 && i==nx+1 % Case (b-2)
                        ey111 = nodeNum(i,j,k,2);
                        ey012 = nodeNum(i-1,j,k+1,2);
                        b(ey111,1) = -s/(2*v0)*(sumE(ey111,1)+sumE(ey012,1));
                    elseif k==nz+1 && i==1 % Case (b-3)
                        ey111 = nodeNum(i,j,k,2);
                        ey210 = nodeNum(i+1,j,k-1,2);
                        b(ey111,1) = -s/(2*v0)*(sumE(ey111,1)+sumE(ey210,1));
                    elseif k==nz+1 && i==nx+1 % Case (b-4)
                        ey111 = nodeNum(i,j,k,2);
                        ey010 = nodeNum(i-1,j,k-1,2);
                        b(ey111,1) = -s/(2*v0)*(sumE(ey111,1)+sumE(ey010,1));
                    end

                end
            end

            % Face

            for i=2:nx
                for k=1:nz:nz+1

                    if k==1 % Case (2-1)
                        ey111 = nodeNum(i,j,k,2);
                        ey112 = nodeNum(i,j,k+1,2);
                        b(ey111,1) = -s/(2*v0)*(sumE(ey111,1)+sumE(ey112,1));
                    elseif k==nz+1 % Case (2-2)
                        ey111 = nodeNum(i,j,k,2);
                        ey110 = nodeNum(i,j,k-1,2);
                        b(ey111,1) = -s/(2*v0)*(sumE(ey111,1)+sumE(ey110,1));
                    end

                end
            end

            for i=1:nx:nx+1
                for k=2:nz

                    if i==1 % Case (2-3)
                        ey111 = nodeNum(i,j,k,2);
                        ey211 = nodeNum(i+1,j,k,2);
                        b(ey111,1) = -s/(2*v0)*(sumE(ey111,1)+sumE(ey211,1));
                    elseif i==nx+1 % Case (2-4)
                        ey111 = nodeNum(i,j,k,2);
                        ey011 = nodeNum(i-1,j,k,2);
                        b(ey111,1) = -s/(2*v0)*(sumE(ey111,1)+sumE(ey011,1));
                    end

                end
            end

        end

        % Outmost ABC boundary for Ez

        for k=1:nz

            % Edge

            for i=1:nx:nx+1
                for j=1:ny:ny+1

                    if i==1 && j==1 % Case (c-1)
                        ez111 = nodeNum(i,j,k,3);
                        ez221 = nodeNum(i+1,j+1,k,3);
                        b(ez111,1) = -s/(2*v0)*(sumE(ez111,1)+sumE(ez221,1));
                    elseif i==1 && j==ny+1 % Case (c-2)
                        ez111 = nodeNum(i,j,k,3);
                        ez201 = nodeNum(i+1,j-1,k,3);
                        b(ez111,1) = -s/(2*v0)*(sumE(ez111,1)+sumE(ez201,1));
                    elseif i==nx+1 && j==1 % Case (c-3)
                        ez111 = nodeNum(i,j,k,3);
                        ez021 = nodeNum(i-1,j+1,k,3);
                        b(ez111,1) = -s/(2*v0)*(sumE(ez111,1)+sumE(ez021,1));
                    elseif i==nx+1 && j==ny+1 % Case (c-4)
                        ez111 = nodeNum(i,j,k,3);
                        ez001 = nodeNum(i-1,j-1,k,3);
                        b(ez111,1) = -s/(2*v0)*(sumE(ez111,1)+sumE(ez001,1));
                    end

                end
            end

            % Face

            for i=1:nx:nx+1
                for j=2:ny

                    if i==1 % Case (3-1)
                        ez111 = nodeNum(i,j,k,3);
                        ez211 = nodeNum(i+1,j,k,3);
                        b(ez111,1) = -s/(2*v0)*(sumE(ez111,1)+sumE(ez211,1));
                    elseif i==nx+1 % Case (3-2)
                        ez111 = nodeNum(i,j,k,3);
                        ez011 = nodeNum(i-1,j,k,3);
                        b(ez111,1) = -s/(2*v0)*(sumE(ez111,1)+sumE(ez011,1));
                    end

                end
            end

            for i=2:nx
                for j=1:ny:ny+1

                    if j==1 % Case (3-3)
                        ez111 = nodeNum(i,j,k,3);
                        ez121 = nodeNum(i,j+1,k,3);
                        b(ez111,1) = -s/(2*v0)*(sumE(ez111,1)+sumE(ez121,1));
                    elseif j==ny+1 % Case (3-4)
                        ez111 = nodeNum(i,j,k,3);
                        ez101 = nodeNum(i,j-1,k,3);
                        b(ez111,1) = -s/(2*v0)*(sumE(ez111,1)+sumE(ez101,1));
                    end

                end
            end 

        end 

    end
    
    %% Solve Ax = b
    

    [x,flag] = bicgstabl(A,b,tol,maxit,D);

    if flag~=0
        fprintf('BiCGSTABL Convergency flag is: %d\n',flag);
        return;
    end

    
    %% Update variables
    
    
    % Update H field
    

    % Update Hx
    for i=1:nx
        for j=1:ny
            for k=1:nz

                ey112 = nodeNum(i,j,k+1,2);
                ey111 = nodeNum(i,j,k,2);
                ez121 = nodeNum(i,j+1,k,3);
                ez111 = nodeNum(i,j,k,3);

                hx(i,j,k) = chz(i,j,k)*(x(ey112,1)-x(ey111,1))-chy(i,j,k)*(x(ez121,1)-x(ez111,1))-2*sumHx(i,j,k)-(2/(s*mu(i,j,k)))*Mx(i,j,k);

            end
        end
    end

    % Update Hy
    for i=1:nx
        for j=1:ny
            for k=1:nz

                ez211 = nodeNum(i+1,j,k,3);
                ez111 = nodeNum(i,j,k,3);
                ex112 = nodeNum(i,j,k+1,1);
                ex111 = nodeNum(i,j,k,1);

                hy(i,j,k) = chx(i,j,k)*(x(ez211,1)-x(ez111,1))-chz(i,j,k)*(x(ex112,1)-x(ex111,1))-2*sumHy(i,j,k)-(2/(s*mu(i,j,k)))*My(i,j,k);

            end
        end
    end

    % Update Hz
    for i=1:nx
        for j=1:ny
            for k=1:nz

                ex121 = nodeNum(i,j+1,k,1);
                ex111 = nodeNum(i,j,k,1);
                ey211 = nodeNum(i+1,j,k,2);
                ey111 = nodeNum(i,j,k,2);

                hz(i,j,k) = chy(i,j,k)*(x(ex121,1)-x(ex111,1))-chx(i,j,k)*(x(ey211,1)-x(ey111,1))-2*sumHz(i,j,k)-(2/(s*mu(i,j,k)))*Mz(i,j,k);

            end
        end
    end


    
    % Update sumHx, sumHy, sumHz

    sumHz = sumHz + hz;
    sumHy = sumHy + hy;
    sumHx = sumHx + hx;

    % Update E field
    
    sumE = sumE+x;
    
    % Print the basis coefficient for the port with the lowest x and y index
    
    if printQ==1
        fprintf('q = %5d;',q);
        for i=1:numProbe
            if probeDirecIndex(1)==1
                fprintf('p%d = %15.5e;',i,x(nodeNum(probeCellIndex(1,i),probeCellIndex(3,i),probeCellIndex(5,i),1)));
            elseif probeDirecIndex(2)==1
                fprintf('p%d = %15.5e;',i,x(nodeNum(probeCellIndex(1,i),probeCellIndex(3,i),probeCellIndex(5,i),2)));
            elseif probeDirecIndex(3)==1
                fprintf('p%d = %15.5e;',i,x(nodeNum(probeCellIndex(1,i),probeCellIndex(3,i),probeCellIndex(5,i),3)));
            else
                fprintf('Probe printing error!\n');
                return;
            end
        end
        fprintf('\n');
    end
                
    % Update probe E field (voltage)
    
    for n=1:numProbe
        
        vtg = zeros(1,tStep);
        
        for i=probeCellIndex(1,n):probeCellIndex(2,n)
            for j=probeCellIndex(3,n):probeCellIndex(4,n)
                for k=probeCellIndex(5,n):probeCellIndex(6,n)
                    
                    if probeDirecIndex(1)==1
                        vtg = vtg+x(nodeNum(i,j,k,1),1)*(lagPoly(3,:).*exp(lagPoly(4,:)))*dx(i);
                    elseif probeDirecIndex(2)==1
                        vtg = vtg+x(nodeNum(i,j,k,2),1)*(lagPoly(3,:).*exp(lagPoly(4,:)))*dy(j);
                    elseif probeDirecIndex(3)==1
                        vtg = vtg+x(nodeNum(i,j,k,3),1)*(lagPoly(3,:).*exp(lagPoly(4,:)))*dz(k);
                    else
                        fprintf('Probe setting error!\n');
                        return;
                    end
                    
                end
            end
        end
        
        vtg = vtg';
        probe(:,n) = probe(:,n)+vtg;
        
    end
    
    % Store field profile for the whole area
    
    if plotFieldEx==1 || plotFieldEy == 1 || plotFieldEz == 1
        recordE = recordE+x*(lagPoly(3,(1:fac:tStep)).*exp(lagPoly(4,(1:fac:tStep))));
        for i = 1:plotStep
            recordHx(:,:,:,i) = recordHx(:,:,:,i) + hx.*(lagPoly(3,(1+(i-1)*fac)).*exp(lagPoly(4,(1+(i-1)*fac))));
            recordHy(:,:,:,i) = recordHy(:,:,:,i) + hy.*(lagPoly(3,(1+(i-1)*fac)).*exp(lagPoly(4,(1+(i-1)*fac))));
            recordHz(:,:,:,i) = recordHz(:,:,:,i) + hz.*(lagPoly(3,(1+(i-1)*fac)).*exp(lagPoly(4,(1+(i-1)*fac))));
            
        end
        
    end
    
    if plot1DField == 1 && TFSF_Update == 1
        for i = 1:plotStep
            recordHy_1D(:,:,:,i) = recordHy_1D(:,:,:,i) + hy_1D.*(lagPoly(3,(1+(i-1)*fac)).*exp(lagPoly(4,(1+(i-1)*fac))));
            recordHy_1D(:,:,:,i) = recordHy_1D(:,:,:,i) + hy_1D.*(lagPoly(3,(1+(i-1)*fac)).*exp(lagPoly(4,(1+(i-1)*fac))));
        end
    end
    
% Update q   
    
q = q+1;  

end

fprintf('5. Calculation done with q = %d\n',qStop);
tFinal = toc;
fprintf('6. Simulation duration is %d seconds\n',round(tFinal));

%% Postprocess

% Plot the field distribution in time using function 'surf'
% Eg. If plotSec==3, plotLayer==4, plotFieldEz==1
% This plot the Ez component of layer 4 in z cross section



% Plot E Field

if plotFieldEx == 1 || plotFieldEy==1 || plotFieldEz==1

    im = LagPos.Plot_E_Field(plotFieldEx,plotFieldEy,plotFieldEz,plotSec,recordE,nodeNum,nx,ny,nz,plotStep,plotLayer);

end

% Plot H Field

if plotFieldHx == 1 || plotFieldHy==1 || plotFieldHz==1

    im = LagPos.Plot_H_Field(plotFieldHx,plotFieldHy,plotFieldHz,plotSec,recordHx,recordHy,recordHz,nx,ny,nz,plotStep,plotLayer);

end
    


%% Write 1D Waveform
