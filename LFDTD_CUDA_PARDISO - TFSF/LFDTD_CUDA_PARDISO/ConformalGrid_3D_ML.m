function [sigma_x,sigma_y,sigma_z] =  ConformalGrid_3D_ML(NX,NY,NZ,dx,dy,nr,sigma,cor)

% NX/NY: Number of grids along x/y direction
% dx/dy: Meshing resolution
% nr: Radius of circle in terms of number of grids, not actual length
% cor: n X 2 matrix, where n is the total number of inhomogeneous material
% [x_cor y_cor] - center cor for material 1
% [x_cor y_cor] - center cor for material 2
% [x_cor y_cor] - center cor for material 3, etc
% cor - offset for the bottom apearture location - [x_cor_start, x_cor_end; y_cor_start, y_cor_end; z_cor_start, z_cor_end]

% |
% |
% |
% |
% |
% |
% |
% |
% | (1,1)
% --------------------------- Coordinate for cor - starts from (1,1)
% (1,1)

nx = NX/2;
ny = NY/2;

[~,dim] = size(nr); % Number of inhomogeneous material

x = 1:nx+1;
y = 1:ny+1;
[x_cor, y_cor] = meshgrid(x,y);
Ex_Q1 = ones(nx,ny+1,dim);
Ey_Q1 = ones(nx+1,ny,dim);
R = nr*dx;

sigma_x_Plane = zeros(NX,NY,dim+2);
sigma_y_Plane = zeros(NX,NY,dim+2);
sigma_z_Hollow = zeros(NX,NY);
sigma_z_Hole = zeros(NX,NY);
sigma_z_Temp = zeros(nx+1,ny+1);

sigma_x = zeros(NX,NY,NZ);
sigma_y = zeros(NX,NY,NZ);
sigma_z = zeros(NX,NY,NZ);

for n = 1:dim
    % Ex length calculation
    for j = 1:nr(n)+1
        theta = asin((y_cor(j,1)-1)*dy/R(n));
        L_x = R(n)*cos(theta);
        Delta_x = ceil(L_x/dx);
        Ex_Q1(Delta_x,j,n) = (dx-abs(Delta_x*dx - L_x))/dx;
        Ex_Q1(Delta_x+1:end,j:end,n) = 0;
    end
    Ex_Q1(:,nr(n)+2:end,n) = 0;

    % Ey length calculation
    for i = 1:nr(n)+1
        theta = asin((x_cor(1,i)-1)*dx/R(n));
        L_y = R(n)*cos(theta);
        Delta_y = ceil(L_y/dy);
        Ey_Q1(i,Delta_y,n) = (dy-abs(Delta_y*dy - L_y))/dy;
        Ey_Q1(i:end,Delta_y+1:end,n) = 0;
    end
    Ey_Q1(nr(n)+2:end,:,n) = 0;

    %% Extend Ex to full circle

    Ex_Q2 = [zeros(nx,ny+1);Ex_Q1(:,:,n)];
    Buffer_Matrix_x = flip(Ex_Q2,1);
    Ex_Q2 = Ex_Q2 + Buffer_Matrix_x;

    [grid_x, grid_y] = size(Ex_Q2);

    Ex = [zeros(grid_x,grid_y-1), Ex_Q2];
    Buffer_Matrix_x = flip(Ex, 2);
    Buffer_Matrix_x(:,ny+1) = 0;
    Ex = Ex + Buffer_Matrix_x;

    %% Extend Ey to full circle

    Ey_Q2 = [zeros(nx,ny); Ey_Q1(:,:,n)];
    Buffer_Matrix_y = flip(Ey_Q2,1);
    Buffer_Matrix_y(nx+1,:) = 0;
    Ey_Q2 = Ey_Q2 + Buffer_Matrix_y;

    [grid_x, grid_y] = size(Ey_Q2);
    Ey = [zeros(grid_x,grid_y), Ey_Q2];
    Buffer_Matrix_y = flip(Ey, 2);
    Ey = Ey + Buffer_Matrix_y;

    %% Assign Sigma_X/Sigma_Y according to Ex/Ey

    if cor(n,1)<0 && cor(n,2)>0

        sigma_x_Plane(1:NX-abs(cor(n,1))+1,cor(n,2):NY,n) = Ex(abs(cor(n,1)):NX,1:NY-cor(n,2)+1);
        sigma_y_Plane(1:NX-abs(cor(n,1))+1,cor(n,2):NY,n) = Ey(abs(cor(n,1)):NX,1:NY-cor(n,2)+1);

    elseif cor(n,1)>0 && cor(n,2)<0

        sigma_x_Plane(cor(n,1):NX,1:NY-abs(cor(n,2))+1,n) = Ex(1:NX-cor(n,1)+1,abs(cor(n,2)):NY);
        sigma_y_Plane(cor(n,1):NX,1:NY-abs(cor(n,2))+1,n) = Ey(1:NX-cor(n,1)+1,abs(cor(n,2)):NY);

    elseif cor(n,1)<0 && cor(n,2)<0

        sigma_x_Plane(1:NX-abs(cor(n,1))+1,1:NY-abs(cor(n,2))+1,n) = Ex(abs(cor(n,1)):NX,abs(cor(n,2)):NY);
        sigma_y_Plane(1:NX-abs(cor(n,1))+1,1:NY-abs(cor(n,2))+1,n) = Ey(abs(cor(n,1)):NX,abs(cor(n,2)):NY);

    else
        sigma_x_Plane(cor(n,1):NX,cor(n,2):NY,n) = Ex(1:NX-cor(n,1)+1,1:NY-cor(n,2)+1);
        sigma_y_Plane(cor(n,1):NX,cor(n,2):NY,n) = Ey(1:NX-cor(n,1)+1,1:NY-cor(n,2)+1);

    end


end

%% Calculate effective conductivity for two material

sigma_x_Plane(:,:,dim+1) = (sigma_x_Plane(:,:,1) - sigma_x_Plane(:,:,2))*sigma(1);
sigma_y_Plane(:,:,dim+1) = (sigma_y_Plane(:,:,1) - sigma_y_Plane(:,:,2))*sigma(1);

sigma_x_Plane(:,:,dim+2) = (sigma_x_Plane(:,:,1) - sigma_x_Plane(:,:,3))*sigma(1);
sigma_y_Plane(:,:,dim+2) = (sigma_y_Plane(:,:,1) - sigma_y_Plane(:,:,3))*sigma(1);

%% Assign Sigma_Z from Sigma_X/Sigma_Y

for i = 1:nx+1
    for j = 1:ny+1
        if (sigma_x_Plane(i,j,1) == 1) || (sigma_y_Plane(i,j,1) == 1)
            sigma_z_Temp(i,j) = sigma(1);
        end
    end
end

Sigma_z_Q2 = flip(sigma_z_Temp,2);
Sigma_z_Q2 = [sigma_z_Temp,Sigma_z_Q2(:,2:end)];

sigma_z_Plane = flip(Sigma_z_Q2);
sigma_z_Plane = [Sigma_z_Q2;sigma_z_Plane(2:end,:)];


for i = 1:NX
    for j = 1:NY
        if (sigma_x_Plane(i,j,dim+1) == sigma(1))
            sigma_z_Hollow(i,j) = sigma(1);
            sigma_z_Hollow(i+1,j) = sigma(1);
        elseif (sigma_x_Plane(i,j,dim+1) >= 0.1) && (sigma_x_Plane(i+1,j,dim+1) >= 0.1)
            sigma_z_Hollow(i+1,j) = sigma(1);
        end
    end
end

for i = 1:NX
    for j = 1:NY
        if (sigma_y_Plane(i,j,dim+1) == sigma(1))
            sigma_z_Hollow(i,j) = sigma(1);
            sigma_z_Hollow(i,j+1) = sigma(1);
        elseif (sigma_y_Plane(i,j,dim+1) >= 0.1) && (sigma_y_Plane(i,j+1,dim+1) >= 0.1)
            sigma_z_Hollow(i,j+1) = sigma(1);
        end
    end
end


for i = 1:NX
    for j = 1:NY
        if (sigma_x_Plane(i,j,dim+2) == sigma(1))
            sigma_z_Hole(i,j) = sigma(1);
            sigma_z_Hole(i+1,j) = sigma(1);
        elseif (sigma_x_Plane(i,j,dim+2) >= 0.1) && (sigma_x_Plane(i+1,j,dim+2) >= 0.1)
            sigma_z_Hole(i+1,j) = sigma(1);
        end
    end
end

for i = 1:NX
    for j = 1:NY
        if (sigma_y_Plane(i,j,dim+2) == sigma(1))
            sigma_z_Hole(i,j) = sigma(1);
            sigma_z_Hole(i,j+1) = sigma(1);
        elseif (sigma_y_Plane(i,j,dim+2) >= 0.1) && (sigma_y_Plane(i,j+1,dim+2) >= 0.1)
            sigma_z_Hole(i,j+1) = sigma(1);
        end
    end
end

%% Extend to 3D Meshing

for i = 9:56
    sigma_x(:,:,i) = sigma_x_Plane(:,:,dim+1);
    sigma_y(:,:,i) = sigma_y_Plane(:,:,dim+1);
    sigma_z(:,:,i) = sigma_z_Hollow(:,:);
end

for i = 10:56
    sigma_x(:,:,i) = sigma_x_Plane(:,:,dim+1);
    sigma_y(:,:,i) = sigma_y_Plane(:,:,dim+1);
end

sigma_x(:,:,8) = sigma_x_Plane(:,:,dim+2);
sigma_y(:,:,8) = sigma_y_Plane(:,:,dim+2);
sigma_x(:,:,9) = sigma_x_Plane(:,:,dim+2);
sigma_y(:,:,9) = sigma_y_Plane(:,:,dim+2);
sigma_z(:,:,8) = sigma_z_Hole(1:NX,1:NY);

sigma_x(:,:,57) = sigma_x_Plane(:,:,1)*sigma(1);
sigma_y(:,:,57) = sigma_y_Plane(:,:,1)*sigma(1);
sigma_x(:,:,58) = sigma_x_Plane(:,:,1)*sigma(1);
sigma_y(:,:,58) = sigma_y_Plane(:,:,1)*sigma(1);
sigma_z(:,:,57) = sigma_z_Plane(1:NX,1:NY);
end