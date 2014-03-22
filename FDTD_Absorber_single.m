%% Application of Finite-Difference-Time-Domain (FDTD) Method in 3D
%% Electromagnetic Wave Propagation: Terahertz Fishnet Metamaterial
%% Energy Absorber

close all; 
clear all; 
clc;
%% Free space constants
epsilon0 = 8.85418782e-12;  % Permittivity of free space
mu0 = 1.25663706e-6;        % Permeability of free space
c = 1.0/sqrt(mu0*epsilon0); % speed of light 

%% Gaussian half-width
f_target = 3e12; 
t_half = 1/f_target;

%% End time
t_end = 100*t_half;

%% Total mesh dimensions and grid cells sizes (without PML)
ratio = 1;
%nx = 50; ny = 50; nz = 30;
ncellx = 50; ncelly = 50;
nf = 18;
nx = nf*ratio+10; ny = nf*ratio+10; nz = 30;
dx = 1e-5; dy = dx; dz = dx;


%% Number of PML layers
PML = 10;

%% Matrix of material's constants
number_of_materials = 3;
% For material of number x = 1,2,3... :
% Material(x,1) - relative permittivity 
% Material(x,2) - relative permeability
% Material(x,3) - specific conductivity

% Free space
Material(1,1) = 1.0;   Material(1,2) = 1.0;   Material(1,3) = 0.0;
% Metal (Copper)
Material(2,1) = 1.0;   Material(2,2) = 1.0;   Material(2,3) = 5.96e+7;
% Substrate material (RT/Duroid 5880)
Material(3,1) = 2.2;   Material(3,2) = 1.0;   Material(3,3) = 0.0;

%% Add PML layers
nx = nx + 2*PML; ny = ny + 2*PML; nz = nz + 2*PML;

%% Calculate dt    
dt = (1.0/c/sqrt( 1.0/(dx^2) + 1.0/(dy^2) + 1.0/(dz^2)))*0.9999;
number_of_iterations = ceil(t_end/dt);

%% 3D array for geometry
Index = ones(nx, ny, nz);

%% Define of terahertz metamaterial energy absorber 
% Bottom electrode: thickness = 1 cell
Index((1+PML):(nx-PML), (1+PML):(ny-PML), PML+1) = 2;

% Fishnet cells (top electrode) : thickness = 1 cell
for iii = 1:1:ratio

    for jjj = 1:1:ratio

        sl = (ncellx/2-7)+18*(iii-1);
        sr = (ncellx/2+7)+18*(iii-1);
        sb = (ncelly/2-7)+18*(jjj-1);
        st = (ncelly/2+7)+18*(jjj-1);
        hl = (ncellx/2-9)+18*(iii-1);
        hr = (ncellx/2+9)+18*(iii-1);
        hb = ncelly/2+18*(jjj-1);
        ht = (ncelly/2+1)+18*(jjj-1);
        vl = (ncellx/2)+18*(iii-1);
        vr = (ncellx/2+1)+18*(iii-1);
        vb =(ncelly/2-9)+18*(jjj-1);
        vt = (ncelly/2+9)+18*(jjj-1);
        Index(sl:sr, sb:st, PML+4) = 2; % square
        Index(hl:hr, hb:ht, PML+4) = 2; %H-rectangle
        Index(vl:vr, vb:vt, PML+4) = 2; %V-rectangle

    end

end

% Dielectric layer between top and bottom electrodes: thickness = 2 cells
Index((1+PML):(nx-PML), (1+PML):(ny-PML), PML+3) = 3;

%% 3D FDTD physical (fields) and additional arrays are defined as 'single' 
%% to increase performance
Ex = zeros(nx, ny+1, nz+1, 'single'); 
Gx = zeros(nx, ny+1, nz+1, 'single'); 
Fx = zeros(nx, ny+1, nz+1, 'single');  
Ey = zeros(nx+1, ny, nz+1, 'single'); 
Gy = zeros(nx+1, ny, nz+1, 'single'); 
Fy = zeros(nx+1, ny, nz+1, 'single');
Ez = zeros(nx+1, ny+1, nz, 'single'); 
Gz = zeros(nx+1, ny+1, nz, 'single'); 
Fz = zeros(nx+1, ny+1, nz, 'single');
Hx = zeros(nx+1, ny, nz, 'single'); 
Bx = zeros(nx+1, ny, nz, 'single'); 
Hy = zeros(nx, ny+1, nz, 'single');
By = zeros(nx, ny+1, nz, 'single'); 
Hz = zeros(nx, ny, nz+1, 'single'); 
Bz = zeros(nx, ny, nz+1, 'single');

%% FDTD PML coefficients arrays. Here they are already filled with values 
%% corresponding to free space 
m = 4; ka_max = 1.0; R_err = 1.0e-16;
eta = sqrt(mu0/epsilon0*Material(1,1)/Material(1,2));
k_Ex_c = ones(nx, ny, nz, 'single')*2.0*epsilon0;  
k_Ex_d = ones(nx, ny, nz, 'single')*(-2.0*epsilon0);
k_Ey_a = ones(nx+1, ny, nz, 'single');
k_Ey_b = ones(nx+1, ny, nz, 'single')/(2.0*epsilon0);
k_Gz_a = ones(nx+1, ny, nz, 'single');
k_Gz_b = ones(nx+1, ny, nz, 'single');
k_Hy_a = ones(nx, ny, nz, 'single'); 
k_Hy_b = ones(nx, ny, nz, 'single')/(2.0*epsilon0);
k_Hx_c = ones(nx+1, ny, nz, 'single')*2.0*epsilon0/mu0;
k_Hx_d = ones(nx+1, ny, nz, 'single')*(-2.0*epsilon0/mu0);
k_Bz_a = ones(nx, ny, nz, 'single');
k_Bz_b = ones(nx, ny, nz, 'single')*dt;
k_Gx_a = ones(nx, ny+1, nz, 'single');
k_Gx_b = ones(nx, ny+1, nz, 'single');
k_Ey_c = ones(nx, ny, nz, 'single')*2.0*epsilon0; 
k_Ey_d = ones(nx, ny, nz, 'single')*(-2.0*epsilon0);
k_Ez_a = ones(nx, ny+1, nz, 'single'); 
k_Ez_b = ones(nx, ny+1, nz, 'single')/(2.0*epsilon0);
k_Bx_a = ones(nx, ny, nz, 'single'); 
k_Bx_b = ones(nx, ny, nz, 'single')*dt;
k_Hy_c = ones(nx, ny+1, nz, 'single')*2.0*epsilon0/mu0;
k_Hy_d = ones(nx, ny+1, nz, 'single')*(-2.0*epsilon0/mu0);
k_Hz_a = ones(nx, ny, nz, 'single');
k_Hz_b = ones(nx, ny, nz, 'single')/(2.0*epsilon0);
k_Ex_a = ones(nx, ny, nz+1, 'single');
k_Ex_b = ones(nx, ny, nz+1, 'single')/(2.0*epsilon0);
k_Gy_a = ones(nx, ny, nz+1, 'single'); 
k_Gy_b = ones(nx, ny, nz+1, 'single');
k_Ez_c = ones(nx, ny, nz, 'single')*2.0*epsilon0;
k_Ez_d = ones(nx, ny, nz, 'single')*(-2.0*epsilon0);
k_Hx_a = ones(nx, ny, nz, 'single');
k_Hx_b = ones(nx, ny, nz, 'single')/(2.0*epsilon0);
k_By_a = ones(nx, ny, nz, 'single');  
k_By_b = ones(nx, ny, nz, 'single')*dt;
k_Hz_c = ones(nx, ny, nz+1, 'single')*2.0*epsilon0/mu0;   
k_Hz_d = ones(nx, ny, nz+1, 'single')*(-2.0*epsilon0/mu0);

%% General FDTD coefficients 
I = 1:number_of_materials;
K_a(I) = (2.0*epsilon0*Material(I,1) - Material(I,3)*dt)./...
(2.0*epsilon0*Material(I,1) + Material(I,3)*dt);
K_b(I) = 2.0*dt./(2.0*epsilon0*Material(I,1) + Material(I,3)*dt);
K_c(I) = Material(I,2);
Ka = single(K_a(Index)); Kb = single(K_b(Index)); Kc = single(K_c(Index));

%% PML coefficients along x-axis
sigma_max = -(m + 1.0)*log(R_err)/(2.0*eta*PML*dx);
for I=0:(PML-1)
    
    sigma_x = sigma_max*((PML - I)/PML)^m;
    ka_x = 1.0 + (ka_max - 1.0)*((PML - I)/PML)^m;
    k_Ey_a(I+1,:,:) = (2.0*epsilon0*ka_x - sigma_x*dt)/...
    (2.0*epsilon0*ka_x + sigma_x*dt);
    k_Ey_a(nx-I,:,:) = k_Ey_a(I+1,:,:);
    k_Ey_b(I+1,:,:) = 1.0/(2.0*epsilon0*ka_x + sigma_x*dt);
    k_Ey_b(nx-I,:,:) = k_Ey_b(I+1,:,:);
    k_Gz_a(I+1,:,:) = (2.0*epsilon0*ka_x - sigma_x*dt)/...
    (2.0*epsilon0*ka_x + sigma_x*dt);
    k_Gz_a(nx-I,:,:) = k_Gz_a(I+1,:,:);
    k_Gz_b(I+1,:,:) = 2.0*epsilon0/(2.0*epsilon0*ka_x + sigma_x*dt);
    k_Gz_b(nx-I,:,:) = k_Gz_b(I+1,:,:);
    k_Hx_c(I+1,:,:) = (2.0*epsilon0*ka_x + sigma_x*dt)/mu0;
    k_Hx_c(nx-I,:,:) = k_Hx_c(I+1,:,:);
    k_Hx_d(I+1,:,:) = -(2.0*epsilon0*ka_x - sigma_x*dt)/mu0;
    k_Hx_d(nx-I,:,:) = k_Hx_d(I+1,:,:);

    sigma_x = sigma_max*((PML - I - 0.5)/PML)^m;
    ka_x = 1.0 + (ka_max - 1.0)*((PML - I - 0.5)/PML)^m;
    k_Ex_c(I+1,:,:) = 2.0*epsilon0*ka_x + sigma_x*dt;
    k_Ex_c(nx-I-1,:,:) = k_Ex_c(I+1,:,:);
    k_Ex_d(I+1,:,:) = -(2.0*epsilon0*ka_x - sigma_x*dt);
    k_Ex_d(nx-I-1,:,:) = k_Ex_d(I+1,:,:);
    k_Hy_a(I+1,:,:) = (2.0*epsilon0*ka_x - sigma_x*dt)/...
    (2.0*epsilon0*ka_x + sigma_x*dt);
    k_Hy_a(nx-I-1,:,:) = k_Hy_a(I+1,:,:);
    k_Hy_b(I+1,:,:) = 1.0/(2.0*epsilon0*ka_x + sigma_x*dt);
    k_Hy_b(nx-I-1,:,:) = k_Hy_b(I+1,:,:);
    k_Bz_a(I+1,:,:) = (2.0*epsilon0*ka_x - sigma_x*dt)/...
    (2.0*epsilon0*ka_x + sigma_x*dt);
    k_Bz_a(nx-I-1,:,:) = k_Bz_a(I+1,:,:);
    k_Bz_b(I+1,:,:) = 2.0*epsilon0*dt/(2.0*epsilon0*ka_x + sigma_x*dt);
    k_Bz_b(nx-I-1,:,:) = k_Bz_b(I+1,:,:);
    
end

%% PML coefficients along y-axis
sigma_max = -(m + 1.0)*log(R_err)/(2.0*eta*PML*dy);
for J=0:(PML-1)
    
    sigma_y = sigma_max*((PML - J)/PML)^m;
    ka_y = 1.0 + (ka_max - 1.0)*((PML - J)/PML)^m;
    k_Gx_a(:,J+1,:) = (2.0*epsilon0*ka_y - sigma_y*dt)/...
    (2.0*epsilon0*ka_y + sigma_y*dt);
    k_Gx_a(:,ny-J,:) = k_Gx_a(:,J+1,:);
    k_Gx_b(:,J+1,:) = 2.0*epsilon0/(2.0*epsilon0*ka_y + sigma_y*dt);
    k_Gx_b(:,ny-J,:) = k_Gx_b(:,J+1,:);
    k_Ez_a(:,J+1,:) = (2.0*epsilon0*ka_y - sigma_y*dt)/...
    (2.0*epsilon0*ka_y + sigma_y*dt);
    k_Ez_a(:,ny-J,:) = k_Ez_a(:,J+1,:);
    k_Ez_b(:,J+1,:) = 1.0/(2.0*epsilon0*ka_y + sigma_y*dt);
    k_Ez_b(:,ny-J,:) = k_Ez_b(:,J+1,:);
    k_Hy_c(:,J+1,:) = (2.0*epsilon0*ka_y + sigma_y*dt)/mu0;
    k_Hy_c(:,ny-J,:) = k_Hy_c(:,J+1,:);
    k_Hy_d(:,J+1,:) = -(2.0*epsilon0*ka_y - sigma_y*dt)/mu0;
    k_Hy_d(:,ny-J,:) = k_Hy_d(:,J+1,:);

    sigma_y = sigma_max*((PML - J - 0.5)/PML)^m;
    ka_y = 1.0 + (ka_max - 1.0)*((PML - J - 0.5)/PML)^m;
    k_Ey_c(:,J+1,:) = 2.0*epsilon0*ka_y+sigma_y*dt;
    k_Ey_c(:,ny-J-1,:) = k_Ey_c(:,J+1,:);
    k_Ey_d(:,J+1,:) = -(2.0*epsilon0*ka_y-sigma_y*dt);
    k_Ey_d(:,ny-J-1,:) = k_Ey_d(:,J+1,:);
    k_Bx_a(:,J+1,:) = (2.0*epsilon0*ka_y-sigma_y*dt)/...
    (2.0*epsilon0*ka_y+sigma_y*dt);
    k_Bx_a(:,ny-J-1,:) = k_Bx_a(:,J+1,:);
    k_Bx_b(:,J+1,:) = 2.0*epsilon0*dt/(2.0*epsilon0*ka_y+sigma_y*dt);
    k_Bx_b(:,ny-J-1,:) = k_Bx_b(:,J+1,:);
    k_Hz_a(:,J+1,:) = (2.0*epsilon0*ka_y-sigma_y*dt)/...
    (2.0*epsilon0*ka_y+sigma_y*dt);
    k_Hz_a(:,ny-J-1,:) = k_Hz_a(:,J+1,:);
    k_Hz_b(:,J+1,:) = 1.0/(2.0*epsilon0*ka_y+sigma_y*dt);
    k_Hz_b(:,ny-J-1,:) = k_Hz_b(:,J+1,:);
    
end

%% PML coefficients along z-axis 
sigma_max = -(m + 1.0)*log(R_err)/(2.0*eta*PML*dz);
for K=0:(PML-1)
    
    sigma_z = sigma_max*((PML - K)/PML)^m;
    ka_z = 1.0 + (ka_max - 1.0)*((PML-K)/PML)^m;
    k_Ex_a(:,:,K+1) = (2.0*epsilon0*ka_z - sigma_z*dt)/...
    (2.0*epsilon0*ka_z + sigma_z*dt);
    k_Ex_a(:,:,nz-K) = k_Ex_a(:,:,K+1);
    k_Ex_b(:,:,K+1) = 1.0/(2.0*epsilon0*ka_z + sigma_z*dt);
    k_Ex_b(:,:,nz-K) = k_Ex_b(:,:,K+1);
    k_Gy_a(:,:,K+1) = (2.0*epsilon0*ka_z - sigma_z*dt)/...
    (2.0*epsilon0*ka_z + sigma_z*dt);
    k_Gy_a(:,:,nz-K) = k_Gy_a(:,:,K+1);
    k_Gy_b(:,:,K+1) = 2.0*epsilon0/(2.0*epsilon0*ka_z + sigma_z*dt);
    k_Gy_b(:,:,nz-K) = k_Gy_b(:,:,K+1);
    k_Hz_c(:,:,K+1) = (2.0*epsilon0*ka_z + sigma_z*dt)/mu0;
    k_Hz_c(:,:,nz-K) = k_Hz_c(:,:,K+1);
    k_Hz_d(:,:,K+1) = -(2.0*epsilon0*ka_z - sigma_z*dt)/mu0;
    k_Hz_d(:,:,nz-K) = k_Hz_d(:,:,K+1);

    sigma_z = sigma_max*((PML - K - 0.5)/PML)^m;
    ka_z = 1.0 + (ka_max - 1.0)*((PML - K - 0.5)/PML)^m;
    k_Ez_c(:,:,K+1) = 2.0*epsilon0*ka_z + sigma_z*dt;
    k_Ez_c(:,:,nz-K-1) = k_Ez_c(:,:,K+1);
    k_Ez_d(:,:,K+1) = -(2.0*epsilon0*ka_z - sigma_z*dt);
    k_Ez_d(:,:,nz-K-1) = k_Ez_d(:,:,K+1);
    k_Hx_a(:,:,K+1) = (2.0*epsilon0*ka_z - sigma_z*dt)/...
    (2.0*epsilon0*ka_z + sigma_z*dt);
    k_Hx_a(:,:,nz-K-1) = k_Hx_a(:,:,K+1);
    k_Hx_b(:,:,K+1) = 1.0/(2.0*epsilon0*ka_z + sigma_z*dt);
    k_Hx_b(:,:,nz-K-1) = k_Hx_b(:,:,K+1);
    k_By_a(:,:,K+1) = (2.0*epsilon0*ka_z - sigma_z*dt)/...
    (2.0*epsilon0*ka_z + sigma_z*dt);
    k_By_a(:,:,nz-K-1) = k_By_a(:,:,K+1);
    k_By_b(:,:,K+1) = 2.0*epsilon0*dt/(2.0*epsilon0*ka_z + sigma_z*dt);
    k_By_b(:,:,nz-K-1) = k_By_b(:,:,K+1);
    
end

%% Main 3D FDTD+UPML routine (operates with 'singles' to increase speed) 
hhh = waitbar(0, 'Calculations in progress...');
tic;
for T=0:(number_of_iterations-1)
    %% Calculate Fx -> Gx -> Ex
    I = 1:nx; J = 2:ny; K = 2:nz;
    Fx_r = Fx(I,J,K);
    Fx(I,J,K) = Ka(I,J,K).*Fx(I,J,K) + Kb(I,J,K).*...
    ((Hz(I,J,K) - Hz(I,J-1,K))/dy - (Hy(I,J,K) - Hy(I,J,K-1))/dz);
    Gx_r = Gx(I,J,K);
    Gx(I,J,K) = k_Gx_a(I,J,K).*Gx(I,J,K) + k_Gx_b(I,J,K).*(Fx(I,J,K) - Fx_r);
    Ex(I,J,K) = k_Ex_a(I,J,K).*Ex(I,J,K) + k_Ex_b(I,J,K).*...
    (k_Ex_c(I,J,K).*Gx(I,J,K) + k_Ex_d(I,J,K).*Gx_r);

    %% Calculate Fy -> Gy -> Ey
    I = 2:nx; J = 1:ny; K = 2:nz;
    Fy_r = Fy(I,J,K);
    Fy(I,J,K) = Ka(I,J,K).*Fy(I,J,K) + Kb(I,J,K).*...
    ((Hx(I,J,K) - Hx(I,J,K-1))/dz - (Hz(I,J,K) - Hz(I-1,J,K))/dx);
    Gy_r = Gy(I,J,K);
    Gy(I,J,K) = k_Gy_a(I,J,K).*Gy(I,J,K) + k_Gy_b(I,J,K).*(Fy(I,J,K) - Fy_r);
    Ey(I,J,K) = k_Ey_a(I,J,K).*Ey(I,J,K) + k_Ey_b(I,J,K).*...
    (k_Ey_c(I,J,K).*Gy(I,J,K) + k_Ey_d(I,J,K).*Gy_r);

    %% Calculate Fz -> Gz -> Ez
    I = 2:nx; J = 2:ny; K = 1:nz;
    Fz_r = Fz(I,J,K);
    Fz(I,J,K) = Ka(I,J,K).*Fz(I,J,K) + Kb(I,J,K).*...
    ((Hy(I,J,K) - Hy(I-1,J,K))/dx - (Hx(I,J,K) - Hx(I,J-1,K))/dy);
    Gz_r = Gz(I,J,K);
    Gz(I,J,K) = k_Gz_a(I,J,K).*Gz(I,J,K) + k_Gz_b(I,J,K).*(Fz(I,J,K) - Fz_r);
    Ez(I,J,K) = k_Ez_a(I,J,K).*Ez(I,J,K) + k_Ez_b(I,J,K).*...
    (k_Ez_c(I,J,K).*Gz(I,J,K) + k_Ez_d(I,J,K).*Gz_r);

    %% Source of vertical electric field Ez applied to the whole face 
    %% plane at the port 1
    if (T*dt<=10.0*t_half)
        Ex(:,:,nz-PML-1) = Source(t_half, 1.0, T*dt);
    end

    Ez_saved(1:(ny-2*PML-1), 1:(nx-2*PML-1), T+1) = ...
    Ez((PML+1):(nx-PML-1), (PML+1):(ny-PML-1), PML+4)';
    Ex_saved(1:(ny-2*PML-1), 1:(nx-2*PML-1), T+1) = ...
    Ex((PML+1):(nx-PML-1), (PML+1):(ny-PML-1), PML+4)';
    Ex_yz(T+1, 1 :(ny-2*PML-1), 1:(nz-2*PML-1)) = ...
    Ex(nx/2, (PML+1):(ny-PML-1), (PML+1):(nz-PML-1));

    %% Calculate Bx -> Hx
    I = 2:nx; J = 1:ny; K = 1:nz;
    Bx_r = Bx(I,J,K);
    Bx(I,J,K) = k_Bx_a(I,J,K).*Bx(I,J,K) + k_Bx_b(I,J,K).*...
    ((Ey(I,J,K+1) - Ey(I,J,K))/dz - (Ez(I,J+1,K) - Ez(I,J,K))/dy);
    Hx(I,J,K) = k_Hx_a(I,J,K).*Hx(I,J,K) + k_Hx_b(I,J,K).*...
    (k_Hx_c(I,J,K).*Bx(I,J,K) + k_Hx_d(I,J,K).*Bx_r)./Kc(I,J,K);

    %% Calculate By -> Hy
    I = 1:nx; J = 2:ny; K = 1:nz;
    By_r = By(I,J,K);
    By(I,J,K) = k_By_a(I,J,K).*By(I,J,K) + k_By_b(I,J,K).*...
    ((Ez(I+1,J,K) - Ez(I,J,K))/dx - (Ex(I,J,K+1) - Ex(I,J,K))/dz);
    Hy(I,J,K) = k_Hy_a(I,J,K).*Hy(I,J,K) + k_Hy_b(I,J,K).*...
    (k_Hy_c(I,J,K).*By(I,J,K) + k_Hy_d(I,J,K).*By_r)./Kc(I,J,K);

    %% Calculate Bz -> Hz
    I = 1:nx; J = 1:ny; K = 2:nz;
    Bz_r = Bz(I,J,K);
    Bz(I,J,K) = k_Bz_a(I,J,K).*Bz(I,J,K) + k_Bz_b(I,J,K).*...
    ((Ex(I,J+1,K) - Ex(I,J,K))/dy - (Ey(I+1,J,K) - Ey(I,J,K))/dx);
    Hz(I,J,K) = k_Hz_a(I,J,K).*Hz(I,J,K) + k_Hz_b(I,J,K).*...
    (k_Hz_c(I,J,K).*Bz(I,J,K) + k_Hz_d(I,J,K).*Bz_r)./Kc(I,J,K);

    %% Progress bar updates at every time percent
    if ( mod(T, ceil(number_of_iterations/100)) == 0 )
        waitbar((T+1)/number_of_iterations, hhh);
    end

end
toc
close(hhh);

%% Plots: signals, scattering parameters and Ez field between plates
figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'doublebuffer', 'on');    

%Electric field strength Ex at the center of fishnet
fig = figure(1); 
object = 'Ex_Ez_metal_surface.avi';
mov = avifile(object);
mov.fps = 5;

for T=40:5:number_of_iterations
    
    subplot(2,2,3);
    contour(1e6*dx*(1:(nx-2*PML-1)), 1e6*dy*(1:(ny-2*PML-1)), ...
    Index((PML+1):(nx-PML-1), (PML+1):(ny-PML-1), PML+4)', 1);
    hold on;
    subplot(2,2,3);
    pcolor(1e6*dx*(1:(nx-2*PML-1)), 1e6*dy*(1:(ny-2*PML-1)), double(Ex_saved(:, :, T)));
    title(['Ex on metal surface @time t = ', num2str(round(T*dt*1e15)),' fs'], 'FontSize', 8);
    xlabel('[um]', 'FontSize', 6);
    ylabel('[um]', 'FontSize', 6);
    shading interp;
    axis image;
    caxis([-1 1]);
    colorbar;
    drawnow;

    subplot(2,2,4);
    contour(1e6*dx*(1:(nx-2*PML-1)), 1e6*dy*(1:(ny-2*PML-1)), ...
    Index((PML+1):(nx-PML-1), (PML+1):(ny-PML-1), PML+4)', 1);
    hold on;
    subplot(2,2,4);
    pcolor(1e6*dx*(1:(nx-2*PML-1)), 1e6*dy*(1:(ny-2*PML-1)), double(Ez_saved(:, :, T)));
    title(['Ez on metal surface @time t = ', num2str(round(T*dt*1e15)),' fs'], 'FontSize', 8);
    xlabel('[um]', 'FontSize', 6);
    ylabel('[um]', 'FontSize', 6);
    shading interp;
    axis image;
    caxis([-1 1]);
    colorbar;
    drawnow;

    subplot(2,2,[1 2]);
    contour(1e6*dy*(1:(ny-2*PML-1)), 1e6*dz*(1:(nz-2*PML-1)), ...
    squeeze(Index(nx/2,(PML+1):(ny-PML-1),(PML+1):(nz-PML-1)))',1);
    hold on;
    subplot(2,2,[1 2]);
    Ex_yzp(:,:) = Ex_yz(T,:,:);
    pcolor(1e6*dy*(1:(ny-2*PML-1)), 1e6*dz*(1:(nz-2*PML-1)), double(Ex_yzp)');
    title(['Ex at central YZ plane @time t = ', num2str(round(T*dt*1e15)),' fs'], 'FontSize', 8);
    xlabel('[um]', 'FontSize', 6);
    ylabel('[um]', 'FontSize', 6);
    shading interp;
    axis image;
    caxis([-1 1]);
    colorbar;
    drawnow;


    pause(0.02);
    F = getframe(fig);
    mov = addframe(mov,F);
    
end

close(fig);
mov = close(mov);

return;
     