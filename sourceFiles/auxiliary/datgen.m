function [W,H0,U0_tot,Y0_tot,Um_tot,Ym_tot,Y0_ind] = datgen(NF,NU,NY,DOF,NM,noiselevel)
% [W,H0,U0_tot,Y0_tot,Um_tot,Ym_tot,Y0_ind] = datgen(NF,NU,NY,DOF,NM,noiselevel)
% GENerates Input/Output DATa.
%   * input arguments [op = optional parameter, m = matrix]
%                      - NF  = number of spectral lines
%                      - NU  = number of input
%                      - NY  = number of outputs
%                      - DOF [op] = degrees of freedom of the system
%                              [Remark: DOF >= max(NU,NY)]
%                                    - default:  DOF = max(NU,NY)
%                      - NM  [op] = number of independent measurements
%                      - noiselevel [op] = the standard deviation of the
%                        i/o noise is proportional to "noiselevel":
%                                    - default:  noiselevel = 1
%                                    - no noise: noiselevel = 0
%   * output arguments [cv = column vector, m = matrix]
%                      - W  [cv] = angular frequencies (max(W) = 2)
%                      - H0 [m] = the "true" frequency response functions
%                                 - size(H0) = (NU*NY) x NF
%                                 - H0k = reshape(H0(:,k),NU,NY)
%                      - Um_tot  [m]  = noisy input Fourier coeffs. 
%                                  - size(U) = NU x (NF*NM)
%                      - Ym_tot  [m]  = noisy output Fourier coeffs. 
%                                  - size(Y) = NY x (NF*NM)
%                      - U0_tot  [m]  = exact input Fourier coeffs. 
%                                  - size(U) = NU x (NF*NM)
%                      - Y0_tot  [m]  = exact output Fourier coeffs. 
%                                  - size(Y) = NY x (NF*NM)
%                      - Y0_ind  [m]  = exact output Fourier coeffs. for SI excitation
%                                  - size(Y) = (NY*NU) x (NF*NM)
%

if nargin<4, DOF = max(NU,NY); end
if nargin<5, NM = 1; end
if nargin<6, noiselevel = 1; end

if DOF<max(NU,NY), error('DOF MUST BE >= max(NU,NY)'), end

randn('seed',0)
rand('seed',0)

NZ = NU + NY;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a mass (M), damping (D), and stiffness (K) matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = eye(DOF);
%
rand('uniform');
D = rand(DOF); D = D*D';
D = D*(0.005/DOF);
%
K = rand(DOF);
K= K*K';
[V,DD]=eig(K);
DD = diag(linspace(0.7,1.9,DOF).^2);
K = V*DD*inv(V);
%K = K*(1/DOF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = linspace(0.5,2,NF)';
Y0_tot = zeros(NY,NM*NF);
U0_tot = zeros(NU,NM*NF);
Ym_tot = zeros(NY,NM*NF);
Um_tot = zeros(NU,NM*NF);
Y0_ind = zeros(NY*NU,NM*NF);

H0 = zeros(NU*NY,NF);
B = zeros(DOF,NU); B(1:NU,:) = eye(NU);

randn('seed',sum(100*clock));
rand('seed',sum(100*clock));

for k = 1:NF
        Sk = sqrt(-1)*W(k);
        Hk = inv(M*Sk^2 + D*Sk + K)*B*Sk^2;
        Hk = Hk(1:NY,:);
        HH = Hk.';
        H0(:,k) = HH(:);

       
        for n = 1:NM

           U0k = randn(NU,1)+sqrt(-1)*randn(NU,1);
           U0k(1) = 2*U0k(1);

           % input noise at frequency W(k)
           NoiseU = (randn(NU,1)+sqrt(-1)*randn(NU,1))/sqrt(2)*0.01*noiselevel;
           % output noise at frequency W(k)
           NoiseY = (randn(NY,1)+sqrt(-1)*randn(NY,1)).*[ones(NY,1)+abs(Hk*U0k)]/sqrt(2)*0.05*noiselevel;

           Umk = U0k + NoiseU;
           Y0k = Hk*U0k;
           Ymk = Hk*U0k + NoiseY;
           Um_tot(:,(n-1)*NF+k) =  Umk;
           Ym_tot(:,(n-1)*NF+k) =  Ymk;
           U0_tot(:,(n-1)*NF+k) =  U0k;
           Y0_tot(:,(n-1)*NF+k) =  Y0k;
           for (nn = 1:NU)
              Yind = Hk*[zeros(1:nn-1,1);U0k(nn);zeros(1:NU-nn,1)];
              Y0_ind(NY*(nn-1)+1:NY*nn,(n-1)*NF+k) =  Yind;
           end;
           
        end

        
end

