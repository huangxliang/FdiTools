function [Cy, G, CvecG, freqn, M] = SpectralAnal(y, u, options);
%
%		Estimate the output noise covariance matrix, the FRF and its
%		uncertainty via the spectral analysis method (cross- and autopower
%		spectra estimation).
%       If no input data is provided only the noise covariance matrix
%       is estimated.
%
%	function [Cy, G, CvecG, freqn, M] = SpectralAnal(y, u, options);
%
%	Output parameters
%
%		Cy          =	struct{'large', 'small'} containing the covariance matrix of the output noise
%                           Cy.large    =   output noise ovariance matrix at the full record frequency resolution
%                                           size ny x ny x N
%                           Cy.small    =   output noise ovariance matrix at the segment frequency resolution
%                                           size ny x ny x Ns
%
%		G           =	struct{'large', 'small'} containing the FRF
%                           G.large     =   FRF at the full record frequency resolution
%                                           size ny x nu x N
%                           G.small     =   FRF at the segment frequency resolution
%                                           size ny x nu x Ns 
%
%		CvecG       =	struct{'large', 'small'} containing cov(vec(FRF))
%                           CvecG.large	=   cov(vec(FRF)) at the full record frequency resolution
%                                           size (ny*nu) x (ny*nu) x N
%                           CvecG.small =   cov(vec(FRF)) at the segment frequency resolution
%                                           size (ny*nu) x (ny*nu) x Ns 
%
%       freqn       =   struct{'large', 'small'} containing the normalised frequency (divided by the sampling frequency)
%                           freqn.large	=   normalised frequency corresponding to the full record
%                                           size 1 x N 
%                           freqn.small	=   normalised frequency corresponding to one segment
%                                           size 1 x Ns 
%
%       M           =   equivalent number of independent experiments = number of block used for calculating the 
%                       (cross)power spectra 
%                       
%
%	Input parameters
%
%		y           =	output signal, size ny x N
%		u           =	input signal, size nu x N (optional)
%       options     =   struct('width', [], 'moment', [], 'window', [])
%                       options.width   =	width frequency window in DFT lines (w.r.t. the N samples) used to average Cy 
%                                           default value width = 1
%                       options.moment  =	moment (optional)
%                                               0 then Cy is of full rank (default) 
%                                               1 then expected value Cy^-1 exists
%                                               2 then second order moments Cy^-1 exist
%                       options.window  =   type of window used
%                                               'diff': difference window (default) 
%                                               'sine': half sine window 
%                                               'rect': rectangular window 
%                                               'hann': hanning window
%
% Rik Pintelon, 2005
% version July 15, 2009
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% verification of the structure options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% does options exist ?
try
    if isempty(options)
        options = struct('width', 1, 'moment', 0, 'window', 'diff');
    end
catch
    options = struct('width', 1, 'moment', 0, 'window', 'diff');
end

% is options a structure ?
if ~isstruct(options)
    options = struct('width', 1, 'moment', 0, 'window', 'diff');
end % if not structure

% does the structure options contain the correct fields ?
if ~isfield(options, 'width')
    options.width = 1;
end

if ~isfield(options, 'moment')
    options.moment = 0;
end

if ~isfield(options, 'window')
    options.window = 'diff';
end

try
    if isempty(u)
        u = [];
    end
catch
    u = [];
end


% internal variables
width = options.width;
mm = options.moment;
TheWindow = lower(options.window);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialisation variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ny, N] = size(y);
nu = size(u, 1);

% minimal number of required segments (blocks) such that
% Cy of full rank: mm = 0
% E{Cy^-1} exists: mm = 1
% E{Cy^-2} exists: mm = 2
M = ny + nu + mm;

% number of points in a segment
Ns = floor(N/M);

% normalised frequencies
freqn = struct('small', [], 'large', []);
freqn.large = [0:1:N-1]/N;                      % normalised DFT frequencies of the full record

switch TheWindow
    
    case 'diff'
        Syy = zeros(ny, ny, Ns-1);
        Suu = zeros(nu, nu, Ns-1);
        Syu = zeros(ny, nu, Ns-1);
        Cy = struct('small', zeros(ny, ny, Ns-1), 'large', zeros(ny, ny, N));
        G = struct('small', zeros(ny, nu, Ns-1), 'large', zeros(ny, nu, N));
        CvecG = struct('small', zeros(ny*nu, ny*nu, Ns-1), 'large', zeros(ny*nu, ny*nu, N));
        freqn.small = [0.5:1:Ns-1.5]/Ns;        % normalised DFT frequencies of one segment after differentiation
        
    case {'sine', 'rect', 'hann'}
        Syy = zeros(ny, ny, Ns);
        Suu = zeros(nu, nu, Ns);
        Syu = zeros(ny, nu, Ns);
        Cy = struct('small', zeros(ny, ny, Ns), 'large', zeros(ny, ny, N));
        G = struct('small', zeros(ny, nu, Ns), 'large', zeros(ny, nu, N));
        CvecG = struct('small', zeros(ny*nu, ny*nu, Ns), 'large', zeros(ny*nu, ny*nu, N));
        freqn.small = [0:1:Ns-1]/Ns;            % normalised DFT frequencies of one segment
       
end % switch TheWindow


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate (cross)-power spectra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reshape signals in M segments of Ns points each
y = y(:,1:M*Ns);
y = reshape(y.',Ns, M, ny);
y = permute(y, [3, 2, 1]);
if nu > 0
    u = u(:,1:M*Ns);
    u = reshape(u.',Ns, M, nu);
    u = permute(u, [3, 2, 1]);
end % if nu > 0

switch TheWindow
    
    case 'diff'
        % take fft of each segment 
        Y = fft(y, [], 3)/sqrt(Ns);		% ny x M x Ns matrix        
        % differentiate to suppress the leakage (plant and noise transients) 
        Y = diff(Y, 1, 3);
        
        if nu > 0
            U = fft(u, [], 3)/sqrt(Ns);
            U = diff(U, 1, 3);
        end % nu > 0
        
    case {'sine', 'rect', 'hann'}
        
        wind = ones(1, 1, Ns);
        switch TheWindow
            case 'sine'
                wind(1, 1, :) = sin(pi*freqn.small);
            case 'hann'
                wind(1, 1, :) = hanning(Ns, 'periodic').';
        end % switch
        
        % fft windowed signals to suppress the leakage (plant and noise transients) 
        y = y.*repmat(wind, [ny, M, 1]);
        Y = fft(y, [], 3)/sqrt(Ns);         % ny x M x Ns matrix
        if nu > 0
            u = u.*repmat(wind, [nu, M, 1]);
            U = fft(u, [], 3)/sqrt(Ns);
        end % if nu > 0
        
end % switch

% estimate output power spectrum = sum true output power + noise power
for ii = 1:ny
	for jj = 1:ii
		Syy(ii, jj, :) = mean(Y(ii, :, :).*conj(Y(jj, :, :)), 2);
		Syy(jj, ii, :) = conj(Syy(ii, jj, :));
	end % jj
end % ii

if nu > 0
    
    % estimate input power spectrum
    for ii = 1:nu
        for jj = 1:ii
            Suu(ii, jj, :) = mean(U(ii, :, :).*conj(U(jj, :, :)), 2);
            Suu(jj, ii, :) = conj(Suu(ii, jj, :));
        end % jj
    end % ii

    % estimate cross power spectrum
    for ii = 1:ny
        for jj = 1:nu
            Syu(ii, jj, :) = mean(Y(ii, :, :).*conj(U(jj, :, :)), 2);
        end % jj
    end % ii

end % nu > 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate output noise power spectrum, FRF, and cov(vec(FRF)) %
% the appropriate scaling is done at the end of the program    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = length(freqn.small);
if nu > 0    
    for kk = 1:F
        invSuu = inv(squeeze(Suu(:,:,kk)));
        G.small(:,:,kk) = Syu(:,:,kk)*invSuu;
        Cy.small(:,:,kk) = Syy(:,:,kk) - G.small(:,:,kk)*Syu(:,:,kk)';
        % remove imaginary part on diagonal
        Cy.small(:,:,kk) = Cy.small(:,:,kk) - diag(sqrt(-1)*imag(diag(Cy.small(:,:,kk))));
        CvecG.small(:,:,kk) = kron(invSuu, Cy.small(:,:,kk));
    end % kk
else
    Cy.small = Syy;
end % if nu > 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% averaging estimated covariance matrices over neighbouring % 
% DFT lines in the coarse frequency resolution              %
% 2*L + 1 = filter width                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = ceil((width/M-1)/2);                                            % width is expressed in the LargeDFT resolution

% noise covariance matrix
[Cysmall, SmallDFTf] = SmoothMatrix(Cy.small, freqn.small.', L);    % SmallDFTf: grid after differentiation and smoothing

% interpolate Cy to grid Ns frequencies if the smoothing is active
if L > 0
    Cy.small = zeros(size(Cy.small));
    for ii = 1:ny
        for jj = 1:ii
            Cy.small(ii, jj, :) = interp1(SmallDFTf, squeeze(Cysmall(ii,jj,:)), freqn.small.', 'nearest', 'extrap');
            Cy.small(jj, ii, :) = conj(Cy.small(ii, jj, :));
        end % jj
    end % ii    
end % if L>0

if nu > 0    
    % covariance matrix vec(G)
    [CvecGsmall, SmallDFTf] = SmoothMatrix(CvecG.small, freqn.small.', L);    % SmallDFTf: grid after differentiation and smoothing

    % interpolate vec(G) to grid Ns frequencies if the smoothing is active
    if L > 0
        CvecG.small = zeros(size(CvecG.small));
        for ii = 1:ny*nu
            for jj = 1:ii
                CvecG.small(ii, jj, :) = interp1(SmallDFTf, squeeze(CvecGsmall(ii,jj,:)), freqn.small.', 'nearest', 'extrap');
                CvecG.small(jj, ii, :) = conj(CvecG.small(ii, jj, :));
            end % jj
        end % ii    
    end % if L>0
end % nu > 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolate to the large frequency resolution %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interpolate Cy to grid N frequencies
for ii = 1:ny
	for jj = 1:ii
		Cy.large(ii, jj, :) = interp1(SmallDFTf, squeeze(Cysmall(ii,jj,:)), freqn.large.', 'nearest', 'extrap');
		Cy.large(jj, ii, :) = conj(Cy.large(ii, jj, :));
	end % jj
end % ii

if nu > 0
    
    % interpolate covariance matrix vec(G) to grid N frequencies
    for ii = 1:ny*nu
        for jj = 1:ii
            CvecG.large(ii, jj, :) = interp1(SmallDFTf, squeeze(CvecGsmall(ii,jj,:)), freqn.large.', 'nearest', 'extrap');
            CvecG.large(jj, ii, :) = conj(CvecG.large(ii, jj, :));
        end % jj
    end % ii

    % interpolate G to grid N frequencies
    for ii = 1:ny
        for jj = 1:nu
            G.large(ii, jj, :) = interp1(freqn.small.', squeeze(G.small(ii,jj,:)), freqn.large.', 'nearest', 'extrap');
        end % jj
    end % ii
    
end % nu > 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correct the covariances for the bias and the diff operation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch TheWindow
    case 'diff'
        TheScaleCy = M/(M-nu)/2;
        TheScaleCvecG = 2*TheScaleCy;
    case 'sine'
        TheScaleCy = 2*M/(M-nu);
        TheScaleCvecG = TheScaleCy/2;
    case 'rect'
        TheScaleCy = M/(M-nu);
        TheScaleCvecG = TheScaleCy;
    case 'hann'
        TheScaleCy = 8/3*M/(M-nu);
        TheScaleCvecG = TheScaleCy*3/8;
end % switch
Cy.small = Cy.small*TheScaleCy;
Cy.large = Cy.large*TheScaleCy;
if nu > 0
    CvecG.small = CvecG.small*TheScaleCvecG/M;
    CvecG.large = CvecG.large*TheScaleCvecG/M;
end % nu > 0
