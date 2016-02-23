% Read me first file: use of the functions for nonparametric characterization of MIMO systems with nu inputs and ny outputs 
% Examples illustrating the use of the functions can be found in the folder "examples". To run these examples one
% also needs the "DesignMultisineExcitations" toolbox. 
%
%%   1. ARBITRARY EXCITATIONS - time signals
%
%           function [CZ, Z, G, CvecG, dof, CL] = ArbLocalPolyAnal(data, method); 
%
%       Starting from (noisy) input and noisy output time signals the routine calculates the ny x nu frequency response matrix (FRM)  
%       G and its covariance matrix CvecG = Cov(vec(G)). In addition the following sample means and sample covariance are calculated: 
%           (i)     CZ.n:            the (ny+nu) x (ny+nu) noise (noise + nonlinear distortions) covariance matrix of the original 
%                                    input-output DFT spectra 
%           (ii)    Z.m, CZ.m:       the (ny+nu) x 1 sample means and the (ny+nu) x (ny+nu) sample covariance of the sample mean 
%                                    without transient removal 
%           (ii)    Z.m_nt, CZ.m_nt: the (ny+nu) x 1 sample means and the (ny+nu) x (ny+nu) sample covariance of the sample mean 
%                                    with transient removal 
%
%       The method can handle data from different experiments measured under similar operational conditions. 
%
%       If besides the input-output data also a known reference signal is provided, then the routine handles automatically the 
%       errors-in-variables case (noisy input and noisy output observations). 
%
%       If the reference signal is not available and the input is noisy then:
%           1.  The FRF G estimate is biased. Similarly to the spectral analysis method, the bias of the local polynomial FRF 
%               estimate is proportional to the square of the input noise-to-signal ratio. 
%           2.  The noise covariance CY.n is then the sum of the output noise covariance and the input noise covariance reflected 
%               to the output (CY.n = Cov(NY-G*NU)).   
%
%       dof stands for the actual degrees of freedom of the covariance estimates, and ± CL is the correlation length over the
%       frequency of the sample mean and sample covariance. 
%
%       The noise and system transients in the data are removed via a local polynomial approximation. The FRM and the covariances 
%       are calculated via a local polynomial approximation of the FRM. No distinction can be made between the noise and the
%       stochastic nonlinear distortions.
%
%       Usage:
%           1.  The sample means and sample covariances Z.m, CZ.m and Z.m_nt, CZ.m_nt are used in the SML cost function 
%               (see MIMO_MaximumLikelihood toolbox). 
%           2.  Z.m, CZ.n or Z.m_nt, CZ.n are used in the Cramér-Rao lower bound routine (see MIMO_MaximumLikelihood toolbox) 
%               for estimating the covariance of the parametric plant model.  
%
%
%       References:
%
%           Pintelon, R., J. Schoukens, G. Vandersteen, and K. Barbé (2010). Estimation of nonparametric noise and FRF models for 
%           multivariable systems - Part I: theory, Mechanical Systems and Signal Processing, vol. 24, no. 3, pp. 573-595.
%
%           Pintelon, R., J. Schoukens, G. Vandersteen, and K. Barbé (2010). Estimation of nonparametric noise and FRF models for 
%           multivariable systems - Part II: extensions, applications, Mechanical Systems and Signal Processing, vol. 24, no. 3, 
%           pp. 596-616.
%
%           Pintelon, R., and J. Schoukens (2012). System Identification: A Frequency Domain Approach, second edition, 
%           IEEE Press-Wiley, Piscataway (USA). 
%
%
%%   1. ARBITRARY EXCITATIONS - DFT spectra
%
%           [CY, Y, TY, G, CvecG, dof, CL] = LocalPolyAnal(data, method);
%
%       Starting from the known input DFT spectrum and the noisy output DFT spectrum the routine calculates the ny x nu frequency   
%       response matrix (FRM) G and its covariance matrix CvecG = Cov(vec(G)). In addition the following sample means and sample 
%       covariance are calculated: 
%           (i)     CY.n:            the ny x ny noise (noise + nonlinear distortions) covariance matrix of the original 
%                                    input-output DFT spectra 
%           (ii)    Y.m, CY.m:       the ny x 1 sample means and the ny x ny sample covariance of the sample mean without transient removal 
%           (ii)    Y.m_nt, CY.m_nt: the ny x 1 sample means and the ny x ny sample covariance of the sample mean with transient removal 
%
%       To handle the errors-in-variables case one should define the reference signal R as the input and Z = [Y; U] as the output. The 
%       FRM from input U to output Y and its covariance is then calculated using 
%
%           [Geiv, CvecGeiv] = FRF_EIV(G, CvecG);
%
%       Usage:
%           1. Handle frequency domain signals.
%           2. Calculation of the generalized sample means and sample covariances as required by the MIMO_MaximumLikelihood toolbox. 
%           2. Preprocessing to generate starting values for Box-Jenkins models using the MIMO_MaximumLikelihood toolbox.
%
%
%       References:
%
%           Pintelon, R., J. Schoukens, G. Vandersteen, and K. Barbé (2010). Estimation of nonparametric noise and FRF models for 
%           multivariable systems - Part I: theory, Mechanical Systems and Signal Processing, vol. 24, no. 3, pp. 573-595.
%
%           Pintelon, R., J. Schoukens, G. Vandersteen, and K. Barbé (2010). Estimation of nonparametric noise and FRF models for 
%           multivariable systems - Part II: extensions, applications, Mechanical Systems and Signal Processing, vol. 24, no. 3, 
%           pp. 596-616.
%
%           Pintelon, R., and J. Schoukens (2012). System Identification: A Frequency Domain Approach, second edition, 
%           IEEE Press-Wiley, Piscataway (USA). 
%
%
%%   2. PERIODIC EXCITATIONS - FAST METHOD - time signals
%
%           [CZ, Z, freq, G, CvecG, dof, CL] = FastLocalPolyAnal(data, method)
%
%       Starting from noisy input, noisy output time signals of P consecutive periods of one MIMO experiment, the ny x nu FRM  
%       G, its noise covariance matrix CvecG.n, its total covariance CvecG.NL, and the covariance of the stochastic nonlinear 
%       distortions are calculated at the excited frequencies. In addition the following sample means and sample covariance of 
%       the sample means of the output-input DFT spectra are calculated: 
%           (i)     mean over the periods with noise (and system) transient removed Z.n (output and and input on top of each other), 
%                   its sample noise covariance CZ.n, its the total covariance, and the covariance of the stochastic nonlinear 
%                   distortions are calculated. 
%       	(ii)    the sample mean of the output-input DFT spectra over the stochastic nonlinear distortions  
%                   Z.m_NL and its sample total covariance CZ.m_NL. 
%
%       The method can handle the transient response to the periodic inputs. 
%
%       If the reference signal is not available and the input is noisy then:
%           1.  The FRM G estimate is biased. Similarly to the spectral analysis method, the bias of the local polynomial FRF 
%               estimate is proportional to the square of the input noise-to-signal ratio. 
%           2.  Since they are based on the FRM estimate G, the mean output-input spectra Z.m_NL over the stochastic nonlinear 
%               distortions and their noise covariance CZ.m_NL are also biased. The total covariance CZ.m_NL is then the sum of
%               the output contributions (NY+YS) and the input contributions (NU+US) reflected to the output: 
%               CZ.m_NL = alpha*(Cov(NY-G*NU)+Cov(YS-G*US)), where alpha < 1 accounts for the covariance reduction of the local
%               polynomial approximation of the FRM. 
%           3.  The mean output-input spectra Z.n over the periods and their noise covariance CZ.n remain unbiased. The total 
%               covariance CZ.NL is then the sum of the output contributions (NY+YS) and the input contributions  (NU+US) 
%               reflected to the output: CZ.NL = Cov(NY-G*NU)+Cov(YS-G*US).   
%
%       dof stands for the actual degrees of freedom of the covariance estimates, and ± CL is the correlation length over the
%       EXCITED frequencies of the sample mean and sample covariance.
%
%       The noise (and system) transients in the data are removed via a local polynomial approximation. The FRM, the total 
%       covariances, and the sample means of the output-input spectra over the stochastic nonlinear distortions are calculated 
%       via a local polynomial approximation of the FRM. 
%
%       Usage:
%           1.  The sample means and sample covariances Z.n, CZ.n (noise only) and Z.m_NL, CZ.m_NL (noise + stoch. NL distort.) 
%               are used in the SML cost function (see MIMO_MaximumLikelihood toolbox). 
%           2.  Z.m_NL, CZ.NL or Z.n, CZ.n are used in the Cramér-Rao lower bound routine (see MIMO_MaximumLikelihood toolbox) 
%               for estimating the covariance of the parametric plant model.   
%
%       References:
%
%           Pintelon, R., K. Barbé,  G. Vandersteen, and J. Schoukens (2011). Improved (non-)parametric identification of dynamic 
%           systems excited by periodic signals, Mechanical Systems and Signal Processing, vol. 25, no. 7, pp. 2683-2704. 
%
%           Pintelon, R., G. Vandersteen, J. Schoukens, and Y. Rolain (2011). Improved (non-)parametric identification of dynamic 
%           of dynamic systems excited by periodic signals - The multivariate case, Mechanical Systems and Signal Processing, 
%           vol. 25, no. 8, pp. 2892-2922.
%
%           Pintelon, R., and J. Schoukens (2012). System Identification: A Frequency Domain Approach, second edition, 
%           IEEE Press-Wiley, Piscataway (USA).
%
%
%%   3. PERIODIC EXCITATIONS - ROBUST METHOD - time signals
%
%           [CZ, Z, freq, G, CvecG, dof, CL] = RobustLocalPolyAnal(data, method) 
%
%       Starting from noisy input, noisy output time signals of M independent realisations of P consecutive periods of nu 
%       independent MIMO experiments, the ny x nu FRM G, its noise covariance matrix CvecG.n, its total covariance CvecG.NL, 
%       and the covarianve of the stochastic nonlinear distortions w.r.t. one realisation are calculated at the excited 
%       frequencies. In addition the sample mean with noise (and system) transient removed Z (output and and input on top of each 
%       other), its  noise covariance CZ.n, its the total covariance, and the covariance of the stochastic nonlinear distortions 
%       are calculated.
%
%       The method can handle the transient response to the periodic inputs. 
%
%       dof stands for the actual degrees of freedom of the covariance estimates, and ± CL is the correlation length over the
%       EXCITED frequencies of the sample mean and sample covariance.
%
%       The noise (and system) transients are removed via a local polynomial approximation. No local polynomial approximation 
%       of the FRM is made. The total covariances are calculated over the M independent realisations. 
%
%       Usage:
%           1.  The sample mean Z and sample covariance CZ.n or CZ.NL are used in the SML cost function 
%               (see MIMO_MaximumLikelihood toolbox). 
%           2.  Z and CZ.n or CZ.NL are used in the Cramér-Rao lower bound routine (see MIMO_MaximumLikelihood toolbox) 
%               for estimating the covariance of the parametric plant model. 
%
%       References:
%
%           Pintelon, R., K. Barbé,  G. Vandersteen, and J. Schoukens (2011). Improved (non-)parametric identification of dynamic 
%           systems excited by periodic signals, Mechanical Systems and Signal Processing, vol. 25, no. 7, pp. 2683-2704. 
%
%           Pintelon, R., G. Vandersteen, J. Schoukens, and Y. Rolain (2011). Improved (non-)parametric identification of dynamic 
%           systems excited by periodic signals - The multivariate case, Mechanical Systems and Signal Processing, vol. 25, no. 8, 
%           pp. 2892-2922.
%
%           Pintelon, R., and J. Schoukens (2012). System Identification: A Frequency Domain Approach, second edition, 
%           IEEE Press-Wiley, Piscataway (USA).
%
%
%%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, 5 November 2009
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% version 20 October 2011
%


