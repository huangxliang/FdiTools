function plotfrfs(FRFs,p,q,freq,ssG,amp_mode,ph_mode,ph_wrap,addphase)
%
% plotfrfs([FRF1 FRF2 ...],p,q,freq,ssG,amp_mode,ph_mode,ph_wrap,addphase)
%
% plots different transfertfunctions FRF1 ('-'), FRF2 ('--'), ...
% (in 2D matrix).  In addition, poles and zeros of ssG can be plotted
%
% [FRF1 FRF2 ...]     = matrix containing the FRFs
%                       nr of rows    = Nfreq = number of frequency lines
%                       nr of columns = (q*p)*(nrofFRFs)
% p        = number of inputs for each FRF
% q        = number of outputs for each FRF
% freq     = frequency axis vector (size: Nfreq * 1)
% ssG      = ss object containing state space model, default: []
% amp_mode = string indicating amplitude plot format: possibilities:
%            'linlin', 'linlog', 'lindb'(default)
%            'loglin', 'loglog', 'logdb'
% ph_mode  = string indicating plot format for phase: possibilities:
%            'log', 'lin'(default)
% ph_wrap  = string indicating phase wrapping: default for unwrapped
%            angle is used.  Enter 'w' to use wrapped angle. 
% addphase = row-vector, (nr of columns = nr of columns FRF), when plotting, addphase(i) times 360 degrees
%            will be add to the phase of each FRF (to compare phase of different FRFs. 
%            default addphase=[0,0,...,0];
%

min_nr_arg = 4;

quote = char(39);

if nargin == min_nr_arg,
   ssG = [];
   amp_mode = 'lindb';
   ph_mode = 'lin';
   ph_wrap = 'u';
   addphase=zeros(1,size(FRFs,2));
elseif nargin == min_nr_arg + 1,
   amp_mode = 'lindb';
   ph_mode = 'lin';
   ph_wrap = 'u';
   addphase=zeros(1,size(FRFs,2));
elseif nargin == min_nr_arg + 2,
   ph_mode = 'lin';
   ph_wrap = 'u';
   addphase=zeros(1,size(FRFs,2));
elseif nargin == min_nr_arg + 3,
   ph_wrap = 'u';
   addphase=zeros(1,size(FRFs,2));
elseif  nargin == min_nr_arg + 4,
    addphase=zeros(1,size(FRFs,2));
end %if

if (size(addphase,1) ~= 1) | (size(addphase,2) ~= size(FRFs,2))
    error('addphase has wrong size, it must be 1*nr of columns of FRFs')
end

if mod(size(FRFs,2),p*q) ~= 0,
   msg = sprintf('Nr of columns of FRFs (%d) should be a multiple of nr.inp*nr. outp (%d)',size(FRFs,2),p*q);
   disp(msg);
   return;
end; %if

NrOfFRFs = size(FRFs,2)/(p*q);
maxnr = 4;
linetypes = '- ---.: ';
colors = 'brmg';

figure;

for n = 1:NrOfFRFs,
   index = mod(n,maxnr);
   if (index == 0),
      index = maxnr;
   end; %if index
   format = cat(2,colors(index),linetypes((index-1)*2+1:2*index));
   
   for k = 1:p*q,
      amp_place = floor((k-1)/p) *p + k;
      nr_in = mod(k,p);
      if nr_in == 0,
         nr_in = p;
      end %if
      nr_out = ceil(k/p);
      column2plot = (n-1)*p*q + k;
      
      subplot(2*q,p,amp_place);
      switch amp_mode,
      case 'linlin',
         plot(freq,abs(FRFs(:,column2plot)),format);
         ylabel('Amplitude');
      case 'linlog',
         semilogy(freq,abs(FRFs(:,column2plot)),format);
         ylabel('Log(Amplitude)');
      case 'lindb',
         plot(freq,20*log10(abs(FRFs(:,column2plot))),format);
         ylabel('Amplitude [dB]')
      case 'loglin',
         semilogx(freq,abs(FRFs(:,column2plot)),format);
         ylabel('Amplitude');
      case 'loglog',
         loglog(freq,abs(FRFs(:,column2plot)),format);
         ylabel('Log(Amplitude)');
      case 'logdb',
         semilogx(freq,20*log10(abs(FRFs(:,column2plot))),format);
         ylabel('Amplitude [dB]');
      otherwise,
         plot(freq,20*log10(abs(FRFs(:,column2plot))),format);
         ylabel('Amplitude [dB]')
      end %switch
            
      %a=sprintf('input%d, output%d',nr_in, nr_out);
      %title(a);
      xlabel('Frequency [Hz]');
      grid on;
      zoom on;
      hold on;
      
      subplot(2*q,p,amp_place+p);
      switch ph_mode,
      case 'log',
         switch ph_wrap,
         case 'w',
            semilogx(freq,(angle(FRFs(:,column2plot)))*180/pi,format);
         otherwise,
            semilogx(freq,unwrap(angle(FRFs(:,column2plot)))*180/pi+360*addphase(1,column2plot),format);
         end %switch ph_wrap
      otherwise,
         switch ph_wrap,
         case 'w',
            plot(freq,(angle(FRFs(:,column2plot)))*180/pi,format);
         otherwise,
            plot(freq,unwrap(angle(FRFs(:,column2plot)))*180/pi,format);
         end %switch ph_wrap
      end %switch ph_mode
      
      xlabel('Frequency [Hz]');
      ylabel('Phase [deg]');
      grid on;
      zoom on;
      hold on;
   end %for k
end; %for n

% Plot 'x' for every pole at corresponding frequency
if ~isempty(ssG),
   subplot(211);
   v = axis;
   [A,B,C,D,ts] = ssdata(ssG);
   [Z,P,K] = ss2zp(A,B,C,D);
   for ZP = 1:2,
      PZ_freqs = [];
      if ZP == 1,  % process poles
         stage = 'poles';
         PZ = P;
      else  % process zeros
         stage = 'zeros';
         PZ = Z;
      end; %if ZP
      if ts ~= 0,  % discrete-time system
         PZ = ln(PZ)/ts;
      end; %if ts
      PZ_freqs = 1/(2*pi)*(sqrt(real(PZ).^2+imag(PZ).^2));
      
      k = 1;
      while k <= length(PZ_freqs),
         aantal = 1;
         if imag(PZ(k)) ~= 0,
            aantal = 2;
            k = k+1;
         end;
         
         if strcmp(stage,'poles'),
            letter = 'bx';
         else
            if real(PZ(k)<0),  % minimum-phase zero
               letter = 'go';
            else
               letter = 'ro';
            end;
         end;  % if strcmp(stage
         
         cmd = sprintf('plot(%f,v(3)+(v(4)-v(3))*0.01,%s%s%s)',PZ_freqs(k),quote,letter,quote);
         eval(cmd);
         if aantal == 2,
            cmd = sprintf('plot(%f,v(3)+(v(4)-v(3))*0.05,%s%s%s)',PZ_freqs(k),quote,letter,quote);
            eval(cmd);
         end; % if aantal
         k = k+1;
      end; %while
   end; % for ZP
   
   hold off

	% re-adjust axis on subplot(212)
	vnew211 = axis;

	subplot(212);
	v212 = axis;
	axis([vnew211(1) vnew211(2) v212(3) v212(4)]);
	hold off;

end; %if isempty   
