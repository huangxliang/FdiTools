function [fqlog, df, cdmax, freqind] = lin2log(freqv, rf)
%LIN2QLOG Quasilogarithmic subset of a linear frequency grid.
%
%       [fqlog, df, cdmax, freqind] = lin2log(freqv, rf)
%
%       The routine starts from a linear frequency grid, given in freqv, and
%       selects a quasi-logarithmic set, providing that the ratio of successive
%       frequencies is about  rf (or larger, if the frequency vector is not
%       dense enough).
%
%       Algorithm: starting from the first non-zero frequency point, the next
%       point f2 will be the one in freqv, closest to f1*rf, and larger than
%       f1. The points between f1 and f2 are deleted. This is repeated until
%       the end of the file: the last point will only be taken if the last
%       frequency is larger than f(n-1)*sqrt(rf).
%
%       Output arguments:
%       fqlog   =   quasi-log frequency vector (column vector)
%       df      =   minimum common divider of the values in freqv
%                   The harmonic numbers in fqlog can be calculated as
%                   harmno = round(fqlog/df);
%       cdmax   =   maximum common divider of the harmonic numbers
%       freqind =   column vector, indices of the selected
%                   frequency points in freqv
%
%       Input arguments:
%       freqv   =   strictly monotonic increasing frequency vector
%                   (usually with constant increments)
%       rf      =   desired frequency ratio
%
%       Usage: [fqlog, df, cdmax, freqind] = lin2qlog(freqv, rf);
%       Example: fqlog = lin2qlog([1:128], sqrt(2));
%
%       See also: LOG2QLOG, LOGSPACE.

%       Copyright (c) I. Kollar and Vrije Universiteit Brussel, ELEC, 1991-2003
%       All rights reserved.
%       $Revision: $
%       Last modified: 30-Dec-2003
%       Line 70 modified by Rik Pintelon: 25 March 2008

error(nargchk(1,2,nargin))
if min(size(freqv))>1, error('freqv is not a vector'), end
freqv=freqv(:); %column vector
if any(imag(freqv)), error('freqv is complex'), end
if any(freqv<0), error('freqv contains negative elements'), end
if any(diff(freqv)<=0), error('freqv is not strictly increasing'), end
if length(rf)~=1, error('rf is not scalar'), end
if rf<=0, error('rf is not positive'), end
%
freqvnz=freqv; if freqvnz(1)==0, freqvnz(1)=[]; end
dfv=diff([0;freqvnz]); df=min(dfv);
remfnegl=max(freqv)/df*eps*max(freqv);
while any(dfv>remfnegl)
  df0=df; df=min([df0;dfv]);
  dfv=sort(rem([df0;dfv],df));
  remfnegl=max(freqv)/df*eps*max(freqv);
  ind=find(dfv<=remfnegl); dfv(ind)=[];
end
%
if freqv(1)~=0
  lastind=1; freqind(1)=1; fqlog(1)=freqv(1);
else
  lastind=2; freqind(1)=2; fqlog(1)=freqv(2);
end
while lastind<length(freqv), %process until end of vector
  [m,nextind]=min(abs(freqv(lastind)*rf-freqv(lastind+1:length(freqv))));
  lastind=lastind+nextind;
%   if freqv(lastind)>fqlog(length(fqlog))*sqrt(rf),
  if freqv(lastind)>fqlog(length(fqlog))*rf^0.7,
    %step is acceptably large
    freqind=[freqind;lastind]; fqlog=[fqlog;freqv(lastind)];
  end
end %while
cdmax=1;
if nargout>1
  if max(freqv)/df>1e6
    disp(sprintf(['WARNING: maximum harmonic index found in ''lin2log''\n',...
            '   is %.0f, with df = %.3g Hz, T = %.3g s'],...
            round(max(freqv)/df),df,1/df))
    disp('Large indexes are often due to inaccurately given frequency values.')
    disp('If this is the case, before invoking msinclip do')
    disp('   freqv = round(freqv*T)/T;')
    disp('where T is the desired period length.')
    disp(' ')
  end
  harmno=round(fqlog/df); if harmno(1)==0, harmno(1)=[]; end
  for ind=2:min(harmno), if all(rem(harmno,ind)==0), cdmax=ind; end, end
end
%%%%%%%%%%%%%%%%%%%%%%%% end of lin2qlog %%%%%%%%%%%%%%%%%%%%%%%%
