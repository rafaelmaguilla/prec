function varargout = randraw(distribName, distribParams, varargin)
%
%   EFFICIENT RANDOM VARIATES GENERATOR
% 
% See alphabetical list of the supported distributions below (over 50 distributions)
% 
% 1)  randraw 
%           presents general help.
% 2)  randraw( distribName ) 
%           presents help for the specific distribution defined 
%           by usage string distribName (see table below).
% 3)  Y = randraw( distribName, distribParams, sampleSize );
%           returns array Y of size = sampleSize of random variates from distribName  
%           distribution with parameters distribParams

%  Version 2.0 - August 2007
%         1) New distributions support: Nakagami and Rician !
%         2) Small typo corrections in comments
%  Version 1.8 - February 2007
%         GIG distribution (thanks to Mr. Demetris Lamnisos)
%           Computational exceptions in the reparameterized GIG generation were fixed
%  Version 1.7 - December 2006
%         GIG distribution (thanks to Dr. Junbin Gao)
%           Computational exceptions in the reparameterized GIG generation were fixed
%  Version 1.6 - September 2006
%         Exception handling: BINOMIAL distribution - special case for n*p~=0 
%         Geometric distibution: additional note in help section 
%  Version 1.5 - December 2005
%        'true' and 'false' functions were replased by ones and zeros to support Matlab releases 
%         below 6.5
%  Version 1.4 - September 2005 -
%      Bugs fix:
%        1) GAMMA distribution (thanks to Earl Lawrence):
%             special case for a<1
%        2) GIG distribution (thanks to Panagiotis Braimakis):
%            typo in help 
%            code adjustment to overcome possible computational overflows
%        3) CHI SQUARE distribution
%            typo in help
%  Version 1.3 - July 2005 -
%      Bug fix:
%         Typo in GIG distribution generation:
%         should be 'out' instead of 'x' in lines 1852 and 1858 
%  Version 1.2 - May 2005  -   
%      Bugs fix: 
%        1) Poisson distribution did not work for lambda < 21.4. Typo ( ti instead of t )
%        2) GIG distribution:  support to chi=0 or psi=0 cases
%        3) Beta distribution: column sampleSize 
%        4) Cauchy distribution: typo in example
%        5) Chi distribution:   typo in example
%        6) Non-central F distribution:  number of input parameters
%        7) INVERSE GAUSSIAN (IG) distribution: typo in example
%
%  Version 1.1 - April 2005 -  Bug fix:   Generation from binomial distribution using only 'binomial'
%                                   usage string was changed to 'binom' ( 'binomial' works too ).
%  Version 1.0 - March 2005 -  Initial version
%  Alex Bar Guy  &  Alexander Podgaetsky
%    alex.barguy@gmail.com

% Copyright (c) 2005, Alex Bar-Guy
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% * Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%                       SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%                       INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% Any comments and suggestions please send to:
%    alex.barguy@gmail.com

% Reference links:
%   1) http://mathworld.wolfram.com/topics/StatisticalDistributions.html
%   2) http://en.wikipedia.org/wiki/Category:Probability_distributions
%   3) http://www.brighton-webs.co.uk/index.asp
%   4) http://www.jstatsoft.org/v11/i03/v11i03.pdf
%   5) http://www.quantlet.com/mdstat/scripts/csa/html/node236.html

funcName = mfilename;

if nargin == 0
     help(funcName);
     return;
elseif nargin == 1
     runMode = 'distribHelp';
elseif nargin == 2
     runMode = 'genRun';
     sampleSize = [1 1];
else
     runMode = 'genRun';
     sampleSize = [varargin{1:end}];
end

distribNameInner = lower( distribName( ~isspace( distribName ) ) );

if strcmp(runMode, 'distribHelp')
     fid = fopen( [ funcName '.m' ], 'r' );
     printHelpFlag = 0;
     while 1
          tline = fgetl( fid );
          if ~ischar( tline )
               fprintf( '\n Unknown distribution name ''%s''.\n', distribName );
               break;
          end
          if ~isempty( strfind( tline, [ 'END ', distribNameInner,' HELP' ] ) )
               printHelpFlag = 0;
               break;
          end
          if printHelpFlag
               startPosition = strfind( tline, ' % ' ) + 3;
               printLine = tline( startPosition : end );
               if ~strcmp( funcName, 'randraw' )
                    indxs = strfind( printLine, 'randraw' );
                    while ~isempty( indxs )
                         headLine = printLine( 1:indxs(1)-1 );
                         tailLine = printLine( indxs(1)+7:end );
                         printLine = [ headLine, funcName, tailLine ];
                         indxs = strfind( printLine, 'randraw' );
                    end
               end
               pause(0.02);
               fprintf( '\n%s', printLine );
          end
          if ~isempty( strfind( tline, [ 'START ', distribNameInner,' HELP' ] ) )
               printHelpFlag = 1;
          end
     end
     fprintf( '\n\n' );
     fclose( fid );
     if nargout > 0
          varargout{1} = [];
     end
     return;
end

if length(sampleSize) == 1
     sampleSize = [ sampleSize, 1 ];
end

if strcmp(runMode, 'genRun')
     runExample = 0;
     plotFlag = 0;

     dbclear if warning;
     out = [];
     if prod(sampleSize) > 0
          switch lower( distribNameInner )
               case {'exp','exponential'}
                    % START exp HELP START exponential HELP
                    % THE EXPONENTIAL DISTRIBUTION
                    %
                    % pdf = lambda * exp( -lambda*y );
                    % cdf = 1 - exp(-lambda*y);
                    %
                    %  Mean = 1/lambda;
                    %  Variance = 1/lambda^2;
                    %  Mode = lambda;
                    %  Median = log(2)/lambda;
                    %  Skewness = 2;
                    %  Kurtosis = 6;
                    %
                    % PARAMETERS:
                    %   lambda - inverse scale or rate (lambda>0)
                    %
                    % SUPPORT:
                    %   y,  y>= 0
                    %
                    % CLASS:
                    %   Continuous skewed distributions
                    %
                    % NOTES:
                    %  The discrete version of the Exponential distribution is 
                    %  the Geometric distribution.
                    %
                    % USAGE:
                    %   randraw('exp', lambda, sampleSize) - generate sampleSize number
                    %         of variates from the Exponential distribution
                    %         with parameter 'lambda';
                    %   randraw('exp') - help for the Exponential distribution;
                    %
                    % EXAMPLES:
                    %  1.   y = randraw('exp', 1, [1 1e5]);
                    %  2.   y = randraw('exp', 1.5, 1, 1e5);
                    %  3.   y = randraw('exp', 2, 1e5 );
                    %  4.   y = randraw('exp', 3, [1e5 1] );
                    %  5.   randraw('exp');
                    %
                    % SEE ALSO:
                    %   GEOMETRIC, GAMMA, POISSON, WEIBULL distributions
                    % END exp HELP END exponential HELP
                    
                    checkParamsNum(funcName, 'Exponential', 'exp', distribParams, [1]);  
                    lambda  = distribParams(1);
                    validateParam(funcName, 'Exponential', 'exp', 'lambda', 'lambda', lambda, {'> 0'});
                    
                    out = -log( rand( sampleSize ) ) / lambda;
               otherwise
                    fprintf('\n RANDRAW: Unknown distribution name: %s \n', distribName);
                    
          end % switch lower( distribNameInner )
          
     end % if prod(sampleSize)>0

     varargout{1} = out;

     return;

end % if strcmp(runMode, 'genRun')

return;


function checkParamsNum(funcName, distribName, runDistribName, distribParams, correctNum)
if ~any( numel(distribParams) == correctNum )
     error('%s Variates Generation:\n %s%s%s%s%s', ...
          distribName, ...
          'Wrong numebr of parameters (run ',...
          funcName, ...
          '(''', ...
          runDistribName, ...
          ''') for help) ');
end
return;


function validateParam(funcName, distribName, runDistribName, distribParamsName, paramName, param, conditionStr)
condLogical = 1;
eqCondStr = [];
for nn = 1:length(conditionStr)
     if nn==1
          eqCondStr = [eqCondStr conditionStr{nn}];
     else
          eqCondStr = [eqCondStr ' and ' conditionStr{nn}];          
     end
     eqCond = conditionStr{nn}(1:2);
     eqCond = eqCond(~isspace(eqCond));
     switch eqCond
          case{'<'}
               condLogical = condLogical & (param<str2num(conditionStr{nn}(3:end)));
          case{'<='}
               condLogical = condLogical & (param<=str2num(conditionStr{nn}(3:end)));               
          case{'>'}
               condLogical = condLogical & (param>str2num(conditionStr{nn}(3:end))); 
          case{'>='}
               condLogical = condLogical & (param>=str2num(conditionStr{nn}(3:end)));
          case{'~='}
               condLogical = condLogical & (param~=str2num(conditionStr{nn}(3:end)));
          case{'=='}
               if strcmp(conditionStr{nn}(3:end),'integer')
                    condLogical = condLogical & (param==floor(param));                    
               else
                    condLogical = condLogical & (param==str2num(conditionStr{nn}(3:end)));
               end
     end
end

if ~condLogical
     error('%s Variates Generation: %s(''%s'',%s, SampleSize);\n Parameter %s should be %s\n (run %s(''%s'') for help)', ...
          distribName, ...
          funcName, ...
          runDistribName, ...
          distribParamsName, ...
          paramName, ...
          eqCondStr, ...
          funcName, ...
          runDistribName);
end
return;

function cdf = normcdf(y)
cdf = 0.5*(1+erf(y/sqrt(2)));
return;

function pdf = normpdf(y)
pdf = 1/sqrt(2*pi) * exp(-1/2*y.^2);
return;

function cdfinv = norminv(y)
cdfinv = sqrt(2) * erfinv(2*y - 1);
return;

function out = randFrom5Tbls( P, offset, sampleSize)
sizeP = length(P);

if sizeP == 0
     out = [];
     return;
end

a = mod(floor([0 P]/16777216), 64);
na = cumsum( a );
b = mod(floor([0 P]/262144), 64);
nb = cumsum( b );
c = mod(floor([0 P]/4096), 64);
nc = cumsum( c );
d = mod(floor([0 P]/64), 64);
nd = cumsum( d );
e =  mod([0 P], 64);
ne = cumsum( e );

AA = zeros(1, na(end));
BB = zeros(1, nb(end));
CC = zeros(1, nc(end));
DD = zeros(1, nd(end));
EE = zeros(1, ne(end));

t1 = na(end)*16777216;
t2 = t1 + nb(end)*262144;
t3 = t2 + nc(end)*4096;
t4 = t3 + nd(end)*64;

k = (1:sizeP)+offset-1;
for ii = 1:sizeP
     AA(na(ii)+(0:a(ii+1))+1) = k(ii);
     BB(nb(ii)+(0:b(ii+1))+1) = k(ii);
     CC(nc(ii)+(0:c(ii+1))+1) = k(ii);
     DD(nd(ii)+(0:d(ii+1))+1) = k(ii);
     EE(ne(ii)+(0:e(ii+1))+1) = k(ii);
end

jj = round(min(sum(P),1073741823) *rand(sampleSize));
out = zeros(sampleSize);
N = prod(sampleSize);
for ii = 1:N
     if jj(ii) < t1
          out(ii) = AA( floor(jj(ii)/16777216)+1 );
     elseif jj(ii) < t2
          out(ii) = BB(floor((jj(ii)-t1)/262144)+1);
     elseif jj(ii) < t3
          out(ii) = CC(floor((jj(ii)-t2)/4096)+1);
     elseif jj(ii) < t4
          out(ii) = DD(floor((jj(ii)-t3)/64)+1);
     else
          out(ii) = EE(floor(jj(ii)-t4) + 1);
     end
end

return;