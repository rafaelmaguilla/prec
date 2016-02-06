function varargout = randomVariableGenerator(lambda, sampleSize, varargin)
%
%   GERADOR DE VARIAVEIS ALEATORIAS COM DISTRIBUICAO EXPONENCIAL
%   Y = randraw(lambda, sampleSize);
%           retorna um array Y de comprimento = sampleSize de variaveis aleatorias 
%           com distribuicao exponencial com parametro lambda

%   Referencias:
%   	1) http://mathworld.wolfram.com/topics/StatisticalDistributions.html
%   	2) http://en.wikipedia.org/wiki/Category:Probability_distributions
%   	3) http://www.brighton-webs.co.uk/index.asp
%   	4) http://www.jstatsoft.org/v11/i03/v11i03.pdf
%   	5) http://www.quantlet.com/mdstat/scripts/csa/html/node236.html

funcName = mfilename;

if length(sampleSize) == 1
     sampleSize = [sampleSize, 1];
end

runExample = 0;
plotFlag = 0;

out = [];
if prod(sampleSize) > 0
    % DISTRIBUICAO EXPONENCIAL
    %
    %  fdp = lambda * exp( -lambda*y );
    %  
    %  Media = 1/lambda;
    %  Variancia = 1/lambda^2;
    %  Moda = lambda;
    %  Mediana = log(2)/lambda;
    %
    % PARAMETROS:
    %   lambda
    %
    % UTILIZACAO:
    %   randraw('exp', lambda, sampleSize) - gera um numero sampleSize
    %         de variaveis com distribuicao exponencial
    %         com parametro 'lambda';
    %
    % EXAMPLOS:
    %  1.   y = randraw('exp', 1, [1 1e5]);
    %  2.   y = randraw('exp', 1.5, 1, 1e5);
    %  3.   y = randraw('exp', 2, 1e5 );
    %  4.   y = randraw('exp', 3, [1e5 1] );
    out = -log(rand(sampleSize))/lambda;
  
end

varargout{1} = out;

return;

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
