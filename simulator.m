%Variaveis globais, obtidas como resultado de resolucao da P3
%considerando uma taxa de transmissao de 56Mbits/s
deta_t = 17.857 %duracao do bit, em ns
T = 66660.181 %duracao da cadeia, em ns
N = 3733 %numero de bits por cadeia
%Item (c) - Gerador de sequencias aleatorias de bits

%Para garantir que o as sequencias serao iguais, basta impor
%a mesma semente para o gerador pseudo-aleatorio do MatLab, 
%por meio do comando rng (utilizarei meu numero USP como semente)

rng(8041626)
rand_seqA = randi([0, 1], 1, N)
rng(8041626)
rand_seqB = randi([0, 1], 1, N)



calculate_dalay(rand_seqA, rand_seqB, N, delta_t)
