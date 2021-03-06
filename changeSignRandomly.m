function y = changeSignRandomly(array)
%Essa funcao tem como objetivo alterar, de maneira aleatoria
%os sinais de um array gerado pela funcao generateRandomVariable,
%de modo a resultar na densidade de probabilidade descrita pela
%equacao (1)

%Parametro:
%	array: vetor coluna, produzido pela funcao generateRandomVariable
%Saida:
%	y: vetor coluna resultante da troca de sinal aleatoria dos elementos de array
randomSignArray = transpose(randi([0, 1], 1, length(array)))
randomSignArray(randomSignArray == 0) = -1

y = array.*randomSignArray

return;
