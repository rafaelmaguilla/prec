%Item (d) - Para determinar o atraso entre os sinais vindos de A e B,
%utilizaremos a funcao circshift do MatLab

%Parametros:
%	seqA, seqB: sequencias a serem comparadas
%Saida
%	delayAB: atraso (em numero de bits) entre seqA e seqB
function delayAB = calculate_delay(seqA, seqB) 
if seqA == seqB
	delayAB = 0 %se as sequencias sao identicas, então nao ha atraso relativo entre elas
else %atrasa a sequencia A um bit de cada vez, ate que ela se iguale a sequencia B
	for t = 1:length(seqA)
		if circshift(seqA, t, 2) == seqB 
			delayAB = t %delayAB recebe o valor do atraso em um numero inteiro t de duracao de bits
			break
		end
	end
	if delayAB >= length(seqA)/2
		delayAB = delayAB - length(seqA) %nesse caso, o atraso eh negativo (isto eh, a sequencia seqB esta adiantada)
	end
end

