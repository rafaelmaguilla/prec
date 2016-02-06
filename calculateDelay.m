%Item (d) - Para determinar o atraso entre os sinais vindos de A e B,
%utilizaremos a funcao circshift do MatLab

%recebe como parametros as sequencias a serem comparadas, o numero N de bits da cadeia e a duracao delta_t de cada bit
function delayAB = calculate_delay(seqA, seqB, N, delta_t) 
if seqA == seqB
	delayAB = 0 %se as sequencias sao identicas, então nao ha atraso relativo entre elas
else %atrasa a sequencia A um bit de cada vez, ate que ela se iguale a sequencia B
	for t = 1:N
		if circshift(seqA, t) == seqB 
			delayAB = t %delayAB recebe o valor do atraso em um numero inteiro t de duracao de bits
			break
		end
	end
	if delayAB >= length(seqA)/2
		delayAB = (delayAB - length(seqA))*delta_t %nesse caso, o atraso eh negativo (isto eh, a sequencia seqB esta adiantada)
	else
		delayAB = delayAB*delta_t %nesse caso, o atraso eh positivo (isto eh, a sequencia seqA esta adiantada)
	end
end

