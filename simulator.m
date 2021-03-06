function s = simulator(x, alpha, n)
	%Duracao de cada bit, com base numa taxa de transmissao de 56Mbits/s
	deltaT = 17.857e-9
	% Para garantir sequencias aleatorias identicas,
	% devemos ter sequencias produzidas a partir da 
	% mesma semente (vamos "setar" a semente para
	% o meu numero USP, 8041626, por meio do comando rng)
	rng(8041626) 
	A = randi([0, 1], 1, 3733)
	rng(8041626)
	B = randi([0, 1], 1, 3733)
	k = round(2*(5000-x)/(1.7857*3))
	B = imposeDelay(B, k)
	positionArray = []
	delayAB = calculateDelay(A, B)
for index = 1:n
	newDelayAB = changeSignRandomly(generateRandomVariable(alpha/deltaT, 1))
	positionArray(index) = calculatePosition(delayAB + newDelayAB)
end	
s = mean(positionArray)
return;
