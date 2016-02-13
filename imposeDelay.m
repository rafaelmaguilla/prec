function newArray = imposeDelay(someArray, k)
	newArray = circshift(someArray, k, 2)
end
