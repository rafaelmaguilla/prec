function y = changeSignRandomly(array)

randomSignArray = transpose(randi([0, 1], 1, length(array)))
randomSignArray(randomSignArray == 0) = -1

y = array.*randomSignArray

return;
