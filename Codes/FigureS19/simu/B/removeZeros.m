function [vecA, vecB, target] = removeZeros(vecA, vecB, target)

zero_inds = (target == 0);
target(zero_inds) = [];
vecA(zero_inds) = [];
vecB(zero_inds) = [];

end