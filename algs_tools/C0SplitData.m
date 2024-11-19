function [Nvec, indjCell] = C0SplitData(D, c0)

Vs = setsLambda(D, c0);       
[Nvec, indjCell] = setsReLambda(D, Vs, c0);