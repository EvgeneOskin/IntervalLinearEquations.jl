using ILESolver;
using ValidatedNumerics;
using Base.Test;

bVector = [
    @interval(1, 2.1); @interval(0, 1); @interval(2, 3);
];
aMatrix = [
    @interval(1, 2) @interval(0, 0) @interval(1, 2);
    @interval(0, 0) @interval(1, 2) @interval(0, 0);
    @interval(0, 0) @interval(0, 0) @interval(1, 2);
];
m = [
    1 0 1;
    -0 -5.1 1;
    -1 -2 10;
];

# @test m*bVector == ILESolver.reverseSTI(ILESolver.constituentMatrix(m)*ILESolver.STI(bVector))

# @test bVector == ILESolver.reverseSTI(ILESolver.STI(bVector))

print(ILESolver.F.solve(aMatrix, bVector, 0.1, 1), "\n")
