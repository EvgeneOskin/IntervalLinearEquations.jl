using ILESolver;
using ValidatedNumerics;
using Base.Test;

bVector = [
    @interval(-2, 2); @interval(-2, 2);
];
aMatrix = [
    @interval(2, 4) @interval(-2, 1);
    @interval(-1, 2) @interval(2, 4);
];
preciseX = [
    @interval(-0.333333); @interval(-0.333333);
]
m = [
    1 -2;
    -0 -5.1;
];

@test m*bVector == ILESolver.reverseSTI(ILESolver.constituentMatrix(m)*ILESolver.STI(bVector))

@test bVector == ILESolver.reverseSTI(ILESolver.STI(bVector))

print(ILESolver.F.solve(aMatrix, bVector, 0.1, 1), "\n")
