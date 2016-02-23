module TestSti
using Base.Test;
using ValidatedNumerics;
using ILESolver, ILESolver.sti;

const bVector = [
    @interval(-2, 2); @interval(-2, 2);
] :: Array{Interval{Float64}};
const m = [
    1 -2;
    -0 -5.1;
];

@test begin
    interval_product = m*bVector
    sti_product = ILESolver.sti.reverseSTI(ILESolver.sti.constituentMatrix(m)*ILESolver.sti.STI(bVector))
    interval_product == sti_product
end

@test begin
    transformed = ILESolver.sti.reverseSTI(ILESolver.sti.STI(bVector))
    bVector == transformed
end

end

