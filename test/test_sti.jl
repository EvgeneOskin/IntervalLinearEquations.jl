module

using Compat

const bVector = [
    @interval(-2, 2); @interval(-2, 2);
] :: Array{Interval{Float64}};
const m = [
    1 -2;
    -0 -5.1;
];

@test begin
    @assert m*bVector == ILESolver.reverseSTI(ILESolver.constituentMatrix(m)*ILESolver.STI(bVector))
end

@test begin
    @assert bVector == ILESolver.reverseSTI(ILESolver.STI(bVector))
end

end

