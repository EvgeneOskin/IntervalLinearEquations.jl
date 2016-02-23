module TestSolver

using Compat
const bVector = [
    @interval(-2, 2); @interval(-2, 2);
] :: Array{Interval{Float64}};
const aMatrix = [
    @interval(2, 4) @interval(-2, 1);
    @interval(-1, 2) @interval(2, 4);
] :: Array{Interval{Float64}};
const preciseX = [
    @interval(-0.333333); @interval(-0.333333);
] :: Array{Interval{Float64}};

@test begin
    @assert ILESolver.G.solve(eye(aMatrix), bVector, 0.1, 1.0) == bVector
end

@test begin
    @assert ILESolver.G.solve(aMatrix, bVector, 0.1, 1.0) == preciseX
end

end
