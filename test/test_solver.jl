module SolverTests
using Base.Test;
using FactCheck;
using ValidatedNumerics;
using IntervalLinearEquations;

function interval_approx_eq(a_interval :: Interval{Float64}, b_interval :: Interval{Float64})
    isapprox(a_interval.lo, b_interval.lo) &&
    isapprox(a_interval.hi, b_interval.hi)
end


facts("Solve 2x2 system with number matrix") do
    bVector = [
        @interval(-2.0, 2.0); @interval(-2.0, 2.0);
    ] :: Array{Interval{Float64}};
    aMatrix = [
        @interval(2, 4) @interval(-2, 1);
        @interval(-1, 2) @interval(2, 4);
    ] :: Array{Interval{Float64}};

    result = IntervalLinearEquations.solve("G", eye(aMatrix), bVector, 0.1, 10, 1.0)
    @fact result --> bVector
end

facts("Solve 2x2 system with random matrix") do
    xVector = [
        @interval(-1.0, 2.0); @interval(-1.0, 2.0);
    ] :: Array{Interval{Float64}};
    aMatrix = [
        @interval(-2, 4) @interval(-3, 0);
        @interval(-1, 0) @interval(2, 5);
    ] :: Array{Interval{Float64}};
    bVector = aMatrix * xVector

    solution = IntervalLinearEquations.solve("G", aMatrix, bVector, 0.1, 10, 1.0)
    bFromSolution = aMatrix * solution
    @fact all(map(interval_approx_eq, bVector, bFromSolution)) --> true
    @fact all(map(interval_approx_eq, solution, xVector)) --> true
end

facts("Solve 2x2 system with interval matrix") do
    bVector = [
        @interval(-2.0, 2.0); @interval(-2.0, 2.0);
    ] :: Array{Interval{Float64}};
    aMatrix = [
        @interval(2, 4) @interval(-2, 1);
        @interval(-1, 2) @interval(2, 4);
    ] :: Array{Interval{Float64}};
    preciseX = [
        @interval(-0.33333333333333337, 0.33333333333333337);
        @interval(-0.33333333333333337, 0.33333333333333337);
    ] :: Array{Interval{Float64}};

    solution = IntervalLinearEquations.solve("G", aMatrix, bVector, 0.1, 10, 1.0)
    bFromSolution = aMatrix * solution
    @fact all(map(interval_approx_eq, bVector, bFromSolution)) --> true
    @fact all(map(interval_approx_eq, solution, preciseX)) --> true
end

end
