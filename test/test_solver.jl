module SolverTests
using Base.Test;
using FactCheck;
using ValidatedNumerics;
using ILESolver;

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

    result = ILESolver.solve("G", eye(aMatrix), bVector, 0.1, 10, 1.0)
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

    solution = ILESolver.solve("G", aMatrix, bVector, 0.1, 10, 1.0)
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

    solution = ILESolver.solve("G", aMatrix, bVector, 0.1, 10, 1.0)
    bFromSolution = aMatrix * solution
    @fact all(map(interval_approx_eq, bVector, bFromSolution)) --> true
    @fact all(map(interval_approx_eq, solution, preciseX)) --> true
end

facts("Solve 7x7 system with interval matrix") do
    bVector = [
        @interval(-347.5265180000002, 399.17321294000016);
        @interval(-230.54111184000013, 310.8722265000002);
        @interval(-891.2806088000004, 888.7951626800003);
        @interval(-119.23151460000004, 173.51377802000005);
        @interval(-191.1484819000001, 394.6258187600002);
        @interval(-103.38426236000002, 147.09842838000006);
        @interval(-320.0635949000001, 749.4182123600002);
    ] :: Array{Interval{Float64}};
    aMatrix = [
        @interval(4, 6) @interval(-9, 0) @interval(0, 12) @interval(2, 3) @interval(5, 9) @interval(-23, -9) @interval(15, 23);
        @interval(0, 1) @interval(6, 10) @interval(-1, 1) @interval(-1, 3) @interval(-5, 1) @interval(1, 15) @interval(-3, -1);
        @interval(0, 3) @interval(-20, -9) @interval(12, 77) @interval(-6, 30) @interval(0, 3) @interval(-18, 1) @interval(0, 1);
        @interval(-4, 1) @interval(-1, 1) @interval(-3, 1) @interval(3, 5) @interval(5, 9) @interval(1, 2) @interval(1, 4);
        @interval(0, 3) @interval(0, 6) @interval(0, 20) @interval(-1, 5) @interval(8, 14) @interval(-6, 1) @interval(10, 17);
        @interval(-7, -2) @interval(1, 2) @interval(7, 14) @interval(-3, 1) @interval(0, 2) @interval(3, 5) @interval(-2, 1);
        @interval(-1, 5) @interval(-3, 2) @interval(0, 8) @interval(1, 11) @interval(-5, 10) @interval(2, 7) @interval(6, 82);
    ] :: Array{Interval{Float64}};
    preciseX = [
        @interval(-1.2247431, 0.5054298);
        @interval(-9.517504, 18.2644433);
        @interval(-0.0281865, 1.6075521);
        @interval(-14.45534, 16.4076957);
        @interval(-1.3435652, 3.9882184);
        @interval(-3.5289385, 4.5434583);
        @interval(-0.674008, 5.43086238);
    ] :: Array{Interval{Float64}};

    solution = ILESolver.solve("G", aMatrix, bVector, 0.1, 10, 1.0)
    bFromSolution = aMatrix * solution
    @fact all(map(interval_approx_eq, bVector, bFromSolution)) --> true
    @fact all(map(interval_approx_eq, solution, preciseX)) --> true
end
end
