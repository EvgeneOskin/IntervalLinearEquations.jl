module TestSolver
using Base.Test;
using ValidatedNumerics;
using ILESolver;


@test begin
    const bVector = [
        @interval(-10, 95);
        @interval(35, 14);
        @interval(-6, 2);
        @interval(30, 7);
        @interval(4, 95);
        @interval(-6, 46);
        @interval(-2, 65);
    ] :: Array{Interval{Float64}};
    const aMatrix = [
        @interval(4, 6) @interval(-9, 0) @interval(0, 12) @interval(2, 3) @interval(5, 9) @interval(-23, -9) @interval(15, 23);
        @interval(0, 1) @interval(6, 10) @interval(-1, 1) @interval(-1, 3) @interval(-5, 1) @interval(1, 15) @interval(-3, -1);
        @interval(0, 3) @interval(-20, -9) @interval(12, 77) @interval(-6, 30) @interval(0, 3) @interval(-18, 1) @interval(0, 1);
        @interval(-4, 1) @interval(-1, 1) @interval(-3, 1) @interval(3, 5) @interval(5, 9) @interval(1, 2) @interval(1, 4);
        @interval(0, 3) @interval(0, 6) @interval(0, 20) @interval(-1, 5) @interval(8, 14) @interval(-6, 1) @interval(10, 17);
        @interval(-7, -2) @interval(1, 2) @interval(7, 14) @interval(-3, 1) @interval(0, 2) @interval(3, 5) @interval(-2, 1);
        @interval(-1, 5) @interval(-3, 2) @interval(0, 8) @interval(1, 11) @interval(-5, 10) @interval(2, 7) @interval(6, 82);
    ] :: Array{Interval{Float64}};
    const preciseX = [
        @interval(-1.2247431, 0.5054298);
        @interval(18.2644433, -9.517504);
        @interval(-0.0281865, 1.6075521);
        @interval(16.4076957, -14.45534);
        @interval(-1.3435652, 3.9882184);
        @interval(-3.5289385, 4.5434583);
        @interval(5.43086238, -0.674008);
    ] :: Array{Interval{Float64}};

    solution = ILESolver.solve("G", aMatrix, bVector, 0.1, 1.0)
    solution == preciseX
end

@test begin
    const bVector = [
        @interval(-2, 2); @interval(-2, 2);
    ] :: Array{Interval{Float64}};
    const aMatrix = [
        @interval(2, 4) @interval(-2, 1);
        @interval(-1, 2) @interval(2, 4);
    ] :: Array{Interval{Float64}};
    const preciseX = [
        @interval(-0.333333, 0.333333);
        @interval(-0.333333, 0.333333);
    ] :: Array{Interval{Float64}};

    solution = ILESolver.solve("G", aMatrix, bVector, 0.1, 1.0)
    solution == preciseX
end

@test begin
    const bVector = [
        @interval(-2, 2); @interval(-2, 2);
    ] :: Array{Interval{Float64}};
    const aMatrix = [
        @interval(2, 4) @interval(-2, 1);
        @interval(-1, 2) @interval(2, 4);
    ] :: Array{Interval{Float64}};

    ILESolver.solve("G", eye(aMatrix), bVector, 0.1, 1.0) == bVector
end

end
