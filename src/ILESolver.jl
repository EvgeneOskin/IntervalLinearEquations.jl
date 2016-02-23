module ILESolver
using ValidatedNumerics;
using ILESolver.F
using ILESolver.G
using ILESolver.System
using ILESolver.sti

function solve{T}(name, A :: Array{Interval{T}}, B :: Array{Interval{T}}, precision, scale)
    @assert ndims(A) == 2
    @assert ndims(B) == 1
    @assert precision > 0
    @assert scale <= 1.0 && scale > 0
    @assert size(A, 1) == size(A, 2)
    @assert size(A, 1) == size(B, 1)
    @assert name == "F" || name == "G"

    if name == "F"
        case_module = ILESolver.F
    else
        case_module = ILESolver.G
    end

    const system = ILESolver.System.System(
        A, B, ILESolver.sti.STI(B), size(A, 1), size(A, 1)*2
    )
    const initial = case_module.initialConditions(system)
    solver = ILESolver.System.Solver(
        initial, initial, initial, system, [initial], 1
    )

    start = false
    while !start || solver.iternation < 5
        solver.before = solver.current
        calculatedSubDifferential = case_module.subDifferential(solver)
        equationValue = case_module.equation(solver)
        solver.current = iterate(solver.before, scale, calculatedSubDifferential, equationValue)
        solver.roots.push(solver.current)
        start = true
        solver.iternation += 1
    end
    solver.current.intervals
end

function iterate(previous, scale, subdiff, equation_value)
    new_sti = previous.sti - scale * inv(subdiff) * equation_value
    IntervalVector(ILESolver.sti.reverseSTI(new_sti), sti)
end


end
