module ILESolver
using ValidatedNumerics;

include("Types.jl")
include("sti.jl")
include("F.jl")
include("G.jl")
using Debug

@debug function solve{T}(name, A :: Array{Interval{T}}, B :: Array{Interval{T}}, precision, scale)
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

    system = Types.Configuration{T}(
        A, Types.IntervalVector(B, ILESolver.sti.STI(B)),
        size(A, 1), size(A, 1)*2
    )
    initial = case_module.initialConditions(system)
    solver = Types.Solver{T}(
        initial,
        initial,
        initial,
        system,
        [initial],
        1
    )

    start = false
    while !start || solver.iternation < 5
        solver.previous = solver.current
        calculatedSubDifferential = case_module.subDifferential(solver)
        equationValue = case_module.equation(solver)
        @bp
        solver.current = iterate(solver.previous, scale, calculatedSubDifferential, equationValue)
        solver.roots.push(solver.current)
        start = true
        solver.iternation += 1
    end
    solver.current.intervals
end

function iterate(previous, scale, subdiff, equation_value)
    new_sti = previous.sti - scale * inv(subdiff) * equation_value
    Types.IntervalVector(ILESolver.sti.reverseSTI(new_sti), sti)
end

end
