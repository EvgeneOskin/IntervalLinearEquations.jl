module ILESolver
using ValidatedNumerics;

include("Types.jl")
include("sti.jl")
include("F.jl")
include("G.jl")

function solve{T}(name :: AbstractString, A :: Array{Interval{T}, 2}, B :: Array{Interval{T}, 1}, precision, scale)
    @assert precision > 0.0
    @assert scale <= 1.0 && scale > 0

    case_module = get_module(name)
    solver = initialize(case_module, A, B)

    is_initial = true
    while is_initial || solver.iternation < 10
        solver.previous = solver.current
        calculatedSubDifferential = case_module.subDifferential(solver)
        equationValue = case_module.equation(solver)
        if det(calculatedSubDifferential) == 0
            break
        end
        solver.current = iterate(solver.previous, scale, calculatedSubDifferential, equationValue)
        solver.roots = vcat(solver.roots, [solver.current])
        is_initial = false
        solver.iternation += 1
    end
    solver.current.intervals
end

function initialize{T}(case_module, A :: Array{Interval{T}, 2}, B :: Array{Interval{T}, 1})
    @assert size(A, 1) == size(A, 2)
    @assert size(A, 1) == size(B, 1)


    system = Types.Configuration{T}(
        A, Types.IntervalVector(B, ILESolver.sti.STI(B)),
        size(A, 1), size(A, 1)*2
    )
    initialInNumbers = case_module.initialConditions(system)
    initial = Types.IntervalVector(sti.reverseSTI(initialInNumbers), initialInNumbers)
    solver = Types.Solver{T}(
        initial,
        initial,
        initial,
        system,
        [initial],
        1
    )
    solver
end

function iterate(previous, scale, subdiff, equation_value)
    new_sti = previous.sti - scale * inv(subdiff) * equation_value
    Types.IntervalVector(ILESolver.sti.reverseSTI(new_sti), new_sti)
end

function get_module(name :: AbstractString)
    @assert name == "F" || name == "G"

    if name == "F"
        F
    else
        G
    end
end

end
