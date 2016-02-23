module G

using ValidatedNumerics;
include("Types.jl")
include("SubDiff.jl")
include("sti.jl")

function equation(solver)
    sti.STI(solver.system.a*solver.previous.intervals) - solver.system.b.sti
end

function initialConditions(system)
    aMatrixMid = map(x -> mid(x), system.a)
    leftMatrix = sti.constituentMatrix(aMatrixMid)
    rightVector = system.b.sti
    initialInNumbers = leftMatrix \ rightVector
    Types.IntervalVector(sti.reverseSTI(initialInNumbers), initialInNumbers)
end

function subDifferential(solver)
    identityMatrix = eye(solver.system.size)

    loSubDiff = map(
        index -> calulateLoSubdifferentialRow(slicedim(solver.system.a, 1, index), solver.previous.sti, slicedim(identityMatrix, 1, index)),
        range(1, solver.system.size)
    )
    hiSubDiff = map(
        index -> calulateHiSubdifferentialRow(slicedim(solver.system.a, 1, index), solver.previous.sti, slicedim(identityMatrix, 1, index)),
        range(1, solver.system.size)
    )
    vcat(loSubDiff..., hiSubDiff...)
end

function calulateLoSubdifferentialRow{T}(aVector :: Array{Interval{T}}, xVector :: Array{T}, identityVector :: Array{T})
    calulateSubdifferentialRowPart(SubDiff.productLoSubDifferential, -1, aVector, xVector, identityVector)
end

function calulateHiSubdifferentialRow{T}(aVector :: Array{Interval{T}}, xVector :: Array{T}, identityVector :: Array{T})
    calulateSubdifferentialRowPart(SubDiff.productHiSubDifferential, 1, aVector, xVector, identityVector)
end

function calulateSubdifferentialRowPart{T}(partFunction, sumFactor, aVector :: Array{Interval{T}}, xVector :: Array{T}, identityVector :: Array{T})
    systemSize = size(xVector, 1)
    systemSizeHalf = int(systemSize/2)
    identityVector = transpose(identityVector)
    calculatedSum = sumFactor * sum(map(
        index -> partFunction(aVector[index], xVector[index], xVector[index + systemSizeHalf]),
        range(1, systemSizeHalf)
    ))
    stiedCalculatedSum = sti.STI(calculatedSum)
    transpose(stiedCalculatedSum)
end
end
