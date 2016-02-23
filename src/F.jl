module ILESolver
module F

using ValidatedNumerics;
using Debug
using ILESolver.System
using ILESolver.sti

function equation{T}(solver)
    intervalMultiplication = aMatrix*solver.before.intervals :: Array{Interval{T}}
    ILESolver.STI(intervalMultiplication) - solver.before.sti + solver.system.b.sti :: Array{T}
end

function initialConditions(system)
    aMatrixMid = map(x -> mid(x), system.a)
    leftMatrix = eye(system.sti_size) - ILESolver.sti.constituentMatrix(aMatrixMid)
    rightVector = system.b.sti
    initialInNumbers = leftMatrix \ rightVector
    ILESolver.System.IntervalVector(ILESolve.sti.reverseSTI(initialInNumbers), initialInNumbers)
end

function subDifferential(solver)
    identityMatrix = eye(solver.system.size)

    loSubDiff = map(
        index -> calulateLoSubdifferentialRow(slicedim(solver.system.a, 1, index), solver.before.sti, slicedim(identityMatrix, 1, index)),
        range(1, solver.system.size)
    )
    hiSubDiff = map(
        index -> calulateHiSubdifferentialRow(slicedim(solver.system.a, 1, index), solver.before.sti, slicedim(identityMatrix, 1, index)),
        range(1, solver.system.size)
    )
    vcat(loSubDiff..., hiSubDiff...)
end

function calulateLoSubdifferentialRow{T}(aVector :: Array{Interval{T}}, xVector :: Array{T}, identityVector :: Array{T})
    calulateSubdifferentialRowPart(ILESolver.productLoSubDifferential, -1, aVector, xVector, identityVector)
end

function calulateHiSubdifferentialRow{T}(aVector :: Array{Interval{T}}, xVector :: Array{T}, identityVector :: Array{T})
    calulateSubdifferentialRowPart(ILESolver.productHiSubDifferential, 1, aVector, xVector, identityVector)
end

function calulateSubdifferentialRowPart{T}(partFunction, sumFactor, aVector :: Array{Interval{T}}, xVector :: Array{T}, identityVector :: Array{T})
    systemSize = size(xVector, 1)
    systemSizeHalf = int(systemSize/2)
    identityVector = transpose(identityVector)
    calculatedSum = sumFactor * sum(map(
        index -> partFunction(aVector[index], xVector[index], xVector[index + systemSizeHalf]),
        range(1, systemSizeHalf)
    ))
    stiedCalculatedSum = ILESolver.STI(calculatedSum)
    transpose(stiedCalculatedSum - identityVector)
end

end
end
