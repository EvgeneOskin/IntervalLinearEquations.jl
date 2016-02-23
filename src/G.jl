module ILESolver

module G

using ValidatedNumerics;

function equation(solver)
    ILESolver.STI(solver.system.a*solver.before.intervals) - solver.system.b.sti
end

function initialConditions(system)
    aMatrixMid = map(x -> mid(x), system.a)
    leftMatrix = ILESolver.constituentMatrix(aMatrixMid)
    rightVector = system.b.sti
    initialInNumbers = leftMatrix \ rightVector
    ILESolver.System.IntervalVector(ILESolve.sti.reverseSTI(initialInNumbers), initialInNumbers)
end

function subDifferential(solver)
    identityMatrix = eye(solver.system.size)

    loSubDiff = map(
        index -> calulateLoSubdifferentialRow(slicedim(solver.system.a, 1, index), solver.before.sti, slicedim(identityMatrix, 1, index)),
        range(1, systemSizeHalf)
    )
    hiSubDiff = map(
        index -> calulateHiSubdifferentialRow(slicedim(solver.system.a, 1, index), solver.before.sti, slicedim(identityMatrix, 1, index)),
        range(1, systemSizeHalf)
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
    transpose(stiedCalculatedSum)
end
end


end
