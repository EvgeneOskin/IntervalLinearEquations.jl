module ILESolver
using ValidatedNumerics;

module F

using ValidatedNumerics;
using ILESolver;
using Debug

function solve{T}(A :: Array{Interval{T}}, B :: Array{Interval{T}}, precision, scale)
    @assert ndims(A) == 2
    @assert ndims(B) == 1
    @assert precision > 0
    @assert scale <= 1.0 && scale > 0
    @assert size(A, 1) == size(A, 2)
    @assert size(A, 1) == size(B, 1)

    xNcurrent = initialConditions(A, B)
    xNbefore = xNcurrent
    start = false
    iter = 1
    while !start || iter < 5
        xNbefore = xNcurrent
        calculatedSubDifferential = subDifferential(A, xNbefore)
        equationValue = equation(xNbefore, A, B)
        xNcurrent = xNbefore - scale * inv(calculatedSubDifferential) * equationValue
        print(iter, "\nbef=\n", xNbefore, "\ncur=\n", xNcurrent, "\n")
        print("subdiff=\n", inv(calculatedSubDifferential), "\n", equationValue)
        start = true
        iter += 1
    end
    ILESolver.reverseSTI(xNcurrent) :: Array{Interval{T}}
end

function equation{T}(xArg :: Array{T}, aMatrix :: Array{Interval{T}}, bVector :: Array{Interval{T}})
    intervalMultiplication = aMatrix*ILESolver.reverseSTI(xArg) :: Array{Interval{T}}
    ILESolver.STI(intervalMultiplication) - xArg + ILESolver.STI(bVector) :: Array{T}
end

function initialConditions(aMatrix, bVector)
    systemSize = size(bVector, 1)*2
    aMatrixMid = map(x -> mid(x), aMatrix)
    leftMatrix = eye(systemSize) - ILESolver.constituentMatrix(aMatrixMid)
    rightVector = ILESolver.STI(bVector)
    initialInNumbers = leftMatrix \ rightVector
    initialInNumbers
end

function subDifferential(aMatrix, xVector)
    systemSize = size(xVector, 1)
    systemSizeHalf = int(systemSize/2)
    identityMatrix = eye(systemSizeHalf)

    loSubDiff = map(
        index -> calulateLoSubdifferentialRow(slicedim(aMatrix, 1, index), xVector, slicedim(identityMatrix, 1, index)),
        range(1, systemSizeHalf)
    )
    hiSubDiff = map(
        index -> calulateHiSubdifferentialRow(slicedim(aMatrix, 1, index), xVector, slicedim(identityMatrix, 1, index)),
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
    transpose(stiedCalculatedSum - identityVector)
end

end

module G

using ValidatedNumerics;
using ILESolver;
using Debug

@debug function solve{T}(A :: Array{Interval{T}}, B :: Array{Interval{T}}, precision, scale)
    @assert ndims(A) == 2
    @assert ndims(B) == 1
    @assert precision > 0
    @assert scale <= 1.0 && scale > 0
    @assert size(A, 1) == size(A, 2)
    @assert size(A, 1) == size(B, 1)

    xNcurrent = initialConditions(A, B)
    xNbefore = xNcurrent
    start = false
    iter = 1
    while !start || iter < 5
        xNbefore = xNcurrent
        calculatedSubDifferential = subDifferential(A, xNbefore)
        equationValue = equation(xNbefore, A, B)
        @bp
        xNcurrent = xNbefore - scale * inv(calculatedSubDifferential) * equationValue
        print(iter, "\nbef=\n", xNbefore, "\ncur=\n", xNcurrent, "\n")
        print("subdiff=\n", inv(calculatedSubDifferential), "\n", equationValue)
        start = true
        iter += 1
    end
    ILESolver.reverseSTI(xNcurrent) :: Array{Interval{T}}
end

function subDifferential(aMatrix, xVector)
    systemSize = size(xVector, 1)
    systemSizeHalf = int(systemSize/2)
    identityMatrix = eye(systemSizeHalf)

    loSubDiff = map(
        index -> calulateLoSubdifferentialRow(slicedim(aMatrix, 1, index), xVector, slicedim(identityMatrix, 1, index)),
        range(1, systemSizeHalf)
    )
    hiSubDiff = map(
        index -> calulateHiSubdifferentialRow(slicedim(aMatrix, 1, index), xVector, slicedim(identityMatrix, 1, index)),
        range(1, systemSizeHalf)
    )
    vcat(loSubDiff..., hiSubDiff...)
end

function equation(xArg, aMatrix, bVector)
    ILESolver.STI(aMatrix*ILESolver.reverseSTI(xArg)) - ILESolver.STI(bVector)
end

function initialConditions(aMatrix, bVector)
    aMatrixMid = map(x -> mid(x), aMatrix)
    leftMatrix = ILESolver.constituentMatrix(aMatrixMid)
    rightVector = ILESolver.STI(bVector)
    leftMatrix \ rightVector
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


function constituentMatrix{T}(matrix :: Array{T, 2})
    @assert ndims(matrix) == 2
    @assert size(matrix, 1) == size(matrix, 2)

    positivMatrix = map(pos, matrix)
    negativMatrix = map(neg, matrix)
    [
        positivMatrix negativMatrix;
        negativMatrix positivMatrix
    ] :: Array{T, 2}
end

function STI{T}(intervalVector :: Array{Interval{T}})
    @assert ndims(intervalVector) == 1

    function stiReducer(reduced, element)
        reducedSize = size(reduced, 1)
        if (reducedSize == 0)
            return [-element.lo; element.hi]
        end
        halfReducedSize = int(reducedSize/2)
        [
            slice(reduced, 1:halfReducedSize)
            [-element.lo]
            slice(reduced, (1 + halfReducedSize):reducedSize)
            [element.hi]
        ] # Could use cat instead
    end

    reduce(stiReducer, [], intervalVector) :: Array{T}
end

pos = x -> x >= zero(x) ? x : zero(x)
neg = x -> x < zero(x) ? -x : zero(x)

function reverseSTI{T}(intervalVector :: Array{T})
    @assert ndims(intervalVector) == 1

    vectorSize = size(intervalVector, 1)
    halfVectorSize = int(vectorSize/2)
    negativLoIntervalPart = slice(intervalVector, 1:halfVectorSize)
    positivHiIntervalPart = slice(intervalVector, (1 + halfVectorSize):vectorSize)

    positivLoIntervalPart = map(-, negativLoIntervalPart)
    zippedIntervalParts = zip(positivLoIntervalPart, positivHiIntervalPart)
    collect(Interval{T}, map(x -> @interval(x[1], x[2]), zippedIntervalParts)) :: Array{Interval{T}}
end

function productLoSubDifferential{T}(aInterval :: Interval{T}, xLo :: T, xHi :: T)
    pos(aInterval.lo) * dXLoNeg(xLo) + neg(aInterval.hi) * dXHiNeg(xHi) - dMaxHiLo(aInterval, pos(aInterval.hi)*pos(xLo), neg(aInterval.lo)*pos(xHi))
end

function productHiSubDifferential{T}(aInterval :: Interval{T}, xLo :: T, xHi :: T)
    - pos(aInterval.lo) * dXHiNeg(xHi) - neg(aInterval.hi) * dXLoNeg(xLo) + dMaxHiHi(aInterval, pos(aInterval.hi)*pos(xHi), neg(aInterval.lo)*pos(xLo))
end

function dXLoPos(x)
    if x < 0
        [@interval(0); @interval(0)]
    elseif x == 0
        [@interval(0, 1); @interval(0)]
    else
        [@interval(1); @interval(0)]
    end
end


function dXLoNeg(x)
    if x < 0
        [@interval(-1); @interval(0)]
    elseif x == 0
        [@interval(-1, 0); @interval(0)]
    else
        [@interval(0); @interval(0)]
    end
end

function dXHiPos(x)
    if x < 0
        [@interval(0); @interval(0)]
    elseif x == 0
        [@interval(0); @interval(0, 1)]
    else
        [@interval(0); @interval(1)]
    end
end

function dXHiNeg(x)
    if x < 0
        [@interval(0); @interval(-1)]
    elseif x == 0
        [@interval(0); @interval(-1, 0)]
    else
        [@interval(0); @interval(0)]
    end
end

function dMaxHiLo(c, left, right)
    if left > right
        [pos(c.hi), 0]
    elseif left == right
        0.5*[pos(c.hi), neg(c.lo)]
    else
        [0, neg(c.lo)]
    end
end

function dMaxHiHi(c, left, right)
    if left > right
        [0, pos(c.hi)]
    elseif left == right
        0.5*[neg(c.lo), pos(c.hi)]
    else
        [neg(c.lo), 0]
    end
end

end
