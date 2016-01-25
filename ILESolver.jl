module ILESolver

using ValidatedNumerics;

module F

using ValidatedNumerics;
using ILESolver;

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
    while !start && iter < 5
        xNbefore = xNcurrent
        calculatedSubDifferential = subDifferential(A, xNbefore)
        print(calculatedSubDifferential, "\n")
        newWithouScaling = equation(xNbefore, A, B) / calculatedSubDifferential
        xNcurrent = xNbefore - scale * newWithouScaling
        start = true
        iter += 1
    end
    xNcurrent :: Array{Interval{T}}
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
    identityMatrix = eye(systemSize)

    loSubDiff = map(
        index -> calulateLoSubdifferentialRow(slicedim(aMatrix, 1, index), xVector, slicedim(identityMatrix, 1, index)),
        range(1, systemSizeHalf)
    )
    hiSubDiff = map(
        index -> calulateHiSubdifferentialRow(slicedim(aMatrix, 1, index), xVector, slicedim(identityMatrix, 1, index)),
        range(1, systemSizeHalf)
    )
    [loSubDiff hiSubDiff]
end

function calulateHiSubdifferentialRow{T}(aVector :: Array{Interval{T}}, xVector :: Array{T}, identityVector :: Array{T})
    systemSize = size(xVector, 1)
    systemSizeHalf = int(systemSize/2)
    calculatedSum = sum(map(
        index -> productHiSubDifferential(aVector[index], xVector[index], xVector[index+systemSizeHalf]),
        range(1, systemSizeHalf)
    ))
    calculatedSum - identityVector
end

function calulateLoSubdifferentialRow{T}(aVector :: Array{Interval{T}}, xVector :: Array{T}, identityVector :: Array{T})
    systemSize = size(xVector, 1)
    systemSizeHalf = int(systemSize/2)
    calculatedSum = - sum(map(
        index -> productLoSubDifferential(aVector[index], xVector[index], xVector[index + systemSizeHalf]),
        range(1, systemSizeHalf)
    ))
    print(calculatedSum, "\n", identityVector, "\n")
    calculatedSum - identityVector
end

function productHiSubDifferential{T}(aInterval :: Interval{T}, xLo :: T, xHi :: T)
    aInterval.lo * dXLoPos(xLo) + aInterval.hi * dXHiNeg(xHi) - dMaxLoHi(aInterval, aInterval.lo*xLo, aInterval.hi*xHi)
end

function productLoSubDifferential{T}(aInterval :: Interval{T}, xLoNegative :: T, xHi :: T)
    - aInterval.lo * dXLoNeg(xHi) - aInterval.hi * dXHiPos(xLoNegative) + dMaxHiLo(aInterval, aInterval.lo*xHi, aInterval.hi*xLoNegative)
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

function dMaxLoHi(c, left, right)
    if left > right
        [c.hi, 0]
    elseif left == right
        0.5*[c.hi, c.lo]
    else
        [0, c.lo]
    end
end

function dMaxHiLo(c, left, right)
    if left > right
        [0, c.hi]
    elseif left == right
        0.5*[c.hi, c.lo]
    else
        [c.lo, 0]
    end
end
end

module G
using ValidatedNumerics;
using ILESolver;

function solve(A, B)
    @assert ndims(A) == 2
    @assert ndims(B) == 1
    @assert size(A, 1) == size(A, 2)
    @assert size(A, 1) == size(B, 1)

    reverseSTI(STI(B))
end

function equation(xArg, aMatrix, bVector)
    STI(aMatrix*reverseSTI(xArg)) - STI(bVector)
end

function initialConditions{T}(aMatrix :: Array{Interval{T}}, bVector :: Array{Interval{T}})
    aMatrixMid = map(x -> mid(x), aMatrix)
    leftMatrix = ILESolver.constituentMatrix(aMatrixMid)
    rightVector = STI(bVector)
    leftMatrix \ transpose(rightVector) :: Array{Interval{T}}
end

end


function constituentMatrix{T}(matrix :: Array{T, 2})
    @assert ndims(matrix) == 2
    @assert size(matrix, 1) == size(matrix, 2)

    positivMatrix = map(x -> x >= zero(x) ? x : zero(x), matrix)
    negativMatrix = map(x -> x < zero(x) ? -x : zero(x), matrix)
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

end
