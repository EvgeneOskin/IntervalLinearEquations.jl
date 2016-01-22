module ILESolver

using ValidatedNumerics;

module F

using ValidatedNumerics;
using ILESolver;

function solve(A, B, precision, scale)
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
        newWithouScaling = equation(xNbefore, A, B) / subDifferential(A, xNbefore)
        xNcurrent = xNbefore - scale * newWithouScaling
        start = true
        iter += 1
    end
    xNcurrent
end

function equation(xArg, aMatrix, bVector)
    print(aMatrix, "\n", ILESolver.reverseSTI(xArg), "\n")
    intervalMultiplication = aMatrix*ILESolver.reverseSTI(xArg)
    print(intervalMultiplication, "\n")
    ILESolver.STI(intervalMultiplication) - xArg + ILESolver.STI(bVector)
end

function initialConditions(aMatrix, bVector)
    systemSize = size(bVector, 1)*2
    aMatrixMid = map(x -> mid(x), aMatrix)
    leftMatrix = eye(systemSize) - ILESolver.constituentMatrix(aMatrixMid)
    rightVector = ILESolver.STI(bVector)
    leftMatrix \ rightVector
end

function subDifferential(aMatrix, xVector)
    systemSize = size(xVector, 1)
    systemSizeHalf = systemSize/2.
    identityMatrix = eye(aMatrix)

    loSubDiff = map(
        index -> calulateLoSubdifferentialRow(aMatrix[index], xVector, identityMatrix[index]),
        1..systemSizeHalf
    )
    hiSubDiff = map(
        index -> calulateHiSubdifferentialRow(aMatrix[index], xVector, identityMatrix[index]),
        (systemSizeHalf + 1)..systemSize
    )
end

function calulateHiSubdifferentialRow(aVector, xVector, identityVector)
    systemSize = size(xVector, 1)
    systemSizeHalf = systemSize/2.
    sum(map(
        index -> productHiSubDifferential(aVector[index], xVector[index], xVector[index+systemSizeHalf]),
        1..systemSizeHalf
    )) - identityVector
end

function calulateLoSubdifferentialRow(aVector, xVector, identityVector)
    systemSize = size(xVector, 1)
    systemSizeHalf = systemSize/2.
    -sum(map(
        index -> productLoSubDifferential(aVector[index], xVector[index], xVector[index+systemSizeHalf]),
        1..systemSizeHalf
    )) - identityVector

end

function productHiSubDifferential(aInterval, xLo, xHi)
    aInterval.lo * dXLoPos(xLo) + aInterval.hi * dXHiNeg(xHi) - dMaxLoHi(aInterval, aInterval.lo*xLo, aInterval.hi*xHi)
end

function productHiSubDifferential(aInterval, xLoNegative, xHi)
    - aInterval.lo * dXLoNeg(xHi) - aInterval.hi * dXHiPos(xLo) + dMaxHiLo(aInterval, aInterval.lo*xHi, aInterval.hi*xLo)
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

function initialConditions(aMatrix, bVector)
    aMatrixMid = map(x -> mid(x), aMatrix)
    leftMatrix = ILESolver.constituentMatrix(aMatrixMid)
    rightVector = STI(bVector)
    leftMatrix \ transpose(rightVector)
end

end


function constituentMatrix(matrix)
    @assert ndims(matrix) == 2
    @assert size(matrix, 1) == size(matrix, 2)

    positivMatrix = map(x -> x >= 0 ? x : 0, matrix)
    negativMatrix = map(x -> x < 0 ? -x : 0, matrix)
    [
        positivMatrix negativMatrix;
        negativMatrix positivMatrix
    ]
end

function STI(intervalVector)
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

    reduce(stiReducer, [], intervalVector)
end

function reverseSTI(intervalVector)
    @assert ndims(intervalVector) == 1

    vectorSize = size(intervalVector, 1)
    halfVectorSize = int(vectorSize/2)
    negativLoIntervalPart = slice(intervalVector, 1:halfVectorSize)
    positivHiIntervalPart = slice(intervalVector, (1 + halfVectorSize):vectorSize)

    positivLoIntervalPart = map(-, negativLoIntervalPart)
    zippedIntervalParts = zip(positivLoIntervalPart, positivHiIntervalPart)
    map(x -> @interval(x[1], x[2]), zippedIntervalParts)
end

end
