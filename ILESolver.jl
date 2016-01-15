module ILESolver

using ValidatedNumerics;

module F

function solve(A, B)
    @assert ndims(A) == 2
    @assert ndims(B) == 1
    @assert size(A, 1) == size(A, 2)
    @assert size(A, 1) == size(B, 1)

    reverseSTI(STI(B))
end

function equation(xArg, aMatrix, bVector)
    STI(aMatrix*reverseSTI(xArg)) - xArg + STI(bVector)
end

function initialConditions(aMatrix, bVector)
    aMatrixMid = map(x -> mid(x), aMatrix)
    leftMatrix = eye(size(bVector, 1)) - constituentMatrix(aMatrixMid)
    rightVector = STI(bVector)
    leftMatrix \ transpose(rightVector)
end

function subDifferential(aMatrix, xVatrox)
end

end

module G

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
    leftMatrix = constituentMatrix(aMatrixMid)
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
