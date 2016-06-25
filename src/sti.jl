module sti
using ValidatedNumerics;

pos = x -> x >= zero(x) ? x : zero(x)
neg = x -> x < zero(x) ? -x : zero(x)

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
        halfReducedSize = Int(reducedSize/2)
        vcat(
            slice(reduced, 1:halfReducedSize),
            [-element.lo],
            slice(reduced, (1 + halfReducedSize):reducedSize),
            [element.hi]
        ) # Could use cat instead
    end

    reduce(stiReducer, [], intervalVector) :: Array{T}
end

function reverseSTI{T}(intervalVector :: Array{T})
    @assert ndims(intervalVector) == 1

    vectorSize = size(intervalVector, 1)
    halfVectorSize = Int(vectorSize/2)
    negativLoIntervalPart = slice(intervalVector, 1:halfVectorSize)
    positivHiIntervalPart = slice(intervalVector, (1 + halfVectorSize):vectorSize)

    positivLoIntervalPart = map(-, negativLoIntervalPart)
    zippedIntervalParts = zip(positivLoIntervalPart, positivHiIntervalPart)
    collect(Interval{T}, map(x -> safeInterval(x[1], x[2]), zippedIntervalParts)) :: Array{Interval{T}}
end

function safeInterval(a, b)
    a < b ? @interval(a, b) : @interval(b, a)
end

end
