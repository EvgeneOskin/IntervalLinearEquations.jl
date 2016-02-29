module SubDiff
using ValidatedNumerics;
include("sti.jl")
using Debug;

function productLoSubDifferential{T}(aInterval :: Interval{T}, xLo :: T, xHi :: T)
    maxFactor = dMaxHiLo(aInterval, sti.pos(aInterval.hi)*sti.pos(xLo), sti.neg(aInterval.lo)*sti.pos(xHi))
    result = sti.pos(aInterval.lo) * dXLoNeg(xLo) + sti.neg(aInterval.hi) * dXHiNeg(xHi) - maxFactor
    result
end

function productHiSubDifferential{T}(aInterval :: Interval{T}, xLo :: T, xHi :: T)
    maxFactor = dMaxHiHi(aInterval, sti.pos(aInterval.hi)*sti.pos(xHi), sti.neg(aInterval.lo)*sti.pos(xLo))
    result = - sti.pos(aInterval.lo) * dXHiNeg(xHi) - sti.neg(aInterval.hi) * dXLoNeg(xLo) + maxFactor
    result
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
        [sti.pos(c.hi); 0]
    elseif left == right
        0.5*[sti.pos(c.hi); sti.neg(c.lo)]
    else
        [0; sti.neg(c.lo)]
    end
end

function dMaxHiHi(c, left, right)
    if left > right
        [0; sti.pos(c.hi)]
    elseif left == right
        0.5*[sti.neg(c.lo); sti.pos(c.hi)]
    else
        [sti.neg(c.lo); 0]
    end
end

end
