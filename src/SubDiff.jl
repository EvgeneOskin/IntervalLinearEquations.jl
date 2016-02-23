module ILESolver

module SubDiff
using ILESolver.sti
using ValidatedNumerics;


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

end
