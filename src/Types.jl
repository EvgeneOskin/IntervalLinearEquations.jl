module Types

using ValidatedNumerics;

immutable IntervalVector{T}
    intervals :: Array{Interval{T}, 1}
    sti :: Array{T, 1}
end

immutable Configuration{T}
  a :: Array{Interval{T}, 2}
  b :: IntervalVector{T}
  size :: Int
  sti_size :: Int
end


type Solver{T}
    current
    previous
    initial
    system :: Configuration{T}
    roots
    iternation :: Int
end

end
