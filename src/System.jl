module ILESolver

module System

immutable System{T}
  a :: Array{Interval{T}}
  b :: IntervalVector{T}
  size :: Int
  sti_size :: Int
end

immutable IntervalVector{T}
    intervals :: Array{Interval{T}}
    sti :: Array{T}
end


type Solver{T}
    current :: IntervalVector{T}
    previous :: IntervalVector{T}
    initial :: IntervalVector{T}
    system :: System{T}
    roots :: Array{T}
    iternation :: Int
end

end
end
