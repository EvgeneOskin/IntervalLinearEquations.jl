# Interval Linear Equations Solver

Julia package to solve [system of linear equations](https://en.wikipedia.org/wiki/System_of_linear_equations)
with [intervals](https://en.wikipedia.org/wiki/Interval_(mathematics)) coefficients instead of [numbers](https://en.wikipedia.org/wiki/Number).

## Installation

```julia
julia> Pkg.add("IntervalLinearEquations")
```

## Usage

```julia
julia> using IntervalLinearEquations;
julia> using ValidatedNumerics
julia> b = [@interval(-2.0, 2.0); @interval(-2.0, 2.0)]
julia> a = [
            @interval(2, 4) @interval(-2, 1);
            @interval(-1, 2) @interval(2, 4);
       ]
julia> solution = IntervalLinearEquations.solve(
        "G",  # G to solve Ax=b, F to solve Ax - x = b
        a,  # Matrix
        b,  # Vector
        0.1,  # precision
        10,  # iterations limit
        1.0  # Scale factor
        )
```

## How do the results compare to using `\`?

The results from this package are better than from the `\`.

```julia
A =
    [
        @interval( 2, 4)    @interval(-2, 1) ;
        @interval(-1, 2)    @interval( 2, 4)
    ]

b =
    [
        @interval(-2, 2),
        @interval(-2, 2)
    ]
x = A \ b
A * x != b
x = IntervalLinearEquations.solve("G", a, b, 0.1, 10, 1.0)
A * x == b
```

So the default implementation of `\` operator return solution that doesn't pass
the substitution check.
