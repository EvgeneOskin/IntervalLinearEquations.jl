module StiTests
using FactCheck

using Base.Test;
using ValidatedNumerics;
using IntervalLinearEquations, IntervalLinearEquations.sti;

const bVector = [
    @interval(-2, 2); @interval(-2, 2);
] :: Array{Interval{Float64}};
const m = [
    1 -2;
    -0 -5.1;
];

facts("Calculate costrituent matrix") do
    interval_product = m*bVector
    sti_product = IntervalLinearEquations.sti.reverseSTI(
        IntervalLinearEquations.sti.constituentMatrix(m)*
        IntervalLinearEquations.sti.STI(bVector)
    )
    @fact interval_product --> sti_product
end

facts("Calculate sti and reverse sti") do
    transformed = IntervalLinearEquations.sti.reverseSTI(
        IntervalLinearEquations.sti.STI(bVector)
    )
    @fact bVector --> transformed
end

end

