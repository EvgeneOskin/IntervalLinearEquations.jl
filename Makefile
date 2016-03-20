.PHONY=test lint

test:
	julia -e 'Pkg.init(); run(`ln -s $(pwd()) $(Pkg.dir("ILESolver"))`); Pkg.pin("ILESolver"); Pkg.resolve()'
	julia --code-coverage test/runtests.jl

lint:
	julia -e 'using Lint; for i in readdir("src") lintfile(string("src/", i)) end'
