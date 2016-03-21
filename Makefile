.PHONY: test lint coverage

test:
	julia -e 'Pkg.init(); run(`ln -sf $$(pwd()) $$(Pkg.dir("ILESolver"))`); Pkg.pin("ILESolver"); Pkg.resolve()'
	julia --code-coverage test/runtests.jl

coverage:
	julia -e 'cd(Pkg.dir("ILESolver")); Pkg.add("ILESolver"); using Coverage; @show get_summary(Coverage.process_folder())'

lint:
	julia -e 'using Lint; for i in readdir("src") lintfile(string("src/", i)) end'
