PROJECT_NAME=IntervalLinearEquations
.PHONY: test lint coverage

test:
	julia -e 'Pkg.init(); run(`ln -sf $$(pwd()) $$(Pkg.dir("${PROJECT_NAME}"))`); Pkg.pin("${PROJECT_NAME}"); Pkg.resolve()'
	julia --code-coverage test/runtests.jl

coverage:
	julia -e 'cd(Pkg.dir("${PROJECT_NAME}")); Pkg.add("${PROJECT_NAME}"); using Coverage; @show get_summary(Coverage.process_folder())'

lint:
	julia -e 'using Lint; for i in readdir("src") if ismatch(r".*\.jl$$", i) lintfile(string("src/", i)) end end'
