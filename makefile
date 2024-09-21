doc:
	julia --project=docs docs/make.jl
	npm --prefix ./docs run docs:dev
