direct collocation: 
	kran_3.jl
	Dependencies: JuMP, Ipopt, PyCall, PyPlot, DelimitedFiles
	change choose_problem in line 10 manually:
		choose_problem = 0: a_min = -a_max
		choose_problem = 1: a_min = 0 (no neg. vel.)

indirect collocation:
	kran_5.m
	Dependencies: kran.dat, kran_a_min_0.dat, kran_a_min_0_dual.dat, kran_dual.dat 
			starting values for bvp4c. Will also be generated when running kran_3.jl
	change choose_problem in line 10 manually:
		choose_problem = 0: a_min = -a_max
		choose_problem = 1: a_min = 0 (no neg. vel.)
