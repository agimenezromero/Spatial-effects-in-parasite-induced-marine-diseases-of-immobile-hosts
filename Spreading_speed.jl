using Distributed

addprocs(10)

@everywhere include("SIRP_IBM_Library.jl")
@everywhere using SharedArrays

@everywhere function R_0_MF(params; ρ_S=1.0)
    
    return  (params[2]/params[1]) / (params[3]/(ρ_S*params[4]) + 1)
    
end

function phase_transition_parallel(lambda_ext, kappa_ext, Ls, M)

    param_arr = [1, lambda_ext, 1, 1, kappa_ext]

    R0 = R_0_MF(param_arr; ρ_S=1.0)
    
	filename = string("Epidemic_velocity_κ_", kappa_ext, "_R0_", R0, ".txt")
    
	f = open(filename, "w")

	println(f, "#R_0\tκ\tL\tt_L")

    for L in Ls
        
        println("L:", L)
        result_k = SharedArray{Float64}(1)

	    @sync @distributed for i in 1 : M

            S = 1 : L^2
            I = []
            R = []

            P0 = 50# * L^2

            P = zeros(L, L)

            P[Int32(floor(L/2)), Int32(floor(L/2))] = P0

            P = vec(P)

		    grid = construct_squared_grid(L, Int32)

		    params = Parameters{Float32}(param_arr, zeros(5))

		    vars = Variables{Int32}(S, I, R, P)

		    t_L = epidemic_velocity(grid, vars, params)

            result_k[1] += t_L
		    
	    end

        t_L_f = result_k[1] / M

	    println(f, R0, "\t", kappa_ext, "\t", L, "\t", t_L_f)

    end

    close(f)
    
end

lambda_ext = parse(Float64, ARGS[1])
kappa_ext = parse(Float64, ARGS[2])

Ls = [10, 20, 30, 40, 50, 60]

M = 1000

phase_transition_parallel(lambda_ext, kappa_ext, Ls, M)

