using Distributed

addprocs(10)

@everywhere include("SIRP_IBM_Library.jl")
@everywhere using SharedArrays

@everywhere function R_0_MF(params; ρ_S=1.0)
    
    return  (params[2]/params[1]) / (params[3]/(ρ_S*params[4]) + 1)
    
end

function phase_transition_parallel(L, M, lambda_ext)
    
	filename = string("PT_κ_", lambda_ext, "_L_", L, ".txt")
    
	f = open(filename, "w")

	println(f, "#R_0\tm_S\tErr_m_S\txi_S\tErr_xi_S\tm_R\tErr_m_R\txi_R\tErr_xi_R")

	result_k = SharedArray{Float64}(6)

    param_arr = [1, lambda_ext, 1, 1, 0.54]

    if lambda_ext > 30
        
        M = 500

    end

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

		S_inf, R_inf = IBM_final_state(grid, vars, params)
		
		result_k[1] += S_inf
		result_k[2] += S_inf^2
		result_k[3] += S_inf^4
		
		result_k[4] += R_inf
		result_k[5] += R_inf^2
		result_k[6] += R_inf^4
		
	end

    PHI =  R_0_MF(param_arr; ρ_S=1.0)

    S0 = L^2

	M_S = result_k[1] / M
	M_squared_S = result_k[2] / M
	M_forth_S = result_k[3] / M

	Xi_S = M_squared_S - M_S^2
	Xi_S_σ = M_forth_S - M_squared_S^2

	Err_M_S = sqrt(abs(Xi_S) / M)
	Err_Xi_S = sqrt(abs(Xi_S_σ) / M)
		
	M_R = result_k[4] / M
	M_squared_R = result_k[5] / M
	M_forth_R = result_k[6] / M

	Xi_R = M_squared_R - M_R^2
	Xi_R_σ = M_forth_R - M_squared_R^2

	Err_M_R = sqrt(abs(Xi_R) / M)
	Err_Xi_R = sqrt(abs(Xi_R_σ) / M)

	println(f, PHI, "\t", M_S/S0, "\t", Err_M_S/S0, "\t", Xi_S/S0, "\t", Err_Xi_S/S0, "\t", M_R/S0, "\t", Err_M_R/S0, "\t", Xi_R/S0, "\t", Err_Xi_R/S0)

	#println(kappa_Ext, "\t", M_S/S0, "\t", Err_M_S/S0, "\t", Xi_S/S0, "\t", Err_Xi_S/S0, "\t", M_R/S0, "\t", Err_M_R/S0, "\t", Xi_R/S0, "\t", Err_Xi_R/S0)
    
    close(f)
    
end

lambda_ext = parse(Float64, ARGS[1])

L = 100

M = 1000

phase_transition_parallel(L, M, lambda_ext)
