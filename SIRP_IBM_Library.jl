using StatsBase
using Plots

struct Grid{T <: Integer}
   
    x :: T
    y :: T
    
    N :: T
    
    r :: Array{T} #Array where r[x] is the new index after moving to the right
    l :: Array{T}
    u :: Array{T}
    d :: Array{T}
    
    id :: Array{T} #Array of indices of the grid
    
    borders :: Array{T}
    
end

struct Parameters{T <: Real}
   
    rates :: Array{T} # γ, λ, μ, β, κ
    
    total_rates :: Array{T} # γ_tot, λ_tot, μ_tot, β_tot, κ_tot
    
end

struct Variables{T <: Integer}
   
    S :: Array{T} #Flat Array of positions of S
    I :: Array{T} #Flat Array of positions of I
    R :: Array{T} #Flat Array of positions of R
    
    P :: Array{T} #Flat Array of P in each cell
    
end

function move_squared(L)
    
    r = zeros(L^2)
    l = zeros(L^2)
    u = zeros(L^2)
    d = zeros(L^2)
    
    for i in 1 : L^2
        
        #Right
        if i % L == 0

            r[i] = i - L + 1

        else

            r[i] = i + 1

        end
        
        #Left
        if (i - 1) % L == 0

            l[i] = i + L - 1

        else

            l[i] = i - 1

        end
        
        #Up
        u[i] = mod1(i + L, L^2)
        
        #Down
        d[i] = mod1(i - L, L^2)
        
        
    end
        
    return r, l, u, d
    
end

function get_system_borders(id, L)
    
    borders = []
    
    M = reshape(id, (L, L))
   
    for item in M[:, 1]
       
        append!(borders, item)
        
    end
    
    for item in M[1, :]
       
        append!(borders, item)
        
    end
    
    for item in M[L, :]
       
        append!(borders, item)
        
    end
    
    for item in M[:, L]
       
        append!(borders, item)
        
    end
    
    return borders
    
end

function construct_squared_grid(L, T)
   
    r, l, u, d = move_squared(L)
    
    id = [i for i in 1 : L^2]
    
    borders = get_system_borders(id, L)
    
    grid = Grid{T}(L, L, L^2, r, l, u, d, id, borders)
    
end

@views function compute_waiting_time(vars, params)

    length_I = length(vars.I)
    sum_P = sum(vars.P)
    
    params.total_rates[1] = length_I * params.rates[1] #N_I * γ
    params.total_rates[2] = length_I * params.rates[2] #N_I * λ
    
    params.total_rates[3] = sum_P * params.rates[3] #N_P * μ
    
    if length(vars.P[vars.S]) != 0
    
        params.total_rates[4] = sum(vars.P[vars.S]) * params.rates[4] #\sum_{i=1}^N_S β*P_{cell}=β\sum_{i=1}^N_S P_{cell}
        
    else
        
        params.total_rates[4] = 0.0
        
    end
        
    params.total_rates[5] = sum_P * params.rates[5]
    
    W = sum(params.total_rates)

    tau = -1/W * log(rand())
    
    return tau, W
    
end

#Microscopic RW
function parasite_mobility_RW(grid, params, vars)
    
    valid_Ps_idx = grid.id[vars.P .> 0] #findall(vars.P .> 0)
    valid_Ps_N = vars.P[valid_Ps_idx]
    
    idx = sample(valid_Ps_idx, Weights(Float64.(valid_Ps_N)))
    
    #Choice movement
    opt = rand([1, 2, 3, 4])
    
    if opt == 1 #up
    
        new_idx = mod1(idx + grid.x, grid.N)
        
    elseif opt == 2 #down
        
        new_idx = mod1(idx - grid.x, grid.N)
        
    elseif opt == 3 #right
    
        if idx % grid.x == 0

            new_idx = idx - grid.x + 1

        else

            new_idx = idx + 1

        end
        
    else #left
    
        if (idx - 1) % grid.x == 0

            new_idx = idx + grid.x - 1

        else

            new_idx = idx - 1

        end
        
    end
    
    vars.P[idx] -= 1
    vars.P[new_idx] += 1    
end

@views function select_random_parasite(grid, params, vars)

    valid_Ps_idx = grid.id[vars.P .> 0] #findall(vars.P .> 0)
    valid_Ps_N = vars.P[valid_Ps_idx]
    
    idx = sample(valid_Ps_idx, Weights(Float64.(valid_Ps_N)))
    
    return idx
    
end

@views function select_random_S(params, vars)
   
    weights = vars.P[vars.S] .* (params.rates[4]/params.total_rates[4]) #Array where each item is the number of P for each Pina cell
    S_ind = [i for i in 1 : length(vars.S)]

    idx = sample(S_ind, Weights(weights)) #Pinna that gets infected
    
    return idx
    
end
    
@views function parasite_mobility(grid, params, vars)

    #Select idx
    idx = select_random_parasite(grid, params, vars)
    
    new_idx = rand([grid.r[idx], grid.l[idx], grid.u[idx], grid.d[idx]])
    
    vars.P[idx] -= 1
    vars.P[new_idx] += 1   
    
end

#Parasite death
@views function parasite_death(grid, params, vars)

    idx = select_random_parasite(grid, params, vars)
    
    vars.P[idx] -= 1
    
end

#Pinna infection
@views function infect_Pinna(grid, params, vars)
   
    idx = select_random_S(params, vars)

    vars.P[vars.S[idx]] -= 1

    append!(vars.I, vars.S[idx]) #Append the position of the S individual in I list

    splice!(vars.S, idx) #Delete the position of S list
    
end

#Pina death
@views function Pinna_death(grid, params, vars)
   
     #All infected pinna have same prob of dying, choose one at random
    idx = rand(1 : length(vars.I))

    append!(vars.R, vars.I[idx])

    splice!(vars.I, idx)
    
end

@views function parasite_production(grid, params, vars)
   
    #All infected pinna have same proba of producing parasite, choose one at random

    vars.P[rand(vars.I)] += 1 #Append the created parasite in the position s
    
end

@views function choose_apply_event(reactions, orders, W, grid, params, vars)
    
    U = rand() * W
    
    for i in 1 : 5
       
        if U < sum(params.total_rates[orders[1:i]])
           
            #Execute corresponding reaction
            reactions[i](grid, params, vars)
            
            if i != 1
            
                aux_r = copy(reactions)
                aux_o = copy(orders)

                reactions[i-1] = aux_r[i]
                reactions[i] = aux_r[i-1]
                
                orders[i-1] = aux_o[i]
                orders[i] = aux_o[i-1]
                
            end
            
            break
            
        end
        
    end
    
    return reactions, orders
        
end 

function get_xy_positions(M, flat_arr)
    
    Cartesian_ind = CartesianIndices(M)[flat_arr]

    y = [idx[1] for idx in Cartesian_ind]
    x = [idx[2] for idx in Cartesian_ind]
    
    return x, y
    
end

function IBM(t_end, grid, vars, params)
    
    S_t = [length(vars.S)]
    I_t = [length(vars.I)]
    R_t = [length(vars.R)]
    P_t = [sum(vars.P)]
    
    time = [0.0]
    
    reactions = [Pinna_death, parasite_production, parasite_death, infect_Pinna, parasite_mobility]

    orders = [1, 2, 3, 4, 5]

    t = 0
    
    while t < t_end
        
        #Compute total rate
        tau, W = compute_waiting_time(vars, params)
        
        if tau < 0 || (length(vars.I) == 0 && sum(vars.P) == 0)
            
            break
            
        end
       
        #Choose event to happen
        reactions, orders = choose_apply_event(reactions, orders, W, grid, params, vars)
        
        t += tau
        
        #Record values of interest
        append!(S_t, length(vars.S))
        append!(I_t, length(vars.I))
        append!(R_t, length(vars.R))
        append!(P_t, sum(vars.P))        
        
        append!(time, t)
        
    end
    
    return S_t, I_t, R_t, P_t, time
    
end

function IBM_fast(t_end, grid, vars, params; filename="record.txt")

    t = 0
    
    reactions = [Pinna_death, parasite_production, parasite_death, infect_Pinna, parasite_mobility]

    orders = [1, 2, 3, 4, 5]
    
    f = open(filename, "w")
    
    println(f, "#t\tS\tI\tR\tP")
    
    println(f, t, "\t", length(vars.S), "\t", length(vars.I), "\t", length(vars.R), "\t", sum(vars.P))
    
    while t < t_end
    
        #Compute total rate
        tau, W = compute_waiting_time(vars, params)
        
        if tau < 0 || (length(vars.I) == 0 && sum(vars.P) == 0)
            
            break
            
        end

        #Choose event to happen
        reactions, orders = choose_apply_event(reactions, orders, W, grid, params, vars)
        
        t += tau
        
        println(f, t, "\t", length(vars.S), "\t", length(vars.I), "\t", length(vars.R), "\t", sum(vars.P))
        
    end
    
    close(f)
    
    return vars.S, vars.I, vars.R, vars.P #Return last config
    
end

function IBM_final_state(grid, vars, params)
    
    reactions = [Pinna_death, parasite_production, parasite_death, infect_Pinna, parasite_mobility]

    orders = [1, 2, 3, 4, 5]
    
    while (length(vars.I) > 0 || sum(vars.P) > 0)
        
        #Compute total rate
        tau, W = compute_waiting_time(vars, params)

        #Choose event to happen
        reactions, orders = choose_apply_event(reactions, orders, W, grid, params, vars)
        
    end
    
    return length(vars.S), length(vars.R)
    
end

function epidemic_velocity(grid, vars, params)
    
    reactions = [Pinna_death, parasite_production, parasite_death, infect_Pinna, parasite_mobility]

    orders = [1, 2, 3, 4, 5]
    
    N_I = 0
    t = 0
    
    while (length(vars.I) > 0 || sum(vars.P) > 0)
        
        #Compute total rate
        tau, W = compute_waiting_time(vars, params)

        t += tau
        
        #Choose event to happen
        reactions, orders = choose_apply_event(reactions, orders, W, grid, params, vars)
        
        length_I = length(vars.I)
        
        if length_I > N_I
            
            #Check if new I in borders
            if vars.I[end] in grid.borders
               
                return t
                
            end
           
            N_I = length_I
            
        end
        
    end
    
    return t
    
end

function IBM_final_vars(t_end, grid, vars, params)
    
    reactions = [Pinna_death, parasite_production, parasite_death, infect_Pinna, parasite_mobility]

    orders = [1, 2, 3, 4, 5]

    t = 0
    
    while t < t_end
        
        #Compute total rate
        tau, W = compute_waiting_time(vars, params)
        
        if tau < 0 || (length(vars.I) == 0 && sum(vars.P) == 0)
            
            break
            
        end
       
        #Choose event to happen
        reactions, orders = choose_apply_event(reactions, orders, W, grid, params, vars)
        
        t += tau
        
    end
    
    return vars.S, vars.I, vars.R, vars.P
    
end

function IBM_animation(t_end, grid, vars, params; every_frame=1, ms=8, clims=(0, 100))
    
    S_t = [length(vars.S)]
    I_t = [length(vars.I)]
    R_t = [length(vars.R)]
    P_t = [sum(vars.P)]
    
    time = [0.0]

    t = 0
    
    #Initialize plots
    M = reshape(vars.P, (grid.x, grid.y))
    
    X_S, Y_S = get_xy_positions(M, vars.S)
    X_I, Y_I = get_xy_positions(M, vars.I)
    X_R, Y_R = get_xy_positions(M, vars.R)

    N_P = sum(vars.P)

    N_S = length(vars.S)
    N_I = length(vars.I)
    N_R = length(vars.R)

    p1 = heatmap(M, c=:blues, clims=clims)
    
    scatter!(X_S, Y_S, xlim=(0.5, x+0.5), ylim=(0.5, y+0.5), color=:green3, m=:v, ms=ms, label="")
    scatter!(X_I, Y_I, xlim=(0.5, x+0.5), ylim=(0.5, y+0.5), color=:red, m=:v, ms=ms, label="")
    scatter!(X_R, Y_R, xlim=(0.5, x+0.5), ylim=(0.5, y+0.5), color=:black, m=:v, ms=ms, label="")
    
    p2 = plot(time, P_t, label="P", color=:deepskyblue)
    
    p3 = plot(time, [S_t I_t R_t], color=[:green3 :red :black], labels=["S" "I" "R"], legend=:topleft)
    
    plot(p1, p2, p3, layout=@layout[a [b ; c]], size=(1200, 500))
    
    reactions = [Pinna_death, parasite_production, parasite_death, infect_Pinna, parasite_mobility]

    orders = [1, 2, 3, 4, 5]
    
    anim = @animate for i in 1 : 1000000 
    
        if t >= t_end
        
            break
        
        end
        
        #Compute total rate
        tau, W = compute_waiting_time(vars, params)
        
        if tau < 0 || (length(vars.I) == 0 && sum(vars.P) == 0)
            
            break
            
        end

        #Choose event to happen
        reactions, orders = choose_apply_event(reactions, orders, W, grid, params, vars)

        t += tau
        
        #Plot
        M = reshape(vars.P, (grid.x, grid.y))
            
        X_S, Y_S = get_xy_positions(M, vars.S)
        X_I, Y_I = get_xy_positions(M, vars.I)
        X_R, Y_R = get_xy_positions(M, vars.R)

        N_P = sum(vars.P)

        N_S = length(vars.S)
        N_I = length(vars.I)
        N_R = length(vars.R)
        
        append!(S_t, N_S)
        append!(I_t, N_I)
        append!(R_t, N_R)
        append!(P_t, N_P)
        
        append!(time, t)
        
        p1 = heatmap(M, c=:blues, clims=clims)
        
        scatter!(X_S, Y_S, xlim=(0.5, x+0.5), ylim=(0.5, y+0.5), color=:green3, m=:v, ms=ms, label="")
        scatter!(X_I, Y_I, xlim=(0.5, x+0.5), ylim=(0.5, y+0.5), color=:red, m=:v, ms=ms, label="")
        scatter!(X_R, Y_R, xlim=(0.5, x+0.5), ylim=(0.5, y+0.5), color=:black, m=:v, ms=ms, label="")
        
        p2 = plot(time, P_t, color=:deepskyblue, lw=3, label="P", legend=:none)
        
        p3 = plot(time, [S_t I_t R_t], color=[:green3 :red :black], lw=3, legend=:none)
        
        plot(p1, p2, p3, layout=@layout[a [b ; c]], size=(1200, 500))
        
    end every every_frame
    
    return anim
    
end

function IBM_animation_fast(t_end, grid, vars, params; every_frame=1, color="default", size=(1000, 800))
    
    t = 0
    
    reactions = [Pinna_death, parasite_production, parasite_death, infect_Pinna, parasite_mobility]

    orders = [1, 2, 3, 4, 5]
    
    #Initialize plots
    M = reshape(vars.P, (grid.x, grid.y))
    
    X_S, Y_S = get_xy_positions(M, vars.S)
    X_I, Y_I = get_xy_positions(M, vars.I)
    X_R, Y_R = get_xy_positions(M, vars.R)

    Mp = zeros(grid.x, grid.y)
    
    for i in 1 : length(X_S)
       
        Mp[X_S[i], Y_S[i]] = 2
        
    end
    
    for i in 1 : length(X_I)

        Mp[X_I[i], Y_I[i]] = 1

    end
    
    for i in 1 : length(X_R)
       
        Mp[X_R[i], Y_R[i]] = -1
        
    end

    p1 = heatmap(Mp, color=color)
    
    plot(p1, size=size)
    
    anim = @animate for i in 1 : 1000000 
    
        if t >= t_end
        
            break
        
        end
        
        #Compute total rate
        tau, W = compute_waiting_time(vars, params)
        
        if tau < 0 || (length(vars.I) == 0 && sum(vars.P) == 0)
            
            break
            
        end

        #Choose event to happen
        reactions, orders = choose_apply_event(reactions, orders, W, grid, params, vars)

        t += tau
        
        #Plot
        M = reshape(vars.P, (grid.x, grid.y))
            
        X_S, Y_S = get_xy_positions(M, vars.S)
        X_I, Y_I = get_xy_positions(M, vars.I)
        X_R, Y_R = get_xy_positions(M, vars.R)

        Mp = zeros(grid.x, grid.y)

        for i in 1 : length(X_S)

            Mp[X_S[i], Y_S[i]] = 2

        end
        
        for i in 1 : length(X_I)

            Mp[X_I[i], Y_I[i]] = 1

        end

        for i in 1 : length(X_R)

            Mp[X_R[i], Y_R[i]] = -1

        end

        p1 = heatmap(Mp, color=color)
        
        plot(p1, size=size)
        
    end every every_frame
    
    return anim
    
end
