
# Modellparameter
K_x = 0.99
l_s = 1.0
g   = 9.81
v_soll = 0.4
v_min = -10.0
a_max  = 0.33
a_min = 0.0
phi_0 = 0
dphi_0 = 0
phi_tf = 0
dphi_tf = 0

N = 3*30

using JuMP, Ipopt, PyCall, PyPlot, DelimitedFiles
mod = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1000, "mumps_mem_percent" => 8000))
@variable(mod, tf >= 0.1, start = 1.5918761117048705)              # Optimierungshorizont [0, tf]
@variable(mod, a_min/K_x <= u[0:N-1] <= a_max/K_x, start = 0)
@variable(mod, x1[0:N])
@variable(mod, x2[0:N])
@variable(mod, x3[0:N])
#@constraint(mod, u[0] == 0)
@constraint(mod, x1[0] == 0)
@constraint(mod, x2[0] == phi_0)
@constraint(mod, x3[0] == dphi_0)
@constraint(mod, x1[N] == v_soll)
@constraint(mod, x2[N] == phi_tf)
@constraint(mod, x3[N] == dphi_tf)
@NLconstraint(mod, x1_collo[i=0:N-1],x1[i+1] - x1[i] == tf / N * u[i])
@NLconstraint(mod, x2_collo[i=0:N-1],x2[i+1] - x2[i] == tf / N * (x3[i]+x3[i+1])/2.0)
@NLconstraint(mod, x3_collo[i=0:N-1],x3[i+1] - x3[i] == tf / N * (- u[i]-g*(x2[i]+x2[i+1])/2.0)/l_s)



@objective(mod, Min, tf)

optimize!(mod)
println("Solver Status: ", termination_status(mod))
if termination_status(mod) == JuMP.MathOptInterface.LOCALLY_SOLVED

    close("all")
    println("Zielfunktional: ", objective_value(mod))
    println("Endzeit: ", value(tf))

    # Ergebnisse auslesen
    t = collect(0:N) * value(tf) / N                                  # Zeit
    x1_sol = value.(x1.data)
    x2_sol = value.(x2.data)
    x3_sol = value.(x3.data)
    u_sol  = value.(u.data)
    for i in 1:N-1
        if abs(u_sol[i]-u_sol[i+1]) > 0.2
            t_i = t[i]
            println(t_i)
        end
    end
    t[N]
    global x1_sol, x2_sol, x3_sol, u_sol

    plt.figure("Direkte Kollokation Kran", figsize=(10, 13))
    clf()
    subplot(4,1,1)
    xlabel("time")
    ylabel("u")
    plot(t[1:N], u_sol)
    subplot(4,1,2)
    plot(t, x1_sol)
    xlabel("time")
    ylabel("dx/dt")
    subplot(4,1,3)
    plot(t,x2_sol)
    xlabel("time")
    ylabel("phi")
    subplot(4,1,4)
    plot(t,x3_sol)
    xlabel("time")
    ylabel("dphi")
    if a_min == 0
        writedlm("kran_a_min_0.dat", [collect(0:N)/N*value(tf) [value.(u.data); NaN] value.(x1.data) value.(x2.data) value.(x3.data)])
        writedlm("kran_a_min_0_dual.dat", [(collect(0:N-1).+0.5)/N*value(tf) dual.(x1_collo.data) dual.(x2_collo.data) dual.(x3_collo.data) ])

    else
        writedlm("kran.dat", [collect(0:N)/N*value(tf) [value.(u.data); NaN] value.(x1.data) value.(x2.data) value.(x3.data)])
        writedlm("kran_dual.dat", [(collect(0:N-1).+0.5)/N*value(tf) dual.(x1_collo.data) dual.(x2_collo.data) dual.(x3_collo.data) ])

    end

end
