# Miniprojekt 5 Zeitoptimale Umsteuerung Brückenkran Numerische Methoden der Optimierung und Optimalen Steuerung
# Exercise 3, direct collocation
# Authors: Laura Kleckner and Marcel Dzubba
# Summer term 2022
# Code inspired by lecture examples!
#
# Excel >> Matlab
#
############# choose case manually!!! #############
choose_problem = 0 # 0: a_min = - a_max; 1: a_min = 0
##################################################
# model parameters

K_x = 0.99
l_s = 1.0
g   = 9.81
v_soll = 0.4
v_min = -10.0
a_max  = 0.33
if choose_problem == 1
    a_min = 0.0
elseif choose_problem == 0
    a_min = -a_max
else
    error("Choose choose_problem to zero or one.")
end
phi_0 = 0
dphi_0 = 0
phi_tf = 0
dphi_tf = 0

N = 3*40 # collocation points

using JuMP, Ipopt, PyCall, PyPlot, DelimitedFiles
# Interior Point solver suits well for collocation problems
mod = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1000, "mumps_mem_percent" => 8000))
@variable(mod, tf >= 0.1, start = 1.5918761117048705)           # time interval [0, tf]
@variable(mod, a_min/K_x <= u[0:N-1] <= a_max/K_x, start = 0)   # initialize Input for k = 0... N-1
@variable(mod, x1[0:N])                                         # -""-       x1 = x for k = 0... N
@variable(mod, x2[0:N])                                         # -""-       x2 = 𝛗
@variable(mod, x3[0:N])                     	                # -""-       x3 = ̇𝛗

# initial values
@constraint(mod, x1[0] == 0)
@constraint(mod, x2[0] == phi_0)
@constraint(mod, x3[0] == dphi_0)
# final values
@constraint(mod, x1[N] == v_soll)
@constraint(mod, x2[N] == phi_tf)
@constraint(mod, x3[N] == dphi_tf)

# direct collocation
@NLconstraint(mod, x1_collo[i=0:N-1],x1[i+1] - x1[i] == tf / N * u[i])
@NLconstraint(mod, x2_collo[i=0:N-1],x2[i+1] - x2[i] == tf / N * (x3[i]+x3[i+1])/2.0)
@NLconstraint(mod, x3_collo[i=0:N-1],x3[i+1] - x3[i] == tf / N * (- u[i]-g*(x2[i]+x2[i+1])/2.0)/l_s)

# julia magic starting here
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
    # print switiching times
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
    # save states + input AND Dual values (= Lagranian multiplicators).
    # Those dual values can be utilized for an adequate Initialization of the bvp4c algorithm
    if choose_problem == 1
        writedlm("kran_a_min_0.dat", [collect(0:N)/N*value(tf) [value.(u.data); NaN] value.(x1.data) value.(x2.data) value.(x3.data)])
        writedlm("kran_a_min_0_dual.dat", [(collect(0:N-1).+0.5)/N*value(tf) dual.(x1_collo.data) dual.(x2_collo.data) dual.(x3_collo.data) ])
    else choose_problem == 0
        writedlm("kran.dat", [collect(0:N)/N*value(tf) [value.(u.data); NaN] value.(x1.data) value.(x2.data) value.(x3.data)])
        writedlm("kran_dual.dat", [(collect(0:N-1).+0.5)/N*value(tf) dual.(x1_collo.data) dual.(x2_collo.data) dual.(x3_collo.data) ])
    end
end
