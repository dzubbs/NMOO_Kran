# does not work!
t1= 0.697
t2= 0.893
tN= 1.592
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
u_max = a_max/K_x
N = 200

u0 = zeros(1,N+1)
for i in 1:N+1
    if tN/N*i <= t1
        u0[i] = u_max
    elseif tN/N*i<= t2
        u0[i] = -u_max
    else
        u0[i] = u_max
    end
end

using JuMP, Ipopt, PyCall, PyPlot, DelimitedFiles
mod = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1000, "mumps_mem_percent" => 8000))
@variable(mod, tf >= 0.1, start = tN)              # Optimierungshorizont [0, tf]
@variable(mod, a_min/K_x <= u[i=0:N] <= a_max/K_x, start = u_max)
@variable(mod, -10.0 <= x1[0:N] <= 10.0)
@variable(mod, x2[0:N], start = 0)
@variable(mod, x3[0:N], start = 0)
@variable(mod, p1[0:N], start = -1.5025)
@variable(mod, p2[0:N], start = 3.4699)
@variable(mod, p3[0:N], start = 1.4894)


#@constraint(mod, u[0] == 0)
@constraint(mod, x1[0] == 0)
@constraint(mod, x2[0] == phi_0)
@constraint(mod, x3[0] == dphi_0)
@constraint(mod, x1[N] == v_soll)
@constraint(mod, x2[N] == phi_tf)
@constraint(mod, x3[N] == dphi_tf)
@constraint(mod, (1.0+ p1[N]*u[N]+p2[N]*x3[N]-p3[N]*(x2[N]*g/l_s+1/l_s*u[N])) == 0.0)
@constraint(mod, u[0] ==  u_max)
@constraint(mod, u[N] ==  u_max)

for i in 0:N-1
    @NLconstraint(mod, x1[i+1] - x1[i] == tf / N * u[i])
    @NLconstraint(mod, x2[i+1] - x2[i] == tf / N * (x3[i]+x3[i+1])/2.0)
    @NLconstraint(mod, x3[i+1] - x3[i] == tf / N * (- u[i]-g*(x2[i]+x2[i+1])/2.0)/l_s)
    @NLconstraint(mod, p1[i+1] - p1[i] == 0.0)
    @NLconstraint(mod, p2[i+1] - p2[i] == tf/N*(g/l_s)*(p3[i]+p3[i+1])/2.0)
    @NLconstraint(mod, p3[i+1] - p3[i] == -tf/N*(p2[i]+p2[i+1])/2.0)
    #@NLconstraint(mod, u[i]*(p1[i]-p3[i]/l_s) <=0)

end

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
    p1_sol = value.(p1.data)
    p2_sol = value.(p2.data)
    p3_sol = value.(p3.data)
    u_sol  = value.(u.data)
#    for i in 1:N-1
#        if abs(u_sol[i]-u_sol[i+1]) > 0.1
#            t_i = t[i]
#            println(t_i)
#        end
#    end
    H = 1 .+ p1_sol.*u_sol.+p2_sol.*x3_sol.-p3_sol.*(x2_sol*g/l_s.+u_sol/l_s)
    global x1_sol, x2_sol, x3_sol, u_sol, p1_sol, p2_sol, p3_sol

    plt.figure("Direkte Kollokation Kran", figsize=(10, 13))
    clf()
    subplot(5,1,1)
    xlabel("time")
    ylabel("u")
    plot(t, u_sol)
    subplot(5,1,2)
    xlabel("time")
    ylabel("u")
    plot(t, p1_sol-p3_sol/l_s)
    subplot(5,1,3)
    plot(t, x1_sol)
    xlabel("time")
    ylabel("dx/dt")
    subplot(5,1,4)
    plot(t,x2_sol)
    xlabel("time")
    ylabel("phi")
    subplot(5,1,5)
    plot(t,x3_sol)
    xlabel("time")
    ylabel("dphi")

    plt.figure("H", figsize=(10, 13))
    clf()
    subplot(4,1,1)
    xlabel("time")
    ylabel("u")
    plot(t, p1_sol)
    subplot(4,1,2)
    xlabel("time")
    ylabel("u")
    plot(t, p2_sol)
    subplot(4,1,3)
    xlabel("time")
    ylabel("u")
    plot(t, p3_sol)
    subplot(4,1,4)
    xlabel("time")
    ylabel("u")
    plot(t, H)
end
