
# this code for the schoemer rattleback is mostly taken from the code provided by Hendrik Ranocha
# which inturn is based on code by Elmar Schömer

function setup_rattleback_schoemer(ω0 = [0.01, -0.02, -2.0])
    p = @SVector [0.5, 0.5, 0.0]
    E = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    parameters = (mE = 0.1, mp = 0.05, g = 9.81,
    p = p, E = E)

    goldenratio = 0.5 * (1 + sqrt(5))
    (; m, d, B0, I0) = inertia(goldenratio, 1, 1,
            parameters)
    parameters = (; m, d, B0, I0, parameters...)

    c = @SVector [0, 0, 1 - d]
    q = Quaternion(1.0, 0, 0, 0)
    w = @SVector [ω0[1], ω0[2], ω0[3]]

    R = to_rotation_matrix(q)
    B = R * B0 * R'
    B_1 = inv(B)

    e3 = SVector(0, 0, 1)
    r = -B_1 * e3 / sqrt(e3' * B_1 * e3)
    r = r + d * R * e3
    v = cross(r, w)


    y = [c[1], c[2], c[3], 
             q[1], q[2], q[3], q[4],
             v[1], v[2], v[3],
             w[1], w[2], w[3]]
    
  
    return y, parameters
end



function inertia(a, b, c, parameters)
    a2 = a^2
    b2 = b^2
    c2 = c^2

    B0 = @SMatrix [inv(a2) 0 0; 0 inv(b2) 0; 0 0 inv(c2)]

    (; mE, mp, p, E) = parameters
    m = mE + 2 * mp

    Ixx = 0.2 * mE * (b2 + c2)
    Iyy = 0.2 * mE * (a2 + c2)
    Izz = 0.2 * mE * (a2 + b2)

    I0 = @SMatrix [Ixx 0 0; 0 Iyy 0; 0 0 Izz]

    P = dot(p, p) * E - p * p'

    d = mE / m * 3 * c/ 8
    cr = @SVector [0, 0, -d]
    CR = dot(cr, cr) * E - cr * cr'

    I0 = I0 + 2 * mp * P - m * CR
    
    return (; m, d, B0, I0)
end

function drdt(w, r, B_1)
    e3 = SVector(0, 0, 1)
    s = -transpose(e3) * r
    result = cross(w, r) + B_1 * cross(w, e3) / s + r * (e3' * cross(w, r) / s)
    return result
end


function rattleback!(dy, y, parameters, t)
    c, q, v, w = unpack_rattleback(y)
    (; m, g, d, B0, I0, E) = parameters # E is the identity matrix
    e3 = SVector(0, 0, 1)

    R = to_rotation_matrix(q)
    B = R * B0 * R'
    B_1 = inv(B)
    I = R * I0 * R'

    r = -B_1 * e3 / sqrt(e3' * B_1 * e3)
    r_0A = r 
    dr = drdt(w, r, B_1)
    r = r + d * R * e3
    dr = dr + d * cross(w, R * e3)
    r2 = dot(r, r)
    I_1 = inv(I - r * r' * m + E * (m * r2))
    u = cross(r, e3) * m * g - cross(w, I * w) - cross(r, cross(w, dr)) * m
    dw = I_1 * u
    dv = -cross(dw, r) - cross(w, dr)
    dc = v
    qw = Quaternion(0, w[1], w[2], w[3])
    dq = 0.5 * qw * q
    pack_rattleback!(dy, dc, dq, dv, dw)
    return nothing
end


function unpack_rattleback(y)
    c = @SVector [y[1], y[2], y[3]]
    q = Quaternion(y[4], y[5], y[6], y[7])
    # Quaternionic.normalize(q)
    v = @SVector [y[8], y[9], y[10]]
    w = @SVector [y[11], y[12], y[13]]
    return c, q, v, w
end

function pack_rattleback!(dy, dc, dq, dv, dw)
    dy[1] = dc[1]
    dy[2] = dc[2]
    dy[3] = dc[3]
    dy[4] = dq[1]
    dy[5] = dq[2]
    dy[6] = dq[3]
    dy[7] = dq[4]
    dy[8] = dv[1]
    dy[9] = dv[2]
    dy[10] = dv[3]
    dy[11] = dw[1]
    dy[12] = dw[2]
    dy[13] = dw[3]
    return dy
end




# reduced version of schoemer rattleback
# basically the same as the original version, but with the velocity v removed
# and dc is now calculate as dc = -cross(w, r) instead of dc = v.
function setup_rattleback_schoemer_reduced(ω0 = [0.01, -0.02, -2.0])
    p = @SVector [0.5, 0.5, 0.0]
    E = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    parameters = (mE = 0.1, mp = 0.05, g = 9.81,
    p = p, E = E)

    goldenratio = 0.5 * (1 + sqrt(5))
    (; m, d, B0, I0) = inertia(goldenratio, 1, 1,
            parameters)
    parameters = (; m, d, B0, I0, parameters...)

    c = @SVector [0, 0, 1 - d]
    q = Quaternion(1.0, 0, 0, 0)
    w = @SVector [ω0[1], ω0[2], ω0[3]]

    R = to_rotation_matrix(q)
    B = R * B0 * R'
    B_1 = inv(B)

    e3 = SVector(0, 0, 1)
    r = -B_1 * e3 / sqrt(e3' * B_1 * e3)
    r = r + d * R * e3



    y = [c[1], c[2], c[3], 
             q[1], q[2], q[3], q[4],
             w[1], w[2], w[3]]
    
  
    return y, parameters
end

function rattleback_reduced!(dy, y, parameters, t)
    c, q, w = unpack_rattleback_reduced(y)
    (; m, g, d, B0, I0, E) = parameters # E is the identity matrix
    e3 = SVector(0, 0, 1)

    R = to_rotation_matrix(q)
    B = R * B0 * R'
    B_1 = inv(B)
    I = R * I0 * R'

    r = -B_1 * e3 / sqrt(e3' * B_1 * e3)
    r_0A = r 
    dr = drdt(w, r, B_1)
    r = r + d * R * e3
    dr = dr + d * cross(w, R * e3)
    r2 = dot(r, r)
    I_1 = inv(I - r * r' * m + E * (m * r2))
    u = cross(r, e3) * m * g - cross(w, I * w) - cross(r, cross(w, dr)) * m
    dw = I_1 * u
    dc = -cross(w, r)
    qw = Quaternion(0, w[1], w[2], w[3])
    dq = 0.5 * qw * q
    pack_rattleback_reduced!(dy, dc, dq, dw)
    return nothing
end


function unpack_rattleback_reduced(y)
    c = @SVector [y[1], y[2], y[3]]
    q = Quaternion(y[4], y[5], y[6], y[7])
    w = @SVector [y[8], y[9], y[10]]
    return c, q, w
end

function pack_rattleback_reduced!(dy, dc, dq, dw)
    dy[1] = dc[1]
    dy[2] = dc[2]
    dy[3] = dc[3]
    dy[4] = dq[1]
    dy[5] = dq[2]
    dy[6] = dq[3]
    dy[7] = dq[4]
    dy[8] = dw[1]
    dy[9] = dw[2]
    dy[10] = dw[3]
    return dy
end


# plot the omega values (in körperfesten komponenten) of the schoemer solution
function plot_schoemer_omega(sol, pos = :topleft)	
    omega1 = zeros(length(sol))
    omega2 = zeros(length(sol))
    omega3 = zeros(length(sol))

    for i in 1:length(sol)
        c, q, v, w = unpack_rattleback(sol[i])
        R = to_rotation_matrix(q)
        ω = R' * w

        omega1[i] = ω[1]
        omega2[i] = ω[2]
        omega3[i] = ω[3]

    end

    lw = 2
    fontsize = 14
    fontsize2 = 18

    p = plot(sol.t, omega1, label=L"\tilde{\omega}_1",
        xguide = L"t",
        yguide = L"\tilde{\omega}_i",
        title = L"\tilde{\omega}" *" "* L"(t)" *" - Wackelstein Schömer", 
        legend = pos,
        size = (700, 500),
        grid = true,
        xlims=(first(sol.t), last(sol.t)),
        dpi = 300,

        lw = lw,
        titlefont = font(fontsize2,"Computer Modern"),
        guidefont = font("Computer Modern"),
        legendfont = font(fontsize,"Computer Modern"),
        xtickfontsize=fontsize,
        ytickfontsize=fontsize,
        ylabelfontsize=fontsize2,
        xlabelfontsize=fontsize2,
    )

        plot!(sol.t, omega2, label=L"\tilde{\omega}_2", lw = lw)
        plot!(sol.t, omega3, label=L"\tilde{\omega}_3", lw = lw)

    return p
end

# calculate the energy of the schoemer rattleback
function calculate_energy_schoemer(sol, parameters)
    (; m, g, d, B0, I0) = parameters
    Eges = zeros(length(sol))
    Ekin = zeros(length(sol))
    Erot = zeros(length(sol))
    Epot = zeros(length(sol))
    for i in 1:length(sol)
        c, q, v, w = unpack_rattleback(sol[i])

        R = to_rotation_matrix(q)
        B = R * B0 * R'
        B_1 = inv(B)
        e3 = SVector(0, 0, 1)

        
        h_S = c[3]

        I_t = R * I0 * R'

        Ekin[i] = 0.5 * m * dot(v, v)
        Erot[i] = 0.5 * w' * I_t * w     
        Epot[i] = m*g*h_S

        Eges[i] = Ekin[i] + Erot[i] + Epot[i]
    end
    return Eges, Ekin, Erot, Epot
end




# calculate the energy of the schoemer_red rattleback
function calculate_energy_schoemer_reduced(sol, parameters)
    (; m, g, d, B0, I0) = parameters
    Eges = zeros(length(sol))
    Ekin = zeros(length(sol))
    Erot = zeros(length(sol))
    Epot = zeros(length(sol))
    for i in 1:length(sol)
        c, q, w = unpack_rattleback_reduced(sol[i])

        R = to_rotation_matrix(q)
        B = R * B0 * R'
        B_1 = inv(B)
        e3 = SVector(0, 0, 1)
        r = -B_1 * e3 / sqrt(e3' * B_1 * e3)
        r = r + d * R * e3
        
        h_S = c[3]
        #h_S = -r[3]
        
        I_t = R * I0 * R'
        v = -cross(w, r)

        Ekin[i] = 0.5 * m * dot(v, v)
        Erot[i] = 0.5 * w' * I_t * w     
        Epot[i] = m*g*h_S

        

        Eges[i] = Ekin[i] + Erot[i] + Epot[i]
    end
    return Eges, Ekin, Erot, Epot
end