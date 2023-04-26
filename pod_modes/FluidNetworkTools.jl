module FluidNetworkTools
using LinearAlgebra
using Plots
using Statistics
# using MATLAB

using PotentialFlow

# -----------------------------------------------------------------------------
# Fluid Network Tools [2D fluid flow]
# Date Created  : 03/30/2022
# Date Modified : --/--/----
# By SBD
# Reference:    A. G. Nair & K. Taira
#               "Network-Theoretic Approach to Sparsified Discrete Vortex
#               Dynamics"; Journal of Fluid Mechanics; 768; 549-571; 2015
#               (doi.org/10.1017/jfm.2015.97)
# -----------------------------------------------------------------------------
# [A] = adjacency_mat(x,y,gamma,alpha,type)
# Function to calculate adjacency matrix given x; y & gamma.
# A         : adjacency matrix using algebraic mean (default)
# x         : x coordinate vector
# y         : y coordinate vector
# gamma     : circulation vector
# alpha     : network direction parameter [0.5=undirected]
# type      : algebraic or geometric mean ("alg' [default] | 'goem")
# -----------------------------------------------------------------------------

#export ADJ, NetworkCentroids, BoundVortexlets

# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

# Bezier functions
function compute_bernstein(i,n; steps=100)
    return [binomial(n,i)*t^i*(1-t)^(n-i) for t in LinRange(0,1,steps)]
end;

function compute_bernstein_poly(px,py; steps=100)
    n = length(px)-1
    bernsteins = [compute_bernstein(i,n) for i=0:n]
    x_vals = [sum(px[k]*bernsteins[k][t] for k=1:n+1) for t=1:steps]
    y_vals = [sum(py[k]*bernsteins[k][t] for k=1:n+1) for t=1:steps]
    return x_vals, y_vals
end;

# Bexier plotting
function plot_with_bernstein(px,py; steps=100, subplot=1)
    x_vals, y_vals = compute_bernstein_poly(px,py; steps=steps)
    plot(x_vals, y_vals, color=:blue, label="", subplot=subplot)
end;

function plot_with_bernstein!(px,py; steps=100, subplot=1)
    x_vals, y_vals = compute_bernstein_poly(px,py; steps=steps)
    plot!(x_vals, y_vals, color=:blue, label="", linealpha = 0.3, linewidth = 5, subplot=subplot)
end;

# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

function plot_fluid_domain(xn, yn, ω)
    scatter(xn, 
    yn, 
    marker_z = ω/maximum(abs.(ω)),
    color = :RdBu,
    legend=false,
    markerstrokewidth=0,
    markersize=3,
    ylims = (-2, 2), # <---- this needs to be handeled programatically!!!
    xlims = (-0.5, 10),
    aspect_ratio=:equal,
    clims = (-0.5, 0.5),
    grid = false,
    axis = ([], false),
    dpi = 250)

end

function plot_vortical_network(centroids, ω̄, x̄, ȳ, xb, yb, itr)

    # Plot vortex-field as scatter
    scatter(x̄, 
            ȳ, 
            marker_z = ω̄/maximum(abs.(ω̄)), 
            color = :RdBu,
            legend=false,
            markerstrokewidth=0,
            markersize=1.5,
            #xlims = (-1,4),
            ylims = (-2,2), # <---- this needs to be handeled programatically!!!
            aspect_ratio=:equal,
            clims=(-0.5, 0.5),
            grid = false,
            axis = ([], false));

    # Plot network
    FluidNetworkTools.plot_network!(centroids);

    # Plot body
    plot!(xb[:,itr],
          yb[:,itr],
          linewidth = 3,
          c = :black,
          xlims = (-1,6),
          ylims = (-2,2),
          aspect_ratio=:equal)

end

function plot_network!(centroids)

    # Plot in- and out-degree edges
    # for point = 1:1:size(centroids, 1)-1
       
    #     # need a programatic way to set the control point value
    #     px = [centroids[point,1], 2, centroids[point+1,1]]
    #     py = [centroids[point,2], -3, centroids[point+1,2]]
        
    #     FluidNetworkTools.plot_with_bernstein!(px, py);
    #     FluidNetworkTools.plot_dual_with_bernstein!(px, py);

    # end

    # Plot centroids
    scatter!(centroids[:,1], 
             centroids[:,2], 
             c = :black, 
             markersize = abs.(centroids[:,3])/maximum(abs.(centroids[:,3]))*10,
             aspect_ratio=:equal);

end;


function plot_network(centroids)

    # open plot to plot into
    plot(legend=false)

    # Plot in- and out-degree edges
    for point = 1:1:size(centroids, 1)-1
       
        # need a programatic way to set the control point value
        px = [centroids[point,1], 2, centroids[point+1,1]]
        py = [centroids[point,2], -3, centroids[point+1,2]]
        
        FluidNetworkTools.plot_with_bernstein!(px, py);
        FluidNetworkTools.plot_dual_with_bernstein!(px, py);

    end

    # Plot centroids
    scatter!(centroids[:,1], 
             centroids[:,2], 
             c = :black, 
             markersize = abs.(centroids[:,3])/maximum(abs.(centroids[:,3]))*10,
             aspect_ratio=:equal);

end;


# plot out-degree
function plot_dual_with_bernstein!(px,py; steps=100, subplot=1)
    hold_val = px[2]
    px[2] = py[2]
    py[2] = hold_val
    x_vals, y_vals = compute_bernstein_poly(px,py; steps=steps)
    plot!(x_vals, y_vals, color=:red, label="", linealpha = 0.3, linewidth = 5, linestyle = :solid, subplot=subplot)
    plot!(x_vals, y_vals, color=:red, label="", linealpha = 0.3, linewidth = 5, linestyle = :dash, subplot=subplot)
end;

# Add degree based opacity to in- and out-degree plots
function plot_with_bernstein_degree!(px,py,ind,outd; steps=100, subplot=1)
    x_vals, y_vals = compute_bernstein_poly(px,py; steps=steps)
    plot!(x_vals, y_vals, color=:blue, label="", linealpha = ind, linewidth = 5, subplot=subplot)
end;

function plot_dual_with_bernstein_degree!(px,py,ind,outd; steps=100, subplot=1)
    hold_val = px[2]
    px[2] = py[2]
    py[2] = hold_val
    x_vals, y_vals = compute_bernstein_poly(px,py; steps=steps)
    plot!(x_vals, y_vals, color=:red, label="", linealpha = outd, linewidth = 5, linestyle = :solid, subplot=subplot)
    plot!(x_vals, y_vals, color=:red, label="", linealpha = outd, linewidth = 5, linestyle = :dash, subplot=subplot)
end;

# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

function Centroids(ω, x, y, labels)
    
    n_clusters = 1 + maximum(labels)
    
    centroids = zeros(Float64, n_clusters, 3);

    for c = 0:1:n_clusters-1;
        centroids[c+1,1] = mean(x[findall(labels .== c)]);
        centroids[c+1,2] = mean(y[findall(labels .== c)]);
        centroids[c+1,3] = sum(ω[:][findall(labels .== c)]);
    end
    
    return centroids
end

function split_vorticity(ω̄, x̄, ȳ)
    
    ω̄₊ = ω̄[ω̄ .> 0];
    x̄₊ = x̄[ω̄ .> 0];
    ȳ₊ = ȳ[ω̄ .> 0];

    ω̄₋ = ω̄[ω̄ .< 0];
    x̄₋ = x̄[ω̄ .< 0];
    ȳ₋ = ȳ[ω̄ .< 0];
    
   return
    
end

# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

function ADJ(x, y, Gamma, alpha)
    # -----------------------------------------------------------------------------
    # Adjacency matrix calculator [2D fluid flow]
    # Date Created  : 03/29/2022
    # Date Modified : --/--/----
    # By SBD
    # Modified from MGM's MATLAB code
    # Reference:    A. G. Nair & K. Taira
    #               "Network-Theoretic Approach to Sparsified Discrete Vortex
    #               Dynamics"; Journal of Fluid Mechanics; 768; 549-571; 2015
    #               (doi.org/10.1017/jfm.2015.97)
    # -----------------------------------------------------------------------------
    # [A] = adjacency_mat(x,y,gamma,alpha,type)
    # Function to calculate adjacency matrix given x; y & gamma.
    # A         : adjacency matrix using algebraic mean (default)
    # x         : x coordinate vector
    # y         : y coordinate vector
    # gamma     : circulation vector
    # alpha     : network direction parameter [0.5=undirected]
    # -----------------------------------------------------------------------------

    pts = length(x)
    
    x = reshape(x,pts,1)
    y = reshape(y,pts,1)
    
    # position vector matrix, r = (xi,yj)
    r = zeros(Float64, pts, 2)
    r[:,1:2] = [x y];               

    println("Calculating adjacency matrix of ", pts, " nodes...")
    
    A = zeros(Float64,pts,pts)
    
    for i = 1:pts   # can use parfor if needed
        
        global u_ij, u_ji = BiotSavart(Gamma, r, i, pts)

        A[i,:] = alpha*u_ij .+ (1 - alpha)*u_ji
   
    end

    return A #u_ij, u_ji, A, dist_mag, dist_ij, r
        
end

function supra(x_struct, y_struct, x_fluid, y_fluid, Gamma_struct, Gamma_fluid, alpha)

    A = zeros(Float64, length(A_struct) + length(A_fluid), length(A_struct) + length(A_fluid));

    Gamma = vcat(Gamma_struct, Gamma_fluid)
    xc = vcat(xc_struct, xc_fluid);
    yc = vcat(yc_struct, yc_fluid);

    r = [xc yc]; # Position vector array

    pts = length(x)

    for i = 1:pts   # can use parfor if needed
        
        global u_ij, u_ji = BiotSavart(Gamma, r, i, pts)

        A[i,:] = alpha*u_ij .+ (1 - alpha)*u_ji
   
    end

    return A

end

function BiotSavart(Gamma, r, i, pts)
    # -----------------------------------------------------------------------------
    # Biot-Savart Law [2D fluid flow]
    # Date Created  : 03/29/2022
    # Date Modified : --/--/----
    # By SBD
    # Modified from MGM's MATLAB code
    # Reference:    A. G. Nair & K. Taira
    #               "Network-Theoretic Approach to Sparsified Discrete Vortex
    #               Dynamics"; Journal of Fluid Mechanics; 768; 549-571; 2015
    #               (doi.org/10.1017/jfm.2015.97)
    # -----------------------------------------------------------------------------
    # u_ij, u_ji = adjacency_mat(x,y,gamma,alpha)
    # Function to calculate adjacency matrix given x; y & gamma.
    # Gamma     : circulation vector
    # r         : 
    # i         : 
    # pts       : 
    # -----------------------------------------------------------------------------
    
    dist_ij = r .- ones(Float64,pts,1)*r[i,:]';      # distance vector, r_[i->j]
    dist_mag = sqrt.(sum(dist_ij.^2, dims = 2));     # ||r_[i->j]||
    dist_ij = 1 ./dist_mag;                          # distance ratio
    dist_ij[isinf.(dist_ij)] .= 0;                   # update r_[i->i] = 0
    u_ij = abs(Gamma[i]) * dist_ij;                  # u_[i->j]
    u_ji = abs.(Gamma) .* dist_ij;                   # u_[j->i]

    # -----------------------------------------------------------------------------
    #    dist_ij = r(1:pts,:) - ones(pts,1)*r(i,:);  % distance vector, r_{i->j}
    #    dist_mag = sqrt(sum(dist_ij.^2,2));         % ||r_{i->j}||
    #    dist_ij = 1./dist_mag;                      % distance ratio
    #    dist_ij(isinf(dist_ij)) = 0;                % update r_{i->i} = 0
    #    u_ij = abs(Gamma(i)) * dist_ij;             % u_{i->j}
    #    u_ji = abs(Gamma) .* dist_ij;               % u_{j->i}
    # -----------------------------------------------------------------------------

    return u_ij, u_ji #, dist_mag, dist_ij
    
end

# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

function BoundVortexlets(wn, xn, yn, vxb, vyb, xb, yb, delta)
    # -----------------------------------------------------------------------------
    # Biot-Savart Law [2D fluid flow]
    # Date Created  : 03/29/2022
    # Date Modified : --/--/----
    # By SBD
    # Modified from MGM's MATLAB code
    # Reference:    A. G. Nair & K. Taira
    #               "Network-Theoretic Approach to Sparsified Discrete Vortex
    #               Dynamics"; Journal of Fluid Mechanics; 768; 549-571; 2015
    #               (doi.org/10.1017/jfm.2015.97)
    # -----------------------------------------------------------------------------
    # u_ij, u_ji = adjacency_mat(x,y,gamma,alpha)
    # Function to calculate adjacency matrix given x; y & gamma.
    # x         : x coordinate vector
    # y         : y coordinate vector
    # gamma     : circulation vector
    # alpha     : network direction parameter [0.5=undirected]
    # u_ij      : out-degree
    # u_ji      : in-degree
    # -----------------------------------------------------------------------------
       
    # delta = 0.25; % smoothing parameter (published 0.0025)
    u = [1 0]; # freestream velocity
    nc = size(vxb,1); # number of control points

    n = unitNormals(xb, yb)
    
    # -------------------------- Normal Velocities -------------------------
    u_plate = [vyb vxb]; # set the normal velocity at each vortexlet equal to the normal velocity of the fluid
      
    vel_control = zeros(Float64, nc)
    
    for i = 1:nc
        vel_control[i] = (u_plate[i,1]^2 + u_plate[i,2]^2)^(1/2);
    end    
        
    # --------------------------- Body Matrix ----------------------------
    # co-locate control points with every vortexlet on the wing surface
    pos_control = hcat(vec(xb'), vec(yb'));  
    
    pos_vortex = pos_control;
    
    Mb = ones(Float64, nc, nc); 

    for i = 1:1:nc
        for j = 1:1:nc
            
            disp = pos_control[i,:] - pos_vortex[j,:]; # displacement from vortexlet to control point
            rsqr = disp[1].^2 + disp[2].^2; # radius squared of the displacement
            Mb[i,j] = sum(conj([disp[2]/(2π*(rsqr+delta^2)), -disp[1]/(2π*(rsqr+delta^2))]).*n[i,:]);
            
        end 
    end

    # -------------------------- Wake Matrix ---------------------------
    pos_vortex = hcat(xn, yn);    
    
    nw = length(wn)
    
    Mw = ones(Float64, nc, nw); 
    
    for i = 1:1:nc
        for j = 1:1:nw
            
            disp = pos_control[i,:] - pos_vortex[j,:]; # displacement from vortexlet to control point
            rsqr = disp[1].^2 + disp[2].^2; # radius squared of the displacement
            Mw[i,j] = sum(conj([disp[2]/(2π*(rsqr+delta^2)), -disp[1]/(2π*(rsqr+delta^2))]).*n[i,:]);
            
        end 
    end
    
    # ------------------ Solve for Body Strengths ---------------------    
    # solve for the strength of vortexlets gamma = inv(M)*-V
    strength_bound = Mb\(vel_control.*n .- Mw*wn); 
    
    circ_at_each_vortex = zeros(nc,2); 
    
    for pos = 1:nc
        circ_at_each_vortex[pos,2] = sum(strength_bound[1:pos]);
    end
        
    wb = circ_at_each_vortex[:,2];
    
    return wb

end

# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

function ADJ_from_edges(edges)

    dim = length(edges)/2 - 1;
    weights = zeros(Float64, dim,dim);

    for (u, v) in edges
        weights[u,v] = weights[v,u] = 1;
    end

    return weights

end

function generate_karate_club1()
    # derived from http://konect.cc/networks/ucidata-zachary/
    
    edges = [
        ( 1,  2), ( 1,  3), ( 2,  3), ( 1,  4), ( 2,  4),
        ( 3,  4), ( 1,  5), ( 1,  6), ( 1,  7), ( 5,  7),
        ( 6,  7), ( 1,  8), ( 2,  8), ( 3,  8), ( 4,  8),
        ( 1,  9), ( 3,  9), ( 3, 10), ( 1, 11), ( 5, 11),
        ( 6, 11), ( 1, 12), ( 1, 13), ( 4, 13), ( 1, 14),
        ( 2, 14), ( 3, 14), ( 4, 14), ( 6, 17), ( 7, 17),
        ( 1, 18), ( 2, 18), ( 1, 20), ( 2, 20), ( 1, 22),
        ( 2, 22), (24, 26), (25, 26), ( 3, 28), (24, 28),
        (25, 28), ( 3, 29), (24, 30), (27, 30), ( 2, 31),
        ( 9, 31), ( 1, 32), (25, 32), (26, 32), (29, 32),
        ( 3, 33), ( 9, 33), (15, 33), (16, 33), (19, 33),
        (21, 33), (23, 33), (24, 33), (30, 33), (31, 33),
        (32, 33), ( 9, 34), (10, 34), (14, 34), (15, 34),
        (16, 34), (19, 34), (20, 34), (21, 34), (23, 34),
        (24, 34), (27, 34), (28, 34), (29, 34), (30, 34),
        (31, 34), (32, 34), (33, 34),
    ]

    return ADJ_from_edges(edges)
end


function generate_karate_club()
    # derived from http://konect.cc/networks/ucidata-zachary/
    weights = zeros(34, 34)
    edges = [
        ( 1,  2), ( 1,  3), ( 2,  3), ( 1,  4), ( 2,  4),
        ( 3,  4), ( 1,  5), ( 1,  6), ( 1,  7), ( 5,  7),
        ( 6,  7), ( 1,  8), ( 2,  8), ( 3,  8), ( 4,  8),
        ( 1,  9), ( 3,  9), ( 3, 10), ( 1, 11), ( 5, 11),
        ( 6, 11), ( 1, 12), ( 1, 13), ( 4, 13), ( 1, 14),
        ( 2, 14), ( 3, 14), ( 4, 14), ( 6, 17), ( 7, 17),
        ( 1, 18), ( 2, 18), ( 1, 20), ( 2, 20), ( 1, 22),
        ( 2, 22), (24, 26), (25, 26), ( 3, 28), (24, 28),
        (25, 28), ( 3, 29), (24, 30), (27, 30), ( 2, 31),
        ( 9, 31), ( 1, 32), (25, 32), (26, 32), (29, 32),
        ( 3, 33), ( 9, 33), (15, 33), (16, 33), (19, 33),
        (21, 33), (23, 33), (24, 33), (30, 33), (31, 33),
        (32, 33), ( 9, 34), (10, 34), (14, 34), (15, 34),
        (16, 34), (19, 34), (20, 34), (21, 34), (23, 34),
        (24, 34), (27, 34), (28, 34), (29, 34), (30, 34),
        (31, 34), (32, 34), (33, 34),
    ]
    for (u, v) in edges
        weights[u,v] = weights[v,u] = 1
    end
    return weights
end

# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------



# function MATLAB_louvain(w, x, y, gamma)
    
#     # open MATLAB session
#     s1 = MSession();
    
#     # calculate Adjacency matrix
#     A = FluidNetworkTools.ADJ(x,y,w,1);
    
#     # pass variable to matlab
#     put_variable(s1, :A, A)
#     put_variable(s1, :gamma, gamma)
    
#     # call Louvain function
#     eval_string(s1, "labels = find_louvain_communities(gamma,A);");
    
#     # convert labels to julia array
#     labels = jarray(get_mvariable(s1, :labels));
    
#     # close MATLAB session
#     close(s1)  
    
#     return Int64.(labels)
# end    



function signed_community_centroids(ω̄, x̄, ȳ)
    
    # Split ± vorticity
    ω̄₊ = ω̄[ω̄ .> 0]
    x̄₊ = x̄[ω̄ .> 0]
    ȳ₊ = ȳ[ω̄ .> 0]

    ω̄₋ = ω̄[ω̄ .< 0]
    x̄₋ = x̄[ω̄ .< 0]
    ȳ₋ = ȳ[ω̄ .< 0];
    
    # run Louvain algorithm on +ω elements
    labels = MATLAB_louvain(ω̄₊, x̄₊, ȳ₊, 0.75)
    centroids = Centroids(ω̄₊, x̄₊, ȳ₊, labels)

    # run Louvain algorithm on -ω elements and concatinate
    labels = MATLAB_louvain(ω̄₋, x̄₋, ȳ₋, 0.75)
    
    return vcat(centroids, Centroids(ω̄₋, x̄₋, ȳ₋, labels));
    
end

# -----------------------------------------------------------------------------------------------------------------------
# POD Stuff
# -----------------------------------------------------------------------------------------------------------------------

struct PODBasis{T}
    coefficients::Matrix{T}
    modes::Matrix{T}
end


function PODsvd(X; subtractmean::Bool = true)

    if subtractmean
        X .-= mean(X,dims=2)
    end

    # Economy sized SVD
    F = svd(X)

    # Mode coefficients
    a = Diagonal(F.S)*F.Vt

    POD = PODBasis(a, F.U)

    return POD, F.S

end


function DMD(X,r)

    # ----------------------------------------------------------------------
    # DMD(X,r)
    # ----------------------------------------------------------------------
    # Preforms the exact Dynamic Mode Decomposition as Tu et al.
    # X - Full data array rows of mesurments, columns in time 
    # r - truncation parameter
    # From Data-Driven Science and Engineering, Brunton 2021 (section 7.2)
    # ----------------------------------------------------------------------



    # Step 1
    # Preform economy sized SVD on X
    # -----------------------------------------
    X .-= mean(X,dims=2) # Subtract mean

    F = svd(X[:,1:end-1])

    U = F.U
    Σ = Diagonal(F.S) 
    V = F.V 

    # Truncate to r
    U = U[:,1:r]
    Σ = Σ[1:r,1:r]
    V = V[:,1:r]

    # Step 2
    # Solve for the best-fit operator of size r 
    # -----------------------------------------
    Ã = U'*X[:,2:end]*V/Σ
   
    # Step 3
    # Preform spectral decomposition
    # -----------------------------------------
    Λ, W = eigen(Ã)

    # Step 4
    # Get DMD modes
    # -----------------------------------------
    Φ = X[:,2:end]*(V/Σ)*W;

    # Step 4+
    # Preform expansion
    # -----------------------------------------
    α₁ = Σ*V[1,:]
    b = (W*Diagonal(Λ))\α₁

    return Φ, Λ, b

end

# function DMDadv(X,r,Δt,t)

#     #----------------------------------------------------------------------
#     # DMDadv(X,r,Δt,t)
#     #----------------------------------------------------------------------
#     # Advance the solution of the DMD system to time: t
#     # X - Full data array: rows of mesurments, columns in time 
#     # r - truncation parameter
#     # Δt - time step size of data array
#     # t - time to return solution of DMD system for
#     # From Data-Driven Science and Engineering, Brunton 2021 (section 7.2)
#     #----------------------------------------------------------------------

#     Φ, Λ, b = DMD(X,r; subtractmean::Bool = true)

#     ω = log.(Λ)./Δt

#     Ω = Diagonal(ω)
    
#     return Φ.*exp(Ω.*t).*b

# end

# -----------------------------------------------------------------------------------------------------------------------
# Force Calculation
# -----------------------------------------------------------------------------------------------------------------------

function unitNormals(xb, yb)
    
    # --------------------- Normal Unit-Vectors ---------------------------
    n = normals(xb, yb)

    return n.*sqrt.(n[:,1].^2 + n[:,2].^2)

end


function normals(xb, yb)
    
    # --------------------- Normal Vectors ---------------------------
    n = zeros(Float64, length(xb), 2);

    n[:, 1] = prepend!(diff(yb),diff(yb)[1]).*length(xb);
    n[:, 2] = prepend!(diff(xb),diff(xb)[1]).*-length(xb);

    return n

end

function impulse(xb,yb,x,y,ẋb,ẏb,omg,delta)

    # $$P = \int_{V_f} x \ \times \ \omega \ dV + \int_{S_b} x \ \times \ (n \ \times \ v) \ dS$$

    Δx = y[1,1] - y[2,1]; # length of grid elements
    Δs = sqrt((xb[1,1] - xb[2,1])^2 + (yb[1,1] - yb[2,1])^2); # length of beam elements

    # Volume integral
    p = [Δx^2*sum(omg.*y'), Δx^2*sum(-omg.*x)]

    # Calculate bound vortexlet strengths
    ωb = BoundVortexlets(vec(omg), vec(x), vec(y), ẋb, ẏb, xb, yb, delta)

    # Calculate body normals
    nx, ny = unitNormals(xb, yb)

    # Surface integral
    z = (ωb./Δs + nx.*ẏb - ny.*ẋb).*Δs # Inner cross product (note ωb is scaled by Δs)
    p += [sum(yb'*z); -sum(xb'*z)] # Outer cross product

    return p

end


# ------------------------------------
# Jeff's code
# ------------------------------------

# function impulse(vm::VortexModel{Nb,Ne}, wphysical::Nodes{Dual}, fphysical::ScalarData) where {Nb,Ne}

#     xg, yg = coordinates(wphysical, vm.g)                     # SAM: get grid coordinates form wphisical, the grid vorticity and the vortex model
#     Δx = cellsize(vm.g)                                       # SAM: get grid cell size
#     _bodypointsvelocity!(vm._bodyvectordata, vm.bodies)       # SAM: ???

#     # Formula 61 (see formula 6.16 in book)               
#     p = [Δx^2*sum(wphysical.*yg'), Δx^2*sum(-wphysical.*xg)]  # SAM: Volume integral

#     for i in 1:Nb
#         r = getrange(vm.bodies,i)                             # SAM: range of points on body i.e 1:50
#         p += _impulsesurfaceintegral(vm.bodies[i].points, fphysical[r], vm._bodyvectordata.u[r], vm._bodyvectordata.v[r])
#     end

#     return p[1], p[2]
# end



# function _impulsesurfaceintegral(body::Body{N,RigidBodyTools.ClosedBody}, f, u, v) where {N}
#     nx,ny = normalmid(body)                                                   # SAM: get body normals
#     Δs = dlengthmid(body)                                                     # SAM: get lenght of body cells
#     return _crossproductsurfaceintegral(body,(f./Δs + nx.*v - ny.*u).*Δs)     # SAM: z = (f./Δs + nx.*v - ny.*u).*Δs
# end



# function _crossproductsurfaceintegral(body::Body, z)
#     rx,ry = collect(body)     # SAM: get xb, yb
#     return [ry'*z, -rx'*z]    # SAM: outer cross porduct
# end

# ------------------------------------
# ------------------------------------






function force(P, ρ, Δt)

    # Calculates force from impulse history
    # $$F = -\rho \frac{dP}{dt}$$

    return -ρ*diff(P)/Δt

end

function force(xb, yb, x, y, ẏb, ẋb, omg, delta, α, Δt, ρ)

    P = impulse(xb, yb, x, y, ẏb, ẋb, omg, delta, α)

    return -ρ*diff(P)/Δt

end


function StrainEnergy(xb,yb,K_B)
    
    # x and y components of flag position
    nb = size(xb,1)
    Δ = 1/(nb - 1)
    ds2 = (2*Δ)^2

    strain_energy = zeros(Float64, nb, size(xb,2))

    # 1st point needs one-sided stencil
    strain_energy[1,:] = K_B.*((1/ds2.*(xb[1,:] .- 2xb[2,:] .+ xb[3,:])).^2 .+ (1/ds2.*(yb[1,:] .- 2yb[2,:] .+ yb[3,:])).^2)

    # Last point needs one-sided stencil
    strain_energy[66,:] = K_B.*((1/ds2.*(xb[nb-2,:] .- 2xb[nb-1,:] .+ xb[nb,:])).^2 .+ (1/ds2.*(yb[nb-2,:] .- 2yb[nb-1,:] .+ yb[nb,:])).^2)

    for j = 2:nb-1

        strain_energy[j,:] = K_B.*((1/ds2.*(xb[j-1,:] .- 2xb[j,:] .+ xb[j+1,:]) ).^2 .+ (1/ds2*(yb[j-1,:] .- 2yb[j,:] .+ yb[j+1,:])).^2)

    end

    return strain_energy
    
end

# -----------------------------------------------------------------------------------------------------------------------
# State Dynamics
# -----------------------------------------------------------------------------------------------------------------------

function Ã(k,l,ξ) 

    # Adjaceny matrix function

    Γ = ξ[:,3] 
    
    return abs(Γ[l])./(2π*abs.(ξ[k,1:2] - ξ[l,1:2]))

end

function g̃(k,l,ξ) 

    # Interaction Function
    
    A = (ξ[k,1:2] - ξ[l,1:2])

    # we need to save ĵ values
    A2 = A[2]

    # k̂ X î = -ĵ 
    A[2] = -A[1]

    # k̂ X ĵ = î
    A[1] = A2

    return A./abs.(ξ[k,1:2] - ξ[l,1:2])
    
end

function stateDynamics(ξ)

    ξ̇ = zeros(Float64, size(ξ, 1), 2)

    for k = 1:1:size(ξ, 1)

        ξ̇ᵏ = zeros(Float64, 1, 2)

        for l = 1:1:size(ξ, 1)

            if l != k

                ξ̇ᵏ += (Ã(k,l,ξ).*g̃(k,l,ξ))'

            end
            
        end

    ξ̇[k,:] = ξ̇ᵏ 

    end

    return ξ̇
    
end

# -----------------------------------------------------------------------------------------------------------------------
# Potential Flow Friends
# -----------------------------------------------------------------------------------------------------------------------


function compute_ẋ!(ẋ, x, t)
    body, ambient_sys = x
    motion = ẋ[1]
    # update the instantaneous motion of the body with the current motion
    motion.ċ, motion.c̈, motion.α̇ = motion.kin(t)
    
    Bodies.enforce_no_flow_through!(body, motion, ambient_sys, t)
    
    # Zero the velocity
    reset_velocity!(ẋ, x)
    
    # Compute the self-induced velocity of the system
    self_induce_velocity!(ẋ, x, t)
    
    # Modify the velocity so that it provides the rate of change in the circle plane.
    Bodies.transform_velocity!(ẋ, ẋ, x, body)
end



# end module
end

# Data Cleaning Routines
# -------------------------------

function dat2csv(dat_path::AbstractString, csv_path::AbstractString)
    open(csv_path, write=true) do io
        for line in eachline(dat_path)
            join(io, split(line), ',')
            println(io)
        end
    end
    return csv_path
end

function dat2csv(dat_path::AbstractString)
    base, ext = splitext(dat_path)
    ext == ".dat" ||
        throw(ArgumentError("file name doesn't end with `.dat`"))
    return dat2csv(dat_path, "$base.csv")
end