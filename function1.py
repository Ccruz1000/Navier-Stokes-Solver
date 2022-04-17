import numpy as np
delta_t_CFL = 1/(np.abs(u[1:num_y-2,1:num_x-2])/dx + np.abs(v[1:num_y-2,1:num_x-2])/dy + np.sqrt(gamma*R*T[1:num_y-2,1:num_x-2])*sqrt(1/dx**2 + 1/dy**2) + 2*v_prime*(1/dx**2 + 1/dy**2))
#################################################Two point scheme##########################

#First order accurate forward difference

def fir_o_for_y(dy,u,N):
    N = u.shape[0]
    dudy = np.zeros_like(u)
    dudy[0] = (u[1]-u[0]) / (dy)

    for i in np.arange(1, N-1):
        dudy[i] = (u[i+1]-u[i]) / (dy)
    return dudy

def fir_o_for_y2d(dy,u,N):
    N = u.shape[0]
    dudy = np.zeros_like(u)
    dudy[0] = (u[1]-u[0]) / (dy)

    for i in np.arange(1, N-1):
        dudy[i] = (u[i+1]-u[i]) / (dy)
    return dudy

#x-direction
def fir_o_for_x(dx,u,N):
    N = u.shape[0]
    dudx = np.zeros_like(u)
    dudx[0] = (u[1]-u[0]) / (dx)

    for i in np.arange(1, N-1):
        dudx[i] = (u[i+1]-u[i]) / (dx)

    return dudx



#First order accurate backward difference

def fir_o_back_y(dy,u,N):
    N = u.shape[0]
    dudy = np.zeros_like(u)
    dudy[0] = (u[1]-u[0]) / (dy)

    for i in np.arange(1, N-1):
        dudy[i] = (u[i]-u[i-1]) / (dy)

    return dudy

#2d
def fir_o_back_y2d(dy, u, N):
    N = u.shape[0]
    dudy1 = np.zeros((len(u), 2))
    dudy1[ :, 0] = (u[ :, 1] - u[ :, 0]) / dy
    dudy1[:,N - 1] = (u[:, -2] - u[:, -1]) / dy


    for i in np.arange(0, 2):
        for j in np.arange(1, N - 1):
            dudy1[i, j] = (u[i, j] - u[i, j-1]) / dy
    return dudy1


#x-direction
def fir_o_back_x(dx,u,N):
    N = u.shape[0]
    dudx = np.zeros_like(u)
    dudx[0] = (u[1]-u[0]) / (dx)

    for i in np.arange(1, N-1):
        dudx[i] = (u[i]-u[i-1]) / (dx)

    return dudx


#2d
def fir_o_back_x2d(dx, u, N):
    N = u.shape[0]
    dudx1 = np.zeros((len(u), 2))
    dudx1[0, :] = (u[1, :] - u[0, :]) / dx
    dudx1[N - 1, :] = (u[-2, :] - u[-1, :]) / dx


    for j in np.arange(0, 2):
        for i in np.arange(1, N - 1):
            dudx1[i, j] = (u[i, j] - u[i-1, j]) / dx
    return dudx1


# First order accurate central difference
def fir_o_cen_y(dy,u,N):
    N = u.shape[0]
    dudy = np.zeros_like(u)
    dudy[0] = (u[2]-u[0]) / (2*dy)

    for i in np.arange(1, N-1):
        dudy[i] = (u[i+1]-u[i-1]) / (2*dy)

    return dudy


# x-direction
def fir_o_cen_x(dx,u,N):
    N = u.shape[0]
    dudx = np.zeros_like(u)
    dudx[0] = (u[2]-u[0]) / (2*dx)

    for i in np.arange(1, N-1):
        dudx[i] = (u[i+1]-u[i-1]) / (2*dx)

    return dudx









##############################################################Three point scheme#######################

#Three point forward difference for first derivative
def three_for_first_y(dy,u,N):
    N = u.shape[0]
    dudy = np.zeros_like(u)
    dudy[0] = (-3*u[0]+4*u[1]-u[2]) / (2*dy)

    for i in np.arange(1, N-1):
        dudy[i] = (-3*u[i]+4*u[i+1]-u[i+2]) / (2*dy)

    return dudy

#x-direction
def three_for_first_x(dx,u,N):
    N = u.shape[0]
    dudx = np.zeros_like(u)
    dudx[0] = (-3*u[0]+4*u[1]-u[2]) / (2*dx)

    for i in np.arange(1, N-1):
        dudx[i] = (-3*u[i]+4*u[i+1]-u[i+2]) / (2*dx)

    return dudx









#Three point backward difference for first derivative

def three_back_first_y(dy,u,N):
    N = u.shape[0]
    dudy = np.zeros_like(u)
    dudy[0] = (u[0]-4*u[1]+3*u[2]) / (2*dy)

    for i in np.arange(1, N-1):
        dudy[i] = (u[i-2]-4*u[i-1]+3*u[i]) / (2*dy)

    return dudy
#x-direction
def three_back_first_x(dx,u,N):
    N = u.shape[0]
    dudx = np.zeros_like(u)
    dudx[0] = (u[0]-4*u[1]+3*u[2]) / (2*dx)

    for i in np.arange(1, N-1):
        dudx[i] = (u[i-2]-4*u[i-1]+3*u[i]) / (2*dx)

    return dudx








##########################################Second derivative##########################################

#Three point backward difference for second derivative

def three_back_sec_y(dy,u,N):
    N = u.shape[0]
    dudy2 = np.zeros_like(u)
    dudy2[0] = (u[0]-2*u[1]+u[2]) / (dy**2)

    for i in np.arange(1, N-1):
        dudy2[i] = (u[i-2]-2*u[i-1]+u[i]) / (dy**2)

    return dudy2
#x-direction
def three_back_sec_x(dx,u,N):
    N = u.shape[0]
    dudx2 = np.zeros_like(u)
    dudx2[0] = (u[0]-2*u[1]+u[2]) / (dx**2)

    for i in np.arange(1, N-1):
        dudx2[i] = (u[i-2]-2*u[i-1]+u[i]) / (dx**2)

    return dudx2








#Three point central difference for second derivative

def three_cen_sec_y(dy,u,N):
    N = u.shape[0]
    dudy2 = np.zeros_like(u)
    dudy2[0] = (u[0]-2*u[1]+u[2]) / (dy**2)

    for i in np.arange(1, N-1):
        dudy2[i] = (u[i-1]-2*u[i]+u[i+1]) / (dy**2)

    return dudy2

#2d
def three_cen_secy2d(dy, u, N):
    N = u.shape[0]
    dudy2 = np.zeros((len(u), 2))
    dudy2[:, 0] = (u[:, 2] + u[:, 0] - 2 * u[:, 1]) / (dy ** 2)
    dudy2[:, N-1] = 0

    for i in np.arange(0, 1):
        for j in np.arange(1, N - 1):
            dudy2[i, j] = (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) / (dy ** 2)
    return dudy2









#x-direction
def three_cen_sec_x(dx,u,N):
    N = u.shape[0]
    dudx2= np.zeros_like(u)
    dudx2[0] = (u[0]-2*u[1]+u[2]) / (dx**2)

    for i in np.arange(1, N-1):
        dudx2[i] = (u[i-1]-2*u[i]+u[i+1]) / (dx**2)

    return dudx2

#2d
def three_cen_sec_x2d(dx, u, N):
    N = u.shape[0]
    dudx2 = np.zeros((len(u), 2))
    dudx2[0, :] = (u[2, :] + u[0, :] - 2 * u[1, :]) / (dx ** 2)
    dudx2[N - 1, :] = 0

    for j in np.arange(0, 1):
        for i in np.arange(1, N - 1):
            dudx2[i, j] = (u[i+1, j] - 2 * u[i, j] + u[i-1, j]) / (dx ** 2)
    return dudx2


# Three point forward difference for second derivative
def three_for_sec_y(dy, u, N):
    N = u.shape[0]
    dudy2 = np.zeros_like(u)
    dudy2[0] =(u[0]-2*u[1]+u[2]) / (dy**2)

    for i in np.arange(1, N-1):
        dudy2[i] = (u[i]-2*u[i+1]+u[i+2]) / (dy**2)

    return dudy2

#x-direction
def three_for_sec_x(dx,u,N):
    N = u.shape[0]
    dudx2 = np.zeros_like(u)
    dudx2[0] =(u[0]-2*u[1]+u[2]) / (dx**2)

    for i in np.arange(1, N-1):
        dudx2[i] = (u[i]-2*u[i+1]+u[i+2]) / (dx**2)

    return dudx2












