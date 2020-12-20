import numpy as np
import scipy.sparse as sparse
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import random

# def flatten(Mu, Mv):
#     #flatten two matrices into one array
#     return np.concatenate((Mu.flatten(), Mv.flatten()))
#     # total_amount = np.shape(Mu)[0]*np.shape(Mu)[1] # number of u_ij, v_ij
#     # flatten_uv = np.zeros(total_amount*2)
#     # flatten_uv[:total_amount] = Mu.flatten()
#     # flatten_uv[total_amount:] = Mv.flatten()
#     # return flatten_uv

def flatten(us):
    return np.concatenate(tuple(u.flatten() for u in us))

def speed_filter(M, axis, num_axis):
    result = M.copy()
    target_slices = tuple(slice(None) if i!=axis else [0,-1] for i in range(num_axis))
    result[target_slices] = 0
    return result

def flattened_idx(idx, sz):
    result = 0
    curr_block_sz = 1
    for i, s in zip(idx[::-1], sz[::-1]):
        result += i * curr_block_sz
        curr_block_sz *= s
    return result

def unflattened_idx(flattened_idx, sz):
    result = []
    for s in sz[::-1]:
        flattened_idx, curr = divmod(flattened_idx, s)
        result.append(curr)
    return result[::-1]

def construct_B(dg, grid_res, num_axis):
    del_x = 1 / dg
    amount = (grid_res + 1) * (grid_res ** (num_axis - 1))
    rows = []
    cols = []
    datas = []
    for i in range(grid_res**num_axis):
        curr_idx = unflattened_idx(i, (grid_res,) * num_axis)
        for j in range(num_axis):
            curr_grid_shape = [grid_res + 1 if k==j else grid_res for k in range(num_axis)]
            rows.append(i)
            datas.append(-del_x)
            curr_flattend_idx = flattened_idx(curr_idx, curr_grid_shape)
            cols.append(curr_flattend_idx + j * amount)
            curr_idx[j] += 1
            rows.append(i)
            datas.append(del_x)
            curr_flattend_idx = flattened_idx(curr_idx, curr_grid_shape)
            cols.append(curr_flattend_idx + j * amount)
            curr_idx[j] -= 1
    return sparse.coo_matrix((datas, (rows, cols)), shape = (grid_res**num_axis, amount * num_axis))
            

def construct_D(dg, grid_res, num_axis):
    del_x = 1 / dg
    amount = (grid_res + 1) * (grid_res ** (num_axis - 1))
    grid_shape = (grid_res,) * num_axis
    rows = []
    cols = []
    datas = []
    for ax in range(num_axis):
        curr_grid_shape = [grid_res + 1 if a==ax else grid_res for a in range(num_axis)]
        start_idx = ax * amount
        for i in range(amount):
            curr_idx = unflattened_idx(i, curr_grid_shape)
            if curr_idx[ax] == 0:
                continue
            rows.append(start_idx + i)
            cols.append(flattened_idx(curr_idx, grid_shape))
            datas.append(del_x)
            curr_idx[ax] -= 1
            rows.append(start_idx + i)
            cols.append(flattened_idx(curr_idx, grid_shape))
            datas.append(-del_x)
    return sparse.coo_matrix((datas, (rows, cols)), shape = (amount * num_axis, grid_res**num_axis))

def conj_grad(A, x, b):
    r = b - np.dot(A, x)
    p = r
    double_r = (r**2).sum()
    
    for _ in range(len(b)):
        Ap = np.dot(A, p)
        alpha = double_r / np.dot(np.transpose(p), Ap)
        x += np.dot(alpha, p)
        r -= np.dot(alpha, Ap)
        rsnew = (r**2).sum()
        if np.sqrt(rsnew) < 1e-10:
            break
        beta  = rsnew / double_r
        p = r + beta*p
        double_r = rsnew
    return x

def avection(q, q_dot, d_t):
    q += q_dot * d_t

def external_force(q_dot, g, d_t):
    q_dot += g * d_t

def linear_interp_weight(grid_coord):
    u = grid_coord[:,0][:,None]
    v = grid_coord[:,1][:,None]

    u = np.concatenate((1 - u, u), axis=1)
    v = np.concatenate((1 - v, v), axis=1)

    w = np.einsum("ij,ik->ijk", u, v)
    return w

def project_to_staggered_grid_ax(relative_coord, qdot, axis, num_axis, grid_res):
    relative_coord = relative_coord + 0.5 # set origin to grid location (-0.5, -0.5)
    relative_coord[:, axis] -= 0.5 # set the origin of aligned axis back to 0
    grid_coord = np.fmod(relative_coord, 1) # get the relative coordinate of points inside cells
    grid_loc = np.int32(relative_coord) # get the index of the cell
    grid_weights = linear_interp_weight(grid_coord)
    shape = np.array((grid_res + 1,) * num_axis)
    shape[axis] -= 1 # set up the temporary result holder for each block in the grid
    curr_blocks = np.zeros(tuple(shape) + (2,) * num_axis)
    curr_counter = np.zeros(shape)
    shape += 1
    result_holder = np.zeros(shape) # used for storing flattened results
    # fill projected values into the block result holder
    for loc, w, u in zip(grid_loc, grid_weights, qdot[:, axis]):
        loc = tuple(loc)
        curr_blocks[loc] += w * u
        curr_counter[loc] += 1
    # fill the results in the block holder into the grid
    for idx in zip(*np.where(curr_counter > 0)):
        idxx = tuple(slice(idx[i], idx[i]+2) for i in range(2))
        result_holder[idxx] += curr_blocks[idx] / curr_counter[idx]
    # remove the out of bound grid points out 
    final_slices = tuple(slice(0 if i==axis else 1, grid_res + 1) for i in range(num_axis))
    return result_holder[final_slices]

def project_to_staggered_grid(q, qdot, num_axis, pix_size, grid_res):
    relative_coord = q / pix_size
    result = []
    for i in range(num_axis):
        result.append(project_to_staggered_grid_ax(relative_coord, qdot, i, num_axis, grid_res))
    return result

def get_speed_updater(D, p, dt, density, num_axis, grid_res, pix_size):
    """
    Should return a list of staggered grid point values
    """
    # D = construct_assem_D(pix_size, pix_size, grid_res)
    dp = np.matmul(D, p)
    ax_dp = dp.reshape(num_axis,-1)
    result = []
    for ax, adp in enumerate(ax_dp):
        adp = adp.reshape(tuple(grid_res + 1 if i==ax else grid_res for i in range(num_axis)))
        result.append(-adp * dt / density)
    return result

def unproject_from_staggered_grid_ax(relative_coord, u, axis, num_axis):
    relative_coord = relative_coord - 0.5 # set origin to grid location (0.5, 0.5)
    relative_coord[:, axis] += 0.5 # set the origin of aligned axis back to 0
    grid_coord = np.fmod(relative_coord + 1, 1) # get the relative coordinate of points inside cells
    grid_loc = np.int32(np.floor(relative_coord)) # get the index of the cell
    grid_weights = linear_interp_weight(grid_coord)
    result = np.empty(grid_loc.shape[0])
    for i in range(grid_coord.shape[0]):
        curr_loc = grid_loc[i]
        grid_slice = []
        weight_slice = []
        for j in range(num_axis):
            if curr_loc[j] < 0:
                grid_slice.append(0)
                weight_slice.append(1)
            elif curr_loc[j] >= u.shape[j] - 1:
                grid_slice.append(curr_loc[j])
                weight_slice.append(0)
            else:
                grid_slice.append(slice(curr_loc[j], curr_loc[j] + 2))
                weight_slice.append(slice(None))
        weight = grid_weights[i][weight_slice]
        value = u[grid_slice]
        result[i] = (weight * value).sum()
    return result        

def unproject_from_staggered_grid(q, us, num_axis, pix_size):
    relative_coord = q / pix_size
    result = np.empty(q.shape)
    for i in range(num_axis):
        result[:,i] = unproject_from_staggered_grid_ax(relative_coord, us[i], i, num_axis)
    return result


# num_particles = 5000
# box_size = 2
# grid_res = 2
# pix_size = box_size / grid_res

# q = np.random.rand((num_particles, 2))
# q_dot = np.zeros((num_particles, 2))

# def construct_assem_B(dx, dy, size):
#     #construct the assemble B martix with grid row and grid col and dx dy
#     del_x = 1/dx
#     del_y = 1/dy
#     amount = (size + 1)*size
#     rows = []
#     cols = []
#     datas = []
#     for i in range(size*size):
#         rows.append(i)
#         cols.append(i)
#         datas.append(-del_x)
#         rows.append(i)
#         cols.append(i + size)
#         datas.append(del_x)
#         rows.append(i)
#         cols.append(i + amount + i // size)
#         datas.append(-del_y)
#         rows.append(i)
#         cols.append(i + amount + i // size + 1)
#         datas.append(del_y)
#     B = sparse.coo_matrix((datas, (rows, cols)), shape = (size*size, amount * 2)).toarray()
#     return B

# def construct_assem_D(dx, dy, size):
#     #construct the assemble D matrix with grid size and dx dy
#     del_x = 1/dx
#     del_y = 1/dy
#     amount = (size+1)*size
#     rows = []
#     cols = []
#     datas = []
#     #construct u
#     for i in range(amount):
#         # if i < size:
#         #     rows.append(i)
#         #     cols.append(i)
#         #     datas.append(del_x)
#         # elif i >= amount - size:
#         #     rows.append(i)
#         #     cols.append(i-size)
#         #     datas.append(-del_x)
#         #     #become zero!
#         #else:
#         if i >= size and i < amount - size:
#             rows.append(i)
#             cols.append(i-size)
#             datas.append(-del_x)
#             rows.append(i)
#             cols.append(i)
#             datas.append(del_x)
#     #construct v
#     for i in range(size):
#         corner_x = amount+i*(size+1)
#         corner_y = size*i
#         # rows.append(corner_x)
#         # cols.append(corner_y)
#         # datas.append(del_y)
#         #become 0?
#         for j in range (size - 1):
#             rows.append(corner_x + 1 + j)
#             cols.append(corner_y + j)
#             datas.append(-del_y)
#             rows.append(corner_x + 1 + j)
#             cols.append(corner_y + j + 1)
#             datas.append(del_y)
#         # rows.append(corner_x + size)
#         # cols.append(corner_y + size -1)
#         # datas.append(-del_y)
#         # become zero
                     
#     D = sparse.coo_matrix((datas, (rows, cols)), shape = (amount *2, size*size)).toarray()
#     return D

# def project_qdot_to_grid(q, q_dot):
#     num_axis = 2
#     relative_coord = q / pix_size + 0.5 # set origin to grid location (-0.5, -0.5)
#     result = []
#     for ax in range(num_axis):
#         ax_relative = relative_coord.copy()
#         ax_relative[:, ax] -= 0.5
#         ax_grid_coord = np.fmod(ax_relative, 1)
#         ax_grid_loc = np.int32(ax_relative)
#         ax_grid_weights = linear_interp_weight(ax_grid_coord)
#         shape = np.array((grid_res + 1,) * num_axis)
#         shape[ax] -= 1
#         curr_blocks = np.zeros(tuple(shape) + (2,) * num_axis)
#         curr_counter = np.zeros(shape)
#         shape += 1
#         result_holder = np.zeros(shape)
#         for loc, w, u in zip(ax_grid_loc, ax_grid_weights, q_dot[:,ax]):
#             loc = tuple(loc)
#             curr_blocks[loc] += w * u
#             curr_counter[loc] += 1
#         for idx in zip(*np.where(curr_counter > 0)):
#             idxx = tuple(slice(idx[i], idx[i]+2) for i in range(num_axis))
#             result_holder[idxx] += curr_blocks[idx] / curr_counter[idx] # TODO: Generalize this line of code
#         starts = [1] * num_axis
#         starts[ax] -= 1
#         result.append(result_holder[starts[0]:grid_res+1, starts[1]:grid_res+1]) # TODO: Generalize this line of code
#     return result
# q = np.array([[0.3,0.4],[1.7,1.3]])
# q_dot = np.ones(q.shape)
# print(project_qdot_to_grid(q, q_dot))


    # u_relative = relative_coord
    # u_relative[:,0] -= 0.5
    # v_relative = relative_coord
    # v_relative[:,1] -= 0.5

    # u_grid_coord = np.fmod(u_relative, 1)
    # u_grid_loc = np.int32(u_relative)
    # v_grid_coord = np.fmod(v_relative, 1)
    # v_grid_loc = np.int32(v_relative)

    # u_grid_weights = linear_interp_weight(u_grid_coord)
    # v_grid_weigths = linear_interp_weight(v_grid_coord)

    # u = np.zeros((grid_res + 2, grid_res + 1))
    # v = np.zerso((grid_res + 1, grid_res + 2))

    # for 


