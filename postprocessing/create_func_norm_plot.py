import matplotlib.pyplot as plt
import random
import numpy as np

def generate_random_piecewise_lin_func(n):
    piece_lin_func_x = range(0, n)
    piece_lin_func_y = []
    for _ in range(n):
        piece_lin_func_y.append(random.randrange(start=0, stop=30))
    return piece_lin_func_x, piece_lin_func_y

def create_grid_new(num_levels, g_x, g_y, grid_sfc, grid_top):
    # area for this function between sfc and top should be num_levels-1
    grid_heights = [] # size is num_levels
    grid_heights.append(grid_sfc) # could also use g_x[0]?
    prev_x_ind = 0
    area_up_to_prev_x = 0
    #for i in range(1, num_levels-1):
    for i in range(1, num_levels-1):
        grid_level_found = False
        while not grid_level_found:
            new_x_area_increment = (g_y[prev_x_ind] + g_y[prev_x_ind+1])/2*(g_x[prev_x_ind+1]-g_x[prev_x_ind])
            desired_area_up_ith_level = i
            if (area_up_to_prev_x + new_x_area_increment >= desired_area_up_ith_level):
                # grid level lies between g_x[prev_x_ind] and g_x[prev_x_ind+1]
                # TODO check what happens if new grid level is exactly g_x[prev_x_ind+1]
                A = desired_area_up_ith_level - area_up_to_prev_x # can be bigger than 1 if more than one grid level lies between the new level and the prev_x level, but A must always be nonnegative
                prev_ind = prev_x_ind
                if abs(g_y[prev_ind+1]-g_y[prev_ind]) < 1e-6:
                    grid_level = A/g_y[prev_ind] + g_x[prev_ind]
                else:
                    m = (g_y[prev_ind+1] - g_y[prev_ind])/(g_x[prev_ind+1] - g_x[prev_ind])
                    b = g_y[prev_ind] - g_x[prev_ind]*m
                    p = 2*b/m
                    q = -1*g_x[prev_ind]**2-2/m*b*g_x[prev_ind]-2/m*A
                    grid_level_1 = -p/2 + np.sqrt((p/2)**2-q)
                    if grid_level_1 > g_x[prev_ind] and grid_level_1 < g_x[prev_ind+1] and grid_level_1 > grid_heights[i-1]:
                        grid_level = grid_level_1
                    else:
                        grid_level_2 = -p/2 - np.sqrt((p/2)**2-q)
                        grid_level = grid_level_2
                        if grid_level_2 < g_x[prev_ind] or grid_level_2 > g_x[prev_ind+1] or grid_level_2 < grid_heights[i-1]:
                            print('something went wrong')
                            # in cases like this, just fall back to equidistant, or leave grid as it is before?
                grid_heights.append(grid_level)
                grid_level_found = True
            else:
                area_up_to_prev_x += new_x_area_increment
                prev_x_ind += 1
    
    grid_heights.append(grid_top) # could also use g_x[-1]?
    return grid_heights

def calc_integral(g_x, g_y):
    sum = 0
    for i in range(len(g_x)-1):
        sum += (g_y[i] + g_y[i+1])/2 * (g_x[i+1] - g_x[i])
    return sum

def get_final_density(g_x, g_y, min_dens, grid_sfc, grid_top, num_levels):
    # grid_sfc and grid_top should be the last and first element in g_x
    h = grid_top - grid_sfc
    integral = calc_integral(g_x, g_y)
    # if integral is zero all g_y were the same and smaller than the min_dens, so after + delta_s - min_dens they are all zero and there is no factor we can apply to change the integral, so in that case just set the g_y to equi_dens
    if integral < 1e-6:
        # set g_y to equi_dens and return
        print('TODO')
    factor = (num_levels-1-h*min_dens)/integral
    for i in range(len(g_y)):
        g_y[i] = factor*(g_y[i]) + min_dens
    return g_x, g_y

def get_factorized_function(g_x, g_y, g_norm_min_x, g_norm_min_y, grid_sfc, grid_top, num_levels):
    # grid_sfc and grid_top should be the last and first element in g_x
    x_new = g_x.copy()
    y_new = g_y.copy()
    x_min_new = g_norm_min_x.copy()
    y_min_new = g_norm_min_y.copy()
    h = grid_top - grid_sfc

    # shift initial grid density function down to zero
    for i in range(len(x_new)):
        y_new[i] -= min(y_new)

    # calculate integrals
    integral = calc_integral(x_new, y_new)
    integral_min = calc_integral(x_min_new, y_min_new)

    # calculate factor
    factor = ((num_levels-1)-integral_min)/integral

    # handle case if integral is zero (if g(z)=0)
    if integral < 1e-6:
        print('TODO')

    # apply factor
    for i in range(len(y_new)):
        y_new[i] = factor*(y_new[i])
    return x_new, y_new

def get_shifted_function(g_x, g_y, y_norm_min_dens):
    x_new = g_x.copy()
    y_new = g_y.copy()
    for i in range(len(y_new)):
        y_new[i] += y_norm_min_dens[i]
    return x_new, y_new

def get_normalized_min_dens(lam, num_levels, x, y):
    integral_min_dens = calc_integral(x,y)
    x_new = x.copy()
    y_new = y.copy()
    for i in range(len(x)):
        y_new[i] *= lam*(num_levels-1)/integral_min_dens
    return x_new, y_new


def get_piecewise_lin_func():
    x = [0, 1, 3, 3.5, 4.5, 5, 5.25, 5.5, 6, 8]
    y = [0, 0.3, 0.5, 1.9, 3.0, 1.6, 0.5, 0.2, 0.1, 0.05]
    return x, y

def get_piecewise_lin_func_min_dens(x):
    y = [3, 1.5, 0.8, 0.7, 0.5, 0.4, 0.3, 0.2, 0.1, 0.1]
    y_final = []
    for i in range(len(x)):
        y_final.append(y[i])
    return x, y

def create_plot(lam, file_out):
    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(1, 6)
    fig.set_figheight(4)
    fig.set_figwidth(9)
    #fig.suptitle('Horizontally stacked subplots')

    ax1.set_title('1. $g_\mathrm{min}^{\mathrm{(nn)}}(z)$')
    ax2.set_title('2. $g_\mathrm{min}(z)$')
    ax3.set_title('3. $g^{\mathrm{(nn)}}(z)$')
    ax4.set_title('4. $\\alpha g^{\mathrm{(nn)}}(z)$')
    ax5.set_title('5.$g(z)$')
    ax6.set_title('Create the grid')

    ax1.set_ylabel('$z$ [km]')
    
    ax2.get_yaxis().set_visible(False)
    ax3.get_yaxis().set_visible(False)
    ax4.get_yaxis().set_visible(False)
    ax5.get_yaxis().set_visible(False)
    ax6.get_yaxis().set_visible(False)

    num_levels = 9
    grid_sfc = 0
    grid_top = 8
    h = grid_top - grid_sfc
    min_dens_integral_dens = lam*(num_levels-1)/(grid_top - grid_sfc)
    #print(min_dens_integral_dens)
    equi_dens = (num_levels - 1)/(h)

    col0 = 'plum'
    col1 = 'teal'
    col2 = 'steelblue'
    col3 = 'darkorange'
    col4 = 'darkgreen'
    grid_col = 'darkred'

    x_init, y_init = get_piecewise_lin_func()
    x_init_min_dens, y_init_min_dens = get_piecewise_lin_func_min_dens(x_init)
    x_norm_min_dens, y_norm_min_dens = get_normalized_min_dens(lam, num_levels, x_init_min_dens, y_init_min_dens)
    x_fact, y_fact = get_factorized_function(x_init, y_init, x_norm_min_dens, y_norm_min_dens, grid_sfc, grid_top, num_levels)
    x_fact_shift, y_fact_shift = get_shifted_function(x_fact, y_fact, y_norm_min_dens)
    grid_heights = create_grid_new(num_levels, x_fact_shift, y_fact_shift, grid_sfc, grid_top)

    max_x_axis = max(*y_init, *y_init_min_dens, *y_norm_min_dens, *y_fact, *y_fact_shift)
    min_x_axis = min(*y_init, *y_init_min_dens, *y_norm_min_dens, *y_fact, *y_fact_shift)
    max_y_axis = max(*x_init, *x_init_min_dens, *x_norm_min_dens, *x_fact, *x_fact_shift)
    min_y_axis = min(*x_init, *x_init_min_dens, *x_norm_min_dens, *x_fact, *x_fact_shift)

    ax1.set_xlim([min_x_axis, max_x_axis])
    ax2.set_xlim([min_x_axis, max_x_axis])
    ax3.set_xlim([min_x_axis, max_x_axis])
    ax4.set_xlim([min_x_axis, max_x_axis])
    ax5.set_xlim([min_x_axis, max_x_axis])
    ax6.set_xlim([min_x_axis, max_x_axis])

    ax1.set_ylim([min_y_axis, max_y_axis])
    ax2.set_ylim([min_y_axis, max_y_axis])
    ax3.set_ylim([min_y_axis, max_y_axis])
    ax4.set_ylim([min_y_axis, max_y_axis])
    ax5.set_ylim([min_y_axis, max_y_axis])
    ax6.set_ylim([min_y_axis, max_y_axis])

    ax1.vlines(min_dens_integral_dens, grid_sfc, grid_top, colors='black', linestyles='solid', label="$\\frac{\lambda (n-1)}{z_\mathrm{top}-z_\mathrm{sfc}}$", zorder=5)
    ax2.vlines(min_dens_integral_dens, grid_sfc, grid_top, colors='black', linestyles='solid', label="$\\frac{\lambda (n-1)}{z_\mathrm{top}-z_\mathrm{sfc}}$", zorder=5)
    ax3.vlines(min_dens_integral_dens, grid_sfc, grid_top, colors='black', linestyles='solid', label="$\\frac{\lambda (n-1)}{z_\mathrm{top}-z_\mathrm{sfc}}$", zorder=5)
    ax4.vlines(min_dens_integral_dens, grid_sfc, grid_top, colors='black', linestyles='solid', label="$\\frac{\lambda (n-1)}{z_\mathrm{top}-z_\mathrm{sfc}}$", zorder=5)
    ax5.vlines(min_dens_integral_dens, grid_sfc, grid_top, colors='black', linestyles='solid', label="$\\frac{\lambda (n-1)}{z_\mathrm{top}-z_\mathrm{sfc}}$", zorder=5)
    ax6.vlines(min_dens_integral_dens, grid_sfc, grid_top, colors='black', linestyles='solid', label="$\\frac{\lambda (n-1)}{z_\mathrm{top}-z_\mathrm{sfc}}$", zorder=5)

    ax1.plot(y_init_min_dens, x_init_min_dens, color=col0, linestyle='solid')
    ax1.fill_betweenx(x_init_min_dens, y_init_min_dens, color=col0, alpha=0.5)
    ax2.plot(y_init_min_dens, x_init_min_dens, color=col0, linestyle='solid', alpha=0.5)
    ax2.fill_betweenx(x_init_min_dens, y_init_min_dens, color=col0, alpha=0.1)
    #ax2.plot(y_norm_min_dens, x_norm_min_dens, color=col1)
    #ax2.fill_betweenx(x_norm_min_dens, y_norm_min_dens, color=col1, alpha=0.5)

    ax2.plot(y_norm_min_dens, x_norm_min_dens, color='black', linestyle='dashed', label="$g_\mathrm{min}(z)$", zorder=5)
    ax3.plot(y_norm_min_dens, x_norm_min_dens, color='black', linestyle='dashed', label="$g_\mathrm{min}(z)$", zorder=5)
    ax4.plot(y_norm_min_dens, x_norm_min_dens, color='black', linestyle='dashed', label="$g_\mathrm{min}(z)$", zorder=5)
    ax5.plot(y_norm_min_dens, x_norm_min_dens, color='black', linestyle='dashed', label="$g_\mathrm{min}(z)$", zorder=5)
    ax6.plot(y_norm_min_dens, x_norm_min_dens, color='black', linestyle='dashed', label="$g_\mathrm{min}(z)$", zorder=5)

    ax3.plot(y_init,x_init, color=col2)
    ax3.fill_betweenx(x_init, y_init, color=col2, alpha=0.5)

    ax4.plot(y_init, x_init, alpha=0.5, color=col2, linestyle='solid')
    ax4.fill_betweenx(x_init, y_init, alpha=0.1, color=col2, linestyle='solid')
    ax4.plot(y_fact, x_fact, color=col3)
    ax4.fill_betweenx(x_fact, y_fact, alpha=0.5, color=col3)

    ax5.plot(y_fact, x_fact, color=col3, alpha=0.5, linestyle='solid')
    ax5.fill_betweenx(x_fact, y_fact, alpha=0.1, color=col3, linestyle='solid')
    ax5.plot(y_fact_shift, x_fact_shift, color=col4)
    ax5.fill_betweenx(x_fact_shift, y_fact_shift, alpha=0.5, color=col4)
    
    ax6.plot(y_fact_shift, x_fact_shift, color=col4, alpha=0.5, linestyle='solid')
    ax6.hlines(grid_heights, min_x_axis, max_x_axis, color=grid_col, linewidths=3.0, zorder=6)
    ax6.fill_betweenx(x_fact_shift, y_fact_shift, alpha=0.5, color=col4)

    ax3.legend()
    
    fig.tight_layout()
    #fig.savefig(f'func_norm_lam_{lam}_new.svg')
    plt.savefig(f'func_norm_lam_{lam}_new.png', dpi=300)
    #plt.show()

def main():
    lam_one_half = 0.5
    file_out_one_half = 'func_norm.png'
    create_plot(lam_one_half, file_out_one_half)
    
    lam_one = 1.0
    file_out_one = 'func_norm_lam_1.png'
    create_plot(lam_one, file_out_one)

    lam_zero = 1e-20
    file_out_zero = 'func_norm_lam_0.png'
    create_plot(lam_zero, file_out_zero)

if __name__ == '__main__':
    main()