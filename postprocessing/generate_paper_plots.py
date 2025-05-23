import subprocess
import os
import shutil
import re

# To get the old plots back, change the VariableGroupds back in the CaseDefintions file
# for each of the three cases

# This is the output directory where the final plot files are stored
OUTPUT_DIR = '/home/carstensen/results_refined_grid_v2/paper_plots'

# The names of the diectories are relevant, since these names are used as names
# in the legend in the plots
# These are the directories where the data is read in from
READ_DIR_HI_RES = '/home/carstensen/results_refined_grid_v2/hi-res'
READ_DIR_NO_ADAPT = '/home/carstensen/results_refined_grid_v2/dycore'
READ_DIR_ADAPT = '/home/carstensen/results_refined_grid_v2/grid adaptation'

#================================================================================================

hi_res_name = READ_DIR_HI_RES.split('/')[-1]
no_adapt_name = READ_DIR_NO_ADAPT.split('/')[-1]
adapt_name = READ_DIR_ADAPT.split('/')[-1]

def clean_path(input_string):
    """
    Removes characters from a path string
    that are not valid for the operating system

    :param input_string: Path string to have characters removed
    :return: A cleaned version of the path
    """
    # string = string.replace('.', '')
    input_string = input_string.replace(',', '')
    input_string = input_string.replace(' ', '_')
    input_string = input_string.replace('*', 'x')
    if 'win' in os.name.lower() or 'nt' in os.name.lower():
        # Replace all ':' except the one specifying the drive
        # (Uses look ahead and look back patterns)
        input_string = re.sub(r'(?<![A-Z]):(?!\\)', '-', input_string)
    else:
        input_string = input_string.replace(':', '-')
    return input_string

def clean_title(string):
    """
    Cleans plot titles so they can be used as part of a file name.
    Matplotlib plot titles may contain latex code enclosed in '$' characters.
    This code should be removed entirely to avoid having invalid characters in file names.

    :param string: Title of a plot which is used to generate a file name
    :return: A cleaned version of the input string
    """
    return re.sub("\$[^$]*\$",'',string)

def removeInvalidFilenameChars(filename):
    """
    Removes characters from a string that are not valid for a filename
    :param filename: Filename string to have characters removed
    :return: A character stripped version of the filename
    """
    filename = filename.replace('/', '')
    filename = clean_path(filename)
    filename = clean_title(filename)
    return filename

hi_res_name = removeInvalidFilenameChars(hi_res_name)
no_adapt_name = removeInvalidFilenameChars(no_adapt_name)
adapt_name = removeInvalidFilenameChars(adapt_name)

# To add variables, add them to VariableGroupPaperPlots.py

# This is only the output directory where the plots are put temporarily
output_dir_tmp = '/tmp/delete_me_plot_generation'
os.mkdir(output_dir_tmp)
os.mkdir(OUTPUT_DIR)

# To get the dycore grid for the grid comparison plot the script automatically reads in the
# dycore.grd file
output_dir_tmp_grid_comp = output_dir_tmp + '/grid_comp'
subprocess.run(['python', 'pyplotgen.py', '--grid-comparison-plot', '-o',
                output_dir_tmp_grid_comp, '-c', READ_DIR_NO_ADAPT], cwd='pyplotgen')
shutil.copy(output_dir_tmp_grid_comp + '/grid_comp.png', OUTPUT_DIR + '/grid_comp.png')

# These are the grid adaptation plots
output_dir_tmp_grid_adapt = output_dir_tmp + '/grid_adapt'
subprocess.run(['python', 'pyplotgen.py', '--grid-adapt-plot', '-o', output_dir_tmp_grid_adapt,
                '-c', READ_DIR_ADAPT], cwd='pyplotgen')
shutil.copy(output_dir_tmp_grid_adapt + '/arm_grid_adapt.png', OUTPUT_DIR + '/arm_grid_adapt.png')
shutil.copy(output_dir_tmp_grid_adapt + '/astex_a209_grid_adapt.png',
            OUTPUT_DIR + '/astex_a209_grid_adapt.png')
shutil.copy(output_dir_tmp_grid_adapt + '/gabls2_grid_adapt.png',
            OUTPUT_DIR + '/gabls2_grid_adapt.png')

# These are the plots that show each refinement criterion term
output_dir_tmp_each_ref_crit = output_dir_tmp + '/each_ref_crit'
subprocess.run(['python', 'pyplotgen.py', '-o', output_dir_tmp_each_ref_crit, '-c',
                READ_DIR_HI_RES, '--hq'], cwd='pyplotgen')
shutil.copy(output_dir_tmp_each_ref_crit + '/arm/each_ref_crit_term_ARM.png',
            OUTPUT_DIR + '/arm_each_ref_crit_term.png')
shutil.copy(output_dir_tmp_each_ref_crit + '/astex_a209/each_ref_crit_term_ASTEX.png',
            OUTPUT_DIR + '/astex_each_ref_crit_term.png')
shutil.copy(output_dir_tmp_each_ref_crit + '/gabls2/each_ref_crit_term_GABLS2.png',
            OUTPUT_DIR + '/gabls2_each_ref_crit_term.png')

# These are the plots that show the normalized (minimum) grid density
output_dir_tmp_norm_grid_dens = output_dir_tmp + '/norm_grid_dens'
subprocess.run(['python', 'pyplotgen.py', '-o', output_dir_tmp_norm_grid_dens, '-c',
                READ_DIR_NO_ADAPT, '--hq'], cwd='pyplotgen')
shutil.copy(output_dir_tmp_norm_grid_dens + '/arm/norm_grid_dens_ARM.png',
            OUTPUT_DIR + '/arm_norm_grid_dens.png')
shutil.copy(output_dir_tmp_norm_grid_dens + '/astex_a209/norm_grid_dens_ASTEX.png',
            OUTPUT_DIR + '/astex_norm_grid_dens.png')
shutil.copy(output_dir_tmp_norm_grid_dens + '/gabls2/norm_grid_dens_GABLS2.png',
            OUTPUT_DIR + '/gabls2_norm_grid_dens.png')

# To make changes to the reduced height, adjust the height in the cases in Case_definitions_red_height.py
# These are the profiles with reduced height (rcm profile)
output_dir_tmp_profiles_red_height = output_dir_tmp + '/profiles_red_height'
subprocess.run(['python', 'pyplotgen.py', '-o', output_dir_tmp_profiles_red_height, '-c',
                READ_DIR_HI_RES, READ_DIR_NO_ADAPT, READ_DIR_ADAPT,
                '--hq', '--red-height'], cwd='pyplotgen')
shutil.copy(output_dir_tmp_profiles_red_height + '/arm/rcm_ARM.png',
            OUTPUT_DIR + '/arm_rcm_profile.png')
shutil.copy(output_dir_tmp_profiles_red_height + '/astex_a209/rcm_ASTEX.png',
            OUTPUT_DIR + '/astex_rcm_profile.png')

# These are the profiles with full height (rtm profile)
output_dir_tmp_profiles_full_height = output_dir_tmp + '/profiles_full_height'
subprocess.run(['python', 'pyplotgen.py', '-o', output_dir_tmp_profiles_full_height, '-c',
                READ_DIR_HI_RES, READ_DIR_NO_ADAPT, READ_DIR_ADAPT, '--hq'], cwd='pyplotgen')
shutil.copy(output_dir_tmp_profiles_full_height + '/arm/rtm_ARM.png',
            OUTPUT_DIR + '/arm_rtm_profile.png')
shutil.copy(output_dir_tmp_profiles_full_height + '/astex_a209/rtm_ASTEX.png',
            OUTPUT_DIR + '/astex_rtm_profile.png')
shutil.copy(output_dir_tmp_profiles_full_height + '/gabls2/rtm_GABLS2.png',
            OUTPUT_DIR + '/gabls2_rtm_profile.png')

# To adjust the minimum and maximum value of the color scales, adjust the constants in VariableGroupPaperPlots (WP2_VAR_MIN_ASTEX etc.)

# These are the hovmoeller plots
output_dir_tmp_hovmoeller = output_dir_tmp + '/hovmoeller'
subprocess.run(['python', 'pyplotgen.py', '-o', output_dir_tmp_hovmoeller, '-c', READ_DIR_HI_RES, READ_DIR_NO_ADAPT, READ_DIR_ADAPT, '--hq', '--red-height', '-t'], cwd='pyplotgen')
shutil.copy(output_dir_tmp_hovmoeller + f'/arm/wpthlp_ARM_{hi_res_name}.png', OUTPUT_DIR + '/arm_wpthlp_hovmoeller_hi_res.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/arm/wpthlp_ARM_{no_adapt_name}.png', OUTPUT_DIR + '/arm_wpthlp_hovmoeller_no_adapt.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/arm/wpthlp_ARM_{adapt_name}.png', OUTPUT_DIR + '/arm_wpthlp_hovmoeller_adapt.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/astex_a209/wpthlp_ASTEX_{hi_res_name}.png', OUTPUT_DIR + '/astex_wpthlp_hovmoeller_hi_res.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/astex_a209/wpthlp_ASTEX_{no_adapt_name}.png', OUTPUT_DIR + '/astex_wpthlp_hovmoeller_no_adapt.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/astex_a209/wpthlp_ASTEX_{adapt_name}.png', OUTPUT_DIR + '/astex_wpthlp_hovmoeller_adapt.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/gabls2/wpthlp_GABLS2_{hi_res_name}.png', OUTPUT_DIR + '/gabls2_wpthlp_hovmoeller_hi_res.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/gabls2/wpthlp_GABLS2_{no_adapt_name}.png', OUTPUT_DIR + '/gabls2_wpthlp_hovmoeller_no_adapt.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/gabls2/wpthlp_GABLS2_{adapt_name}.png', OUTPUT_DIR + '/gabls2_wpthlp_hovmoeller_adapt.png')

shutil.copy(output_dir_tmp_hovmoeller + f'/arm/cloud_frac_ARM_{hi_res_name}.png', OUTPUT_DIR + '/arm_cloud_frac_hovmoeller_hi_res.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/arm/cloud_frac_ARM_{no_adapt_name}.png', OUTPUT_DIR + '/arm_cloud_frac_hovmoeller_no_adapt.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/arm/cloud_frac_ARM_{adapt_name}.png', OUTPUT_DIR + '/arm_cloud_frac_hovmoeller_adapt.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/astex_a209/cloud_frac_ASTEX_{hi_res_name}.png', OUTPUT_DIR + '/astex_cloud_frac_hovmoeller_hi_res.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/astex_a209/cloud_frac_ASTEX_{no_adapt_name}.png', OUTPUT_DIR + '/astex_cloud_frac_hovmoeller_no_adapt.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/astex_a209/cloud_frac_ASTEX_{adapt_name}.png', OUTPUT_DIR + '/astex_cloud_frac_hovmoeller_adapt.png')

shutil.copy(output_dir_tmp_hovmoeller + f'/arm/wp2_ARM_{hi_res_name}.png', OUTPUT_DIR + '/arm_wp2_hovmoeller_hi_res.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/arm/wp2_ARM_{no_adapt_name}.png', OUTPUT_DIR + '/arm_wp2_hovmoeller_no_adapt.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/arm/wp2_ARM_{adapt_name}.png', OUTPUT_DIR + '/arm_wp2_hovmoeller_adapt.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/astex_a209/wp2_ASTEX_{hi_res_name}.png', OUTPUT_DIR + '/astex_wp2_hovmoeller_hi_res.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/astex_a209/wp2_ASTEX_{no_adapt_name}.png', OUTPUT_DIR + '/astex_wp2_hovmoeller_no_adapt.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/astex_a209/wp2_ASTEX_{adapt_name}.png', OUTPUT_DIR + '/astex_wp2_hovmoeller_adapt.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/gabls2/wp2_GABLS2_{hi_res_name}.png', OUTPUT_DIR + '/gabls2_wp2_hovmoeller_hi_res.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/gabls2/wp2_GABLS2_{no_adapt_name}.png', OUTPUT_DIR + '/gabls2_wp2_hovmoeller_no_adapt.png')
shutil.copy(output_dir_tmp_hovmoeller + f'/gabls2/wp2_GABLS2_{adapt_name}.png', OUTPUT_DIR + '/gabls2_wp2_hovmoeller_adapt.png')

shutil.rmtree(output_dir_tmp)


#================================================================================================

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
    plt.savefig(file_out, dpi=300)
    #plt.show()

def main():
    lam_one_half = 0.5
    file_out_one_half = OUTPUT_DIR + '/func_norm.png'
    create_plot(lam_one_half, file_out_one_half)
    
    lam_one = 1.0
    file_out_one = OUTPUT_DIR + '/func_norm_lam_1.png'
    create_plot(lam_one, file_out_one)

    lam_zero = 1e-20
    file_out_zero = OUTPUT_DIR + '/func_norm_lam_0.png'
    create_plot(lam_zero, file_out_zero)

main()