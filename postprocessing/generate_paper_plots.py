import subprocess
import os
import shutil
import re


# This is the output directory where the final plot files are stored
OUTPUT_DIR = '/home/carstensen/results_refined_grid/plots'

# The names of the diectories are relevant, since these names are used as names
# in the legend in the plots
# These are the directories where the data is read in from
READ_DIR_HI_RES = '/home/carstensen/results_refined_grid/hi-res'
READ_DIR_NO_ADAPT = '/home/carstensen/results_refined_grid/dycore'
READ_DIR_ADAPT = '/home/carstensen/results_refined_grid/grid adaptation'

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
                output_dir_tmp_grid_comp], cwd='pyplotgen')
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
            OUTPUT_DIR + '/arm_grid_adapt.png')
shutil.copy(output_dir_tmp_each_ref_crit + '/astex_a209/each_ref_crit_term_ASTEX.png',
            OUTPUT_DIR + '/astex_grid_adapt.png')
shutil.copy(output_dir_tmp_each_ref_crit + '/gabls2/each_ref_crit_term_GABLS2.png',
            OUTPUT_DIR + '/gabls2_grid_adapt.png')

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