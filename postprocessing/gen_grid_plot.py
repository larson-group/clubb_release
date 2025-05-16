import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

hi_res = [i for i in range(0,9000,20)]
hi_res_grid_spacings = []
hi_res_grid_spacings_heights = []
for i in range(len(hi_res)-1):
    hi_res_grid_spacings.append(hi_res[i+1]-hi_res[i])
    hi_res_grid_spacings_heights.append((hi_res[i+1]+hi_res[i])/2)

dycore = [0.0, 100.0, 325.0, 625.0, 950.0, 1300.0, 1675.0, 2065.0, 2475.0, 2915.0, 3400.0, 3890.0, 4390.0, 4915.0, 5450.0, 6010.0, 6670.0, 7490.0, 8490.0]
dycore_grid_spacings = []
dycore_grid_spacings_heights = []
for i in range(len(dycore)-1):
    dycore_grid_spacings.append(dycore[i+1]-dycore[i])
    dycore_grid_spacings_heights.append((dycore[i+1]+dycore[i])/2)

plt.scatter(hi_res_grid_spacings, hi_res_grid_spacings_heights, s=0)
for i, (x, y) in enumerate(zip(hi_res_grid_spacings, hi_res_grid_spacings_heights)):
    plt.text(x, y, str(i), fontsize=12, ha='center', va='center', color='black')

plt.scatter(dycore_grid_spacings, dycore_grid_spacings_heights, s=0)
for i, (x, y) in enumerate(zip(dycore_grid_spacings, dycore_grid_spacings_heights)):
    plt.text(x, y, str(i), fontsize=12, ha='center', va='center', color='red')

plt.title('Grid comparison')
plt.xlabel('$\Delta z \quad [\mathrm{m}]$')
plt.ylabel('$z \quad [\mathrm{m}]$')

# Create a custom legend entry
legend_patch_hi_res = mpatches.Patch(color='black', label="hi-res")
legend_patch_dycore = mpatches.Patch(color='red', label="dycore")

# Add the legend
plt.legend(handles=[legend_patch_hi_res, legend_patch_dycore])

#plt.plot(hi_res)
#plt.plot(dycore)
#plt.savefig('test.png')
plt.show()