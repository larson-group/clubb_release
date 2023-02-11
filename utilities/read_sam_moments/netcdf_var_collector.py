import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

gravity = 9.81

# Change netcdf_file path to the file you want to use.
netcdf_file = "/home/pub/les_and_clubb_benchmark_runs/sam_benchmark_runs/JULY_2017/DYCOMS_RF01_96x96x320/DYCOMS_RF01_96x96x320.nc"

# "/home/pub/les_and_clubb_benchmark_runs/sam_benchmark_runs/JULY_2017/DYCOMS_RF01_96x96x320/DYCOMS_RF01_96x96x320.nc"

# "/home/pub/les_and_clubb_benchmark_runs/sam_benchmark_runs/BOMEX_64x64x75/BOMEX_64x64x75_100m_40m_1s.nc"


def get_var_from_netcdf(nc_file, var_name, df):
    '''
    This function takes in a variable name and a netcdf file and returns a dataframe with the variable name as the column name and the values of the variable as the column values.
    '''
    nc = netCDF4.Dataset(nc_file)
    var = np.array(nc.variables[var_name][:])
    #flat = var.flatten()
    # df.to_csv('test.csv', index=False)
    #print(var.flatten())
    df[var_name] = var.flatten()
    nc.close()

def make_index(df, ncdf):
    '''
    This function creates a multi-index for the dataframe that matches how netcdf stores data. We have each of the 75 levels for each time step. 
    This function creates a list of the time step index for each level and then sets the index of the dataframe to that list.
    This index is equivalent to getting the time step index for each level in the netcdf file.
    '''
    index = []
    for i in range(len(ncdf.dimensions['time'])):
        for j in range(len(ncdf.dimensions['z'])):
            index.append(ncdf.variables['time'][i])

    df["Timestep"] = index
    #df.index = index

def make_z_col(ncdf):
    '''
    This function creates a column in the dataframe that contains the z level for each row. This is equivalent to the z level in the netcdf file.
    '''
    # z = []
    # for i in range(len(ncdf.dimensions['time'])):
    #     for j in range(len(ncdf.dimensions['z'])):
    #         z.append(ncdf.variables['z'][i])
    z =  ncdf.variables['z'][:]
    reps = len(ncdf.dimensions['time'])


    # print(len(ncdf.variables['z']))
    # print(len(ncdf.dimensions['time']))


    # print(len(ncdf.dimensions['time']))

    # file = open("test.txt", "w")
    # file.write())
    # file.close()

    z = np.tile(z, reps)

    # print(len(z))

    return np.array(z), len(ncdf.dimensions['z'])


def create_df_1040_scalar_dissipation(df, ncdf):
    '''
    This function creates the dataframe for clubb issue 1040. It will have 75 row that represent each z level. There will be 360 of these for each timestep in the clubb case.
    Each variable will take up 1 column. an index will be listed as well for reference in the output csv file.
    Predictand: Q2DIFTR/QT2.

    Predictors: W2/(U2 + V2), sqrt(g/THETAL)*dTHETAL/dz, sqrt(U2 + V2)*d(sqrt(U2 + V2))/dz
    Other Variables: CLD, WSKEW, THETAL, U, V, RTP2, U2, V2, W2, QT2, Q2DIFTR
    '''
    # create index
    make_index(df, ncdf)
    z, levels = make_z_col(ncdf)


    # Q2DIFTR = np.array(ncdf.variables['Q2DIFTR'][:]).flatten('F')
    # QT2 = np.array(ncdf.variables['QT2'][:]).flatten('F')
    # W2 = np.array(ncdf.variables['W2'][:]).flatten('F')
    # U2 = np.array(ncdf.variables['U2'][:]).flatten('F')
    # V2 = np.array(ncdf.variables['V2'][:]).flatten('F')
    # CLD = np.array(ncdf.variables['CLD'][:]).flatten('F')
    # WSKEW = np.array(ncdf.variables['WSKEW'][:]).flatten('F')
    # THETAL = np.array(ncdf.variables['THETAL'][:]).flatten('F')
    # U = np.array(ncdf.variables['U'][:]).flatten('F')
    # V = np.array(ncdf.variables['V'][:]).flatten('F')
    # print(np.array(ncdf.variables['Q2DIFTR'][:]).shape)

    Q2DIFTR = np.array(ncdf.variables['Q2DIFTR'][:]).flatten()
    QT2 = (1.e-6)*np.array(ncdf.variables['QT2'][:]).flatten()
    W2 = np.array(ncdf.variables['W2'][:]).flatten()
    U2 = np.array(ncdf.variables['U2'][:]).flatten()
    V2 = np.array(ncdf.variables['V2'][:]).flatten()
    CLD = np.array(ncdf.variables['CLD'][:]).flatten()
    WSKEW = np.array(ncdf.variables['WSKEW'][:]).flatten()
    THETAL = np.array(ncdf.variables['THETAL'][:]).flatten()
    U = np.array(ncdf.variables['U'][:]).flatten()
    V = np.array(ncdf.variables['V'][:]).flatten()
    #RTP2 = np.array(ncdf.variables['RTP2'][:]).flatten()

    Q2DIFTR_QT2 = Q2DIFTR/np.maximum((5.e-8),QT2)
    W2_U2_V2 = W2/(U2 + V2)
    
    sqrt_g_THETAL_dTHETAL_dz = np.sqrt((gravity/THETAL)*np.maximum(0,original_gradient(THETAL, z)))
    d_sqrt_U2_V2_dz = original_gradient(np.sqrt(U**2 + V**2), z)
    d_sqrt_U2_plus_V2_plus_W2_dz = original_gradient(np.sqrt(U2 + V2 + W2), z)


    df["z"] = z
    df['Q2DIFTR_on_QT2'] = Q2DIFTR_QT2 # goal value | dimesions: 1/time
    
    # dimensionless
    df['W2_on_U2_plus_V2'] = W2_U2_V2
    df['CLD'] = CLD
    df['WSKEW'] = WSKEW

    #1/time
    df['sqrt_g_on_THETAL_dTHETAL_dz'] = sqrt_g_THETAL_dTHETAL_dz
    df['d_sqrt_U2_plus_V2_dz'] = d_sqrt_U2_V2_dz
    df['d_sqrt_U2_plus_V2_plus_W2_dz'] = d_sqrt_U2_plus_V2_plus_W2_dz 

    # dimensional

    df['W2'] = W2
    df['U2'] = U2
    df['V2'] = V2
    df['THETAL'] = THETAL
    df['QT2'] = QT2
    df['Q2DIFTR'] = Q2DIFTR

    # height = np.array([i for i in range(len(THETAL))])


    # plt.plot(THETAL, height)
    # plt.savefig('test.png')

    # np.array([ncdf.variables['THETAL']]).flatten('F').tofile('thetal.csv', sep=',\n')
    # test_df = pd.DataFrame()
    # test_df['THETAL'] = np.array([ncdf.variables['THETAL']]).flatten('F')
    # test_df.to_csv('thetal_pandas.csv', index=False)

    return levels




def gradient(var, z):
    '''
    This function takes in a variable and returns the gradient of that variable. The gradient is calculated using an approximation to avoid divide by zero errors.

    dTHETAL/dz(n) = ( THETAL(n+1)-THETAL(n) ) / ( z(n+1) - z(n) )
    This function is for a different file format, and is not used in the final model. Use original_gradient instead.
    '''
    # print(len(var)-1)
    # print(len(var)-361)
    # print(len(var)-360)
    gradient = []
    for i in range(len(var)-1):
        if i > (len(var)-361):
            gradient.append(gradient[-359])
        else:
            gradient.append((var[i+360]-var[i])/(z[i+360]-z[i]))

    

    return np.array(gradient)

def original_gradient(var, z):
    '''
    This function takes in a variable and returns the gradient of that variable. The gradient is calculated using an approximation to avoid divide by zero errors.

    dTHETAL/dz(n) = ( THETAL(n+1)-THETAL(n) ) / ( z(n+1) - z(n) )
    '''
    # print(len(var)-1)
    # print(len(var)-361)
    # print(len(var)-360)
    gradient = []
    for i in range(len(var)):
        if i == (len(var)-1):
            gradient.append(gradient[-1])
        else:
            gradient.append((var[i+1]-var[i])/(z[i+1]-z[i]))

        

    

    return np.array(gradient)

def graph_last_3_hours(df, variable, filename, levels):
    '''
    This function graphs the average of the last 3 hours of the data
    This is done by getting the average of each height level for the last 3 hours
    3 hours is 180 minutes. The number of rows is the variable levels, which represents the rows per timestep (min). So 180*levels = x rows
    '''
    vars = []

    offset = levels*180

    var_last_3_hours = np.array(df[variable][-offset:])
    z_last_3_hours = np.array(df['z'][-offset:])

    # print('amde it here')
    # print(len(var_last_3_hours))
    # print(len(z_last_3_hours))

    for i in range(levels):
        vars.append([var_last_3_hours[i+j] for j in range(0,(13500-74), levels)])

    vars = np.array(vars)

    vars = np.mean(vars, axis=1)

    plt.plot(vars, z_last_3_hours[-levels:])
    plt.xlabel(variable)
    plt.ylabel('height')
    plt.savefig(filename + '.png')
    plt.clf()


    

def graphing_testing(df, levels):

    variables = ['Q2DIFTR_on_QT2', 'W2_on_U2_plus_V2', 'sqrt_g_on_THETAL_dTHETAL_dz', 'd_sqrt_U2_plus_V2_dz', 'd_sqrt_U2_plus_V2_plus_W2_dz', 'W2', 'U2', 'V2', 'THETAL', 'QT2', 'Q2DIFTR', 'CLD', 'WSKEW']

    for variable in variables:
        graph_last_3_hours(df, variable, variable + '_last_3_hours', levels)



if __name__ == '__main__':
    df = pd.DataFrame()

    ncdf = netCDF4.Dataset(netcdf_file)

    # make_index(df, ncdf)
    # make_z_col(df, ncdf)
    levels = create_df_1040_scalar_dissipation(df, ncdf)
    
    graphing_testing(df, levels)

    # df["z"].to_csv('z.csv', index=False)

    df.to_csv('dycoms_rf01.csv', index=False)

    
