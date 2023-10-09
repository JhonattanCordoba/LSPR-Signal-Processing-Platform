# Universidade Federal de Minas Gerais - UFMG

# LSPR signal processing platform

# Authors:
# Felipe Moreira Fernandes Teixeira,
# Felipe Aragão Nogueira de Freitas,
# Gabriel Lopes Machado, and
# Júlia de Backer Pacífico

# Standardization of data needed for this code: please check the README.txt file

# Imports the relevant libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.misc import derivative
import seaborn as sns
import os
import pyts
from pyts import image
from pyts.multivariate.image import JointRecurrencePlot
import nolds

# The folder_and_inputs function gets the global parameters and creates folder to save the results.
# Before running the code, please change the path for the file analysis_parameters.csv
# at the variable inputs on line 31 of the code
def folder_and_inputs():
  # Reads relevant input parameters from spreadsheet
  inputs = pd.read_csv('E:/Python 3.11.5/TCC test/analysis_parameters.csv', sep = ';')
  inputs = inputs.iloc[0]

  min_intensity = inputs["min_intensity"]
  max_intensity = inputs["max_intensity"]

  min_wavelength = inputs["min_wavelength"]
  max_wavelength = inputs["max_wavelength"]

  n_col_data_str = inputs["n_col_data"]
  n_col_data_str = n_col_data_str.replace("[", "")
  n_col_data_str = n_col_data_str.replace("]", "")
  n_col_data_str = n_col_data_str.replace(" ", "")
  n_col_data = n_col_data_str.split(",")
  n_col_data = list(filter(None, n_col_data))
  n_col_data = map(int, n_col_data)
  n_col_data = list(n_col_data)

  groups_str = inputs["groups"]
  groups_str = groups_str.replace("[", "")
  groups_str = groups_str.replace("]", "")
  groups_str = groups_str.replace(" ", "")
  groups_str = groups_str.replace("'", "")
  groups = groups_str.split(",")

  width = inputs["fig_width"]
  height = inputs["fig_height"]
  resolution = inputs["fig_dpi"]
  font = inputs["fig_font"]
  legend_font = inputs["fig_legend_font"]
  legend_loc = inputs["legend_loc"] # options: upper left, upper right, lower left, lower right
  data_name = inputs["data_name"]

  # Percentage of the threshold for the recurrence plot
  recurrence_percentage = inputs["recurrence_percentage"] # 10 generates the best graph

  # Polinomial size
  polinomial_size = inputs["polinomial_size"] #choose polinomial size for interpolation

  # Reference for the shift calculation
  # A value of 0 means first data column after the wavelength column
  reference = inputs["reference"]

  # Gets path for the input data and the name of the .xlsx file
  file_path = inputs["file_path"]
  file_name = inputs["file_name"]

  # Creates folder to save figures and spreadsheets in
  if os.path.isdir(file_path+'results') == 0:
    os.makedirs(file_path+'results')
  results_path = 'results/'

  full_path = file_path+results_path
  full_file_path = file_path+file_name

  return [full_path,full_file_path,min_intensity,max_intensity,min_wavelength,max_wavelength,n_col_data_str,n_col_data,groups_str,groups,width,height,resolution,font,legend_font,legend_loc,recurrence_percentage,polinomial_size,reference,data_name]

[full_path,full_file_path,min_intensity,max_intensity,min_wavelength,max_wavelength,n_col_data_str,n_col_data,groups_str,groups,width,height,resolution,font,legend_font,legend_loc,recurrence_percentage,polinomial_size,reference,data_name] = folder_and_inputs()

# The functions interpolator_dec and segmentation were first created on Matlab by Gabriel Lopes Machado
# and translated to Python by Felipe Moreira Fernandes Teixeira and Felipe Aragão Nogueira de Freitas
# The comments on these two functions were written by Gabriel L. Machado and translated by Felipe M. F. Teixeira
# Both functions work to prepare the data for the boxplot
# Calculates the distante between lambda values for the same intensity
def interpolator_dec(s_prev, segmentation, i, intensity):

    pre_ref = 0
    pre_data = 0

    for Lambda in range(segmentation[0][2], segmentation[0][3]):
        if(s_prev[0][Lambda]<=intensity):
            pre_ref = Lambda
            break

    for Lambda in range(segmentation[i][2], segmentation[i][3]):
        if(s_prev[i][Lambda]<=intensity):
            pre_data = Lambda
            break

    if((pre_data==0) or (pre_ref==0)):
        return -1
    else:
        dif1_ref = abs(s_prev[1][pre_ref] - intensity)
        dif2_ref = abs(s_prev[1][pre_ref-1] - intensity)
        dif_ref = dif1_ref + dif2_ref

        pre_ref = pre_ref + (wavelength[0]-1)

        dif1_ref_aux = dif1_ref

        dif1_ref = dif2_ref/dif_ref
        dif2_ref = dif1_ref_aux/dif_ref

        delta_lambd1 = dif1_ref*pre_ref + dif2_ref*(pre_ref-1)


        dif1_data = abs(s_prev[1][pre_data] - intensity)
        dif2_data = abs(s_prev[1][pre_data-1] - intensity)
        dif_data = dif1_data + dif2_data

        pre_data = pre_data + (wavelength[0]-1)

        dif1_data_aux = dif1_data

        dif1_data = dif2_data/dif_data
        dif2_data = dif1_data_aux/dif_data

        delta_lambd2 = dif1_data*pre_data + dif2_data*(pre_data-1)

        delta_lambd = (delta_lambd2 - delta_lambd1)

    return delta_lambd

def segmentation(sample, reference_signal, max_intensity, min_intensity):
    n_data = 2
    data = []

    data.append(reference_signal)
    data.append(sample)

    s_prev = data.copy() # auxiliar variable to store in use data

    vec = []
    # This section of the code subtracts the minimum value to stablish the lower limit
    # of the data at 0, then divides the data by the maximum to define the superior limit at 1,
    # meaning it normalizes in the interval between 0 and 1
    for t in range (n_data): # Or n_data + 1
        s_prev[t] = (data[t]-data[t].min()) / (data[t].max() - data[t].min())


    # Getting the 2 regions of analysis, taking the first minimum point, maximum point and last minimum point after the global maximum
    segmentation = []

    for t in range(n_data):
        seg_aux = []
        for t in range(4): # 4 means four points to define the 2 interpolation regions
            seg_aux.append(0)
        segmentation.append(seg_aux)

    for t in range(n_data):
        indices = [i for i, x in np.ndenumerate(s_prev[t]) if x == 0.] # first 0 point
        segmentation[t][0] = indices[-1][0]

        indices = [i for i, x in np.ndenumerate(s_prev[t]) if x == 1.] # first maximum  point
        segmentation[t][1] = indices[0][0]

        indices = [i for i, x in np.ndenumerate(s_prev[t]) if x == 1.] # second maximum point wich is expected to be the same as the first 1 point
        segmentation[t][2] = indices[-1][0]

        minimum = np.amin(s_prev[t][segmentation[t][1]:])
        indices = [i for i, x in np.ndenumerate(s_prev[t]) if x == minimum] # second minimum point
        segmentation[t][3] = indices[-1][0]

    vec = []

    for intensity in list(np.arange(min_intensity, max_intensity, 0.01)):

        lbd = interpolator_dec(s_prev, segmentation, 1, intensity)

        if lbd !=-1:
            vec.append(lbd)
        else:
            print("Error on lambda calculation.")
            break

    return vec

# This section of the code defines a data frame from the input data
df = pd.read_excel(full_file_path)

# The original_data function treats the input data in order to separate the wavelength vector
# from the rest of the data and excludes cells that are empty or contain values that are not numbers
# It also plots the unmodified and the normalized input data
def original_data(df,data_name):

  # Gets initial size of file
  N_row_df = df.shape[0]
  N_col_df = df.shape[1]

  # Finds empty columnns (code interprets empty columns in excel as unnamed columns)
  cols = list()
  for column in df.columns:
      if "Unnamed" in column:
          cols.append(column)

  # Removes last row that indicates end of file
  df.drop([0, N_row_df-1], inplace = True)

  # Removes empty columns and replaces remaining NaN with 0
  df.drop(cols, axis = 1, inplace = True)
  df = df.fillna(0)

  # Gets wavelength column and removes it from the data frame
  wavelength = df["Wavelength"].values
  df.drop(["Wavelength"], axis = 1, inplace = True)

  header = df.columns.tolist()

  # Gets size of file after removing unecessary data
  N_row_df = df.shape[0]
  N_col_df = df.shape[1]

  # Checks if the reference value is within the number of columns in the dataframe
  if reference > N_col_df:
    print("Error: Column chosen as reference does not exist")

  # Normalization of the signals
  # Normalized data between 0 and 1
  normalized_df = (df-df.min()) / (df.max() - df.min())

  # Gets size of the normalized data
  N_row_normalized = normalized_df.shape[0]
  N_col_normalized = normalized_df.shape[1]

  # plots data
  plt.figure(figsize = (width,height), dpi = resolution)
  plt.xticks(fontsize = font)
  plt.yticks(fontsize = font)
  plt.xlabel('Wavelength (nm)', fontsize = font)
  plt.ylabel('Intensity', fontsize = font)
  plt.title('Light Absorption', fontsize = font)
  plt.plot(wavelength, df, lw=2, label = header)
  plt.legend(fontsize=legend_font,loc=legend_loc)
  fig_name = data_name + "_Unmodified_Data"
  plt.savefig(full_path + fig_name, dpi=resolution)
  plt.show()

  # plots normalized data
  plt.figure(figsize = (width,height), dpi = resolution)
  plt.xticks(fontsize = font)
  plt.yticks(fontsize = font)
  plt.xlabel('Wavelength (nm)', fontsize = font)
  plt.ylabel('Normalized Intensity', fontsize = font)
  plt.title('Normalized Light Absorption', fontsize = font)
  plt.plot(wavelength, normalized_df, lw=2, label = header)
  plt.legend(fontsize=legend_font,loc=legend_loc)
  fig_name = data_name + "_Normalized_Data"
  plt.savefig(full_path + fig_name, dpi=resolution)
  plt.show()

  return [wavelength,df,normalized_df,N_row_df,N_col_df]

[wavelength,df,normalized_df,N_row_df,N_col_df] = original_data(df,data_name)

# The interpolation function interpolates the input data with a polinomial
# whose size is set by the analysis_parameters.csv file
# This interpolation is done in order to remove noise in the data and
# prevent inflection points caused by noise, which would have an impact
# when calculating the derivative
# This function also plots the interpolated and normalized interpolated data
def interpolation(df,N_row_df,N_col_df,data_name):
  df_interpol = []
  polynoms = []
  normalized_interpol = np.empty((N_row_df, N_col_df))

  cont = 0
  for column in df:
    x = wavelength
    y = df[column].values
    s = np.polyfit(x, y, polinomial_size)
    p = np.poly1d(s)
    polynoms.append(p)

    y_new = np.polyval(s, x)
    df_interpol.append(y_new)
    normalized_interpol[:,cont] = (y_new-y_new.min()) / (y_new.max() - y_new.min())
    cont = cont + 1

  plt.figure(figsize = (width,height), dpi = resolution)
  cont = 0
  for column in df:
    plt.plot(x, df_interpol[cont], label=column,lw=2)
    cont = cont + 1

  plt.xticks(fontsize = font)
  plt.yticks(fontsize = font)
  plt.xlabel('Wavelength (nm)', fontsize = font)
  plt.ylabel('Intensity', fontsize = font)
  plt.title('Interpolated Light Absorption', fontsize = font)
  plt.legend(fontsize=legend_font,loc=legend_loc)
  fig_name = data_name + "_Interpolated_Data"
  plt.savefig(full_path + fig_name, dpi=resolution)
  plt.show()

  plt.figure(figsize = (width,height), dpi = resolution)
  cont = 0
  for column in df:
    plt.plot(x, normalized_interpol[:,cont], label=column,lw=2)
    cont = cont + 1

  plt.xticks(fontsize = font)
  plt.yticks(fontsize = font)
  plt.xlabel('Wavelength (nm)', fontsize = font)
  plt.ylabel('Normalized Intensity', fontsize = font)
  plt.title('Normalized Interpolated Light Absorption', fontsize = font)
  plt.legend(fontsize=legend_font,loc=legend_loc)
  fig_name = data_name + "_Normalized_Interpolated_Data"
  plt.savefig(full_path + fig_name, dpi=resolution)
  plt.show()

  return [df_interpol,normalized_interpol,polynoms]

[df_interpol,normalized_interpol,polynoms] = interpolation(df,N_row_df,N_col_df,data_name)

# The derivative function calculates the derivative of the input data
# Additionally, this function plots the derivative data and calculates the
# distance of the signals from the reference at the inflection points
def derivative(df,polynoms):
  plt.figure(figsize = (width,height), dpi = resolution)

  df_derivative = []
  zeros = []
  i = 0
  legends = df.columns
  for column in df_interpol:

      x = wavelength
      d = polynoms[i].deriv()
      zeros.append(d.r)

      derivative_df = d(x)
      y_new_d = (derivative_df-derivative_df.mean()) / derivative_df.std()
      df_derivative.append(y_new_d)
      plt.plot(x, y_new_d, label=legends[i],lw=2)
      i+=1

  plt.plot([min(wavelength),min_wavelength,max_wavelength,max(wavelength)],[0,0,0,0], label ='zero_line')
  plt.xticks(fontsize = font)
  plt.yticks(fontsize = font)
  plt.xlabel('Wavelength (nm)', fontsize = font)
  plt.ylabel('Normalized Intensity', fontsize = font)
  plt.title('Derivative of the Interpolated Light Absorption', fontsize = font)
  plt.legend(fontsize=legend_font,loc=legend_loc)
  fig_name = data_name + "_Derivative_of_Interpolated_Data"
  plt.savefig(full_path + fig_name, dpi=resolution)
  plt.grid()
  plt.show()

  all_zeros = np.round(np.real(zeros)).astype(int)

  for zeros in all_zeros:
      for i in range(len(zeros)):
          if zeros[i] <= min_wavelength or zeros[i] >= max_wavelength:
              zeros[i] = 0

  base_values = all_zeros[reference]
  zero_distances = []
  for zeros in all_zeros:
      zero_distances.append(zeros - base_values)

  return [df_derivative,zero_distances]

[df_derivative,zero_distances] = derivative(df,polynoms)

# The boxplot function performs an statistical analysis of the data
# illustrating the spread and differences of the chosen groups of data.
# This function provides informations such as the mean and standard deviation of the data.
# Additionally, this function creates the boxplot
def boxplot(df_interpol):
  conj_vec = []
  ref = df_interpol[reference] 
  index = 0
  counter = 0

  for i in n_col_data:
      vec = []
      while index < (counter + i):
          sample = df_interpol[index] 
          vec.extend(segmentation(sample, ref, max_intensity, min_intensity))
          index += 1

      counter = index
      conj_vec.append(vec)

  boxplot_data = []

  for group_idx in range(len(groups)):
      for i in range(len(conj_vec[group_idx])):
          aux_dic = {}
          aux_dic['Type'] = groups[group_idx]
          aux_dic['Lambda'] = conj_vec[group_idx][i]
          boxplot_data.append(aux_dic)

  data_df = pd.DataFrame(boxplot_data)

  sns.set(rc={'figure.figsize':(width,height),'figure.dpi':resolution})
  sns.set_theme(style="whitegrid")
  ax = sns.boxplot(data_df,x="Type", y="Lambda", width=0.5)
  ax.set_xlabel("Type",fontsize=font)
  ax.set_ylabel("Delta lambda",fontsize=font)
  plt.xticks(fontsize = font)
  plt.yticks(fontsize = font)
  plt.title('Boxplot '+ data_name, fontsize = font)
  fig_name = data_name + "_Boxplot"
  plt.savefig(full_path + fig_name, dpi=resolution)

  pivot_data = {}
  i = 0
  for Type in groups:
      pivot_data[Type] = conj_vec[i]
      i+=1

  pivot_data_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in pivot_data.items() ]))
  print(pivot_data_df.describe())

  return pivot_data_df

pivot_data_df = boxplot(df_interpol)

# Maximum Lyapunov Exponent (MLE) and Detrended Fluctuation Analysis (DFA)

# Source for the MLE analysis
# Jiaru Yang, et al., Maximum Lyapunov exponent‑based multiple chaotic slime mold algorithm for real‑world optimization
# DOI: 10.1038/s41598-023-40080-1
# As defined by Jiaru Yang, et al., the MLE shows if a system is performing chaotic or periodic motion

# Source for the DFA analysis:
# R. M. Bryce and K. B. Sprague, Revisiting detrended fluctuation analysis
# DOI: 10.1038/srep00315
# As defined by R. M. Bryce and K. B. Sprague, the  DFA enables the estimation 
# of the power law scaling (Hurst exponent) of a systems’ signal
# in the presence of (extrinsic) nonstationaries while eliminating 
# spurious detection of long-range dependence 

# The MLE_DFA function calculates and saves the ME and DFA for the input data
# Additionaly, it provides the meaning found for each value
def MLE_DFA(df_interpol):
  df_for_MLE = pd.DataFrame(np.transpose(np.array(df_interpol)))
  MLE = pd.DataFrame(np.zeros((1,len(df_interpol))))
  DFA = pd.DataFrame(np.zeros((1,len(df_interpol))))
  MLE_meaning = list([])
  DFA_meaning = list([])

  for i in range(len(df_interpol)):
    MLE[i] = nolds.lyap_r(df_for_MLE[i],lag=27, min_tsep=137)
    DFA[i] = nolds.dfa(df_for_MLE[i])

    if np.array(MLE[i]) >= 0: MLE_meaning.append("Positive MLE: exponential divergence, chaotic motion")
    elif np.array(MLE[i]) < 0: MLE_meaning.append("Negative MLE: periodic motion")

    if np.array(DFA[i]) < 0.5: DFA_meaning.append("anti-correlated")
    elif np.array(DFA[i]) == 0.5: DFA_meaning.append("uncorrelated")
    elif 1 > np.array(DFA[i]) > 0.5: DFA_meaning.append("correlated")
    elif np.array(DFA[i]) > 1: DFA_meaning.append("detrending was not performed")

  MLE_DFA_interpretation = np.array(np.zeros((4,len(df_interpol)), dtype=object))

  for i in range(len(df_interpol)):
    MLE_DFA_interpretation[0,i] = np.array(MLE[i])
    MLE_DFA_interpretation[1,i] = MLE_meaning[i]
    MLE_DFA_interpretation[2,i] = np.array(DFA[i])
    MLE_DFA_interpretation[3,i] = DFA_meaning[i]

  MLE_DFA_interpretation = pd.DataFrame(MLE_DFA_interpretation)
  MLE_DFA_interpretation.to_csv(full_path +'MLE_DFA_interpretation.csv',index=False, sep=',')

MLE_DFA(df_interpol)

# The recurrence function aids the analysis of the data, indicating the self-similarity
# of the data and distinguishing chaos from randomness
# This function creates the recurrence data and plots it
def recurrence(df_interpol):
  rp = pyts.image.RecurrencePlot(dimension=1, threshold='point', percentage=recurrence_percentage)
  X_rp = rp.fit_transform(df_interpol)

  plt.figure(figsize = (width,height), dpi = resolution)
  plt.imshow(X_rp[0], cmap='binary', origin='lower')
  plt.xlabel('Wavelength (nm)', fontsize = font)
  plt.ylabel('Wavelength (nm)y', fontsize = font)
  plt.title('Recurrence Plot '+ data_name, fontsize = font)
  fig_name = data_name + "_Recurrence"
  plt.savefig(full_path + fig_name, dpi=resolution)
  plt.show()

recurrence(df_interpol)

# The save_data function prepares the relevant dataframes and saves them to .csv spreadsheets
def save_data(normalized_df,df_derivative,zero_distances,pivot_data_df,normalized_interpol):
  derivative_df = {}
  zero_distances_df = {}
  columns = normalized_df.columns

  for i in range(len(columns)):
      derivative_df[columns[i]] = df_derivative[i]
      zero_distances_df[columns[i]] = zero_distances[i]

  derivative_df = pd.DataFrame(derivative_df)
  zero_distances_df = pd.DataFrame(zero_distances_df)
  zero_distances_df = zero_distances_df.loc[~(zero_distances_df==0).all(axis=1)]
  normalized_interpol = pd.DataFrame(normalized_interpol)

  # save data to .csv
  normalized_df.to_csv(full_path + data_name + '_normalized_data.csv',index=False, sep=',')
  derivative_df.to_csv(full_path + data_name + '_derivative.csv',index=False, sep=',')
  zero_distances_df.to_csv(full_path + data_name + '_shift.csv',index=False, sep=',')
  pivot_data_df.to_csv(full_path + data_name + '_boxplot.csv',index=False, sep=',')
  normalized_interpol.to_csv(full_path + data_name + '_normalized_interpolated_data.csv',index=False, sep=',')

save_data(normalized_df,df_derivative,zero_distances,pivot_data_df,normalized_interpol)
