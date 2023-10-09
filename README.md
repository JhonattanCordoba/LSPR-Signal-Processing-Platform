**Universidade Federal de Minas Gerais - UFMG**

**LSPR Signal Processing Platform**

Felipe Moreira Fernandes Teixeira

This document has instructions and brief descriptions of the functions present in the LSPR signal processing platform. First a description of the input data standardization is given, then the functionality of the code is explained and a list of the libraries utilized is presented. Please be mindful of the file format of the inputs and the expected layout of the variables within them. Additionally, it is important to note that the code was tested using Python version 3.11.5.

This platform requires two input files. The first one is the excel spreadsheet containing the biosensing data, which should be presented in a .xlsx format. Additionally, the following requirements must be met:

1) Excel files must only contain a single sheet, because the code is only prepared to read the first sheet of a file with multiple sheets
2) The first column must be the wavelength value, used for the x axis when plotting graphs
3) All columns must be the same length, meaning you must have a data value for each presented wavelength value
4) The first value of each column must be the variable's name, used for the legends on graphs

The second input file is named analysis_parameters.csv and contains key information for the code to function as expected. Each column must contain one of the variables and the first row of the file must contain the variable names, which must not be changed. The information needed is as follows:

1) n_col_data:
	
	 a) For each group of input data, the variable n_col_data informs the number of columns each group presents. For example, when using biosensing data, positive data for a certain virus count as one group while negative data counts as another
		
	 b) If one data does not fit inside a group, it is counted as a group of its own and its n_col_data number is 1
		
	 c) The order of the values of n_col_data must correspond to the order the data is presented in the .xlsx input data. The wavelength column does not count for this definition
		
	 d) The values of n_col_data must be presented within square brackets, [ ], and separeted with commas
		
	 e) Example: if the data has 1 reference signal, 5 signals seropositive for Virus A, and 3 signals seronegative for Virus A, the n_col_data value is [1,5,3]

2) groups:
	
	 a) As with n_col_data, the variable groups refers to the names of each group and similarly should be presented within square brackets, [ ], and separated with commas
		
	 b) Since this is a variable that contains text, each group name must be presented within simple quotation marks, ' '
		
	 c) Example: ['Virus A Reference','Virus A Positive','Virus A Negative']

3) min_intensity and max_intensity:
	
	 a) These variables are used by the segmentation function when defining the range where the peak of the normalized data within your region of interest is. Generally, when normalizing between 0 and 1, the peak of the function is going to be set to 1, meaning if the peak of the function is within the region of interest, the value of max_intensity should always be 1.
		
	 b) These variables accept float values with dots as the decimal marker

4) min_wavelength and max_wavelength:
	
	 a) These variables define your region of interest. For example, if you know the peak of the resonance of your data exists between 500nm and 800nm, min_wavelength = 500 and max_wavelength = 800
		
	 b) These variables should only receive values that are present in the wavelength column of the .xlsx data spreadsheet

5) fig_width, fig_height, fig_dpi, fig_font, fig_legend_font, legend_loc:
	
	a) These variables are all necessary to format the output graphs
		
	b) fig_width and fig_height:
	
		i) These variables determine the dimensions of the plotted graphs
		
		ii) The library matplotlib natively utilizes values in inches, meaning fig_width and fig_height must be in inches
		
		iii) Example: fig_width = 15, fig_height = 10
		
	 c) fig_dpi:
	
		i) This variable determines the resolution of the plotted graphs and should be changed depending on the desired use of the images
	
		ii) Example: fig_dpi = 300
		
	 d) fig_font and fig_legend_font:
		
		i) These variables determine the font size for the tile and axis as well as the legend on the figures produced by the platform
		
		ii) Example: fig_font = 12 and legend_font = 10
		
	 e) legend_loc:
	
		i) This variable determines the position of the legend on the graphs
		
		ii) Its value must be one of the following options: upper left, lower left, upper right, lower right
		
		iii) Example: legen_loc = upper left

6) recurrence_percentage, polinomial_size and reference:
	
	 a) These variables are necessary to determine the behavior of certain functions to analyze the data
		
	 b) recurrence_percentage:
	
		i) When generating the recurrence plot, a threshold percentage is needed to plot only the values that are above this threshold
		
		ii) This variable expects an integer value that corresponds to the percentage
		
		iii) Example: for 10%, recurrence_percentage = 10
		
	 c) polinomial_size: 
	
		i) When interpolating the data, the size for the polinomial generated must be chosen
		
		ii) This variable expects an integer value
		
		iii) Example: polinomial_size = 3
		
	 d) reference:
	
		i) When calculating the wavelength shift of the data as well as when performing the statistical analysis for the boxplot, it is necessary to have a value as the reference
		
		ii) The reference value must be an integer corresponding to the column of the reference data. Here two considerations are important: the wavelength column is not counted for the reference value determination and Python begins counting array indexes at 0
		
		iii) Example: with the wavelength column as the first column of the data input and the reference data as the second, since the wavelength isn't counted and the data used as reference is the first presented, reference = 0. If the reference was the second column of relevant data, again not counting the wavelength column, reference value would be 1 

7) data_name, file_path and file_name
		
	 a) data_name:
		
			i) This variable corresponds to the name of the data input and is used for the name and title of the figures
			
			ii) Example: data_name = Virus A
	
	 b) file_path:
		
		i) This variable corresponds to the path in which the libraries, the .py code, the data .xlsx file and the analysis_parameters.csv file are located

		ii) Its value must always end with a forward slash mark, /
		
		iii) Example: E:/Folder1/Folder2/

		iv) Note: the path for the analysis_parameters.csv file must be changed manually within the .py code at the variable "inputs" within the "folder_and_inputs" function, which can be found below the libraries imports
		
	 c) file_name: 
	
		i) This variable corresponds to the name of the .xlsx file
		
		ii) Its value must include the .xlsx file format at the end
		
		iii) Example: Virus_A_Data.xlsx

This code was created in order to allow fast and precise analysis of LSPR-based biosensing data. In favor of accomplishing this goal, the analysis platform performs normalization of the entry data, interpolation, derivative, calculates the horizontal shift at inflection points, segments the data and performs statistical analysis with boxplot and recurrence plot, and provides information about the system presented with the calculation of the Maximum Lyapunov Exponent and of the Detrended Fluctuation Analysis coefficient. The libraries utilized were pandas, numpy, matplotlib, scipy, seaborn, os, pyts and nolds, as well as their respective dependencies.

The functions implemented on the code are described below, with their respective inputs and outputs.

1) folder_and_inputs():

	 a) Description: The folder_and_inputs function gets the global parameters and creates folder to save the results.
	
	 b) Input: 
	
		i) none
	
	 c) Output: 
	
		i) all of the variables present at the analysis_parameters.csv file

2) interpolator_dec(s_prev, segmentation, i, intensity):

	 a) Description: The interpolator_dec function finds the distance between segments at the same intensity
	
	 b) Input: 
	
		i) s_prev: the previous segment
	
		ii) segmentation: the section of the data
	
		iii) i: always equal to 1
	
		iv) intensity: a value between min_intensity and max_intensity
	
	 c) Output: 

		i) delta_lambd: the distance between the current and the previous segments in nanometers

3) segmentation(sample, reference_signal, max_intensity, min_intensity)

	 a) Description: The segmentation function creates segments of normalized data within the interest region
	
	 b) Input: 
	
		i) sample: a sample of the data from the boxplot function
	
		ii) reference_signal: the signal used as reference for shift calculation
	
		iii) max_intensity: the chosen max intensity, usually 1, defined in analysis_parameters.csv
	
		iv) min_intensity: the chosen min intensity, defined in analysis_parameters.csv
	
	 c) Output: 
	
		i) vec: a vector with the calculated distance between segments at the same intensity

4) original_data(df,data_name)

	 a) Description: The original_data function reads and prepares the input data, and plots its unmodified and normalized versions
	
	 b) Input: 
	
		i) df: a dataframe created from the .xlsx file
	
		ii) data_name: the name of the data analyzed, defined in analysis_parameters.csv
	
	 c) Output: 
	
		i) wavelength: the wavelength vector
	
		ii) df: the modified dataframe, ready to be used in other functions
	
		iii) normalized_df: the dataframe that contains the normalized data
	
		iv) N_row_df and N_col_df: the dimensions of the dataframe, i.e., the number of rows and columns

5) interpolation(df,N_row_df,N_col_df,data_name)

	 a) Description: The interpolation function interpolates the input data using a polinomial of chosen size
	
	 b) Input: 
	
		i) df: the prepared dataframe
	
		ii) N_row_df and N_col_df: the dimensions of the dataframe, i.e., the number of rows and columns
	
		iii) data_name: the name of the data analyzed, defined in analysis_parameters.csv
	
	 c) Output: 
	
		i) df_interpol: the dataframe that contains the interpolated data
	
		ii) normalized_interpol: the dataframe that contains the normalized interpolated data
	
		iii) polynoms: the polinomial generated for interpolating

6) derivative(df,polynoms)

	 a) Description: The derivative function calculates the derivative of the input data and the wavelength shifts from the reference data at inflection points
	
	 b) Input: 
	
		i) df: the prepared dataframe
	
		ii) polynoms: the polinomial generated for interpolating
		
	 c) Output: 
	
		i) df_derivative: the dataframe that contains the derivative data
		
		ii) zero_distances: the wavelength shifts from the reference data at inflection points 

7) boxplot(df_interpol)

	 a) Description: Performs statistical analysis of the data and plots the boxplot. This function might produce a warning if the seaborn library installed isn't version 0.13.0 or later
	
	 b) Input: 
	
		i) df_interpol: the dataframe that contains the interpolated data
	
	 c) Output: 
	
		i) pivot_data_df: the dataframe that contains the boxplot data

8) MLE_DFA(df_interpol)

	 a) Description: The MLE_DFA function calculates the Maximum Lyapunov Exponent (MLE) and the Detrended Fluctuation Analysis (DFA) coefficient in order to characterize the analyzed system
	
	 b) Input: 
	
		i) df_interpol: the dataframe that contains the interpolated data
	
	 c) Output: 
		
		i) none

9) recurrence(df_interpol)

	 a) Description: The recurrence function creates and plots the recurrence of the input data
	
	 b) Input: 
		
		i) df_interpol: the dataframe that contains the interpolated data
	
	 c) Output: 
		
		i) none

10) save_data(normalized_df,df_derivative,zero_distances,pivot_data_df,normalized_interpol)

	 a) Description: The save_data function prepares and saves the data as .csv files
	
	 b) Input: 
		
		i) normalized_df: the dataframe that contains the normalized data 
					
		ii) df_derivative: the dataframe that contains the derivative data
					
		iii) zero_distances: the wavelength shifts from the reference data at inflection points
						
		iv) pivot_data_df: the dataframe that contains the boxplot data
						
		v) normalized_interpol: the dataframe that contains the normalized interpolated data
	
	 c) Output: 

		i) none
