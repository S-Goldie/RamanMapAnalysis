__author__ = "Stuart Goldie"
__copyright__ = "Copyright 2018, Stuart Goldie"
__date__ = "Tue Mar 19 15:02:46 2019"
__credits__ = ["Stuart Goldie", "Tim Callow", "Matteo Degiacomi"]
__version__ = "3.3.1"
__email__ = "stuart.j.goldie@durham.ac.uk"

import csv    #csv to read and write csv files
import os     #os to check before overwriting existing files
import numpy as np     #numpy for various mathematical operations in data processing
import matplotlib.pyplot as plt    #matplot to produce histograms and display spectra
from timeit import default_timer as timer    #timer to allow an estimated fitting time to be displayed
import sys
from lmfit.models import PolynomialModel, LorentzianModel    #lmfit does the actual least squares minimisation fitting
try:
    from tkinter import Tk    #Tk simply gui to allow file location by explorer window
    from tkinter.filedialog import askopenfilename
except:
    pass

###############################################################################
    #open file, extract data and confirm user inputs#
###############################################################################
    
file_check = False
#user inputs file name, file_check loop repeats until valid file name is entered
while file_check == False:
    try:                #try Tk module to open file, if module incompatable with os then defult text based input will appear
        Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
        input_name = askopenfilename() # show an "Open" dialog box and return the path to the selected file
        try:
            with open(input_name, 'r') as fin:
                cin = csv.reader(fin)
                datalist = [row for row in cin] #using a list comprehension
                shortname = input_name[:(input_name.find('.'))]     #used as name without the file extension for creating results files
                longname = input_name                               #used for importing file that requires the .csv file extension
                data = np.array(datalist, dtype = float)            #read csv file containing data and store in array called data.
                file_check = True
                print('File opened')
        except:         #if unable to open, try again
            print('File not valid')
    except:             #if the Tk module returns an error, a backup text based input will operate
        input_name = input("Open File: ")
        if (os.path.isfile(str(input_name))) == True:
            shortname = input_name[:(input_name.find('.'))]     #used as name without the file extension for creating results files
            longname = input_name                               #used for importing file that requires the .csv file extension
            file_check = True
        elif (os.path.isfile(str(input_name + '.csv'))) == True:
            shortname = input_name
            longname = str(input_name + '.csv')
            file_check = True
        else:
            print('File not found')
            input("Press enter to exit")
            sys.exit()

go_on = False
while go_on == False:       #only when user confirms all inputs will go_on become true and fitting starts
    #user input to control fitting of D+G peak
    DGinput_flag = False
    while DGinput_flag == False:
        DGinput = input('Attempt to fit D+G peak? [Y/N]')
        if str(DGinput) == "Y" or str(DGinput) == "y" or str(DGinput) == "yes" or str(DGinput) == "Yes":
            DGflag = True
            DGinput_flag = True
        elif str(DGinput) == "N" or str(DGinput) == "n" or str(DGinput) == "no" or str(DGinput) == "No":
            DGflag = False
            DGinput_flag = True
        else:
            print("Input error")
            DGinput_flag = False
    #user input to control fitting of D prime peak
    Dpinput_flag = False
    while Dpinput_flag == False:
        Dpinput = input('Attempt to fit D\' peak? [Y/N]')
        if str(Dpinput) == "Y" or str(Dpinput) == "y" or str(Dpinput) == "yes" or str(Dpinput) == "Yes":
            Dpflag = True
            Dpinput_flag = True
        elif str(Dpinput) == "N" or str(Dpinput) == "n" or str(Dpinput) == "no" or str(Dpinput) == "No":
            Dpflag = False
            Dpinput_flag = True
        else:
            print("Input error")
            Dpinput_flag = False
    #printflag is used for users to specify additional details to be printed to monitor fit quality
    printflag = input('Print options. \n[Press Enter for Default]')
    print('________')
    print('Fitting ' + str(shortname))
    print('D+G fit: ' + str(DGflag))
    print('D\' fit: ' + str(Dpflag))
    go_oninput = input('Confirm parameters for fitting. [Y/N]  \n[Type EXIT to quit]')
    if str(go_oninput) == "Y" or str(go_oninput) == "y" or str(go_oninput) == "yes" or str(go_oninput) == "Yes":
        go_on = True
    elif str(go_oninput) == "EXIT":
        sys.exit()
    else:
        go_on = False

#define empty lists to hold the integral of the three peaks for use later as an initial guess for the amplitude of the fitting
#each peak has a different list that contains all the parameters for the different spectra
null_spectra_count = 0      #stores the number of points that have no valid G peak, therefore considered null
null_spectra_list = []      #stores the specific rows that had no valid G peak, returned later for reference
null_fit_count = 0          #stores the number of spectra that returned a ValueError during analysis (usually because errors were not analysed)
null_fit_list = []          #stores the specific rows that had a ValueError
results = []                #used later to store a list of dictionaries that holds the fitted parameters
mean_time = []              #empty list to store the time taken for each spectrum, used to generate estimated time

#cycle through the top row of the data, containing the raman shift values, to find the data array coordinates that hold the data relating to the D, G, 2D peaks 
#and background region at the center. xb1 etc. are placeholders that are increased until the data is found, the placeholders then acting as coordinates
xd1 = 0
xg1 = int(xd1)
xg3 = int(xg1)
xg2 = int(xg3)
xp1 = int(xg2)
xp2 = int(xp1)
xe1 = int(xp2)
xe2 = int(xe1)
xf2 = int(xe2)
try:
    while data[0,xd1] < 1200:       #lower limit of D peak
        xd1 += 1
    while data[0,xg1] < 1500:       #lower limit of G peak and crude upper limit of D peak
        xg1 += 1
    while data[0,xg3] < 1600:       #middle of G peak for D' guess
        xg3 += 1
    while data[0,xg2] < 1700:        #upper limit of G peak
        xg2 += 1
    while data[0,xp1] < 2000:       #lower limit of 2D & D+G only fitting
        xp1 += 1
    while data[0,xp2] < 2250:       #upper limit of D & G only fitting
        xp2 += 1
    while data[0,xe1] < 2550:       #lower limit of 2D peak
        xe1 += 1
    while data[0,xe2] < 2800:       #upper limit of 2D peak and lower limit of D+G
        xe2 += 1
    if DGflag == True:
        while data[0,xf2] < 3050:       #upper limit of D+G Peak
            xf2 += 1
except:                      
    print('Data range incomplete')
    input("Press enter to close")
    sys.exit()

#generate an array of x values along the baselile only, excluding all possible peak data
bkgdx = list()
for value in data[0,:xd1]:
    bkgdx.append(value)
for value in data[0,xg2:xe1]:
    bkgdx.append(value)
for value in data[0,xf2:]:
    bkgdx.append(value)
x_background = np.array(bkgdx, dtype = float)

rank = np.linalg.matrix_rank(data)  #the number of rows in the data array.
print('Found ' + str(rank-1) + 'rows')
x = data[0,:]                                   #Raman x coordinates stored on first row of data

###############################################################################
    #error estimation for weighting of chi squared and signal to noise#
###############################################################################

#generate an estimate of the error of each Raman count value using the standard deviation of 11 sequential data points
#This assumes that an average 11 points are close enough together that a perfect spectrometer would measure them as the same background value. 
#These deviations are averaged and a single value is used for the error.
error_list = []     #create an empty list for holding the average error data for each row 
stop = False        #flag used for terminating the cycle that simultanously counts rows and estimates error
error_counter_row = 1   #row=0 contains the x axis data, not included in this step
print('Estimating errors...')
while error_counter_row < rank:   #increment the 'counter row' until all data is analysed, this also provides a measure of the number of rows of data
   yplace_error = 0        #this increments over the x-range moving along the data
   row_deviations = []     #empty list for holding the standard deviation estimate of each point
   while (yplace_error+11 < xd1):     #whilst yplace_error is less than the 'x' location of the d peak, the Raman counts are analysed
       row_deviations.append(np.std(data[error_counter_row,yplace_error:(yplace_error+11)])) #this finds the standard deviation of 11 points and stores then in row_deviations
       yplace_error += 1       #this increments over the x-range
   yplace_error = xg2          #once the data before the D peak is analysed, skip ahead to the 'background' region where no peaks are observed
   while (yplace_error+11 < xe1):     #repeat the above, until the 2D peak region
       row_deviations.append(np.std(data[error_counter_row,yplace_error:(yplace_error+11)]))
       yplace_error += 1
   yplace_error = xf2
   while (yplace_error+10 < x.size):   #terminate before the end so the 11 values exist for standard deviation
       row_deviations.append(np.std(data[error_counter_row,yplace_error:(yplace_error+11)]))
       yplace_error += 1
   error_est = np.array(row_deviations)    #this converts the list of standard deviations into an array for processing
   row_errors = np.array(np.full(shape=x.size,fill_value=(np.mean(error_est)))).tolist()   #find the mean of all the standard deviations
   error_list.append(row_errors)           #add a row, the same size as the data set containing only the mean deviation value. 
                                                #i.e the single value used as an error the same number of times at there are x,y coordinates as lmfit requires an entry for error and data 
   if (error_counter_row % 10) ==0 :
       print('Current row: '+ str(error_counter_row))
   error_counter_row += 1                        #increment to the next row

error_array = np.array(error_list)

###############################################################################
    #define functions that extract data from fitting outputs#
###############################################################################

#define a function to search through the parameter string produced by the fitting module and return the values as a list of dictionaries
def search(output,parameter):
    string = str(output.params[parameter])
    index_initial = string.find('value=')+6     #defines the first character of the fitted value in the output string
    index_mid = string.find('+/-')              #defines the middle of the output string containing '+/-'
    index_final = string.find('bounds')-2       #defines the end of the output string
    value = float(string[index_initial:index_mid-1])        #returns the fitted parameter as a float
    error = float(string[index_mid+4:index_final])          #returns the error in the fitted parameter as a float
    return {parameter:value,parameter+str('_error'):error}

def Goutput(destination, source):               #define a function to extract key parameters of G peak from an fit output (source) and update the destination dictionary
    destination.update(search(source,'G_amplitude'))
    destination.update(search(source,'G_center'))
    destination.update(search(source,'G_fwhm'))
    destination.update(search(source,'G_height'))
    destination.update({'SignalNoise':((destination['G_height']-np.mean(y_background))/errors[1])})

def Doutput(destination, source):               #define a function to extract fitted parameters from D peak
    destination.update(search(source,'D_amplitude'))
    destination.update(search(source,'D_center'))
    destination.update(search(source,'D_fwhm'))
    destination.update(search(source,'D_height'))
    destination.update({'ID/IG':(destination['D_height']/destination['G_height'])})
    destination.update({'ID/IG_error':(destination['ID/IG']*np.sqrt(((destination['D_height_error']/destination['D_height'])**2)+((destination['G_height_error']/destination['G_height'])**2)))})

def DDoutput(destination, source):               #define a function to extract fitted parameters from 2D peak
    destination.update(search(source,'DD_amplitude'))
    destination.update(search(source,'DD_center'))
    destination.update(search(source,'DD_fwhm'))
    destination.update(search(source,'DD_height'))
    destination.update({'IDD/IG':(destination['DD_height']/destination['G_height'])})   #calculate ratios with physical meaning from parameter data
    destination.update({'IDD/IG_error':(destination['IDD/IG']*np.sqrt(((destination['DD_height_error']/destination['DD_height'])**2)+((destination['G_height_error']/destination['G_height'])**2)))})

def Dpoutput(destination, source):               #define a function to extract fitted parameters from D' peak
    destination.update(search(source,'Dp_amplitude'))
    destination.update(search(source,'Dp_center'))
    destination.update(search(source,'Dp_fwhm'))
    destination.update(search(source,'Dp_height'))

def DGoutput(destination, source):               #define a function to extract fitted parameters from D+G peak
    destination.update(search(source,'DG_amplitude'))
    destination.update(search(source,'DG_center'))
    destination.update(search(source,'DG_fwhm'))
    destination.update(search(source,'DG_height'))

#define blank dictionary entries for outputs with no valid peak
blankD = {'D_center': 0.0,'D_center_error': 0.0,'D_height': 0.0,'D_height_error': 0.0,'D_fwhm': 0.0,'D_fwhm_error': 0.0,'D_amplitude': 0.0,'D_amplitude_error': 0.0,'ID/IG': 0.0, 'ID/IG_error': 0.0}
blankDD = {'DD_center': 0.0,'DD_center_error': 0.0,'DD_height': 0.0,'DD_height_error': 0.0,'DD_fwhm': 0.0,'DD_fwhm_error': 0.0,'DD_amplitude': 0.0,'DD_amplitude_error': 0.0,'IDD/IG': 0.0, 'IDD/IG_error': 0.0}
blankDp = {'Dp_center': 0.0,'Dp_center_error': 0.0,'Dp_height': 0.0,'Dp_height_error': 0.0,'Dp_fwhm': 0.0,'Dp_fwhm_error': 0.0,'Dp_amplitude': 0.0,'Dp_amplitude_error': 0.0}
blankDG = {'DG_center': 0.0,'DG_center_error': 0.0,'DG_height': 0.0,'DG_height_error': 0.0,'DG_fwhm': 0.0,'DG_fwhm_error': 0.0,'DG_amplitude': 0.0,'DG_amplitude_error': 0.0}

#define a function to extract a single value from one fitting output. intended to pass fit output into parameters of other models
def par_extract(output,parameter):
    string = str(output.params[parameter])
    parameter_found = False
    try:
        index_initial = string.find('value=')+6     #defines the first character of the fitted value in the output string
        index_mid = string.find('+/-')              #defines the middle of the output string containing '+/-'
        value = float(string[index_initial:index_mid-1])        #returns the fitted parameter as a float
    except ValueError:
        parameter_found = False
        index_start = 12
        while parameter_found == False:
            index_mid = string.find('bounds')-2
            try:
                value = float(string[index_start:index_mid])
                parameter_found = True
            except ValueError:
                index_start += 1
    return float(value)

#function to close a row fitting and print relevant outputs
def row_end(row_counter,start_time):
    print('Row ' + str(row_counter) + ' of ' +str(rank-1))
    end_time = timer()                                                          #records the end time (start time defined at start of fitting)
    mean_time.append(float(end_time - start_time))                              #estimate of time remaining
    min_estimate = ((np.mean(mean_time))*(rank - 1 - row_counter) // 60)
    sec_estimate = ((np.mean(mean_time))*(rank - 1 - row_counter) % 60)
    print('Estimated time remaining: %.0d m %.0d s' % (min_estimate, sec_estimate))

###############################################################################
    #specify most common models to be used during fitting#
###############################################################################

G_peak = LorentzianModel(nan_policy='propagate',prefix='G_')
D_peak = LorentzianModel(nan_policy='propagate',prefix='D_')
DD_peak = LorentzianModel(nan_policy='propagate',prefix='DD_')  #2D peak
DG_peak = LorentzianModel(nan_policy='propagate',prefix='DG_')  #D+G peak
Dp_peak = LorentzianModel(nan_policy='propagate',prefix='Dp_')  #D' peak
poly6 = PolynomialModel(6,nan_policy='propagate', prefix='Background_')
poly2 = PolynomialModel(2,nan_policy='propagate', prefix='Background_')

#main models, additional specific models will be generated during fitting process as required 
model_p = poly6
model_g = poly6 + G_peak
model_d = poly6 + G_peak + D_peak
model_dd = poly6 + G_peak + DD_peak
model_main = poly6 + G_peak + D_peak + DD_peak
model_full = poly6 + G_peak + D_peak + DD_peak + DG_peak

###############################################################################
    #fitting cycle and output decisions#
###############################################################################

#complete the fitting cycle then increase the row counter by one to iterate over each row until all data has been fitted
#row counter used to enable progress monitoring and time estimates
row_counter = 1

while row_counter < rank:
    start_time = timer()
    y = data[row_counter,:]
    errors = error_array[row_counter-1,:]
    print('--------')
    
    bkgdy = list()                              #generate new array containing y data of baseline only excluding peak data, for use in poly guess
    for value in data[row_counter,:xd1]:
        bkgdy.append(value)
    for value in data[row_counter,xg2:xe1]:
        bkgdy.append(value)
    for value in data[row_counter,xf2:]:
        bkgdy.append(value)
    y_background = np.array(bkgdy, dtype = float)
    
    pars_poly6 = poly6.guess(y_background, x=x_background)              #the fit parameters are guessed from the data
    out_poly = model_p.fit(y, pars_poly6, x=x, weights=(1/errors))    #lmfit module completes least squared regression to fit polynomial model
    
#create new parameters and fit the background + G peak
    pars_G = poly6.guess(y_background, x=x_background)
    pars_G.update(G_peak.guess(y[xg1:xg2], x=x[xg1:xg2]))
    pars_G['G_center'].set(1580, min=1550, max=1620)        #parameters for position of 'G' peak limited by definition of peak
    pars_G['G_amplitude'].set(min=1.e-09)   #parameter amplitude set to initial guess based on integral of peak area
    pars_G['G_sigma'].set(max=100)
    out_G = model_g.fit(y, pars_G, x=x, weights=(1/errors))
     
    if ((out_G.redchi)/(out_poly.redchi)) < 0.975:           #This conditional progression checks the G peak is valid
        pars_D = poly6.guess(y_background, x=x_background)
        pars_D.update(G_peak.guess(y[xg1:xg2], x=x[xg1:xg2]))                   #guess parameters for peak found within G range
        pars_D['G_center'].set(1580, min=1550, max=1620)                              #parameters for position of 'G' peak limited by definition of peak
        pars_D['G_amplitude'].set(min=1.e-09)                                   #peak must be positive, but can be very small
        pars_D['G_sigma'].set(max=100)                                          #peak width limited to only real signals
        pars_D.update(D_peak.guess(y[xd1:xg1], x=x[xd1:xg1]))                   #guess parameters for peak found within D range
        pars_D['D_center'].set(1350, min=1300, max=1380)                              #parameters for position of 'G' peak limited by definition of peak
        pars_D['D_amplitude'].set(min=1.e-09)                                   #peak must be positive, but can be very small
        pars_D['D_sigma'].set(50, max=250)                                          #peak width limited to only real signals
        out_D = model_d.fit(y, pars_D, x=x, weights=(1/errors))
        
        #define a function to automatically extract the fitted outcomes from out_D, and use them to create new input parameters for another fit model
        def Dparameters(pars):
            pars['G_center'].set(par_extract(out_D,'G_center'), min=1550, max=1620)
            pars['G_amplitude'].set(par_extract(out_D,'G_amplitude'), min=1.e-09)
            pars['G_sigma'].set(par_extract(out_D,'G_sigma'), max=100)
            pars['D_center'].set(par_extract(out_D,'D_center'), min=1300, max=1380)
            pars['D_amplitude'].set(par_extract(out_D,'D_amplitude'), min=1.e-09)
            pars['D_sigma'].set(par_extract(out_D,'D_sigma'), max=250)
        
        #proceed to fit the 2D peak after establishing output from D & G model
        pars_DD = pars_G
        pars_DD.update(DD_peak.guess(y[xe1:xe2], x=x[xe1:xe2]))
        pars_DD['DD_center'].set(min=2660, max=2730)
        pars_DD['DD_amplitude'].set(min=1.e-09)
        pars_DD['DD_sigma'].set(max=150)
        out_DD = model_dd.fit(y, pars_DD, x=x, weights=(1/errors))
        
        if ((out_D.redchi)/(out_G.redchi)) < 0.975 or ((out_DD.redchi)/(out_G.redchi)) < 0.975:      #if either the D or 2D peaks made a difference, check for full model.
            pars_main = pars_D
            pars_main.update(DD_peak.guess(y[xe1:xe2], x=x[xe1:xe2]))
            pars_main['DD_center'].set(2700, min=2660, max=2730)
            pars_main['DD_amplitude'].set(min=1.e-09)
            pars_main['DD_sigma'].set(max=150)
            out_main = model_main.fit(y, pars_main, x=x, weights=(1/errors))
            
            #define a function to automatically extract the fitted outcomes in the D & G range from output, and use them to create new input parameters for another fit model
            def fullDparameters(pars, output, bkgd, fixed):
                pars['G_center'].set(par_extract(output,'G_center'), min=1550, max=1620, vary=fixed)
                pars['G_amplitude'].set(par_extract(output,'G_amplitude'), min=1.e-09, vary=fixed)
                pars['G_sigma'].set(par_extract(output,'G_sigma'), max=100, vary=fixed)
                pars['D_center'].set(par_extract(output,'D_center'), min=1300, max=1380, vary=fixed)
                pars['D_amplitude'].set(par_extract(output,'D_amplitude'), min=1.e-09, vary=fixed)
                pars['D_sigma'].set(par_extract(output,'D_sigma'), max=250, vary=fixed)
                pars['Background_c0'].set(par_extract(bkgd,'Background_c0'), vary=False)
                pars['Background_c1'].set(par_extract(bkgd,'Background_c1'), vary=False)
                pars['Background_c2'].set(par_extract(bkgd,'Background_c2'), vary=False)
                pars['Background_c3'].set(par_extract(bkgd,'Background_c3'), vary=False)
                pars['Background_c4'].set(par_extract(bkgd,'Background_c4'), vary=False)
                pars['Background_c5'].set(par_extract(bkgd,'Background_c5'), vary=False)
                pars['Background_c6'].set(par_extract(bkgd,'Background_c6'), vary=False)
            
            if Dpflag == True:
                def Dpfixed(pars, output):
                    pars['Dp_center'].set(par_extract(output,'Dp_center'), vary = False)
                    pars['Dp_amplitude'].set(par_extract(output,'Dp_amplitude'), vary = False)
                    pars['Dp_sigma'].set(par_extract(output,'Dp_sigma'), vary = False)
                model_DDp = poly6 + G_peak + D_peak + Dp_peak
                model_mainDp = poly6 + G_peak + D_peak + DD_peak + Dp_peak
                model_fullDp = poly6 + G_peak + D_peak + DD_peak + DG_peak + Dp_peak
        
            if ((out_main.redchi)/(out_D.redchi)) < 0.975 and ((out_main.redchi)/(out_DD.redchi)) < 0.975 and DGflag == True:
                smodel_dg = poly2 + DD_peak + DG_peak                           #generate a model to fit only in the region of 2D and D+G
                pars_sDD = smodel_dg.make_params()
                pars_sDD.update(poly2.guess(y[xp1:], x=x[xp1:]))                #guess parameters for the subrange above 2000 cm-1
                pars_sDD.update(DD_peak.guess(y[xe1:xe2], x=x[xe1:xe2]))
                pars_sDD.update(DG_peak.guess(y[xe2:xf2], x=x[xe2:xf2]))
                out_sDD = smodel_dg.fit(y[xp1:], pars_sDD, x=x[xp1:])           #fit the subrange with a reduced model. This fitting is much faster than the full model.
                #use previous outputs as inputs to full model to significantly increase the speed of fitting
                pars_full = model_full.make_params()
                pars_full.update(poly6.guess(y_background, x=x_background))
                Dparameters(pars_full)
                pars_full['DD_center'].set(par_extract(out_sDD,'DD_center'), min=2660, max=2730)
                pars_full['DD_amplitude'].set(par_extract(out_sDD,'DD_amplitude'), min=1.e-09)
                pars_full['DD_sigma'].set(par_extract(out_sDD,'DD_sigma'), max=150)
                pars_full['DG_center'].set(par_extract(out_sDD,'DG_center'), min=2850, max=2950)
                pars_full['DG_amplitude'].set(par_extract(out_sDD,'DG_amplitude'), min=1.e-09)
                pars_full['DG_sigma'].set(par_extract(out_sDD,'DG_sigma'), max=100)
                out_full = model_full.fit(y, pars_full, x=x, weights=(1/errors))
                
                if ((out_full.redchi)/(out_main.redchi)) < 0.975 and Dpflag == True:
                    #fit D' peak, using parameters from previous fitting in the reduced range x:2250 with fixed polynomial coefficients so 2D and D+G will not be significantly affected
                    pars_Dp = model_DDp.make_params()
                    fullDparameters(pars_Dp, out_full, out_full, True)                  #fullDparameters extarcted the D and G outputs to use as inputs and also extracted the polynomial coefficients and locks them
                    pars_Dp.update(Dp_peak.guess(y[xg3:xg2], x=x[xg3:xg2]))
                    pars_Dp['Dp_center'].set(1620, min=1600, max=1640)
                    pars_Dp['Dp_amplitude'].set(min=1.e-09)
                    pars_Dp['Dp_sigma'].set(max=75)
                    out_Dp = model_DDp.fit(y[:xp2], pars_Dp, x=x[:xp2])
                    
                    pars_fullDp = model_fullDp.make_params()
                    fullDparameters(pars_fullDp, out_Dp, out_full, False)
                    Dpfixed(pars_fullDp, out_Dp)
                    pars_fullDp['DD_center'].set(par_extract(out_full,'DD_center'), vary = False)
                    pars_fullDp['DD_amplitude'].set(par_extract(out_full,'DD_amplitude'), vary = False)
                    pars_fullDp['DD_sigma'].set(par_extract(out_full,'DD_sigma'), vary = False)
                    pars_fullDp['DG_center'].set(par_extract(out_full,'DG_center'), vary = False)
                    pars_fullDp['DG_amplitude'].set(par_extract(out_full,'DG_amplitude'), vary = False)
                    pars_fullDp['DG_sigma'].set(par_extract(out_full,'DG_sigma'), vary = False)
                    out_fullDp = model_fullDp.fit(y, pars_fullDp, x=x, weights=(1/errors))
                    Dpfull_redchi = float(out_fullDp.chisqr / (out_fullDp.ndata - 20))                    

                    if ((Dpfull_redchi)/(out_full.redchi)) < 0.975:       #D' shown to be a valid peak, output from Dp in relevant range and full for 2D and D+G peaks
                        full_result = {}
                        if 'reports' in printflag:
                            print(out_Dp.fit_report(show_correl=False))
                            print(out_full.fit_report(show_correl=False))
                        try:
                            Goutput(full_result, out_Dp)
                            Doutput(full_result, out_Dp)
                            Dpoutput(full_result, out_Dp)
                            DDoutput(full_result, out_full)
                            DGoutput(full_result, out_full)
                            full_result.update({'Red Chi2':(Dpfull_redchi)})
                            results.append(full_result)
                        except ValueError:
                            null_fit_count += 1
                            null_fit_list.append(row_counter)  
                        row_end(row_counter,start_time)
                        row_counter += 1
                    
                    else:                               #D' fit was no better than just D & G, so return only full output
                        full_result = {}
                        full_result.update(blankDp)
                        if 'reports' in printflag:
                            print(out_full.fit_report(show_correl=False))
                        try:
                            Goutput(full_result, out_full)
                            Doutput(full_result, out_full)
                            DDoutput(full_result, out_full)
                            DGoutput(full_result, out_full)
                            full_result.update({'Red Chi2':(out_full.redchi)})
                            results.append(full_result)
                        except ValueError:
                            null_fit_count += 1
                            null_fit_list.append(row_counter)
                        row_end(row_counter,start_time)
                        row_counter += 1
                    
                elif ((out_full.redchi)/(out_main.redchi)) < 0.975 and Dpflag == False:
                    full_result = {}
                    if 'reports' in printflag:
                        print(out_full.fit_report(show_correl=False))
                    try:
                        Goutput(full_result, out_full)
                        Doutput(full_result, out_full)
                        DDoutput(full_result, out_full)
                        DGoutput(full_result, out_full)
                        full_result.update({'Red Chi2':(out_full.redchi)})
                        results.append(full_result)
                    except ValueError:
                        null_fit_count += 1
                        null_fit_list.append(row_counter)
                    row_end(row_counter,start_time)
                    row_counter += 1
                    
                elif ((out_full.redchi)/(out_main.redchi)) > 0.975 and Dpflag == True:
                    #DG peak was not valid, attempt to fit D' to D & G only region using main output as starting point for 2D parameters
                    pars_Dp = model_DDp.make_params()
                    fullDparameters(pars_Dp, out_main, out_main, True)
                    pars_Dp.update(Dp_peak.guess(y[xg3:xg2], x=x[xg3:xg2]))
                    pars_Dp['Dp_center'].set(1620, min=1600, max=1640)
                    pars_Dp['Dp_amplitude'].set(min=1.e-09)
                    pars_Dp['Dp_sigma'].set(max=75)
                    out_Dp = model_DDp.fit(y[:xp2], pars_Dp, x=x[:xp2])
                    
                    pars_mainDp = model_mainDp.make_params()
                    fullDparameters(pars_mainDp, out_Dp, out_main, False)
                    Dpfixed(pars_mainDp, out_Dp)
                    pars_mainDp['DD_center'].set(par_extract(out_main,'DD_center'), vary = False)
                    pars_mainDp['DD_amplitude'].set(par_extract(out_main,'DD_amplitude'), vary = False)
                    pars_mainDp['DD_sigma'].set(par_extract(out_main,'DD_sigma'), vary = False)
                    out_Dpmain = model_mainDp.fit(y, pars_mainDp, x=x, weights=(1/errors))
                    Dpmain_redchi = float(out_Dpmain.chisqr / (out_Dpmain.ndata - 17))                    
                    
                    if ((Dpmain_redchi)/(out_main.redchi)) < 0.975:       #D' shown to be a valid peak, output from Dp in relevant range and full for 2D and D+G peaks
                        main_result = {}
                        main_result.update(blankDG)
                        if 'reports' in printflag:
                            print(out_Dp.fit_report(show_correl=False))
                            print(out_main.fit_report(show_correl=False))
                        try:
                            Goutput(main_result, out_Dp)
                            Doutput(main_result, out_Dp)
                            Dpoutput(main_result, out_Dp)
                            DDoutput(main_result, out_main)
                            main_result.update({'Red Chi2':(Dpmain_redchi)}) 
                            results.append(main_result)
                        except ValueError:
                            null_fit_count += 1
                            null_fit_list.append(row_counter)  
                        row_end(row_counter,start_time)
                        row_counter += 1
                    
                    else:                               #D' fit was no better than just D & G, so return only full output
                        main_result = {}
                        main_result.update(blankDG)
                        main_result.update(blankDp)
                        if 'reports' in printflag:
                            print(out_main.fit_report(show_correl=False))
                        try:
                            Goutput(main_result, out_main)
                            Doutput(main_result, out_main)
                            DDoutput(main_result, out_main)
                            main_result.update({'Red Chi2':(out_main.redchi)})
                            results.append(main_result)
                        except ValueError:
                            null_fit_count += 1
                            null_fit_list.append(row_counter)
                        row_end(row_counter,start_time)
                        row_counter += 1
                                    
                else:                                   #D+G peak was not valid, user specified not D' so output the main model only with no D' fitting
                    main_result = {}
                    main_result.update(blankDG)
                    if 'reports' in printflag:
                        print(out_main.fit_report(show_correl=False))
                    try:
                        Goutput(main_result, out_main)
                        Doutput(main_result, out_main)
                        DDoutput(main_result, out_main)
                        main_result.update({'Red Chi2':(out_main.redchi)})
                        results.append(main_result)
                    except ValueError:
                        null_fit_count += 1
                        null_fit_list.append(row_counter)
                    row_end(row_counter,start_time)
                    row_counter += 1
            
            elif ((out_main.redchi)/(out_D.redchi)) < 0.975 and ((out_main.redchi)/(out_DD.redchi)) < 0.975 and DGflag == False and Dpflag == True:
                #main model valid and D' needs to be fitted, use same split method as above
                pars_Dp = model_DDp.make_params()
                fullDparameters(pars_Dp, out_main, out_main,  True)
                pars_Dp.update(Dp_peak.guess(y[xg3:xg2], x=x[xg3:xg2]))
                pars_Dp['Dp_center'].set(1620, min=1600, max=1640)
                pars_Dp['Dp_amplitude'].set(min=1.e-09)
                pars_Dp['Dp_sigma'].set(max=75)
                out_Dp = model_DDp.fit(y[:xp2], pars_Dp, x=x[:xp2])
                
                pars_mainDp = model_mainDp.make_params()
                fullDparameters(pars_mainDp, out_Dp, out_main, False)
                Dpfixed(pars_mainDp, out_Dp)
                pars_mainDp['DD_center'].set(par_extract(out_main,'DD_center'), vary = False)
                pars_mainDp['DD_amplitude'].set(par_extract(out_main,'DD_amplitude'), vary = False)
                pars_mainDp['DD_sigma'].set(par_extract(out_main,'DD_sigma'), vary = False)
                out_Dpmain = model_mainDp.fit(y, pars_mainDp, x=x, weights=(1/errors))
                Dpmain_redchi = float(out_Dpmain.chisqr / (out_Dpmain.ndata - 17))                    
                        
                if ((Dpmain_redchi)/(out_main.redchi)) < 0.975:       #D' shown to be a valid peak, output from Dp in relevant range and full for 2D and D+G peaks
                    main_result = {}
                    if 'reports' in printflag:
                        print(out_Dp.fit_report(show_correl=False))
                        print(out_main.fit_report(show_correl=False))
                    try:
                        Goutput(main_result, out_Dp)
                        Doutput(main_result, out_Dp)
                        Dpoutput(main_result, out_Dp)
                        DDoutput(main_result, out_main)
                        main_result.update({'Red Chi2':(Dpmain_redchi)})
                        results.append(main_result)
                    except ValueError:
                        null_fit_count += 1
                        null_fit_list.append(row_counter)  
                    row_end(row_counter,start_time)
                    row_counter += 1
                    
                else:                               #D' fit was no better than just D & G, so return only full output
                    main_result = {}
                    main_result.update(blankDp)
                    if 'reports' in printflag:
                        print(out_main.fit_report(show_correl=False))
                    try:
                        Goutput(main_result, out_main)
                        Doutput(main_result, out_main)
                        DDoutput(main_result, out_main)
                        main_result.update({'Red Chi2':(out_main.redchi)})
                        results.append(main_result)
                    except ValueError:
                        null_fit_count += 1
                        null_fit_list.append(row_counter)
                    row_end(row_counter,start_time)
                    row_counter += 1
            
            elif ((out_main.redchi)/(out_D.redchi)) < 0.975 and ((out_main.redchi)/(out_DD.redchi)) < 0.975 and DGflag == False and Dpflag == False:
                main_result = {}
                if 'reports' in printflag:
                    print(out_main.fit_report(show_correl=False))
                try:
                    Goutput(main_result, out_main)
                    Doutput(main_result, out_main)
                    DDoutput(main_result, out_main)
                    main_result.update({'Red Chi2':(out_main.redchi)})
                    results.append(main_result)
                except ValueError:
                    null_fit_count += 1
                    null_fit_list.append(row_counter)
                row_end(row_counter,start_time)
                row_counter += 1                                                            #moves the fitting onto the next row
    
            elif ((out_D.redchi)/(out_G.redchi)) < 0.975 and DGflag == True:                #if only the D peak is real, but the user DG flag is set then double check a full fit
                smodel_DG = poly2 + DD_peak + DG_peak                           #generate a model to fit only in the region of 2D and D+G
                pars_sDD = smodel_DG.make_params()
                pars_sDD.update(poly2.guess(y_background[:(xd1+(xp2-xg2))], x=x_background[:(xd1+(xp2-xg2))]))                #guess parameters for the subrange above 2000 cm-1
                pars_sDD.update(DD_peak.guess(y[xe1:xe2], x=x[xe1:xe2]))
                pars_sDD.update(DG_peak.guess(y[xe2:xf2], x=x[xe2:xf2]))
                out_sDD = smodel_DG.fit(y[xp1:], pars_sDD, x=x[xp1:])           #fit the subrange with a reduced model. This fitting is much faster than the full model.
                #use previous outputs as inputs to full model to significantly increase the speed of fitting
                pars_full = model_full.make_params()
                pars_full.update(poly6.guess(y_background, x=x_background))
                Dparameters(pars_full)
                pars_full['DD_center'].set(par_extract(out_sDD,'DD_center'), min=2660, max=2730)
                pars_full['DD_amplitude'].set(par_extract(out_sDD,'DD_amplitude'), min=1.e-09)
                pars_full['DD_sigma'].set(par_extract(out_sDD,'DD_sigma'), max=150)
                pars_full['DG_center'].set(par_extract(out_sDD,'DG_center'), min=2850, max=2950)
                pars_full['DG_amplitude'].set(par_extract(out_sDD,'DG_amplitude'), min=1.e-09)
                pars_full['DG_sigma'].set(par_extract(out_sDD,'DG_sigma'), max=100)
                out_full = model_full.fit(y, pars_full, x=x, weights=(1/errors))
        
                if ((out_full.redchi)/(out_main.redchi)) < 0.975 and ((out_full.redchi)/(out_D.redchi)) < 0.975 and Dpflag == True:     #if forcing the D+G improves the fit of the 2D so the full model is an imporovment then continue with that
                    #fit D' peak, using parameters from previous fitting in the reduced range x:2250 with fixed polynomial coefficients so 2D and D+G will not be significantly affected
                    pars_Dp = model_DDp.make_params()
                    fullDparameters(pars_Dp, out_full, out_full, True)                  #fullDparameters extarcted the D and G outputs to use as inputs and also extracted the polynomial coefficients and locks them
                    pars_Dp.update(Dp_peak.guess(y[xg3:xg2], x=x[xg3:xg2]))
                    pars_Dp['Dp_center'].set(1620, min=1600, max=1640)
                    pars_Dp['Dp_amplitude'].set(min=1.e-09)
                    pars_Dp['Dp_sigma'].set(max=75)
                    out_Dp = model_DDp.fit(y[:xp2], pars_Dp, x=x[:xp2])
                    
                    pars_fullDp = model_fullDp.make_params()
                    fullDparameters(pars_fullDp, out_Dp, out_full, False)
                    Dpfixed(pars_fullDp, out_Dp)
                    pars_fullDp['DD_center'].set(par_extract(out_full,'DD_center'), vary = False)
                    pars_fullDp['DD_amplitude'].set(par_extract(out_full,'DD_amplitude'), vary = False)
                    pars_fullDp['DD_sigma'].set(par_extract(out_full,'DD_sigma'), vary = False)
                    pars_fullDp['DG_center'].set(par_extract(out_full,'DG_center'), vary = False)
                    pars_fullDp['DG_amplitude'].set(par_extract(out_full,'DG_amplitude'), vary = False)
                    pars_fullDp['DG_sigma'].set(par_extract(out_full,'DG_sigma'), vary = False)
                    out_Dpfull = model_fullDp.fit(y, pars_fullDp, x=x, weights=(1/errors))
                    Dpfull_redchi = float(out_Dpfull.chisqr / (out_Dpfull.ndata - 20))                    

                    if ((Dpfull_redchi)/(out_full.redchi)) < 0.975:       #D' shown to be a valid peak, output from Dp in relevant range and full for 2D and D+G peaks
                        full_result = {}
                        if 'reports' in printflag:
                            print(out_Dp.fit_report(show_correl=False))
                            print(out_full.fit_report(show_correl=False))
                        try:
                            Goutput(full_result, out_Dp)
                            Doutput(full_result, out_Dp)
                            Dpoutput(full_result, out_Dp)
                            DDoutput(full_result, out_full)
                            DGoutput(full_result, out_full)
                            full_result.update({'Red Chi2':(Dpfull_redchi)})
                            results.append(full_result)
                        except ValueError:
                            null_fit_count += 1
                            null_fit_list.append(row_counter)  
                        row_end(row_counter,start_time)
                        row_counter += 1
                    
                    else:                               #D' fit was no better than just D & G, so return only full output
                        full_result = {}
                        full_result.update(blankDp)
                        if 'reports' in printflag:
                            print(out_full.fit_report(show_correl=False))
                        try:
                            Goutput(full_result, out_full)
                            Doutput(full_result, out_full)
                            DDoutput(full_result, out_full)
                            DGoutput(full_result, out_full)
                            full_result.update({'Red Chi2':(out_full.redchi)})
                            results.append(full_result)
                        except ValueError:
                            null_fit_count += 1
                            null_fit_list.append(row_counter)
                        row_end(row_counter,start_time)
                        row_counter += 1
                        
                elif ((out_full.redchi)/(out_main.redchi)) < 0.975 and ((out_full.redchi)/(out_D.redchi)) < 0.975 and Dpflag == False:  #if forcing the D+G produced a good fit on the 2D as well then output that. No D' fitting
                    full_result = {}
                    if 'reports' in printflag:
                        print(out_full.fit_report(show_correl=False))
                    try:
                        Goutput(full_result, out_full)
                        Doutput(full_result, out_full)
                        DDoutput(full_result, out_full)
                        DGoutput(full_result, out_full)
                        full_result.update({'Red Chi2':(out_full.redchi)})
                        results.append(full_result)
                    except ValueError:
                        null_fit_count += 1
                        null_fit_list.append(row_counter)
                    row_end(row_counter,start_time)
                    row_counter += 1
                
                elif (((out_full.redchi)/(out_main.redchi)) > 0.975 or ((out_full.redchi)/(out_D.redchi)) > 0.975) and Dpflag == True:
                    pars_Dp = model_DDp.make_params()
                    Dparameters(pars_Dp)
                    pars_Dp.update(poly6.guess(y_background, x=x_background))
                    pars_Dp.update(Dp_peak.guess(y[xg3:xg2], x=x[xg3:xg2]))
                    pars_Dp['Dp_center'].set(1620, min=1600, max=1640)
                    pars_Dp['Dp_amplitude'].set(min=1.e-09)
                    pars_Dp['Dp_sigma'].set(max=75)
                    out_Dp = model_DDp.fit(y[:xp2], pars_Dp, x=x[:xp2])
                
                    if ((out_Dp.redchi)/(out_D.redchi)) < 0.975:
                        Dp_result = {}
                        Dp_result.update(blankDD)
                        Dp_result.update(blankDG)
                        if 'reports' in printflag:
                            print(out_Dp.fit_report(show_correl=False))
                        try:
                            Goutput(Dp_result, out_Dp)
                            Doutput(Dp_result, out_Dp)
                            Dpoutput(Dp_result, out_Dp)
                            Dp_result.update({'Red Chi2':(out_Dp.redchi)})
                            results.append(Dp_result)
                        except ValueError:
                            null_fit_count += 1
                            null_fit_list.append(row_counter)
                        row_end(row_counter, start_time)
                        row_counter += 1
                
                    else:               #2D not valid, but D' was requested so return output from small range with n=2 background. Consistent with D' fitting output
                        Dp_result = {}
                        Dp_result.update(blankDD)
                        Dp_result.update(blankDp)
                        Dp_result.update(blankDG)
                        if 'reports' in printflag:
                            print(out_D.fit_report(show_correl=False))
                        try:
                            Goutput(Dp_result, out_D)
                            Doutput(Dp_result, out_D)
                            Dp_result.update({'Red Chi2':(out_D.redchi)})
                            results.append(Dp_result)
                        except ValueError:
                            null_fit_count += 1
                            null_fit_list.append(row_counter)
                        row_end(row_counter, start_time)
                        row_counter += 1
                
                else:                                   #if the full fit is not better and the Dp flag is not set then just return the D parameters with 2D and DG set to no
                    D_result = {}
                    D_result.update(blankDD)
                    D_result.update(blankDG)
                    if 'reports' in printflag:
                        print(out_D.fit_report(show_correl=False))
                    try:
                        Goutput(D_result, out_D)
                        Doutput(D_result, out_D)
                        D_result.update({'Red Chi2':(out_D.redchi)})
                        results.append(D_result)
                    except ValueError:
                        null_fit_count += 1
                        null_fit_list.append(row_counter)
                    row_end(row_counter,start_time)
                    row_counter += 1
            
            elif ((out_D.redchi)/(out_G.redchi)) < 0.975 and DGflag == False and Dpflag == True:
                pars_Dp = model_DDp.make_params()
                Dparameters(pars_Dp)
                pars_Dp.update(poly6.guess(y_background, x=x_background))
                pars_Dp.update(Dp_peak.guess(y[xg3:xg2], x=x[xg3:xg2]))
                pars_Dp['Dp_center'].set(1620, min=1600, max=1640)
                pars_Dp['Dp_amplitude'].set(min=1.e-09)
                pars_Dp['Dp_sigma'].set(max=75)
                out_Dp = model_DDp.fit(y[:xp2], pars_Dp, x=x[:xp2])
                
                if ((out_Dp.redchi)/(out_D.redchi)) < 0.975:
                    Dp_result = {}
                    Dp_result.update(blankDD)
                    if 'reports' in printflag:
                        print(out_Dp.fit_report(show_correl=False))
                    try:
                        Goutput(Dp_result, out_Dp)
                        Doutput(Dp_result, out_Dp)
                        Dpoutput(Dp_result, out_Dp)
                        Dp_result.update({'Red Chi2':(out_Dp.redchi)})
                        results.append(Dp_result)
                    except ValueError:
                        null_fit_count += 1
                        null_fit_list.append(row_counter)
                    row_end(row_counter, start_time)
                    row_counter += 1
                
                else:               #2D not valid, but D' was requested so return output from small range with n=2 background. Consistent with D' fitting output
                    Dp_result = {}
                    Dp_result.update(blankDD)
                    Dp_result.update(blankDp)
                    if 'reports' in printflag:
                        print(out_D.fit_report(show_correl=False))
                    try:
                        Goutput(Dp_result, out_D)
                        Doutput(Dp_result, out_D)
                        Dp_result.update({'Red Chi2':(out_D.redchi)})
                        results.append(Dp_result)
                    except ValueError:
                        null_fit_count += 1
                        null_fit_list.append(row_counter)
                    row_end(row_counter, start_time)
                    row_counter += 1
            
            elif ((out_D.redchi)/(out_G.redchi)) < 0.975 and DGflag == False and Dpflag == False:
                D_result = {}
                D_result.update(blankDD)
                if 'reports' in printflag:
                    print(out_D.fit_report(show_correl=False))
                try:
                    Goutput(D_result, out_D)
                    Doutput(D_result, out_D)
                    D_result.update({'Red Chi2':(out_D.redchi)})
                    results.append(D_result)
                except ValueError:
                    null_fit_count += 1
                    null_fit_list.append(row_counter)
                row_end(row_counter, start_time)
                row_counter += 1
          
            elif ((out_DD.redchi)/(out_G.redchi)) < 0.975 and DGflag == True and Dpflag == True:
                DD_result = {}
                DD_result.update(blankD)
                DD_result.update(blankDp)
                DD_result.update(blankDG)
                if 'reports' in printflag:
                    print(out_DD.fit_report(show_correl=False))
                try:
                    Goutput(DD_result, out_DD)
                    DDoutput(DD_result, out_DD)
                    DD_result.update({'Red Chi2':(out_DD.redchi)})
                    results.append(DD_result)
                except:
                    null_fit_count += 1
                    null_fit_list.append(row_counter)
                row_end(row_counter, start_time)
                row_counter += 1
            
            elif ((out_DD.redchi)/(out_G.redchi)) < 0.975 and DGflag == True and Dpflag == False:
                DD_result = {}
                DD_result.update(blankD)
                DD_result.update(blankDG)
                if 'reports' in printflag:
                    print(out_DD.fit_report(show_correl=False))
                try:
                    Goutput(DD_result, out_DD)
                    DDoutput(DD_result, out_DD)
                    DD_result.update({'Red Chi2':(out_DD.redchi)})
                    results.append(DD_result)
                except:
                    null_fit_count += 1
                    null_fit_list.append(row_counter)
                row_end(row_counter, start_time)
                row_counter += 1
                
            elif ((out_DD.redchi)/(out_G.redchi)) < 0.975 and DGflag == False and Dpflag == True:
                DD_result = {}
                DD_result.update(blankD)
                DD_result.update(blankDp)
                if 'reports' in printflag:
                    print(out_DD.fit_report(show_correl=False))
                try:
                    Goutput(DD_result, out_DD)
                    DDoutput(DD_result, out_DD)
                    DD_result.update({'Red Chi2':(out_DD.redchi)})
                    results.append(DD_result)
                except:
                    null_fit_count += 1
                    null_fit_list.append(row_counter)
                row_end(row_counter, start_time)
                row_counter += 1
                
            elif ((out_DD.redchi)/(out_G.redchi)) < 0.975 and DGflag == False and Dpflag == False:
                DD_result = {}
                DD_result.update(blankD)
                if 'reports' in printflag:
                    print(out_DD.fit_report(show_correl=False))
                try:
                    Goutput(DD_result, out_DD)
                    DDoutput(DD_result, out_DD)
                    DD_result.update({'Red Chi2':(out_DD.redchi)})
                    results.append(DD_result)
                except:
                    null_fit_count += 1
                    null_fit_list.append(row_counter)
                row_end(row_counter, start_time)
                row_counter += 1
                
            else:
                print('Something unexpected has rather got in the way')
                row_end(row_counter, start_time)
                row_counter += 1
                
        #neither the 2D or D peaks improved on the G, hence return only the G peak with blank depending on user input
        elif DGflag == True and Dpflag == True:         #if all peaks were requested, return all blank lists
            g_result = {}
            g_result.update(blankD)
            g_result.update(blankDD)
            g_result.update(blankDp)
            g_result.update(blankDG)
            if 'reports' in printflag:
                print(out_G.fit_report(show_correl=False))
            try:
                Goutput(g_result, out_G)
                g_result.update({'Red Chi2':(out_G.redchi)})
                results.append(g_result)
            except ValueError:
                null_fit_count += 1
                null_fit_list.append(row_counter)
            row_end(row_counter,start_time)
            row_counter += 1                                                            #moves the fitting onto the next row
        
        elif DGflag == True and Dpflag == False:        #if the D+G peak was requested, return this blank
            g_result = {}
            g_result.update(blankD)
            g_result.update(blankDD)
            g_result.update(blankDG)
            if 'reports' in printflag:
                print(out_G.fit_report(show_correl=False))
            try:
                Goutput(g_result, out_G)
                g_result.update({'Red Chi2':(out_G.redchi)})
                results.append(g_result)
            except ValueError:
                null_fit_count += 1
                null_fit_list.append(row_counter)
            row_end(row_counter,start_time)
            row_counter += 1                                                            #moves the fitting onto the next row
        
        elif DGflag == False and Dpflag == True:
            g_result = {}
            g_result.update(blankD)
            g_result.update(blankDD)
            g_result.update(blankDp)
            if 'reports' in printflag:
                print(out_G.fit_report(show_correl=False))
            try:
                Goutput(g_result, out_G)
                g_result.update({'Red Chi2':(out_G.redchi)})
                results.append(g_result)
            except ValueError:
                null_fit_count += 1
                null_fit_list.append(row_counter)
            row_end(row_counter,start_time)
            row_counter += 1                                                            #moves the fitting onto the next row
       
        elif DGflag == False and Dpflag == False:       #no extra peaks were requested by user, only return blank defult D and 2D
            g_result = {}
            g_result.update(blankD)
            g_result.update(blankDD)
            if 'reports' in printflag:
                print(out_G.fit_report(show_correl=False))
            try:
                Goutput(g_result, out_G)
                g_result.update({'Red Chi2':(out_G.redchi)})
                results.append(g_result)
            except ValueError:
                null_fit_count += 1
                null_fit_list.append(row_counter)
            row_end(row_counter,start_time)
            row_counter += 1                                                            #moves the fitting onto the next row
        
        else:
            print('Some other error has occured')
            row_end(row_counter,start_time)
            row_counter += 1                                                            #moves the fitting onto the next row
    
        
    else:                       #if G > poly there is no data
        null_spectra_count += 1
        null_spectra_list.append(row_counter)   #Append the row with no G peak to list to be returned at end of fitting
        if 'reports' in printflag:
            print(out_poly.fit_report(show_correl=False))
        if 'null_plots' in printflag:           #if user entered null_plots a plot will be generated for every null spectrum. This will halt the fitting until it is closed
            plt.plot(x, y, 'b', label='Data')
            plt.plot(x, out_poly.best_fit, 'r', label='Fit')
            plt.xlabel("Raman Shift / cm" + u'\u207b\u00b9', fontsize=14)
            plt.ylabel("Counts", fontsize=14)
            plt.legend()
            plt.show()
        row_end(row_counter,start_time)
        row_counter += 1                                                            #moves the fitting onto the next row

print('Rows processed: ' + str(row_counter-1))      #Confirm the size of data array processed

if null_spectra_count == 0:                         #If all spectrum were valid, confirm no null spectra
    print('Percentage null spectra: 0%')
else:                                               #If null spectra were present, confirm the number and rows
    print('Null Spectra: ' + str(null_spectra_count))
    print(null_spectra_list)

if null_fit_count == 0:                             #If no errors is reading fit output strings, confirm good fits
    print('Percentage bad fit: 0%')
else:                                               #If errors were returned when reading output string, confirm number and rows
    print('Bad fit: ' + str(null_fit_count))
    print(null_fit_list)
    
print('******************************************************')

numberadded = 1                                                 #establish number to add to file name in event of file already existing
file_check = False                                              #flag set to False, loop to continue in event of file name existing 
param_name = str(shortname + '-Parameters')                            #produce histogram output file name
while file_check == False:                                      #False means the file name already exists and a number will be added to prevent overwriting
    if (os.path.isfile(param_name + '.csv')) == True:
        param_name = str(shortname + str(numberadded) + '-Parameters')
        numberadded += 1
    else:                                                       #if the os.file check does not return True then the file name isn't already used
        file_check = True

with open(param_name + '.csv', 'w', newline='') as fout:    #create the results file as a .csv
    if DGflag == True and Dpflag == True:
        cout = csv.DictWriter(fout, ['G_center','G_center_error','G_height','G_height_error','G_fwhm','G_fwhm_error','G_amplitude','G_amplitude_error','D_center','D_center_error','D_height','D_height_error','D_fwhm','D_fwhm_error','D_amplitude','D_amplitude_error','ID/IG','ID/IG_error','DD_center','DD_center_error','DD_height','DD_height_error','DD_fwhm','DD_fwhm_error','DD_amplitude','DD_amplitude_error','IDD/IG','IDD/IG_error','Dp_center','Dp_center_error','Dp_height','Dp_height_error','Dp_fwhm','Dp_fwhm_error','Dp_amplitude','Dp_amplitude_error','DG_center','DG_center_error','DG_height','DG_height_error','DG_fwhm','DG_fwhm_error','DG_amplitude','DG_amplitude_error','Red Chi2','SignalNoise'])
    elif DGflag == True and Dpflag == False:
        cout = csv.DictWriter(fout, ['G_center','G_center_error','G_height','G_height_error','G_fwhm','G_fwhm_error','G_amplitude','G_amplitude_error','D_center','D_center_error','D_height','D_height_error','D_fwhm','D_fwhm_error','D_amplitude','D_amplitude_error','ID/IG','ID/IG_error','DD_center','DD_center_error','DD_height','DD_height_error','DD_fwhm','DD_fwhm_error','DD_amplitude','DD_amplitude_error','IDD/IG','IDD/IG_error','DG_center','DG_center_error','DG_height','DG_height_error','DG_fwhm','DG_fwhm_error','DG_amplitude','DG_amplitude_error','Red Chi2','SignalNoise'])
    elif DGflag == False and Dpflag == True:
        cout = csv.DictWriter(fout, ['G_center','G_center_error','G_height','G_height_error','G_fwhm','G_fwhm_error','G_amplitude','G_amplitude_error','D_center','D_center_error','D_height','D_height_error','D_fwhm','D_fwhm_error','D_amplitude','D_amplitude_error','ID/IG','ID/IG_error','DD_center','DD_center_error','DD_height','DD_height_error','DD_fwhm','DD_fwhm_error','DD_amplitude','DD_amplitude_error','IDD/IG','IDD/IG_error','Dp_center','Dp_center_error','Dp_height','Dp_height_error','Dp_fwhm','Dp_fwhm_error','Dp_amplitude','Dp_amplitude_error','Red Chi2','SignalNoise'])
    else:
        cout = csv.DictWriter(fout, ['G_center','G_center_error','G_height','G_height_error','G_fwhm','G_fwhm_error','G_amplitude','G_amplitude_error','D_center','D_center_error','D_height','D_height_error','D_fwhm','D_fwhm_error','D_amplitude','D_amplitude_error','ID/IG','ID/IG_error','DD_center','DD_center_error','DD_height','DD_height_error','DD_fwhm','DD_fwhm_error','DD_amplitude','DD_amplitude_error','IDD/IG','IDD/IG_error','Red Chi2','SignalNoise'])
    cout.writeheader()          #convert the list on the previous line into headers for the results table
    cout.writerows(results)     #Because the headers and dictionary keys match, each list of dictionaries has the values converted into a line in the results table with every value in the correct column

#further analysis to plot fitted data as a histogram
# hist_x = D/G ratios from fitting output. Since they are stored as a list of dictionaries the for loop iterates 
#through the list, extracting the values with the dictionary key 'ID/IG'
hist_x = []
for entry in results:
    hist_x.append(float(entry['ID/IG']))

#y = 2D/G ratio
hist_y = []
for entry in results:
    hist_y.append(float(entry['IDD/IG']))

#Generate an initial plot of the 2D histogram using automatic axis ranges. Quick visual check for user.
plt.figure(num=1, frameon=False)
plt.hist2d(hist_x, hist_y, bins=(rank//20), cmin=1)  # cmin=1 specifies all ratios with zero population as white
plt.colorbar()                                  #add the colour bar key
plt.xlabel(r'$I_D / I_G$', fontsize=16)    #specify the x axis title
plt.ylabel(r'$I_{2D} / I_G$', fontsize=16)   #specify the y axis title
#plt.rc('xtick', labelsize='large')            #fontsize of the numbers of the tick labels
#plt.rc('ytick', labelsize='large')            #fontsize of the numbers of the tick labels
plt.show()

plotgoon = input('Parameters saved. Plot and save 3D heatmap histogram? \n[Y/N]')
if str(plotgoon)=='Y' or str(plotgoon)=='y' or str(plotgoon)=='yes' or str(plotgoon)=='YES' or str(plotgoon)=='Yes':
    pass
else:
    sys.exit()

#specify list to hold user input later (mutable list needed rather than immutable string for user input)
saveflag = ['entry']

Dmin = min(hist_x)
Dmax = max(hist_x)
Dbin = (rank//20)
DDmin = min(hist_y)
DDmax = max(hist_y)
DDbin = (rank//20)

#if user enters y histogram plotting is complete and program will export data and terminate. Otherwise re-plotting continues.
while ('y' in saveflag) == False:
    print('-Previous-')
    print('2D/G Min = ' + str(DDmin))
    print('2D/G Max = ' + str(DDmax))
    print('2D/G Bins = ' + str(DDbin))
    print('D/G Min = ' + str(Dmin))
    print('D/G Max = ' + str(Dmax))
    print('D/G Bins = ' + str(Dbin))
    print('-Input-')
    #establish false flags for use in error checking of user inputs.
    DDmin_flag = False
    DDmax_flag = False
    Dmin_flag = False
    Dmax_flag = False
    DDbin_flag = False
    Dbin_flag = False
    
    while DDmin_flag != True:
        try:    #the loop will continue until a valid input is entered and the flag is changed to true
            DDmin = float(input("Specify minimum 2D/G axis range: "))   #select the lowest value for the y axis
            DDmin_flag = True
        except ValueError:  #error raised if non-numerical error entered and loop allows re-entry without crashing
            DDmin_flag = False
            print("Non-numerical value entered. Please re-enter range.")
    while DDmax_flag != True:
        try:
            DDmax = float(input("Specify maximum 2D/G axis range: "))   #select the highest value for the y axis
            if DDmax > DDmin:
                DDmax_flag = True
            else:
                Dmax_flag = False
                print("Ensure 2D max is higher than 2D min.\nEnter and re-plot if needed")
        except ValueError:
            DDmax_flag = False
            print("Non-numerical value entered. Please re-enter range.")
    while DDbin_flag != True:
        try:
            DDbin = int(input("Specify number of 2D/G bins:"))
            DDbin_flag = True
        except ValueError:
            print("Non-numerical value entered. Please re-enter range.")
    while Dmin_flag != True:
        try:
            Dmin = float(input("Specify minimum D/G axis range: "))     # select the lowest value for the x axis
            Dmin_flag = True
        except ValueError:
            Dmin_flag = False
            print("Non-numerical value entered. Please re-enter range.")
    while Dmax_flag != True:
        try:
            Dmax = float(input("Specify maximum D/G axis range: "))     # select the highest value for the x axis
            if Dmax > Dmin:
                Dmax_flag = True
            else:
                Dmax_flag = False
                print("Ensure D max is higher than 2D min.\nEnter and re-plot if needed")
        except ValueError:
            Dmax_flag = False
            print("Non-numerical value entered. Please re-enter range.")
    while Dbin_flag != True:
        try:
            Dbin = int(input("Specify number of D/G bins: "))     # select the highest value for the x axis
            Dbin_flag = True
        except ValueError:
            print("Non-numerical value entered. Please re-enter range.")

    #check the axis values entered include all the data, if not a warning message is displayed. 
    for entry in hist_y:    #iterate over all data
        if entry < DDmin:   #if a data point outside the range is found the 'for' loop is cancelled  
            print("Caution: 2D/G min too high, data below axis")    #print message highlighting error
            break
    for entry in hist_y:
        if entry > DDmax:
            print("Caution: 2D/G max too low, data above axis")
            break
    for entry in hist_x:
        if entry < Dmin:
            print("Caution: D/G min too high, data below axis")
            break
    for entry in hist_x:
        if entry > Dmax:
            print("Caution: D/G max too low, data above axis")
            break
    
    plt.figure(frameon=True)    #plot a new 3D histogram using the axis limits entered by the user above
    plt.hist2d(hist_x, hist_y, bins=[Dbin,DDbin], range=[[Dmin,Dmax],[DDmin,DDmax]], cmin=1)  #range now defined by user input
    plt.colorbar()
    plt.xlabel(r'$I_D / I_G$', fontsize=16)
    plt.ylabel(r'$I_{2D} / I_G$', fontsize=16)
    #plt.rc('xtick', labelsize='large')
    #plt.rc('ytick', labelsize='large')
    plt.tight_layout()
    plt.show()
    
    saveflag[0] = str(input("To save close plot window and type y.\nType n to replot data."))     # user inputs to save image with consistent size
    
    if ('y' in saveflag) == True:   #if user has input save above then the plot is saved and program closed
        plt.figure(figsize=(3.3,2.6), frameon=True, dpi=600)    #when printing a plot the plot is closed, need to re-plot to save image
        plt.hist2d(hist_x, hist_y, bins=[Dbin,DDbin], range=[[Dmin,Dmax],[DDmin,DDmax]], cmin=1)  #range now defined by user input
        plt.colorbar()
        plt.xlabel(r'$I_D / I_G$', fontsize=11)
        plt.ylabel(r'$I_{2D} / I_G$', fontsize=11)
        #plt.rc('xtick', labelsize='large')
        #plt.rc('ytick', labelsize='large')
        file_check = False                                              #flag set to False, loop to continue in event of file name existing 
        image_name = str(shortname + '-3D Heatmap')                            #produce histogram output file name
        numberadded = 1
        while file_check == False:                                      #False means the file name already exists and a number will be added to prevent overwriting
            if (os.path.isfile(image_name + '.png')) == True:
                image_name = str(shortname + str(numberadded) + '-3D Heatmap')
                numberadded += 1
            else:                                                       #if the os.file check does not return True then the file name isn't already used
                file_check = True
        plt.tight_layout()
        plt.savefig(image_name + '.png')
