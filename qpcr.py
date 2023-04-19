import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression
import openpyxl
import sys


"""
standard_curve_plot : function

parameters : avg_CT list of avg CT per triplet 

creates a pyplot of log_copies vs avg CT values that are linearized
through a scikit LinearRegression model

creates an excel sheet and appends the pyplot along with slope/intercept values

returns slope, intercept, r2 of CT vs log copies plot
"""
def standard_curve_plot(avg_CT):
    CT = [i for i in avg_CT if avg_CT.index(i) % 4 == 0] #extrapolating std's averages
    del CT[7] #remove 8th value, since std goes from 1-7
    CT.reverse()
    log_copies = [4,5,6,7,8,9,10]

    x = np.array(log_copies).reshape((-1,1))
    y = np.array(CT)
    model = LinearRegression().fit(x, y) #creates line of best fit for CT and log values
    CT_lin = model.predict(x) #obtains linear CT y values
    slope = float(model.coef_)
    intercept = model.intercept_
    r2 = model.score(x,y)
    
    plt.figure(figsize=(4,4)) #creating plot with pyplot
    plt.plot(log_copies,CT_lin, marker = 'o',ls = '--')
    plt.title('Viral Titer Standard Curve')
    plt.xlabel('Log copies/ml') ; plt.ylabel('CT')
    plt.savefig('std_c.png',dpi=150)

    wb = openpyxl.Workbook()
    sh = wb.create_sheet('graph')
    img = openpyxl.drawing.image.Image('std_c.png')
    sh.add_image(img, 'A3')
    plot_col = ['slope = ', str(slope), 'intercept = ', str(intercept), 'r^2 = ', str(r2)]
    sh.append(plot_col)
    del wb['Sheet']
    wb.save(sys.argv[2])
    wb.close()
    return slope, intercept, r2 #returns slope and intercept values
    

"""
create_frame : function

parameters : table dataframe (excel sheet 1) raw data dataframe (excel sheet 2)

while creating columns for frame, calls standard_curve_plot function to 
retrieve slope and intercept information necessary for further calculations

returns frame : a dataframe with all calculations and results
"""
def create_frame(table,data):
    well96 = [] #list of all tests done in loading table 
    well_pos = [] #creating list for 'Well Position' column
    del_list = [] # blank's index values to be removed from CT list later
    for i in range(8): #i from 0-7 : A-H
        row = [] 
        for j in range(0,4): # j from 1-4 : '1-3' - '10-12'
            index = table.iloc[i,j+1] #name of test at [Letter(A-H),triplet]
            if pd.notna(index): #doesn't append blanks
                row.append(index)
                for q in range(1,4):
                    well_pos.append('{}{}'.format(chr(65+i),q+3*j)) #adding 8 4 triplets of well positions
            elif pd.isna(index): 
                del_list.append((i+1)*(3*(j+1))+-3) #add blank's index to del_list
        well96.append(row)
    frame = pd.DataFrame(well_pos,columns=['Well Position']) #creating dataframe with first column as well position
    
    CT = data.iloc[43:,8].tolist() #splicing CT data from raw data excel
    for i in del_list:
       for j in range(3): #delete the triplet from CT list
            del CT[i]
    frame['CT'] = CT #adding CT column to dataframe

    #change CT from list to list of lists of triplets to make average and removal easier
    CT_t = []
    for i in range(int(len(CT)/3)):
        triplet = []
        triplet.append(CT[3*i])
        triplet.append(CT[3*i+1])
        triplet.append(CT[3*i+2])
        CT_t.append(triplet)

    avg_CT = [] #creating list of average CT for each triplet 
    for i in CT_t:
        avg = sum(i)/3
        for j in i:
            if j > (avg + 0.5) or j < (avg - 0.5): 
                #then the value is an outlier and we remove it and find new average
                i.remove(j)
        avg_CT.append(sum(i)/len(i))
        avg_CT.append('')
        avg_CT.append('')
    frame['Average CT'] = avg_CT #adding average CT column
    
    slope, intercept, r2 = standard_curve_plot(avg_CT) #create viral titer standard curve plot and get slope/intercept values
    
    log_copies = [(x-intercept)/slope for x in CT] #creates list of log copies by applying formula on each CT value
    frame['log copies/ ml'] = log_copies
    
    copies_ml = [np.power(10,x) for x in log_copies] #creates list of copies/ml by applying formula on each log copy
    frame['copies/ml'] = copies_ml

    qpcr_titer = [x*4520*2 for x in copies_ml]
    frame['qpcr titer'] = qpcr_titer

    qpcr_t = []  #creating list of triplets for qpcr_titer
    for i in range(int(len(qpcr_titer)/3)):
        triplet = []
        triplet.append(qpcr_titer[3*i])
        triplet.append(qpcr_titer[3*i+1])
        triplet.append(qpcr_titer[3*i+2])
        qpcr_t.append(triplet)
     
    avg_qpcr = [] #creating list of average qpcr for each triplet 
    sd_qpcr = [] #creating list of standard deviations for each triplet
    for i in qpcr_t:
        avg = sum(i)/3
        sd_qpcr.append(np.std(i))
        for j in i:
            if j > (avg + 0.5) or j < (avg - 0.5): 
                #then the value is an outlier and we remove it and find new average
                i.remove(j)
        avg_qpcr.append(sum(i)/len(i))
        avg_qpcr.append('')
        avg_qpcr.append('')
        sd_qpcr.append('')
        sd_qpcr.append('')
    frame['Average qpcr titer'] = avg_qpcr #adding average CT column
    frame['qpcr SD'] = sd_qpcr #adding standard deviation column

    return frame

#main handles command line inputs and adds frame dataframe to excel with pyplot
if __name__ == "__main__":
    if (len(sys.argv) != 3):
        print("No file input, argument format: './qpcr <input_file_name.(xlsx, xls, xlsm)> <output_file_name.(xlsx, xls, xlsm)>'")
    filename = sys.argv[1]
    data = pd.read_excel(filename, sheet_name = 0)
    loading_table = pd.read_excel(filename, sheet_name = 1)
    frame = create_frame(loading_table,data)
    with pd.ExcelWriter(sys.argv[2], mode='a') as writer:  
        frame.to_excel(writer, sheet_name='results')