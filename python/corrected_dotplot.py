import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
def f(value):
    return((-1)*math.log(value,2))

def one_d_dot(input):
    # open the excelfile
    xlsx = pd.ExcelFile(input)
# get the first sheet=2_masterproteins as an object
    masterproteins = xlsx.parse(2)
# only the rows with defined values of Abundance Ratio: (Heavy) / (Light) are kept
    df = masterproteins[pd.notnull(masterproteins['Abundance Ratio: (Heavy) / (Light)'])]
    df=df.sort_values(['Abundance Ratio: (Heavy) / (Light)'], ascending = False)
    foldchange = np.array(df.loc[:,'Abundance Ratio: (Heavy) / (Light)'])
    f1 = np.vectorize(f, otypes=[np.float])
    x = f1(foldchange)
    df['LOG2 (L/H)']=x
    


#plot the dots with 8 colors   
    fh, ax = plt.subplots(1,1)
    
    
    p = np.linspace(1,0,len(x))
    for i in range(0,len(x)):
        ax.scatter(p[i],x[i],facecolors='none', edgecolors='grey')
        
 #add annotations for proteins in dic(indicators)       
    indicators = {'P50747':'Biotin ligase','Q99873':'PRMT1','P17181':'IFNAR1','P42224':'STAT1'}
    for item in indicators:
        for row_num in range(len(foldchange)):
            if item == df.loc[:,'Accession'][row_num]:
                ax.annotate(indicators[item],(p[row_num],df.loc[:,'LOG2 (L/H)'][row_num]),arrowprops=dict(facecolor='black', shrink=0.05))
                break
                               
    #ax.axes.get_yaxis().set_visible(False)
    plt.xlabel('Fraction')
    plt.ylabel('LOG2 (L/H)')
    plt.show()
      
one_d_dot("20180319_04_Qp1_Nyberg_beads.xlsx")
