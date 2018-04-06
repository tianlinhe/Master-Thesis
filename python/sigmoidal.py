#a "sigmoidal" curve for log2 (l/H) in SILAC of Mar
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
def f(value):
    return((-1)*math.log(value,2))

def one_d_dot(input):
    # open the excelfile
    xlsx = pd.ExcelFile(input)
# get the first sheet=master proteins as an object
    masterproteins = xlsx.parse(1)
# only the rows with defined values of Abundance Ratio: (Heavy) / (Light) are kept
    df = masterproteins[pd.notnull(masterproteins['Abundance Ratio: (Heavy) / (Light)'])]
    df=df.sort_values(['Abundance Ratio: (Heavy) / (Light)'], ascending = False)
    foldchange = np.array(df.loc[:,'Abundance Ratio: (Heavy) / (Light)'])
    f1 = np.vectorize(f, otypes=[np.float])
    x = f1(foldchange)
    df['LOG2 (L/H)']=x
    
#split dots in 8 color categories according to their x-value
#value of p makes it a sigmoidal shape. 
    p = np.linspace(1,0,len(x))
    k=p
    x_color = []
    k_color = []
    separator = np.array([-6,-4,-2,-0,2,4,6,8])
    
    for number in separator:
        count =0
        
        for item in x:
            if item < number:
                count +=1
                
        alist.append(count)
        x2=np.split(x,[count,len(x)])
        x_color.append(x2[0])
        x=x2[1]
        k2=np.split(k,[count,len(k)])
        k_color.append(k2[0])
        k=k2[1]       

#plot the dots with 8 colors   
    fh, ax = plt.subplots(1,1)
    colors = ['darkblue','cornflowerblue','lightskyblue','powderblue','mistyrose','lightcoral','firebrick','maroon']
    
    for i in range(0,len(x_color)):
        ax.scatter(x_color[i],k_color[i],color=colors[i])
        
 #add annotations for proteins in dic(indicators)       
    indicators = {'P50747':'Biotin ligase','Q99873':'PRMT1','P17181':'IFNAR1','P42224':'STAT1'}
    for item in indicators:
        for row_num in range(len(foldchange)):
            if item == df.loc[:,'Accession'][row_num]:
                ax.annotate(indicators[item],(df.loc[:,'LOG2 (L/H)'][row_num],p[row_num]),arrowprops=dict(facecolor='black', shrink=0.05))
                break
                               
    #ax.axes.get_yaxis().set_visible(False)
    plt.ylabel('Fraction')
    plt.xlabel('LOG2 (L/H)')
    plt.show()
      
one_d_dot("20180319_04_Qp1_Nyberg_beads.xlsx")
