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
    
    p = np.linspace(1,0,384)
    k=p
    
   
    
    separator = np.array([-6,-4,-2,-0,2,4,6,8])
    x_color = []
    k_color = []
    
    #x_color = np.array(x,[],[],[],[],[],[],[])
    alist = []
    blist = [0]
    blist = blist + alist
    
    
    #print (separator)
    for number in separator:
        print (number)
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
        #print (x)
        #print (x_color)
        
        #print (x2)
    #print (alist)
    #print (k_color[4])
    
    #print (len(x_color))
    



    
    colors = ['darkblue','cornflowerblue','lightskyblue','powderblue','mistyrose','lightcoral','firebrick','maroon']
    y = np.random.random(size=len(foldchange))
    #print (y)
    #k = np.linspace(1,0,384)
    #print (k)
    fh, ax = plt.subplots(1,1)
    for i in range(0,len(x_color)):
        #k = np.random.random(size=len(x_color[i]))
        ax.scatter(x_color[i],k_color[i],color=colors[i])
    indicators = {'P50747':'Biotin ligase','Q99873':'PRMT1','P17181':'IFNAR1','P42224':'STAT1'}
    #x_coor = []
    #y_coor = []
    print (len(y))
    for item in indicators:
        for row_num in range(len(foldchange)):
            if item == df.loc[:,'Accession'][row_num]:
                print (item)
                #y_coor.append(row_num)
                #x_coor.append(df.loc[:,'LOG2 (L/H)'][row_num])
                print (y[row_num])
                ax.annotate(indicators[item],(df.loc[:,'LOG2 (L/H)'][row_num],p[row_num]),arrowprops=dict(facecolor='black', shrink=0.05))
                break
    #print (x_coor)
    #print (y_coor)
                
                
    #ax.axes.get_yaxis().set_visible(False)
    plt.ylabel('Fraction')
    plt.xlabel('LOG2 (L/H)')
    plt.show()
    
    
    
one_d_dot("20180319_04_Qp1_Nyberg_beads.xlsx")
