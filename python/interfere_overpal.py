#it shows the uncorrected Log2 (L/H)
#it matches MS-hits with Interferome database
#input 1)sheet=masterproteins from raw file
    #  2) Interferom
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
def f(value):
    return((-1)*math.log(value,2))

def interferome_overlap(input1, input2):
    # open the rwo file, sheet masterproteins
    xlsx = pd.ExcelFile(input1)
# get the first sheet=master proteins as an object
    masterproteins = xlsx.parse(1)
# only the rows with defined values of Abundance Ratio: (Heavy) / (Light) are kept
    df = masterproteins[pd.notnull(masterproteins['Abundance Ratio: (Heavy) / (Light)'])]
    df=df.sort_values(['Abundance Ratio: (Heavy) / (Light)'], ascending = False)
    foldchange = np.array(df.loc[:,'Abundance Ratio: (Heavy) / (Light)'])
    f1 = np.vectorize(f, otypes=[np.float])
    x = f1(foldchange)
    df['LOG2 (L/H)']=x
    genelist = np.array(df.loc[:,'Description'])
    
  # open the rwo file, sheet masterproteins
    xlsx = pd.ExcelFile(input2)
#df2=foldchange Up2 and Down10000 in interferome database -> up-regulation
    df2 = xlsx.parse(1)
    df2.drop(df2.index[0:18], inplace=True)
    df2.columns = df2.iloc[0]
    df2.drop(df2.index[0],inplace = True)
#df3=foldchange Up10000 and Down2 in interferome datbase -> down-regulation  
    df3 = xlsx.parse(2)
    df3.drop(df3.index[0:18], inplace=True)
    df3.columns = df3.iloc[0]
    df3.drop(df3.index[0],inplace = True)
 
    fh, (ax1, ax) = plt.subplots(1, 2, sharey=True,sharex=True,figsize=(15,15))
    p = np.linspace(1,0,len(x))
    ax1.scatter(p,x,facecolors='none', edgecolors='grey',label='significant hits')
   
   
    count = 0
    for index, row in df2.iterrows():
        
        for row_num in range(len(foldchange)):
            
            if row['Gene Name'] == df.loc[:,'Description'][row_num][0:-1]:
                up=ax.scatter(p[row_num],df.loc[:,'LOG2 (L/H)'][row_num],facecolors='none', edgecolors='red')
                count += 1
                break          
    print ("Up-regulated IRGs: ",  count)        
   
    count = 0
    for index, row in df3.iterrows():
        
        for row_num in range(len(foldchange)):
            
            if row['Gene Name'] == df.loc[:,'Description'][row_num][0:-1]:
                down=ax.scatter(p[row_num],df.loc[:,'LOG2 (L/H)'][row_num],facecolors='none', edgecolors='black')
                count += 1
                break          
    print ("Down-regulated IRGs: ", count)        
           
    ax.legend((up,down),('up-regulated IRGs','down-regulated IRGs'),loc='upper right')
    ax1.legend('hits',loc='upper right')
    plt.xlabel('Fraction')
    plt.ylabel('LOG2 (L/H)')   
    plt.show()
    
interferome_overlap("20180319_04_Qp1_Nyberg_beads.xlsx","Interferome/Interferome_genelist.xlsx")
