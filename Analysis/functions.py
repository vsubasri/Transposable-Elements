import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


""" 
Bar Graph Normalized
"""
def BarGraphNormalized(label1: str, label2:str, df: pd.DataFrame, col1: str, col2: str,
                        xTitle: str, yTitle: str, divisor1: float, divisor2: float,
                        labels: list):
    #'kics', 'lfs', invdf, 'kchrom', 'lchrom', 'Chromosomes','Frequency','Graph2-INV',kicsSVnum, lfsSVnum

    title = label1 + " " + label2
    count1 = df[col1].value_counts()
    count2 = df[col2].value_counts()
    
    x_axis = np.arange(len(labels))
        
    for i in labels:
        if (not(i in count1)):
            add = pd.Series([0], index=[i])
            count1 = count1.append(add)
        if (not(i in count2)):
            add = pd.Series([0], index=[i])
            count2 = count2.append(add)
    
    plt.bar(x_axis - 0.2, [count1[a] for a in labels]/divisor1, 0.4, label = label1)
    plt.bar(x_axis + 0.2, [count2[a] for a in labels]/divisor2, 0.4, label = label2)

    plt.xticks(x_axis, labels)
    plt.xlabel(xTitle)
    plt.ylabel(yTitle)
    plt.title(title)
    plt.legend()
    plt.show()
    
    
"""
Filter to find largest in GDAF
"""
def maxGDAFFilter(colName: str, df: pd.DataFrame):
    temp = df[colName].str.contains(pat=',')
    Index = list(np.where(temp == True)[0])
    temp = df[colName].iloc[Index]
    tempIndex = temp.index

    for i in tempIndex:
        tempList = temp[i].split(',')
        tempList = [float(a) for a in tempList]
        maxTemp = max(tempList)
        df[colName][i] = maxTemp

    df[colName] = df[colName].astype('float64')
    
"""
Making Bar Graphs
"""

def BoxGraphMulti(df: pd.DataFrame, xCol, yCol, compCol):
    sns.boxplot(data=df, x=xCol, y=yCol, hue=compCol, palette='spring')
    #plt.show()

def boxplotPoints(title:list, column:list, df: pd.DataFrame , sizeH=20.50, sizeV=17.50, col='red', trans=0.25):
    plt.rcParams["figure.figsize"] = [sizeH, sizeV]
    plt.rcParams["figure.autolayout"] = True
    data = pd.DataFrame({
        t: df[c] for t,c in zip(title, column)})

    plotTest = data
    plotTest = plotTest.stack().to_frame().reset_index().rename(columns={'level_1': 'iden', 0: 'value'}).drop('level_0', axis='columns')
    plotTest['cat'] = 1

    sns.boxplot(data=plotTest, x='iden', y='value')
    
    #plot with boxplot
    sns.stripplot(x = 'iden',
              y = 'value',
                  color = 'red',
                  alpha = 0.25,
              data = plotTest)
    
    for i, d in enumerate(data):
        print(i)
        y = data[d]
        x = np.random.normal(i, 0.04, len(y))
        plt.scatter(x, y, color = col, alpha = trans)
    

"""
Data prep for boxplot
"""

def makeDataSet(identifiers, series: pd.Series, secondIndex: str)->list:
    
    dataList = []
    
    for i in identifiers:
        try:
            dataList.append(series[(i, secondIndex)])
        except:
            #pass
            dataList.append(0)

    return(dataList)

        
def makeUnequalDF(list1: list, list2: list) -> pd.DataFrame: #only for kics vs lfs D:
    tempDict = dict(kics = list1, lfs = list2)
    df = pd.DataFrame(dict([(k, pd.Series(v)) for k,v in tempDict.items()]))
    return df


"""

Normalized and Regular Boxplot, multiple side by side

"""

def formatDataFrame(df: pd.DataFrame, groupByList: list, normalizeList: list, 
                    iterateList: list, labels: list, typeSV: str, ) -> tuple:
    
    #kicsMore, grouping, numBPChrom, uniqueK, uniqueLabels, 'DEL' 

    dfGrouped = df.groupby(groupByList).size().unstack(fill_value=0)
    dfReg = pd.DataFrame()
    dfNorm = pd.DataFrame()
    
    for i in iterateList:
        tempList = []
        
        try:
            series = dfGrouped.loc[(i,typeSV)]
            
        except:
            series = pd.Series(0, index=labels)
            
        for j,div in zip(labels, normalizeList):
            tempList.append(series[j]/div)
            
        d = {'labels':labels, 'normalized':tempList}
        tempdf = pd.DataFrame(d)

        dfNorm = pd.concat([dfNorm, tempdf], axis=0)
        dfReg = pd.concat([dfReg, series], axis=0)
    
    return(dfReg, dfNorm)

    
def BoxGraphMulti(df: pd.DataFrame, xCol, yCol, compCol):
    sns.boxplot(data=df, x=xCol, y=yCol, hue=compCol, palette='spring')
    plt.show()