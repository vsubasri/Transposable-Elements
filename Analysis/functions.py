import os
import re
import pandas as pd
import numpy as np
import seaborn as sns
from statannotations.Annotator import Annotator
from matplotlib import pyplot as plt
from scipy.stats import mannwhitneyu


#in the GDAF column of SV files, 
#there are rows with multiple AF values separated by commas
#this function replaces the entry with the maximum value from the original entry
def maxGDAFFilter(colName: str, df: pd.DataFrame)->pd.DataFrame:
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
    
    return df


# finds outliers of a column and returns them in a list
def findOutliers(df, std,yCol)->list:
    #return a list of sample_id that are outliers
    #give me somethig that is for one box in the box plot
    q1 = df.quantile(0.25).values[0]
    q3 = df.quantile(0.75).values[0]
    
    IQR = q3-q1
    
    upperBound = (q1-std*IQR)
    lowerBound = (q3+std*IQR)
    
    outliers = df.loc[(df[yCol]< upperBound) | (df[yCol]>lowerBound)]
    print(outliers)
    
    try:
        return list(outliers['Samples_ID'].unique())
    except:
        try:
            return list(outliers['sample_id'].unique())
        except:
            pass


# breaks df down into digestable columns so that findOutliers can be applied
def outliersInDf(df, xCol,yCol, std = 1.5):
    #unique list 
    uniqueList = list(df[xCol].unique())
    #get all types of boxes the stuff on the x axis categories
    for i in uniqueList:
        tempdf = df[df[xCol]==i]
        list_ids = findOutliers(tempdf,std,yCol)
        print("the outliers for "+i+" are:")
        print(list_ids)
        print()
        
        
# graph box graphs
def BoxGraphMulti(df: pd.DataFrame, xCol, yCol, compCol, orderList):    
    
    ax = sns.boxplot(data=df, x=xCol, y=yCol, hue=compCol, medianprops={"linewidth": 4, 'color':'black'}, showfliers=True)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, ha="right")
    
    lfs = df[df[compCol]=='lfs']
    kics = df[df[compCol]=='kics']
    
    print('LFS')
    outliersInDf(lfs, xCol,yCol)
    print('KICS')
    outliersInDf(kics, xCol,yCol)
    
    excludeList = list(set(lfs[xCol]).symmetric_difference(set(kics[xCol])))

    #add statistical test here
    listPairing = []
    #pairings for multi correction
    for i in orderList:
        if i not in excludeList:
            p1 = (i, 'kics')
            p2 = (i, 'lfs')
            p = (p1,p2)
            listPairing.append(p)
            
    if xCol == 'age6':
        p1 = ('>=6','kics')
        p2 = ('<6','kics')
        p = (p1,p2)
        listPairing.append(p)
        p1 = ('>=6','lfs')
        p2 = ('<6','lfs')
        p = (p1,p2)
        listPairing.append(p)
        p1 = ('>=6','lfs')
        p2 = ('Unaffected','lfs')
        p = (p1,p2)
        listPairing.append(p)
        p1 = ('Unaffected','lfs')
        p2 = ('<6','lfs')
        p = (p1,p2)
        listPairing.append(p)
        
        
    annot = Annotator(ax, listPairing, data=df, x=xCol, y=yCol, hue=compCol)
    annot.configure(test='Mann-Whitney', 
                    text_format='star', loc='outside', verbose=2)
    annot.apply_and_annotate()

    plt.show()
    
    
#plot bar graphs and normalize them by a factor
def BarGraphNormalized(label1: str, label2:str, df: pd.DataFrame, col1: str, col2: str,
                        xTitle: str, yTitle: str, divisor1: float, divisor2: float,
                        labels: list):

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
            
    plt.bar(x_axis - 0.2, [count1[a]/divisor1 for a in labels], 0.4, label = label1)
    plt.bar(x_axis + 0.2, [count2[a]/divisor2 for a in labels], 0.4, label = label2)
    
    plt.xticks(x_axis, labels, ha='right', rotation = 45)
    plt.xlabel(xTitle)
    plt.ylabel(yTitle)
    plt.title(title)
    plt.legend()
    plt.show()
    

#kics file matching (SV <> Clinical file)
def kIdAbbv(clinicDf, svDf, clinicCol, svCol):
    
    clinicDf.drop(clinicDf[clinicDf[clinicCol]=='Not applicable'].index, inplace = True)
    
    svDf[svCol] = svDf[svCol].astype(str)
    clinicDf[clinicCol] = clinicDf[clinicCol].astype(str)
    svDf[svCol] =  svDf[svCol].str.replace(" ", "")
    clinicDf[clinicCol] = clinicDf[clinicCol].str.replace(" ", "")
    
    list1 = svDf[svCol].unique()
    list2 = clinicDf[clinicCol].unique()
    
    svDict = dict()
    clinicDict = dict()
    
    for i in range (len(list1)):
        tempList = re.findall('\d{3,}', list1[i])

        #if empty, take the second one
        if not tempList:
            tempList = re.findall('\d+', list1[i])
            svDict[list1[i]] = tempList[1]

        #if not empty 
        else:
            svDict[list1[i]] = tempList[0]


    for i in range (len(list2)):    
        tempList = re.findall('\d{3,}', list2[i])
    
        #if empty, take the second one
        if not tempList:
            tempList = re.findall('\d+', list2[i])
            clinicDict[list2[i]] = tempList[1]

        #if not empty 
        else:
            clinicDict[list2[i]] = tempList[0]
    

    clinicDf['abbv_id'] = clinicDf[clinicCol].map(clinicDict)
    svDf['abbv_id'] = svDf[svCol].map(svDict)
    
    print(len(svDf[svCol].unique()))
        
    try:
        df = pd.merge(svDf[[svCol,'SV_type','abbv_id','ALT']],                   
                         clinicDf,
                         left_on = 'abbv_id', 
                         right_on = 'abbv_id', 
                         how='inner')
    except:
        try:
            print('try2')
            df = pd.merge(svDf[[svCol,'SV type','abbv_id']],                 
                         clinicDf,
                         left_on = 'abbv_id', 
                         right_on = 'abbv_id', 
                         how='inner')
        except:
            print('not good')
            df = pd.DataFrame()
    
    return df


#filters categories that have less than 5 samples/points
def filterBig(df, xcol):
    uniquelist = list(df[xcol].unique())
    n=5
    sampleIdCol = df.columns
    sampleIdCol = sampleIdCol[0]
    
    #if there is less than 5 (<5) unique ids then add to remove list
    removelist = []
    for i in uniquelist:
        dog = df.loc[df[xcol]==i]
        num = dog[sampleIdCol].nunique()        
        if num < n:
            removelist.append(i)
                
    df = df[~df[xcol].isin(removelist)]
        
    return df


#graphs the Tissue and Cancer diagnosis box plots
def graphSVTissue(kdf, ldf, isTissue):
    if isTissue:
        var = 'tissue_type'

    else:
        var = 'cancer_diagnosis'
    
    try:
        kdf = kdf.drop(['abbv_id', 'SV type','ageofonset'], axis=1)
        kdf = kdf.groupby(['sample_id']).value_counts()
    except:
        try:
            kdf = kdf.drop(['ALT'], axis=1)
            kdf = kdf.groupby(['Samples_ID']).value_counts()
        except:
            pass
        
    kdf = kdf.to_frame().reset_index()
    kdf.rename(columns = {0:'freq'}, inplace = True)
    
    kdf = filterBig(kdf,var)
    kdf['dataset'] = 'kics'

    try:
        ldf = ldf.drop(['SV chrom', 'GD_AF',
                            'SV type','ageofonset'], axis=1)
        ldf = ldf.groupby(['sample_id']).value_counts()
    except:
        try:
            ldf = ldf.drop(['ALT'], axis=1)
            ldf = ldf.groupby(['Samples_ID']).value_counts()
        except:
            pass
        
    ldf = ldf.to_frame().reset_index()
    ldf.rename(columns = {0:'freq'}, inplace = True)

    ldf = filterBig(ldf,var)
    ldf['dataset'] = 'lfs'

    mergedDf = pd.concat([ldf, kdf])
    
    #excluding all that don't have a pair         
    excludeList = list(set(ldf[var]).symmetric_difference(set(kdf[var])))
    mergedDf = mergedDf[~mergedDf[var].isin(excludeList)]
        
    uniqueList = list(mergedDf[var].unique())
    
    BoxGraphMulti(mergedDf, var, 'freq', 'dataset', uniqueList) 
    

#plots kics vs lfs and only one category (used for Overall graphs)
def boxplotPoints(title:list, column:list, df: pd.DataFrame , sizeH=20.50, sizeV=17.50, col='red', trans=0.25):    
    plt.rcParams["figure.figsize"] = [sizeH, sizeV]
    plt.rcParams["figure.autolayout"] = True
    data = pd.DataFrame({
        t: df[c] for t,c in zip(title, column)})
    
    plotTest = data
    plotTest = plotTest.stack().to_frame().reset_index().rename(columns={'level_1': 'iden', 0: 'value'}).drop('level_0', axis='columns')

    sns.boxplot(data=plotTest, x='iden', y='value')
        
    #plot with boxplot
    sns.stripplot(x = 'iden',
              y = 'value',
                  color = 'red',
                  alpha = 0.25,
              data = plotTest)
    
    list1 = df[column[0]].dropna()
    list2 = df[column[1]].dropna()
    print(mannwhitneyu(list1,list2))
    plt.show()
    
    
    
#plots overall tissue or cancer graphs, one on big category with comparisons between all
def graphBoxGen(xColumnName, df, ycol = 0):
    
    df = filterBig(df, xColumnName)
    outliersInDf(df,xColumnName,ycol)

    plt.rcParams["figure.figsize"] = [15, 10]
    ax = sns.boxplot(data=df, x=xColumnName, y=ycol, medianprops={"linewidth": 4, 'color':'black'})
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, ha="right")

    #add statistical test here

    # initializing list 
    uniqueTissue = list(df[xColumnName].unique()) 

    #does all possible pairings given a list
    listPairing = [(a,b) for x, a in enumerate(uniqueTissue) for b in uniqueTissue[x+1:]]
    #print(listPairing)

    annot = Annotator(ax, listPairing, data=df, x=xColumnName, y=ycol)
    annot.configure(test='Mann-Whitney',
                    text_format='star', loc='outside', verbose=2)
    annot.apply_and_annotate()
    
    
#functions used to make a new column for age of onset graphs
def ageOfOnsetKics(row):
    if row['ageofonset']<2190:
        return '<6'
    if row['ageofonset']>=2190:
        return '>=6'
    
def ageOfOnsetLfs(row):
    if row['ageofonset']<72:
        return '<6'
    if row['ageofonset']>=72:
        return '>=6' 
    
#implementing and cleaning up the df before sending to be plotted
def cleaning(df):
    try:
        df = df.drop(['SV type','ageofonset'],axis=1)
        df = df.groupby(['sample_id']).value_counts()
    except:
        try:
            df = df.drop(['ALT','ageofonset'],axis=1)
            df = df.groupby(['Samples_ID']).value_counts()
        except:
            pass
    df = df.to_frame().reset_index()
    df.rename(columns = {0: 'freq'},inplace = True)
    return df

def mergeDFs(lfs1, lfs2, kics1):
    lfs1 = cleaning(lfs1)
    lfs2 = cleaning(lfs2)
    kics1 = cleaning(kics1)
    lfsMerged = pd.concat([lfs1,lfs2])
    lfsMerged['dataset'] = 'lfs'
    kics1['dataset'] = 'kics'
    mergedDf = pd.concat([lfsMerged, kics1])
    
    BoxGraphMulti(mergedDf, 'age6', 'freq', 'dataset', ['<6','>=6','Unaffected'])

    
#plotting RMS subtype, needed to create a column for the RMS subtype
def rmsType(row):
    if ('lveolar' in row['type']) and ('habdomyosarcoma' in row['type']):
        return 'ARMS'
    if ('mbryonal' in row['type']) and ('habdomyosarcoma' in row['type']):
        return 'ERMS'
    if('ERMS' in row['type']) and ('ARMS' not in row['type']):
        return 'ERMS'
    return row['type']

#getting df into usable format
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


#formatting dataframe to get a normalized and original df
def formatDataFrame(df: pd.DataFrame, groupByList: list, normalizeList: list, 
                    iterateList: list, labels: list, typeSV: str ) -> tuple:
    
    dfGrouped = df.groupby(groupByList).size().unstack(fill_value=0)
    dfReg = pd.DataFrame()
    dfNorm = pd.DataFrame()
    df1WC = pd.DataFrame()
    df2WC = pd.DataFrame()
    
    for i in iterateList:
        tempList = []
        
        try:
            series = dfGrouped.loc[(i,typeSV)]
            
        except:
            series = pd.Series(0, index=labels)
            
        for j,div in zip(labels, normalizeList):
            tempList.append(series[j]/div)
            
        d = {'chrom':labels, 'normalized':tempList}
        tempdf = pd.DataFrame(d)

        dfNorm = pd.concat([dfNorm, tempdf], axis=0)
        dfReg = pd.concat([dfReg, series], axis=0) 
    
    return(dfReg, dfNorm)