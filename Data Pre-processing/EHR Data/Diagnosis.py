#!/usr/bin/env python
# coding: utf-8

from glob import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def LoadAllDiagnosis():
    
    date_parser=lambda x: pd.to_datetime(x,utc=True,format='%Y/%m/%d')
    files=[z for z in glob(f"/rawdata/*diagnoses_all*.csv")]
    rawindf=pd.concat([pd.read_csv(f,low_memory=False,parse_dates=['DIAGNOSIS DATE'],date_parser=date_parser) for f in files])
    rawindf['diag']=rawindf['DIAGNOSIS DATE'].dt.strftime('%Y%m')
    rawindf.sort_values(by=['ID','DIAGNOSIS DATE'],ascending=True,inplace=True)
    rawindf.rename(columns={'ICD10 CODE':'ICD10'},inplace=True)
    ## keep one record for same date, same ICD code
    indf=rawindf.drop_duplicates(subset=['ID','ICD10'],keep='first')
    
    return indf

# (1) Load Diagnosis results
DiagDF=LoadAllDiagnosis() 
## Remove records with missing ICD-10 fields 
print("-> remove missing in ICD10 field")
DiagDF.dropna(subset=['ICD10'],inplace=True)
DiagDF['ICD code']=DiagDF['ICD10'].apply(lambda x : x.split('.')[0])
## Count the occurrences of each ICD-10 code by patient
VisitCountDF=pd.crosstab(DiagDF['ID'],DiagDF['ICD code'])
VisitDF=VisitCountDF.copy()
VisitDF=VisitDF.applymap(lambda x : 1 if x>0 else 0)

# (2) Identify the ICD-10 code with the highest number of unique patients
TopN=20
SampleSizeDF=pd.DataFrame(VisitDF.sum(axis=0).to_dict().items(),columns=['ICD code','#samples'])
SampleSizeDF=SampleSizeDF[~SampleSizeDF['ICD code'].str.startswith('Z')]
SampleSizeDF.sort_values(by=['#samples'],ascending=False,inplace=True)
TopICD10=list(SampleSizeDF['ICD code'].head(TopN))
TopSampleSizeDF=SampleSizeDF.head(TopN).copy()

# (3) Load sample list and associated metadata
usecols=['MID','Birth','Age','Gender_EMR']
SampleInfoFile=f'{WorkPath}/SampleInfo/TpmiSample.tsv'
SampleDF=pd.read_csv(SampleInfoFile,sep='\t',usecols=usecols)
print(f"Number of samples : {SampleDF.shape[0]}")
BirthDict=dict((u,v) for u,v in zip(SampleDF['MID'],SampleDF['Birth']))

SampleTopICD=SampleDF.set_index('MID').join(VisitDF[TopICD10])

# (4) Calculate gender distribution for top 20 ICD-10 codes
TopSampleSizeDF.set_index('ICD code',inplace=True)
n_Female=SampleTopICD['Gender_EMR'].value_counts()['F']
n_Male=SampleTopICD['Gender_EMR'].value_counts()['M']
n_Total=n_Female+n_Male

for item in TopICD10:
    TopSampleSizeDF.loc[item,'Total(%)']=round(TopSampleSizeDF.loc[item,'#samples']/n_Total*100,2)
    TopSampleSizeDF.loc[item,'Male']=pd.crosstab(SampleTopICD['Gender_EMR'],SampleTopICD[item]).to_dict()[1]['M']
    TopSampleSizeDF.loc[item,'Female']=pd.crosstab(SampleTopICD['Gender_EMR'],SampleTopICD[item]).to_dict()[1]['F']
    TopSampleSizeDF.loc[item,'Male(%)']=round(TopSampleSizeDF.loc[item,'Male']/n_Male*100,2)
    TopSampleSizeDF.loc[item,'Female(%)']=round(TopSampleSizeDF.loc[item,'Female']/n_Female*100,2)
TopSampleSizeDF.reset_index(inplace=True)

# (5) Onset Age
TopDiagDF=DiagDF[DiagDF['ICD code'].isin(TopICD10)]
TopDiagDF['BIRTH']=TopDiagDF['ID'].apply(lambda x : BirthDict.get(x))
TopDiagDF['OnsetAge']=((pd.to_datetime(TopDiagDF['diag'],format="%Y%m")-pd.to_datetime(TopDiagDF['BIRTH'],format="%Y%m"))/365).astype('timedelta64[D]')

for i in range(len(TopICD10)):
    tempdf=TopDiagDF[TopDiagDF['ICD code']==TopICD10[i]]
    OnSetAgeDict=dict((u,v) for u,v in zip(tempdf['ID'],tempdf['OnsetAge']))
    SampleDF[TopICD10[i]]=SampleDF['MID'].apply(lambda x : OnSetAgeDict.get(x))

OnsetAgeDF=SampleDF[['Gender_EMR']+TopICD10].copy()
StaticsDict=OnsetAgeDF.groupby('Gender_EMR').describe().to_dict()


# (6) plot
StatDF=OnsetAgeDF.describe().T.sort_values(by=['mean'],ascending=False)
StatDF['order']=np.arange(StatDF.shape[0])+1 
PlotSortedItems=list(StatDF.index)
ItemYaxis=[[v,u+1] for u,v in enumerate(PlotSortedItems)]

def ExtractData(Dicts,inTrait,inSex):
    stats = {
        "label": inTrait,  
        "med": Dicts[(inTrait,'50%')][inSex],
        "mean":Dicts[(inTrait,'mean')][inSex],
        "q1": Dicts[(inTrait,'25%')][inSex],
        "q3": Dicts[(inTrait,'75%')][inSex],
        "whislo": Dicts[(inTrait,'min')][inSex],  
        "whishi": Dicts[(inTrait,'max')][inSex],
        "fliers": []
        }
    return stats


fs = 14 
ls = 2 
medianprops={'linewidth':ls,'color':'black'}
meanprops={'marker':'o','markerfacecolor':'black','markeredgecolor':'black'}
capprops={'linewidth':ls,'color':'black'} 
whiskerprops={'linewidth':ls,'color':'black'}
ColorDict={'F':'pink','M':'lightblue'}
BoxWidths=0.3
GapWidths=0.4

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 10))
plt.subplots_adjust(wspace=0)

for i in range(len(ItemYaxis)):
    item=ItemYaxis[i][0]
    POS=[ItemYaxis[i][1]]
    Sex='M'
    ax.bxp([ExtractData(StaticsDict,item,Sex)],widths=BoxWidths,positions=[z-(GapWidths/2) for z in POS],vert=False,showmeans=False,patch_artist=True,boxprops={'facecolor':ColorDict[Sex],'linewidth':ls},meanprops=meanprops,medianprops=medianprops,capprops=capprops,whiskerprops=whiskerprops)
    Sex='F'
    ax.bxp([ExtractData(StaticsDict,item,Sex)],widths=BoxWidths,positions=[z+(GapWidths/2) for z in POS],vert=False,showmeans=False,patch_artist=True,boxprops={'facecolor':ColorDict[Sex],'linewidth':ls},meanprops=meanprops,medianprops=medianprops,capprops=capprops,whiskerprops=whiskerprops)

ax.set_yticks([z[1] for z in ItemYaxis])
ax.set_yticklabels([z[0] for z in ItemYaxis], fontsize=14, weight='bold')
ax.set_xticks([z for z in np.arange(0,110,20)])
ax.set_xticklabels([z for z in np.arange(0,110,20)], fontsize=14, weight='bold')
ax.set_xlabel('Onset Age', fontsize=fs, weight='bold')
ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)


