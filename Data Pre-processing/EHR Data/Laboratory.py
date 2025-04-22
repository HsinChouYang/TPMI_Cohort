#!/usr/bin/env python
# coding: utf-8

from glob import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def LabFiles(AllFiles):
    infiles=[z for z in AllFiles if 'lab_tests' in z]
    print(f'\n... {len(infiles)} lab test files found ...')
    print(f" ---> {infiles}")
    indf=pd.concat([pd.read_csv(f,encoding="utf-8-sig",parse_dates=['Sampling Date'],low_memory=False) for f in infiles])
    print(f"\nNumber of original records : {indf.shape[0]}")
    print(f"Number of samples : {len(indf['ID'].unique())}")
    return indf


# (1) Load sample info
SampleInfoFile=f'TpmiSample.tsv'
SampleDF=pd.read_csv(SampleInfoFile,sep='\t')
SampleRecordsDF=SampleDF.copy()
SampleTrackYearsDF=SampleDF.copy()

# (2) Load Lab results
Files=glob(f"/rawData/*")
LabDF=LabFiles(Files)
LabDF.rename(columns={'ID':'MID'},inplace=True)

## Data preprocessing for code book:
## 1. Load and standardize column names in the code book file
## 2. Filter items that contain specific suffixes (_UNIT, _REF, _TYPE)
## 3. Remove secondary test result items from the filtered list
RenameDict={'VARIABLE Name':'NAME',
        'VARIABLE Description(English)':'Description_ENG',
        'VARIABLE Description (Chinese)':'Description_CHN',
        'NHI Code':'NHI'
       }
CodeBookFile=[z for z in Files if 'code_book' in z][0]
CodeBookDF=pd.read_csv(CodeBookFile)
CodeBookDF.rename(columns=RenameDict,inplace=True)
Items=[z for z in list(CodeBookDF['NAME'].unique()) if '_UNIT' not in z and '_REF' not in z and '_TYPE' not in z]
for a in ['ANTI_DELTA_2','ANTI_HBC_2','ANTI_HBE_2','ANTI_HBS_2','HBSAG_2','HCV_AB_2']:
    Items.remove(a)

def Statistics(item):
    global SampleSizeDF
    
    ItemDF=LabDF[['MID','Sampling Date',item]].dropna(subset=[item]).copy() ## Remove records with missing values
    ItemDF=ItemDF[ItemDF['Sampling Date'].dt.year>=2014] ## Filter out EMR records before 2014
    ItemDF.sort_values(by=['MID','Sampling Date'],inplace=True) ## Sort records by sample ID and Sampling Date
    ItemDF.drop_duplicates(subset=['MID','Sampling Date'],keep='last',inplace=True) ## Remove same-day duplicates for each patient
    
    RecordsDict=ItemDF['MID'].value_counts().to_dict() ##  Number of records per patient
    tempdf=ItemDF.drop_duplicates(subset=['MID'],keep='first').copy()
    StartDict=dict((u,v) for u,v in zip(tempdf['MID'],tempdf['Sampling Date'])) ## First sampling date per patient
    tempdf=ItemDF.drop_duplicates(subset=['MID'],keep='last').copy()
    EndDict=dict((u,v) for u,v in zip(tempdf['MID'],tempdf['Sampling Date'])) ## ## Last sampling date per patient

    SampleCheck=SampleDF.loc[:,'MID':'Generation'].copy() ## Create a temporary dataframe for sample information
    ## Add columns
    SampleCheck['Records']=SampleCheck['MID'].apply(lambda x : RecordsDict.get(x))
    SampleCheck['Start']=SampleCheck['MID'].apply(lambda x : StartDict.get(x))
    SampleCheck['End']=SampleCheck['MID'].apply(lambda x : EndDict.get(x))
    SampleCheck['TrackYears']=round((SampleCheck['End']-SampleCheck['Start']).dt.days/365.2425,2)
    SampleCheck.dropna(subset=['TrackYears'],inplace=True)
    TrackYearsDict=dict((u,v) for u,v in zip(SampleCheck['MID'],SampleCheck['TrackYears'])) ##  tracking years per patient
    
    ## Calculate descriptive statistics
    Stats_Records_Dict=SampleCheck['Records'].describe().to_dict()
    Stats_TrackYear_Dict=SampleCheck['TrackYears'].describe().to_dict()
    Stats_BySex_Dict=SampleCheck.groupby('Gender_EMR').describe().to_dict()
   
    ## Format statistical results for reporting
    ReturnSampleSize_ALL=f"{Stats_Records_Dict['count']} ({round(Stats_Records_Dict['count']/SampleDF.shape[0]*100 ,2)}%)"
    ReturnSampleSize_M=f"{Stats_BySex_Dict[('Records', 'count')]['M']} ({round(Stats_BySex_Dict[('Records', 'count')]['M']/SampleDF.shape[0]*100 ,2)}%)"
    ReturnSampleSize_F=f"{Stats_BySex_Dict[('Records', 'count')]['F']} ({round(Stats_BySex_Dict[('Records', 'count')]['F']/SampleDF.shape[0]*100 ,2)}%)"
    ReturnRecords_ALL=f"{round(Stats_Records_Dict['mean'],2)} ({round(Stats_Records_Dict['std'],2)})"
    ReturnTrackYear_ALL=f"{round(Stats_TrackYear_Dict['mean'],2)} ({round(Stats_TrackYear_Dict['std'],2)})"
    ReturnRecords_M=f"{round(Stats_BySex_Dict[('Records', 'mean')]['M'],2)} ({round(Stats_BySex_Dict[('Records', 'std')]['M'],2)})"
    ReturnRecords_F=f"{round(Stats_BySex_Dict[('Records', 'mean')]['F'],2)} ({round(Stats_BySex_Dict[('Records', 'std')]['F'],2)})"
    ReturnTrackYear_M=f"{round(Stats_BySex_Dict[('TrackYears', 'mean')]['M'],2)} ({round(Stats_BySex_Dict[('TrackYears', 'std')]['M'],2)})"
    ReturnTrackYear_F=f"{round(Stats_BySex_Dict[('TrackYears', 'mean')]['F'],2)} ({round(Stats_BySex_Dict[('TrackYears', 'std')]['F'],2)})"
    
    SampleSizeDF.loc[item,'Sample Size']=Stats_Records_Dict['count']
    SampleSizeDF.loc[item,'Percentage']=round(Stats_Records_Dict['count']/SampleDF.shape[0]*100 ,2)
    
    return RecordsDict,TrackYearsDict,ReturnSampleSize_ALL,ReturnSampleSize_M,ReturnSampleSize_F,ReturnRecords_ALL,ReturnRecords_M,ReturnRecords_F,ReturnTrackYear_ALL,ReturnTrackYear_M,ReturnTrackYear_F


TraitDF=CodeBookDF[CodeBookDF['NAME'].isin(Items)].copy()
TraitDF.set_index('NAME',inplace=True)
SampleSizeDF=CodeBookDF[CodeBookDF['NAME'].isin(Items)].copy()
SampleSizeDF.set_index('NAME',inplace=True)

for item in Items:
    [RecordsDict,TrackYearsDict,NALL,NM,NF,RALL,RM,RF,TALL,TM,TF]=Statistics(item)
    SampleRecordsDF[item]=SampleDF['MID'].apply(lambda x : RecordsDict.get(x))
    SampleTrackYearsDF[item]=SampleDF['MID'].apply(lambda x : TrackYearsDict.get(x))
    TraitDF.loc[item,'Sample Size(ALL)']=NALL
    TraitDF.loc[item,'Sample Size(Male)']=NM
    TraitDF.loc[item,'Sample Size(Female)']=NF
    TraitDF.loc[item,'Records(ALL)']=RALL
    TraitDF.loc[item,'Records(Male)']=RM
    TraitDF.loc[item,'Records(Female)']=RF
    TraitDF.loc[item,'Track Years(ALL)']=TALL
    TraitDF.loc[item,'Track Years(Male)']=TM
    TraitDF.loc[item,'Track Years(Female)']=TF

## Top20 items
Total_N=SampleDF.shape[0]
TraitIDs=list(TraitDF.index)
SummaryStatDict=SampleTrackYearsDF.describe().to_dict()
SampleSizeDict=dict()
for item in TraitIDs:
    SampleSizeDict[item]=SummaryStatDict[item]['count']
SampleSizeDF=pd.DataFrame(SampleSizeDict.items(),columns=['Trait','Sample size'])
SampleSizeDF['percentage']=SampleSizeDF['Sample size']/Total_N
SampleSizeDF.sort_values(by=['Sample size'],ascending=False,inplace=True)
PlotItems=list(SampleSizeDF['Trait'][:20])
    
## Apply Winsorization at the 95th percentile to the data volume
Threshold=95
RecordCount_Top20=SampleRecordsDF[['Gender_EMR']+PlotItems].copy()
for i in range(len(PlotItems)):
    item=PlotItems[i]
    CutValue=np.nanpercentile(RecordCount_Top20[item],Threshold)
    RecordCount_Top20[item]=np.where(RecordCount_Top20[item]>CutValue,CutValue,RecordCount_Top20[item])
StaticsDict_RC=RecordCount_Top20.groupby('Gender_EMR').describe().to_dict()

TrackYear_Top20=SampleTrackYearsDF[['Gender_EMR']+PlotItems].copy()
StaticsDict_TY=TrackYear_Top20.groupby('Gender_EMR').describe().to_dict()

## sort by medium
TopStatDF=RecordCount_Top20.describe().T.sort_values(by=['count'],ascending=True)
TopStatDF['percentage']=TopStatDF['count'].apply(lambda x : round(x/Total_N*100,1))
TopStatDF['order']=np.arange(TopStatDF.shape[0])+1
PlotSortedItems=list(TopStatDF.index)
ItemYaxis=[[v,u+1] for u,v in enumerate(PlotSortedItems)]

## Plot
def ExtractData(Dicts,inTrait,inSex):
    stats = {
        "label": inTrait,  # not required
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

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(7, 10), sharey=True)
plt.subplots_adjust(wspace=0)

##### record count #####
for i in range(len(ItemYaxis)):
    item=ItemYaxis[i][0]
    POS=[ItemYaxis[i][1]]
    Sex='M'
    axes[0].bxp([ExtractData(StaticsDict_RC,item,Sex)],widths=BoxWidths,positions=[z-(GapWidths/2) for z in POS],vert=False,showmeans=False,patch_artist=True,boxprops={'facecolor':ColorDict[Sex],'linewidth':ls},meanprops=meanprops,medianprops=medianprops,capprops=capprops,whiskerprops=whiskerprops)
    Sex='F'
    axes[0].bxp([ExtractData(StaticsDict_RC,item,Sex)],widths=BoxWidths,positions=[z+(GapWidths/2) for z in POS],vert=False,showmeans=False,patch_artist=True,boxprops={'facecolor':ColorDict[Sex],'linewidth':ls},meanprops=meanprops,medianprops=medianprops,capprops=capprops,whiskerprops=whiskerprops)

axes[0].set_yticks([z[1] for z in ItemYaxis])
axes[0].set_yticklabels([z[0] for z in ItemYaxis], fontsize=14, weight='bold')
axes[0].set_xticks([z for z in np.arange(0,50,10)])
axes[0].set_xticklabels([z for z in np.arange(0,50,10)], fontsize=14, weight='bold')
axes[0].set_xlabel('# of records', fontsize=fs, weight='bold')
axes[0].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

##### sample size ######
axes[1].barh(TopStatDF['order'],TopStatDF['percentage'], color='lightgrey',height=0.6)
axes[1].set_xticks([z for z in np.arange(0,100,20)])
axes[1].set_xticklabels([z for z in np.arange(0,100,20)], fontsize=14, weight='bold')
axes[1].set_xlabel('Proportion (%)', fontsize=fs, weight='bold')
axes[1].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
axes[1].set_xlim(left=0, right=100)

for index, row in TopStatDF.iterrows():
    axes[1].text(row['percentage'], row['order'], str(row['percentage']), fontsize=fs, weight='bold')

##### track year #####
for i in range(len(ItemYaxis)):
    item=ItemYaxis[i][0]
    POS=[ItemYaxis[i][1]]
    Sex='M'
    axes[2].bxp([ExtractData(StaticsDict_TY,item,Sex)],widths=BoxWidths,positions=[z-(GapWidths/2) for z in POS],vert=False,showmeans=False,patch_artist=True,boxprops={'facecolor':ColorDict[Sex],'linewidth':ls},meanprops=meanprops,medianprops=medianprops,capprops=capprops,whiskerprops=whiskerprops)
    Sex='F'
    axes[2].bxp([ExtractData(StaticsDict_TY,item,Sex)],widths=BoxWidths,positions=[z+(GapWidths/2) for z in POS],vert=False,showmeans=False,patch_artist=True,boxprops={'facecolor':ColorDict[Sex],'linewidth':ls},meanprops=meanprops,medianprops=medianprops,capprops=capprops,whiskerprops=whiskerprops)

axes[2].set_yticks([z[1] for z in ItemYaxis])
axes[2].set_yticklabels([z[0] for z in ItemYaxis], fontsize=14, weight='bold')
axes[2].set_xticks([z for z in np.arange(0,12,2)])
axes[2].set_xticklabels([z for z in np.arange(0,12,2)], fontsize=14, weight='bold')
axes[2].set_xlabel('# of track years', fontsize=fs, weight='bold')
axes[2].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)


