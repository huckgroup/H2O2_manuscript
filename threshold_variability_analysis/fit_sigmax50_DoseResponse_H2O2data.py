# -*- coding: utf-8 -*-
"""

All content of this script refers to "Dynamic hydrogen peroxide levels reveal a 
rate-dependent sensitivity in B-cell lymphoma signaling" by Witmond et al.
For details, refer to the Materials and Methods of the manuscript.

The enclosed code can be used to estimate the level of cell-to-cell variability
in the response threshold to predict the existence of bimodality in the population
distribution described by Hill-type kinetics as derived by Dobrzynski et al. 
(J R Soc Interface 2014).  

"""

# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.optimize as opt
from sklearn.neighbors import KernelDensity

# Functions
def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .1, label, fontsize=28, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)
    
def DRcurve(s,n,b,Rmax,EC50):
    '''
    Parameters
    ----------
    s : stimulus concentration
    n : hill coefficient
    b : basal response
    Rmax : maximal response (NOTE: not max. response amplitude!!)
    EC50 : response threshold

    Returns
    -------
    response concentration

    '''
    Rmax2=Rmax-b
    return b+Rmax2*s**n/(s**n + EC50**n)

def gcut(x,mx50,sigx50,scale):
    gx = scale*(1/(mx50*sigx50*np.sqrt(2*np.pi)))*np.exp(-(sigx50**2/2)-(1/(2*sigx50**2))*np.log(mx50/x)**2)
    return gx

def alpha(H,sigx50):
    aplus=np.nan
    aminus=np.nan
    if H**2*sigx50**2>=2:
        aplus=sigx50*np.sqrt(H**2*sigx50**2-2) + (1/H)*np.log(H**2*sigx50**2-1-H*sigx50*np.sqrt(H**2*sigx50**2-2))
        aminus=-sigx50*np.sqrt(H**2*sigx50**2-2) + (1/H)*np.log(H**2*sigx50**2-1+H*sigx50*np.sqrt(H**2*sigx50**2-2))
    return [aminus,aplus]

#%% Load data
df=pd.read_csv('C:\\Users\\Emma\\Documents\\Python Scripts\\H2O2_data\\DynSign.095_annotated\\annotated_data_tidy_DS095.csv')
df.dropna(inplace=True)
df.reset_index(drop=True,inplace=True)


#%% Continue with 1 timepoint
allproteins=['pCD79a (Y182)', 'pPLCy2 (Y759)','pSYK (Y525/Y526)']
allconc=np.unique(df['H2O2_conc_mM'])
alltime=[5]#np.unique(df['time_min'])
allreps=np.unique(df['replicate'])

df2=df[['H2O2_conc_mM','protein','time_min','fluorescence','replicate']].copy()
df2=df2[df2.protein.isin(allproteins)]
df2=df2[df2.replicate.isin(allreps)]
df2=df2[df2.time_min.isin(alltime)]
df2.reset_index(drop=True,inplace=True)


#%% Compute fold change w.r.t. concentration 0
mfiControl=pd.DataFrame({'protein': [], 'H2O2_conc_mM': [], 'time_min':[], 'replicate': [], 'mfiControl': []})
for prot in allproteins:
    for rep in allreps:
        data=df2[(df2['H2O2_conc_mM']==0)&(df2['protein']==prot)&(df2['replicate']==rep)]
        mRep=data['fluorescence'].median()
        mfiControl = pd.concat([mfiControl, pd.DataFrame.from_records([{'protein': prot, 'H2O2_conc_mM': 0,'time_min': alltime[0], 'replicate': rep, 'mfiControl': mRep}])], ignore_index=True)
g1=mfiControl.groupby(['protein', 'time_min','replicate'])
g2=df2.groupby(['protein', 'time_min','replicate'])

df2['mfiControl']=np.nan
for prot in allproteins:
    for time in alltime:    
        for rep in allreps:
            key=(prot,time,rep)
            gr1=g1.get_group(key)
            gr2=g2.get_group(key)
            mfi=gr1['mfiControl'].iloc[0]
            df2['mfiControl'].iloc[gr2.index]=mfi

df2['medFC']=np.log10(df2['fluorescence']/df2['mfiControl'])


#%% Compute modes of density
data1=df2[(df2.H2O2_conc_mM==0)].round(2)
data2=df2[(df2.H2O2_conc_mM==25)].round(2)
ymin=data1.groupby(['protein','replicate'])['medFC'].agg(lambda x: pd.Series.mode(x)[0])
ymax=data2.groupby(['protein','replicate'])['medFC'].agg(lambda x: pd.Series.mode(x)[0])
ystar=(ymin+ymax)/2


#%% Extract kde
sigx50dens=pd.DataFrame({'protein': [], 'H2O2_conc_mM': [], 'replicate': [], 'density': []})
xtest = np.linspace(-1,3,256)

fig,axes=plt.subplots(1,3,figsize=(10,3),sharey=True)
for ip,prot in enumerate(allproteins):
    ax=axes[ip]
    for conc in allconc:
        for rep in allreps:
            data=df2[(df2.protein==prot)&(df2.H2O2_conc_mM==conc)&(df2.replicate==rep)].medFC
            kde = KernelDensity(kernel='gaussian', bandwidth=0.05).fit(data.values.reshape(-1,1))
            log_dens=kde.score_samples(xtest.reshape(-1,1))
            ax.set(title=prot)
            ax.plot(xtest,np.exp(log_dens))
            ax.axvline(ystar[prot,rep])
            # find nearest value:
            idx=(np.abs(xtest-ystar[prot,rep])).argmin()
            ax.scatter(xtest[idx],np.exp(log_dens[idx]))
            sigx50dens = pd.concat([sigx50dens, pd.DataFrame.from_records([{'protein': prot, 'H2O2_conc_mM': conc, 'replicate': rep, 'density': np.exp(log_dens[idx])}])], ignore_index=True)
            prob=np.exp(log_dens)
            scale=np.trapz(prob.ravel(),xtest.ravel())


fitData=[]
groups=sigx50dens.groupby(['protein','replicate'])
for name,group in groups:
    scale=np.trapz(group.density.ravel(), group.H2O2_conc_mM.ravel())
    fitCoefs, covMatrix = opt.curve_fit(gcut, group.H2O2_conc_mM, group.density,bounds=((0, 0, 0), (np.inf,np.inf, np.inf)))
    resids = group.density-group.H2O2_conc_mM.apply(lambda x: gcut(x,*fitCoefs))
    curFit = dict(zip(['mx50','sigx50','scale'],fitCoefs))
    curFit['protein']=name[0]
    curFit['replicate']=name[1]
    curFit['residuals']=sum(resids**2)
    fitData.append(curFit)
fitTable = pd.DataFrame(fitData)
fitTable = fitTable.groupby('protein').mean().reset_index()

protlabels=['pCD79a','pPLC$\gamma$2','pSYK']

x2=np.linspace(0.01,100,200)
x2=np.logspace(-3,3,200)
fig,ax=plt.subplots(1,1,figsize=(5,4),sharey=True)
fig.suptitle('Time = %.0f min'%alltime[0])
for ip,prot in enumerate(allproteins):
    scaledens=(1/fitTable[fitTable.protein==prot].scale.iloc[0])*gcut(x2,fitTable[fitTable.protein==prot].mx50.iloc[0],fitTable[fitTable.protein==prot].sigx50.iloc[0],fitTable[fitTable.protein==prot].scale.iloc[0])
    ax.plot(x2,scaledens,label=prot)
for ip,prot in enumerate(allproteins):
    ax.scatter(sigx50dens[sigx50dens.protein==prot].H2O2_conc_mM,(1/fitTable[fitTable.protein==prot].scale.iloc[0])*sigx50dens[sigx50dens.protein==prot].density,label=None,marker='x')
ax.set_xlabel('$H_2O_2$ conc. (mM)')
ax.set_ylabel('Probability density')
ax.set_xscale('symlog')
ax.legend(protlabels,facecolor='w', framealpha=1, loc='upper right',fancybox=True, shadow=True, ncol=1)


#%% fit DR curve
compoundData = df2.groupby(['protein','time_min','replicate'])
fitDataDR = []
for name,group in compoundData:
    fitCoefs, covMatrix = opt.curve_fit(DRcurve, group.H2O2_conc_mM, group.medFC,bounds=((0, 0, 0, 0), (10,10,100,np.inf)))#(6, 1, 1, 1000)))
    resids = group.medFC-group.H2O2_conc_mM.apply(lambda x: DRcurve(x,*fitCoefs))
    curFit = dict(zip(['n','b','Rmax','EC50'],fitCoefs))
    curFit['protein']=name[0]
    curFit['time_min']=name[1]
    curFit['replicate']=name[2]
    curFit['residuals']=sum(resids**2)
    fitDataDR.append(curFit)
fitCompound = [ item['time_min'] for item in fitDataDR]

fitTableDR = pd.DataFrame(fitDataDR)#.set_index('time_min')
fitTableDR = fitTableDR.groupby('protein').mean().reset_index()

stim=np.linspace(0,100,200)

fig,ax=plt.subplots(1,1,figsize=(5,4),sharex=True,sharey=True)
for ip,prot in enumerate(allproteins):
    ax.set_xscale('symlog')
    p=fitTableDR[fitTable.protein==prot]
    ax.plot(stim,DRcurve(stim,p.n.mean(),p.b.mean(),p.Rmax.mean(),p.EC50.mean()))
    ax.set_xlim([stim[0], stim[-1]])
    ax.set_xlabel('$H_2O_2$ (mM)')
    ax.set_ylabel('Log$_{10}$(Median FC)')
    ax.set_xticks(allconc)
    ax.legend(allproteins,facecolor='w', framealpha=1)
fig.suptitle('Time = %.0f min'%alltime[0])

#%% Test condition for bimodality
for ip,prot in enumerate(allproteins):
    H=fitTableDR.iloc[ip].n
    sigx50=fitTable.iloc[ip].sigx50 
    a=alpha(H,sigx50)
    print([prot, a, H**2*sigx50**2])

print('\n')

conc2=np.linspace(1,10,3001)
for i in range(1,len(allproteins)):
    bimodalconc=[];
    for conc in conc2:#allconc:
        logf=np.log(conc/fitTable.iloc[i].mx50)
        H=fitTableDR.iloc[i].n
        sigx50=fitTable.iloc[i].sigx50 
        a=alpha(H,sigx50)
        if (logf>a[0] and logf<a[1]):
            #print(allproteins[i], conc, a, logf)
            bimodalconc.append(conc)
    print(allproteins[i],np.min(bimodalconc), np.max(bimodalconc))
        


