"""Export manuscript figures from the supplied prediction arrays."""
from pathlib import Path
import json, shutil
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, precision_recall_curve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from .features import ROOT, DATA
from .evaluation import valid_coordinates
BLUE='#0072B2'; ORANGE='#D55E00'; GREY='#93999F'; DARK='#30363D'
RESULTS=ROOT/'reference_results'
def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,'axes.titlesize':11,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.unicode_minus':False,'axes.linewidth':.8,'hatch.linewidth':.5})
def clean(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.set_axisbelow(True);ax.tick_params(length=3,width=.7)
def functional_evidence():
    scores=np.load(RESULTS/'full_scores.npz');d=pd.read_csv(DATA/'fivecell.csv.gz');mask=valid_coordinates(d)
    fig=plt.figure(figsize=(6.5,5.0))
    a=fig.add_axes([.10,.63,.38,.26]);b=fig.add_axes([.67,.63,.30,.26])
    c=fig.add_axes([.10,.13,.38,.29]);e=fig.add_axes([.67,.13,.30,.29])
    rows=[];curves=[]
    for name,label,col,ls in [('full','Atlas-E2G',BLUE,'-'),('published','rE2G',ORANGE,'--'),('matched','rE2G-L2',GREY,':')]:
        p=scores['fivecell/'+name].copy();p[~d.mapped.to_numpy()]=0
        ap=average_precision_score(d.truth[mask],p[mask],sample_weight=d.weight[mask])
        precision,recall,_=precision_recall_curve(d.truth[mask],p[mask],sample_weight=d.weight[mask])
        a.plot(recall,precision,color=col,ls=ls,lw=1.5,label=f'{label}  {ap:.3f}')
        ix=np.flatnonzero(mask);top=ix[np.argsort(-p[mask],kind='stable')[:100]]
        rows.append(dict(analysis='fivecell',method=name,AP=ap,top100=int(d.truth.iloc[top].sum()),n=int(mask.sum())))
        curves.append(pd.DataFrame(dict(method=name,recall=recall,precision=precision)))
    assert [r['top100'] for r in rows[:2]]==[77,65]
    a.set(xlim=(0,1),ylim=(0,1.02),xlabel='Recall',ylabel='Precision',xticks=[0,.5,1],yticks=[0,.5,1])
    a.legend(loc='lower left',frameon=False,handlelength=1.2,handletextpad=.4,labelspacing=.25);clean(a)
    for i,r in enumerate(rows[:2]):
        b.bar(i,r['top100'],.58,color=[BLUE,ORANGE][i],edgecolor=DARK,lw=.5,hatch=['','//'][i])
        b.text(i,r['top100']+3,str(r['top100']),ha='center',weight='bold')
    b.set(ylim=(0,100),xlim=(-.65,1.65),xticks=[0,1],xticklabels=['Atlas-E2G','rE2G'],ylabel='Positives in top 100',yticks=[0,25,50,75,100])
    b.tick_params(axis='x',length=0);b.grid(axis='y',color='#E3E6E8',lw=.6);clean(b)
    q=pd.read_csv(DATA/'reporter.csv.gz');x=scores['reporter/full'];y=q.response.to_numpy()
    for lib,marker,col in [('500P','o',BLUE),('78P','^',ORANGE)]:
        sel=q.library.eq(lib).to_numpy()
        c.scatter(x[sel],y[sel],s=15,alpha=.65,marker=marker,color=col,lw=.3,edgecolors='white',label=lib)
    c.set(xlabel='Promoter component',ylabel='Measured responsiveness',xlim=(-.88,.65),ylim=(-.22,4.4),xticks=[-.75,-.25,.25],yticks=[0,1,2,3,4])
    c.legend(loc='upper left',frameon=False,handletextpad=.2,labelspacing=.2,borderaxespad=.2);c.text(.98,.04,'ρ = 0.654',transform=c.transAxes,ha='right');clean(c)
    ci=pd.read_csv(RESULTS/'figure_intervals.csv').set_index('method')
    for i,(name,col,hatch) in enumerate(zip(['full','shared','without_identity'],[BLUE,GREY,'#D4D8DB'],['','//','..'])):
        r=ci.loc[name];assert abs(spearmanr(scores['reporter/'+name],y).statistic-r.value)<1e-12
        e.bar(i,r.value,.58,color=col,edgecolor=DARK,lw=.55,hatch=hatch)
        e.errorbar(i,r.value,yerr=[[r.value-r.low],[r.high-r.value]],fmt='none',ecolor=DARK,capsize=3,lw=1)
        e.text(i,r.high+.025,f'{r.value:.3f}',ha='center')
        rows.append(dict(analysis='reporter',method=name,rho=r.value,low=r.low,high=r.high,n=len(q)))
    e.set(ylim=(0,.86),xlim=(-.6,2.6),yticks=[0,.2,.4,.6,.8],xticks=range(3),xticklabels=['Separate\nroles','Shared\nroles','Identity\nremoved'],ylabel='Spearman correlation')
    e.tick_params(axis='x',length=0,pad=6);e.grid(axis='y',color='#E3E6E8',lw=.6);clean(e)
    for xx,yy,t in [(.10,.945,'a  Functional-pair ranking'),(.62,.945,'b  Observed top-100 hits'),(.10,.48,'c  Reporter responsiveness'),(.62,.48,'d  Regional contributions')]:fig.text(xx,yy,t,weight='bold',fontsize=11)
    return fig,pd.DataFrame(rows),pd.concat(curves,ignore_index=True)
def public_prediction():
    public=np.load(RESULTS/'public_scores.npz');full=np.load(RESULTS/'full_scores.npz')
    cases=[('fivecell','non_K562'),('fivecell','non_K562_distal'),('astrocyte','covered')];rows=[]
    methods=['baseline','K562_for_all_cells','complete','full_input']
    for kind,subset in cases:
        d=pd.read_csv(DATA/f'{kind}.csv.gz');mask=valid_coordinates(d)&d.mapped.to_numpy()
        if subset.startswith('non_K562'):mask&=d.cell.ne('K562').to_numpy()
        if subset.endswith('distal'):mask&=d.benchmark_distance.gt(100000).to_numpy()
        groups=[list(ix) for ix in d[mask].groupby(['cell','study','gene']).groups.values() if len(ix)>=4]
        for method in methods:
            p=full[kind+'/full'] if method=='full_input' else public['activity/'+kind+'/'+method]
            hits=sum(int(d.truth.iloc[sorted(ix,key=lambda i:(-p[i],i))[:3]].sum()) for ix in groups)
            rows.append(dict(dataset=kind,subset=subset,method=method,AP=average_precision_score(d.truth[mask],p[mask],sample_weight=d.weight[mask] if kind=='fivecell' else None),hits=hits,slots=3*len(groups),ceiling=sum(min(3,int(d.truth.iloc[ix].sum())) for ix in groups),n=int(mask.sum())))
    values=pd.DataFrame(rows);matrix=values.AP.to_numpy().reshape(3,4)
    fig=plt.figure(figsize=(6.5,4.5));a=fig.add_axes([.32,.55,.65,.24])
    cmap=LinearSegmentedColormap.from_list('public_ap',['#F2F6F9','#AFCEDD','#0072B2'])
    a.pcolormesh(np.arange(5)-.5,np.arange(4)-.5,matrix,cmap=cmap,vmin=0,vmax=1,edgecolors='white',linewidth=2)
    a.set_ylim(2.5,-.5)
    for i in range(3):
        for j in range(4):a.text(j,i,f'{matrix[i,j]:.3f}',ha='center',va='center',color=DARK,fontsize=11,weight='bold')
    a.set(yticks=range(3),yticklabels=['Non-K562','Non-K562 >100 kb','Astrocytes'],xticks=range(4),xticklabels=['Gene\nonly','K562 for\nall cells','Matched\npublic cell','Full\ninputs'])
    a.tick_params(axis='both',length=0,pad=7);a.xaxis.tick_top()
    for sp in a.spines.values():sp.set_visible(False)
    a.set_xticks(np.arange(-.5,4,1),minor=True);a.set_yticks(np.arange(-.5,3,1),minor=True)
    a.grid(which='minor',color='white',lw=3);a.tick_params(which='minor',length=0)
    fig.text(.32,.49,'Average precision · identical covered pairs',fontsize=10,color=DARK)
    q=json.loads((RESULTS/'public_metrics.json').read_text())
    ci={r['method']:r for r in q['control_reporter'] if r['subset']=='all' and r['metric']=='Spearman'}
    reporter=[('Unconstrained refit',ci['cell/pooled_hard']),('Output distillation',ci['cell/pooled_distilled']),('Retained promoter',ci['anchored_promoter'])]
    assert abs(spearmanr(public['transfer/reporter/standard'],pd.read_csv(DATA/'reporter.csv.gz').response).statistic-ci['anchored_promoter']['value'])<1e-12
    b=fig.add_axes([.36,.12,.61,.23])
    for i,(name,r) in enumerate(reporter):
        b.barh(i,r['value'],.60,color=BLUE if i==2 else GREY,edgecolor=DARK,lw=.5,hatch=['//','..',''][i])
        b.errorbar(r['value'],i,xerr=[[r['value']-r['low']],[r['high']-r['value']]],fmt='none',ecolor=DARK,capsize=3,lw=1)
        b.text(r['high']+.025,i,f"{r['value']:.3f}",va='center')
    b.set(xlim=(0,.85),ylim=(-.6,2.6),yticks=range(3),yticklabels=[r[0] for r in reporter],xticks=[0,.2,.4,.6,.8],xlabel='Spearman correlation with reporter responsiveness')
    b.tick_params(axis='y',length=0,pad=7);b.grid(axis='x',color='#E3E6E8',lw=.6);clean(b)
    fig.text(.04,.96,'a  Public references supply missing activity information',weight='bold',fontsize=11)
    fig.text(.04,.405,'b  Component preservation retains the regional association',weight='bold',fontsize=11)
    transfer=pd.DataFrame([dict(label=name,**r) for name,r in reporter])
    return fig,values,transfer
def run(output,formats=('png','pdf','svg')):
    style();out=Path(output);out.mkdir(parents=True,exist_ok=True)
    f2,v2,curves=functional_evidence();f3,v3,transfer=public_prediction()
    v2.to_csv(out/'figure2_values.csv',index=False);curves.to_csv(out/'figure2_pr_curves.csv',index=False)
    v3.to_csv(out/'figure3_values.csv',index=False);transfer.to_csv(out/'figure3_transfer.csv',index=False)
    shutil.copyfile(ROOT/'docs/overview.svg',out/'figure1.svg')
    for i,fig in enumerate([f2,f3],2):
        for ext in formats:fig.savefig(out/f'figure{i}.{ext}',dpi=600)
    return [f2,f3]

