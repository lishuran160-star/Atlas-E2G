"""Coordinate eligibility and covariate-adjusted reporter correlation."""
import numpy as np
from scipy.stats import rankdata,pearsonr


def partial(q,values,amount):
    controls=q[['GC','CpG','length','ubiquitous','basal']].to_numpy()
    controls=np.c_[controls,q.library.eq('78P').astype(int),amount]
    Z=np.c_[np.ones(len(q)),np.apply_along_axis(rankdata,0,controls)]
    x=rankdata(values);y=rankdata(q.response)
    x-=Z@np.linalg.lstsq(Z,x,rcond=None)[0];y-=Z@np.linalg.lstsq(Z,y,rcond=None)[0]
    return float(pearsonr(x,y).statistic)

def valid_coordinates(d):
    canonical=[f'chr{i}' for i in range(1,23)]+['chrX','chrY']
    return (d.chrom.isin(canonical)&d.pchrom.eq(d.chrom)&d.start.ge(0)&d.end.gt(d.start)&d.tss.ge(1000)).to_numpy()
