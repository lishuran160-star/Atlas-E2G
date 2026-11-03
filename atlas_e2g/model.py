"""Regularized logistic fitting and grouped validation splits."""
import warnings
from scipy.special import expit
from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import ConvergenceWarning


def splits(d):
    mask=d.mapped.to_numpy()
    for chrom in sorted(d.chrom.unique()):
        yield 'chromosome',chrom,mask&d.chrom.eq(chrom).to_numpy(),mask&d.chrom.ne(chrom).to_numpy()
    for study in sorted(d.study.unique()):
        test=mask&d.study.eq(study).to_numpy()
        train=(mask&d.study.ne(study).to_numpy()&~d.gene.isin(d.loc[test,'gene']).to_numpy()
               &~d.enhancer_uid.isin(d.loc[test,'enhancer_uid']).to_numpy())
        yield 'study_gene_element_disjoint',study,test,train

def fit_logistic(x,y,unpenalized=False):
    kwargs={'penalty':None,'tol':1e-8,'max_iter':20000} if unpenalized else {'C':.1,'max_iter':10000}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        m=LogisticRegression(random_state=42,**kwargs).fit(x,y)
    if any(issubclass(w.category,ConvergenceWarning) for w in caught):
        raise RuntimeError('Logistic fit did not converge')
    return {'coef':m.coef_[0],'intercept':float(m.intercept_[0])}


def score(m,x):return expit(x@m['coef']+m['intercept'])

