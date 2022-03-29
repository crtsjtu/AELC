#coding:ISO-8859-1
from __future__ import division
import scipy.stats
import scipy.stats as stats
import sys
import numpy as np
import pandas as pd
import os
import csv
import random
random.seed(9001)
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
import networkx as nx
import math

def variation_mat(frame):
    """
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/***
    Return the variation matrix of frame.
    Element i,j is the variance of the log ratio of components i and j."""
    x = 1.*np.asarray(frame)
    n, m = x.shape
    xx = np.tile(x.reshape((n, m, 1)), (1, 1, m))
    xx_t = xx.transpose(0, 2, 1)
    l = np.log(1.*xx/xx_t)
    #v_mat = l.var(axis=0, ddof=1)
    #col_mean = np.nanmean(l, axis=0)
    #inds = np.where(np.isnan(l))
    #l[inds] = np.take(col_mean, inds[1])
    v_mat = np.nanvar(l,axis=0,ddof=1)
    #print('v_mat',v_mat)
    return(v_mat)
def basis_var(var_mat, m, v_min=1e-10):
    """
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/***
    Estimate the variances of the basis of the compositional data x.
    Assumes that the correlations are sparse (mean correlation is small).
    The element of var_mat are refered to as t_ij in the SparCC paper.
    """
    # compute basis variances
    m_inv = np.linalg.inv(m)
    v_vec = var_mat.sum(axis=1)  # elements are t_i's of SparCC paper
    v_base = np.dot(m_inv, v_vec)  # basis variances.
    v_base[v_base <= 0] = 0.1#min(v_base[v_base > 0])#variance 
    # if any variances are <0 set them to V_min
    #v_base[v_base <= 0]= 1e-10#v_min#
    #print('v_base',v_base)
    return(v_base)
def per_sparcc(frame,samplename):
    #print('frame',frame)
    #sample_table=frame.fillna(0.000000001)
    #sample_table=clr(frame.fillna(0.01))
    sample_table=frame.loc[frame.index!=samplename,]#.fillna(0.0001)
    sample_table2=frame
    #sample_table=clr(frame)

    sample_table.to_csv(sys.argv[1]+'sample_table.csv')
    x = 1.*np.asarray(sample_table)
    sampleSize, nodeD = x.shape
    xx = np.tile(x.reshape((sampleSize,nodeD, 1)), (1, 1, nodeD))
    #print('xx',xx)
    xx_t = xx.transpose(0, 2, 1)
    #print('xx_t',xx_t)
    llog = np.log(1.*xx/xx_t)#logratio matrix [[[D]*D]*S]
    #print('xx/xx_t',xx/xx_t)
    #col_mean = np.nanmean(llog, axis=0)
    #inds = np.where(np.isnan(llog))
    #llog[inds] = np.take(col_mean, inds[1])
    #print('llog',llog)
    X_mean = np.nanmean(llog,axis=0)#[D]*D
    #print('X_mean',X_mean)
    #X_mean = llog.mean(axis=0)
    logdif=(llog-X_mean)**2#response to logratio matrix [[[D]*D]*S]
    #print('logdif',logdif)
    var_mat = variation_mat(sample_table)#variance of each pair (log ratio) [D]*D
    var_mat_temp = var_mat.copy()
    # Make matrix from eqs. 13 of SparCC paper such that: t_i = M * Basis_Varainces
    d = var_mat.shape[0]  # number of OTUs
    m = np.ones((d, d)) + np.diag([d-2]*d)
    # get approx. basis variances and from them basis covariances/correlations
    v_base = basis_var(var_mat_temp, m)
    v_base[v_base <= 0] = min(v_base[v_base > 0])#
    vbase_multi=np.array([v_base.tolist()]*nodeD)
    vbase_multiS =np.array(vbase_multi)
    vbase_multiT =np.array(vbase_multi.T)#**2
    vbase_multiST=np.sqrt(vbase_multiS)*np.sqrt(vbase_multiT)
    #import ipdb;ipdb.set_trace()

    sample_table.corr(method='pearson').fillna(0).to_csv(sys.argv[1]+'_pear.csv')
    cor_form=(vbase_multiT+vbase_multiS-var_mat)/(2*vbase_multiST)
    cor_form=pd.DataFrame(cor_form,columns=sample_table.columns.tolist(),index=sample_table.columns.tolist())
    cor_form.to_csv(sys.argv[1]+'_form.csv')
    
    vbase_multi=np.array([[v_base.tolist()]*nodeD]*sampleSize)
    ##w1w2 #D
    
    w1w22=np.array([np.dot(v_base[:,None],v_base[None,:]).tolist()]*sampleSize)
    w1w22_2=np.array([np.dot(v_base[:,None],v_base[None,:]).tolist()]*(sampleSize+1))
    #print('w1w22_2',w1w22_2)
    ##sumS_multi 
    sumS=np.nansum(logdif,axis=0)##sum across all s size ratio,????????? no use in SSN
    sumS_multi=np.array([sumS.tolist()]*sampleSize)
    ##sumNodeik_multi#sum of all i/j
    sumNodeik=np.nansum(logdif,axis=1)#S*S ???i????????.addition of 1 sample, add 1 dimension
    sumNodeik_zero=np.array([np.zeros([nodeD,nodeD]).tolist()]*sampleSize)
    sumNodeik_zero.transpose()[0,:]=sumNodeik.transpose()
    ones=np.array([np.ones([nodeD,nodeD]).tolist()]*sampleSize)
    sumNodeik_multi=np.array([np.dot(sumNodeik_zero[i],ones[i]).tolist() for i in range(len(ones))])

    sumNodek=np.nansum(sumNodeik,axis=1)#i k
    #print('sumNodeik',sumNodeik)
    #print('sumNodek',sumNodek)
    sumNodek_multi=np.dot(np.diag(sumNodek),np.array([np.ones([sampleSize,nodeD]).tolist()]*nodeD))

    ##sumNodejk_multi;sumNodeik?difference
    sumNodejk_multi=np.array([matr.transpose() for matr in sumNodeik_multi])
    a=2*nodeD-3
    b=1/(2*(nodeD-1)*(nodeD-2))
    cor_content=1/(2*(sampleSize-1)*np.sqrt(w1w22))*(b*(a+1)*sumNodeik_multi-2*b*sumNodek_multi+b*(a+1)*sumNodejk_multi-logdif)
    #add 1 sample: problem, the addition of one sample changes the summuation
    #solution: add 1 to sample size in [ sumNodeik_multis1,sumNodek_multis1,sumNodejk_multis1,logdifs1]
    #print('cor_content',cor_content)
    sampleSize2=sampleSize+1
    x2 = 1.*np.asarray(sample_table2)
    xx2 = np.tile(x2.reshape((sampleSize2,nodeD, 1)), (1, 1, nodeD))
    xx_t2 = xx2.transpose(0, 2, 1)
    llog2 = np.log(1.*xx2/xx_t2)#logratio matrix [[[D]*D]*S]
    logdif2=(llog2-X_mean)**2#response to logratio matrix [[[D]*D]*S]
    ##sumS_multi 
    sumS2=np.nansum(logdif2,axis=0)##sum across all s size ratio,?????????
    sumS_multi2=np.array([sumS2.tolist()]*sampleSize2)
    ##sumNodeik_multi#sum of all i/j
    sumNodeik2=np.nansum(logdif2,axis=1)#S*S ???i????????.addition of 1 sample, add 1 dimension
    sumNodeik_zero2=np.array([np.zeros([nodeD,nodeD]).tolist()]*sampleSize2)
    sumNodeik_zero2.transpose()[0,:]=sumNodeik2.transpose()
    ones=np.array([np.ones([nodeD,nodeD]).tolist()]*sampleSize2)
    sumNodeik_multi2=np.array([np.dot(sumNodeik_zero2[i],ones[i]).tolist() for i in range(len(ones))])

    sumNodek2=np.nansum(sumNodeik2,axis=1)#i k
    sumNodek_multi2=np.dot(np.diag(sumNodek2),np.array([np.ones([sampleSize2,nodeD]).tolist()]*nodeD))

    ##sumNodejk_multi;sumNodeik?difference
    sumNodejk_multi2=np.array([matr.transpose() for matr in sumNodeik_multi2])
    a=2*nodeD-3
    b=1/(2*(nodeD-1)*(nodeD-2))
    cor_content2=1/(2*(sampleSize2-1)*np.sqrt(w1w22_2))*(b*(a+1)*sumNodeik_multi2-2*b*sumNodek_multi2+b*(a+1)*sumNodejk_multi2-logdif2)
    #replace subarray
    #print('llog2',llog2)
    #print('logdif2',logdif2)
    #print('w1w22_2',w1w22_2)
    #print('sumNodeik_multi2',sumNodeik_multi2)
    #print('sumNodek_multi2',sumNodek_multi2)
    #print('sumNodejk_multi2',sumNodejk_multi2)
    #import ipdb;ipdb.set_trace()
    cor_content2=np.append(cor_content,[cor_content2[sampleSize2-1]],axis=0)
    cor_sum2=pd.DataFrame(np.nansum(cor_content2,axis=0),index=sample_table.columns,columns=sample_table.columns)
    cor_sum2.to_csv(sys.argv[1]+'_sum2.csv')
    (sample_table2.corr(method='pearson')-cor_sum2).to_csv(sys.argv[1]+'_minus.csv')
    cor_matr2=cor_sum2.fillna(0)#nx.to_numpy_array(G1, nodelist=None, weight='weight')
    #
    cor_matr_gml2=cor_matr2.mask(np.eye(len(cor_matr2),dtype=bool),0).values
    #Crates graph using the data of the correlation matrix
    G2 = nx.from_numpy_matrix(cor_matr_gml2)
    #relabels the nodes to match the  stocks names
    G2 = nx.relabel_nodes(G2,lambda x: sample_table2.columns[x])
    #shows the edges with their corresponding weights
    G2.edges(data=True)
    nx.write_gml(G2, sys.argv[1]+"_cor2.gml")
    
    #import ipdb;ipdb.set_trace()
    cor_content2.shape = -1, nodeD # reshape it
    onetri2=np.triu(np.ones([nodeD,nodeD]),k=1)>0
    tri_frame2=pd.DataFrame(onetri2.tolist()*sampleSize2,index=sample_table2.columns.tolist()*sampleSize2,columns=sample_table2.columns)
    cor_frame2=pd.DataFrame(cor_content2,index=sample_table2.columns.tolist()*sampleSize2,columns=sample_table2.columns).fillna(0)
    cor_frame2=cor_frame2[tri_frame2]
    #print('cor_frame2',cor_frame2)
    cor_frame2.to_csv(sys.argv[1]+'cor_frame2.csv')
    return(cor_frame2.values)
def sampling_sparcc(frame,samplename, iters=1):
    cor_list = []  # list of cor matrices from different random fractions
    for i in range(iters):
        cor_sparse = per_sparcc(frame,samplename)
        cor_list.append(cor_sparse)
    cor_array = np.array(cor_list)
    cor_med = np.nanmedian(cor_array.astype(float), axis=0)  # median correlations
    cor = pd.DataFrame(cor_med)
    return(cor)
def permute_w_replacement(frame):
    from numpy.random import randint
    s = frame.shape[0]
    fun = lambda x: x.values[randint(0,s,(1,s))][0]
    perm = frame.apply(fun, axis=0)
    return(perm)
    
def bootstrap_correlations(df, cor,samplename, bootstraps=10, procs=1):
    #import ipdb; ipdb.set_trace()
    # take absolute value of all values in cor for calculating two-sided p-value
    #abs_cor = np.abs(squareform(cor, checks=False))#distance vecter between point
    abs_cor = np.abs(cor)
    # create an empty array of significant value counts in same shape as abs_cor
    n_sig = pd.DataFrame(np.zeros(abs_cor.shape))
    
    for i in range(bootstraps):#bootstraps 100次
        bootstrap = permute_w_replacement(df)#permuted df
        perm_cor = sampling_sparcc(bootstrap,samplename)
        #import ipdb; ipdb.set_trace()
        print i
        n_sig[np.abs(perm_cor) >= np.abs(cor)]+= 1#f

    # get p_values out
    p_vals = 1.*n_sig/bootstraps#
    #import ipdb; ipdb.set_trace()
    #p_vals.values[np.diag_indices_from(p_vals.values)] = 1 
    return(p_vals)
toy_g=pd.read_csv(sys.argv[1], sep="\t", header=0,index_col=0,encoding='utf8').T
samplename=sys.argv[1].split('_')[2]
# if contain Na
#toy_g[toy_g==0]=None#for relative abundance
#toy_g = toy_g.fillna(1)/toy_g.fillna(1).sum()#toy_g=toy_g.div(toy_g.sum(axis=1), axis=0)#relative abundance
#toy_g=toy_g.fillna(toy_g.mean())
#toy_g=toy_g.fillna(0)
#toy_g=toy_g.loc[:,toy_g.std()>0].T

Framecor=sampling_sparcc(toy_g,samplename)
Framecor.to_csv(sys.argv[1]+'Framecor10.csv')
perm_pval=bootstrap_correlations(toy_g, Framecor,samplename,bootstraps=1000)
perm_pval.to_csv(sys.argv[1]+'perm_pval10.csv')
#import ipdb;ipdb.set_trace()
sparcc_corfilter=Framecor[perm_pval<=0.05]
sparcc_corfilter.columns=toy_g.columns
sampleSize, nodeD = toy_g.shape
sparcc_corfilter.index=toy_g.columns.tolist()*sampleSize
sparcc_corfilter['samples']=np.repeat(toy_g.index.tolist(),nodeD)
#sparcc_corfilter.to_csv(sys.argv[1]+'meancorP.csv')
#Framecor.index=toy_g.columns.tolist()*sampleSize
#Framecor['samples']=np.repeat(toy_g.index.tolist(),nodeD)
#Framecor.to_csv(sys.argv[1]+'_var.csv',index=0)
sparcc_cor_list=pd.melt(sparcc_corfilter.reset_index(),id_vars=['samples','index'])
sparcc_cor_list=sparcc_cor_list[sparcc_cor_list['value'].notnull()]
sparcc_cor_list.columns=['samples','index','variable','value']
sparcc_cor_list.to_csv(sys.argv[1]+'_comselP.csv',index=0)
#full comsel_P
sparcc_corfilter=Framecor#[perm_pval<=0.05]
sparcc_corfilter.columns=toy_g.columns
sparcc_corfilter.index=toy_g.columns.tolist()*sampleSize
sparcc_corfilter['samples']=np.repeat(toy_g.index.tolist(),nodeD)
sparcc_cor_list=pd.melt(sparcc_corfilter.reset_index(),id_vars=['samples','index'])
sparcc_cor_list=sparcc_cor_list[sparcc_cor_list['value'].notnull()]
sparcc_cor_list.columns=['samples','index','variable','value']
sparcc_cor_list.to_csv(sys.argv[1]+'_comsel.csv',index=0)


