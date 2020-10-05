#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import sys,os
import copy
import numpy as np
import pandas as pd
import pylab as py
import matplotlib
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
from  matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

def smooth(X,Y):
    x=[X[j] for j in range(X.size) if np.isnan(Y[j])==False]
    y=[Y[j] for j in range(Y.size) if np.isnan(Y[j])==False]
    f = interp1d(x, y,fill_value="extrapolate")
    for i in range(1,X.size-1):
      if np.isnan(Y[i]): Y[i]=f(X[i])

    for i in range(X.size):
        x=[X[j] for j in range(X.size) if j!=i]
        y=[Y[j] for j in range(Y.size) if j!=i]
        f = interp1d(x, y,fill_value="extrapolate")
        rel_err=abs(f(X[i])-Y[i])/f(X[i])
        if rel_err>0.05: Y[i]=f(X[i])
    return Y

def plot_compass17():

  data=pd.read_excel('expdata/compass17.xlsx')


  bins={}
  bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  bins[(0,7)]="17<Q2 and 0.2<x "
  bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  bins[(1,5)]="6<Q2 and Q2<17 and 0.06<x  and x<0.09 "
  bins[(1,6)]="6<Q2 and Q2<17 and 0.09<x  and x<0.2 "
  bins[(1,7)]="6<Q2 and Q2<17 and 0.2<x "
  bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  bins[(2,3)]="3<Q2 and Q2<6 and 0.02<x  and x<0.03 "
  bins[(2,4)]="3<Q2 and Q2<6 and 0.03<x  and x<0.06 "
  bins[(2,5)]="3<Q2 and Q2<6 and 0.06<x  and x<0.09 "
  bins[(2,6)]="3<Q2 and Q2<6 and 0.09<x  and x<0.2 "
  bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  bins[(3,1)]="1.5<Q2 and Q2<3.0 and 0.009<x and x<0.012 "
  bins[(3,2)]="1.5<Q2 and Q2<3.0 and 0.012<x and x<0.02 "
  bins[(3,3)]="1.5<Q2 and Q2<3.0 and 0.02<x  and x<0.03 "
  bins[(3,4)]="1.5<Q2 and Q2<3.0 and 0.03<x  and x<0.06 "
  bins[(3,5)]="1.5<Q2 and Q2<3.0 and 0.06<x  and x<0.09 "
  bins[(4,0)]="Q2<1.5 and x<0.009 "
  bins[(4,1)]="Q2<1.5 and 0.009<x and x<0.012 "
  bins[(4,2)]="Q2<1.5 and 0.012<x and x<0.02 "
  bins[(4,3)]="Q2<1.5 and 0.02<x  and x<0.03 "
  bins[(4,4)]="Q2<1.5 and 0.03<x  and x<0.06 "

  Q2bins={}
  Q2bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  Q2bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  Q2bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  Q2bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  Q2bins[(4,0)]="Q2<1.5 and x<0.009 "
  Q2b=[]
  for k in Q2bins:
      d=data.query('%s and  had=="h+"'%Q2bins[k])
      Q2b.append(d.Q2.mean())
  Q2b=np.sort(np.unique(Q2b))

  xbins={}
  xbins[(4,0)]="x<0.009 "
  xbins[(4,1)]="0.009<x and x<0.012 "
  xbins[(4,2)]="0.012<x and x<0.02 "
  xbins[(4,3)]="0.02<x  and x<0.03 "
  xbins[(4,4)]="0.03<x  and x<0.06 "
  xbins[(2,3)]="0.02<x  and x<0.03 "
  xbins[(2,4)]="0.03<x  and x<0.06 "
  xbins[(2,5)]="0.06<x  and x<0.09 "
  xbins[(2,6)]="0.09<x  and x<0.2 "
  xbins[(1,7)]="0.2<x "
  
  xb=[]
  for k in xbins:
      d=data.query('%s and  had=="h+"'%xbins[k])
      xb.append(d.x.mean())
  xb=np.sort(np.unique(xb))


  zbins=[[0.24,0.3],[0.3,0.4],[0.4,0.5],[0.65,0.7]]

  nrows,ncols=5,8
  fig = py.figure(figsize=(ncols*3,nrows*3))
  gs = gridspec.GridSpec(5,8)
  gs.update(wspace=0.,hspace=0,left=0.12, right=0.96,bottom=0.13,top=0.96)
  AX={}
  for k in bins:
      ir,ic=k
      ax = py.subplot(gs[ir,ic])
      d=data.query('%s and  had=="h+"'%bins[k])
      for i in range(len(zbins)):
          dd=d.query('%f<z and z<%f'%(zbins[i][0],zbins[i][1]))
          if dd.index.size==0: continue
        
          #if   i==0: c='r'
          #elif i==1: c='b'
          #elif i==2: c='g'
          #elif i==3: c='m'

          e=ax.errorbar(dd.qT,dd.value,dd.alpha,fmt='.',label=r'$%0.2f<z<%0.2f$'%(zbins[i][0],zbins[i][1]))
          c=e[0].get_color()

          ax.fill_between(dd.qT,dd.LOmin/dd.idisNLO,dd.LOmax/dd.idisNLO,color=c,alpha=0.1)



          ax.plot(dd.qT,dd.DDSLO/dd.idisNLO,color=c,ls=':')

          tot=dd.DDSLO+dd.DDSNLOdelta+dd.DDSNLOplus+dd.DDSNLOrest
          tot=tot.values
          #for itot in range(tot.size):
          #  if np.isnan(tot[itot]): 
          #    tot[itot]=tot[itot+1]*5.5
          tot=smooth(dd.qT.values,tot)
          ax.plot(dd.qT,tot/dd.idisNLO,color=c,ls='-',label='')
          #ax.plot(dd.qT,dd.LO/dd.idisNLO,'%s--'%c)
          #ax.plot(dd.qT,tot/dd.idisNLO,'%s-'%c)
  
          #ax.plot(dd.qT,tot/dd.idisNLO,'%s:'%c,label='')
      #ax.axvline(d.Q.values[0])
      ax.semilogy()
      ax.set_ylim(1e-4,10)
      ax.set_xlim(0,8)
      #ax.set_yticklabels([])
  
      ax.tick_params(axis='both', which='major', labelsize=30, direction='in')
      if all(k!=_ for _ in [(4,0),(4,1),(4,2),(4,3),(4,4),(3,5),(2,6),(1,7)]): 
          ax.set_xticklabels([])
      else:
          ax.set_xticks([2,4,6])
          #ax.set_xlabel(r'$q_{\rm T}$',size=35)
          #ax.xaxis.set_label_coords(0.95, -0.02)
      if all(k!=_ for _ in [(4,0),(3,0),(2,2),(1,4),(0,6)]): 
          ax.set_yticklabels([])
          ax.set_yticks([1e-3,1e-2,1e-1,1])
      else:
          ax.set_yticks([1e-3,1e-2,1e-1,1])
      
      if k==(4,0):
          ax.annotate('', xy=(-0.35, 5.2), xycoords='axes fraction', xytext=(-0.35, -0.1),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))
          ax.annotate('', xy=(8.2,-0.3), xycoords='axes fraction', xytext=(-0.1, -0.3),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))        
          ax.annotate(r'$Q^2~({\rm GeV}^2)$', xy=(-1,3),xycoords='axes fraction',size=50,rotation=90)
          ax.annotate(r'$x_{\rm bj}$', xy=(3.9,-0.7),xycoords='axes fraction',size=50)
                  
          for i in range(8):
              if xb[i]<2e-2: msg=r'$%0.3f$'%xb[i]
              else:msg=r'$%0.2f$'%xb[i]
              ax.text(0.5+i,-0.5,msg,transform=ax.transAxes,size=30,ha="center")
              ax.annotate('',xy=(i,-0.35),xycoords='axes fraction',xytext=(i+1, -0.35), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
          for i in range(5):
              ax.text(-0.6,0.5+i,r'$%0.1f$'%Q2b[i],transform=ax.transAxes,size=30,rotation=90,va="center")
              ax.annotate('',xy=(-0.4,i),xycoords='axes fraction',xytext=(-0.4,i+1), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
      ax.plot([d.Q.values[0],8],[2e-4,2e-4],c='y',lw=10,alpha=0.5)
      if k==(2,2): 
          LO1b = mpatches.Rectangle((0,0), 0.05, 0.1, ec="none",color='k',alpha=0.1)
          LO1l=py.Line2D([0,2], [0,0], color='k',ls=':')
          LO2=py.Line2D([0], [0], color='k',ls='-')
          qTrange = mpatches.Rectangle((0,0), 0, 0, ec="none",color='y',alpha=0.5)
          #qTrange=py.Line2D([0], [0], color='y',ls='-',lw=10,alpha=0.5)
          ax.legend([(LO1b,LO1l),LO2,qTrange],[r'$\rm DDS~(LO)$',r'$\rm DDS~(NLO)$',r'$q_{\rm T}>Q$']\
                    ,bbox_to_anchor=[-1.2, 1.]\
                    ,loc='center',fontsize=40,frameon=0)
          msg=r'${\rm COMPASS~17}~h^+$'
          ax.text(-2,2.8,msg,transform=ax.transAxes,size=50)
          msg =r'$\frac{d\sigma}{dx_{\rm bj}dQ^2dzdP^2_{\rm T}}/\frac{d\sigma}{dx_{\rm bj}dQ^2}'
          msg+=r'({\tiny {\rm GeV}^{-2}})'
          msg+=r'~{\rm vs.}~q_{\rm T}{~\rm (GeV)}$'
          ax.text(-2,2.2,msg,transform=ax.transAxes,size=50)
      if k==(1,4): 
          ax.legend(bbox_to_anchor=[3, -2.5], loc='center',fontsize=40,frameon=0\
                   ,markerscale=2,handletextpad=0.1)
      #ax.text(0.6,0.8,'(%d,%d)'%(ir,ic),transform=ax.transAxes,size=20)
  
          
  py.savefig('gallery/compass17.pdf')

def plot_compass17_K():

  data=pd.read_excel('expdata/compass17.xlsx')


  bins={}
  bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  bins[(0,7)]="17<Q2 and 0.2<x "
  bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  bins[(1,5)]="6<Q2 and Q2<17 and 0.06<x  and x<0.09 "
  bins[(1,6)]="6<Q2 and Q2<17 and 0.09<x  and x<0.2 "
  bins[(1,7)]="6<Q2 and Q2<17 and 0.2<x "
  bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  bins[(2,3)]="3<Q2 and Q2<6 and 0.02<x  and x<0.03 "
  bins[(2,4)]="3<Q2 and Q2<6 and 0.03<x  and x<0.06 "
  bins[(2,5)]="3<Q2 and Q2<6 and 0.06<x  and x<0.09 "
  bins[(2,6)]="3<Q2 and Q2<6 and 0.09<x  and x<0.2 "
  bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  bins[(3,1)]="1.5<Q2 and Q2<3.0 and 0.009<x and x<0.012 "
  bins[(3,2)]="1.5<Q2 and Q2<3.0 and 0.012<x and x<0.02 "
  bins[(3,3)]="1.5<Q2 and Q2<3.0 and 0.02<x  and x<0.03 "
  bins[(3,4)]="1.5<Q2 and Q2<3.0 and 0.03<x  and x<0.06 "
  bins[(3,5)]="1.5<Q2 and Q2<3.0 and 0.06<x  and x<0.09 "
  bins[(4,0)]="Q2<1.5 and x<0.009 "
  bins[(4,1)]="Q2<1.5 and 0.009<x and x<0.012 "
  bins[(4,2)]="Q2<1.5 and 0.012<x and x<0.02 "
  bins[(4,3)]="Q2<1.5 and 0.02<x  and x<0.03 "
  bins[(4,4)]="Q2<1.5 and 0.03<x  and x<0.06 "

  Q2bins={}
  Q2bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  Q2bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  Q2bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  Q2bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  Q2bins[(4,0)]="Q2<1.5 and x<0.009 "
  Q2b=[]
  for k in Q2bins:
      d=data.query('%s and  had=="h+"'%Q2bins[k])
      Q2b.append(d.Q2.mean())
  Q2b=np.sort(np.unique(Q2b))

  xbins={}
  xbins[(4,0)]="x<0.009 "
  xbins[(4,1)]="0.009<x and x<0.012 "
  xbins[(4,2)]="0.012<x and x<0.02 "
  xbins[(4,3)]="0.02<x  and x<0.03 "
  xbins[(4,4)]="0.03<x  and x<0.06 "
  xbins[(2,3)]="0.02<x  and x<0.03 "
  xbins[(2,4)]="0.03<x  and x<0.06 "
  xbins[(2,5)]="0.06<x  and x<0.09 "
  xbins[(2,6)]="0.09<x  and x<0.2 "
  xbins[(1,7)]="0.2<x "
  
  xb=[]
  for k in xbins:
      d=data.query('%s and  had=="h+"'%xbins[k])
      xb.append(d.x.mean())
  xb=np.sort(np.unique(xb))


  zbins=[[0.24,0.3],[0.3,0.4],[0.4,0.5],[0.65,0.7]]

  nrows,ncols=5,8
  fig = py.figure(figsize=(ncols*3,nrows*3))
  gs = gridspec.GridSpec(5,8)
  gs.update(wspace=0.,hspace=0,left=0.12, right=0.96,bottom=0.13,top=0.96)
  AX={}
  for k in bins:
      ir,ic=k
      ax = py.subplot(gs[ir,ic])
      d=data.query('%s and  had=="h+"'%bins[k])
      for i in range(len(zbins)):
          dd=d.query('%f<z and z<%f'%(zbins[i][0],zbins[i][1]))
          if dd.index.size==0: continue
        
          if   i==0: c='r'
          elif i==1: c='b'
          elif i==2: c='g'
          elif i==3: c='m'

          #e=ax.errorbar(dd.qT,dd.value,dd.alpha,fmt='.',label=r'$%0.2f<z<%0.2f$'%(zbins[i][0],zbins[i][1]))
          #c=e[0].get_color()

          #ax.fill_between(dd.qT,dd.LOmin/dd.idisNLO,dd.LOmax/dd.idisNLO,color=c,alpha=0.1)

          #ax.plot(dd.qT,dd.DDSLO/dd.idisNLO,color=c,ls=':')

          tot=dd.DDSLO+dd.DDSNLOdelta+dd.DDSNLOplus+dd.DDSNLOrest
          tot=tot.values
          #for itot in range(tot.size):
          #  if np.isnan(tot[itot]): 
          #    tot[itot]=tot[itot+1]*5.5
          K=smooth(dd.qT.values,tot/dd.DDSLO.values)
          #ax.plot(dd.qT,tot/dd.idisNLO,color=c,ls='-',label='')
          #ax.plot(dd.qT,dd.LO/dd.idisNLO,'%s--'%c)
          ax.plot(dd.qT,K,'%so-'%c,label=r'$%0.2f<z<%0.2f$'%(zbins[i][0],zbins[i][1]))
  
          #ax.plot(dd.qT,tot/dd.idisNLO,'%s:'%c,label='')
      #ax.axvline(d.Q.values[0])
      #ax.semilogy()
      #ax.set_ylim(1e-4,10)
      #ax.set_xlim(0,8)
      ax.set_xlim(0,1.8)
      #ax.set_yticklabels([])
  
      ax.tick_params(axis='both', which='major', labelsize=30, direction='in')
      if all(k!=_ for _ in [(4,0),(4,1),(4,2),(4,3),(4,4),(3,5),(2,6),(1,7)]): 
          ax.set_xticklabels([])
      else:
          pass
          #ax.set_xticks([2,4,6])
          #ax.set_xlabel(r'$q_{\rm T}$',size=35)
          #ax.xaxis.set_label_coords(0.95, -0.02)
      if all(k!=_ for _ in [(4,0),(3,0),(2,2),(1,4),(0,6)]): 
          ax.set_yticklabels([])
          ax.set_yticks([1e-3,1e-2,1e-1,1])
      else:
          pass
          #ax.set_yticks([1e-3,1e-2,1e-1,1])
          ax.set_ylim(1,2.5)     
 
      if k==(4,0):
          ax.annotate('', xy=(-0.35, 5.2), xycoords='axes fraction', xytext=(-0.35, -0.1),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))
          ax.annotate('', xy=(8.2,-0.3), xycoords='axes fraction', xytext=(-0.1, -0.3),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))        
          ax.annotate(r'$Q^2~({\rm GeV}^2)$', xy=(-1,3),xycoords='axes fraction',size=50,rotation=90)
          ax.annotate(r'$x_{\rm bj}$', xy=(3.9,-0.7),xycoords='axes fraction',size=50)
                  
          for i in range(8):
              if xb[i]<2e-2: msg=r'$%0.3f$'%xb[i]
              else:msg=r'$%0.2f$'%xb[i]
              ax.text(0.5+i,-0.5,msg,transform=ax.transAxes,size=30,ha="center")
              ax.annotate('',xy=(i,-0.35),xycoords='axes fraction',xytext=(i+1, -0.35), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
          for i in range(5):
              ax.text(-0.6,0.5+i,r'$%0.1f$'%Q2b[i],transform=ax.transAxes,size=30,rotation=90,va="center")
              ax.annotate('',xy=(-0.4,i),xycoords='axes fraction',xytext=(-0.4,i+1), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
      #ax.plot([d.Q.values[0],8],[2e-4,2e-4],c='y',lw=10,alpha=0.5)
      #if k==(2,2): 
      #    LO1b = mpatches.Rectangle((0,0), 0.05, 0.1, ec="none",color='k',alpha=0.1)
      #    LO1l=py.Line2D([0,2], [0,0], color='k',ls=':')
      #    LO2=py.Line2D([0], [0], color='k',ls='-')
      #    qTrange = mpatches.Rectangle((0,0), 0, 0, ec="none",color='y',alpha=0.5)
      #    #qTrange=py.Line2D([0], [0], color='y',ls='-',lw=10,alpha=0.5)
      #    ax.legend([(LO1b,LO1l),LO2,qTrange],[r'$\rm DDS~(LO)$',r'$\rm DDS~(NLO)$',r'$q_{\rm T}>Q$']\
      #              ,bbox_to_anchor=[-1.2, 1.]\
      #              ,loc='center',fontsize=40,frameon=0)
      #    msg=r'${\rm COMPASS~17}~h^+$'
      #    ax.text(-2,2.8,msg,transform=ax.transAxes,size=50)
      #    msg =r'$\frac{d\sigma}{dx_{\rm bj}dQ^2dzdP^2_{\rm T}}/\frac{d\sigma}{dx_{\rm bj}dQ^2}'
      #    msg+=r'({\tiny {\rm GeV}^{-2}})'
      #    msg+=r'~{\rm vs.}~q_{\rm T}{~\rm (GeV)}$'
      #    ax.text(-2,2.2,msg,transform=ax.transAxes,size=50)
      if k==(1,4): 
          ax.legend(bbox_to_anchor=[3, -2.5], loc='center',fontsize=40,frameon=0\
                   ,markerscale=2,handletextpad=0.1)
      #ax.text(0.6,0.8,'(%d,%d)'%(ir,ic),transform=ax.transAxes,size=20)
  
          
  py.savefig('gallery/compass17_K.pdf')

def plot_compass17_rat():

  data=pd.read_excel('expdata/compass17.xlsx')

  bins={}
  bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  bins[(0,7)]="17<Q2 and 0.2<x "#
  bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  bins[(1,5)]="6<Q2 and Q2<17 and 0.06<x  and x<0.09 "#
  bins[(1,6)]="6<Q2 and Q2<17 and 0.09<x  and x<0.2 " #
  bins[(1,7)]="6<Q2 and Q2<17 and 0.2<x "
  bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  bins[(2,3)]="3<Q2 and Q2<6 and 0.02<x  and x<0.03 "
  bins[(2,4)]="3<Q2 and Q2<6 and 0.03<x  and x<0.06 "
  bins[(2,5)]="3<Q2 and Q2<6 and 0.06<x  and x<0.09 "
  bins[(2,6)]="3<Q2 and Q2<6 and 0.09<x  and x<0.2 " #
  bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  bins[(3,1)]="1.5<Q2 and Q2<3.0 and 0.009<x and x<0.012 "
  bins[(3,2)]="1.5<Q2 and Q2<3.0 and 0.012<x and x<0.02 "
  bins[(3,3)]="1.5<Q2 and Q2<3.0 and 0.02<x  and x<0.03 "
  bins[(3,4)]="1.5<Q2 and Q2<3.0 and 0.03<x  and x<0.06 "
  bins[(3,5)]="1.5<Q2 and Q2<3.0 and 0.06<x  and x<0.09 "
  bins[(4,0)]="Q2<1.5 and x<0.009 "
  bins[(4,1)]="Q2<1.5 and 0.009<x and x<0.012 "
  bins[(4,2)]="Q2<1.5 and 0.012<x and x<0.02 "
  bins[(4,3)]="Q2<1.5 and 0.02<x  and x<0.03 "
  bins[(4,4)]="Q2<1.5 and 0.03<x  and x<0.06 "

  Q2bins={}
  Q2bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  Q2bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  Q2bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  Q2bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  Q2bins[(4,0)]="Q2<1.5 and x<0.009 "
  Q2b=[]
  for k in Q2bins:
      d=data.query('%s and  had=="h+"'%Q2bins[k])
      Q2b.append(d.Q2.mean())
  Q2b=np.sort(np.unique(Q2b))

  xbins={}
  xbins[(4,0)]="x<0.009 "
  xbins[(4,1)]="0.009<x and x<0.012 "
  xbins[(4,2)]="0.012<x and x<0.02 "
  xbins[(4,3)]="0.02<x  and x<0.03 "
  xbins[(4,4)]="0.03<x  and x<0.06 "
  xbins[(2,3)]="0.02<x  and x<0.03 "
  xbins[(2,4)]="0.03<x  and x<0.06 "
  xbins[(2,5)]="0.06<x  and x<0.09 "
  xbins[(2,6)]="0.09<x  and x<0.2 "
  xbins[(1,7)]="0.2<x "
  
  xb=[]
  for k in xbins:
      d=data.query('%s and  had=="h+"'%xbins[k])
      xb.append(d.x.mean())
  xb=np.sort(np.unique(xb))


  zbins=[[0.24,0.3],[0.3,0.4],[0.4,0.5],[0.65,0.7]]

  nrows,ncols=2,3
  fig = py.figure(figsize=(ncols*3,nrows*3))

  cnt=0
  for k in [(2,6),(1,6),(0,7)]:
      cnt+=1
      ir,ic=k
      ax = py.subplot(nrows,ncols,cnt)
      d=data.query('%s and  had=="h+"'%bins[k])
      for i in range(len(zbins)):
          dd=d.query('%f<z and z<%f'%(zbins[i][0],zbins[i][1]))
          #if   i==0: c='r'
          #elif i==1: c='b'
          #elif i==2: c='g'
          #elif i==3: c='m'

          thy=(dd.DDSLO)/dd.idisNLO
          ax.errorbar(dd.qT**2,dd.value/thy,dd.alpha/thy,fmt='.')

          #ax.plot(dd.qT,(dd.DDSLO+dd.DDSNLO)/dd.idisNLO,'%s-'%c)
          #ax.plot(dd.qT,dd.LO/dd.idisNLO,'%s--'%c)
          #tot=dd.daleoLO+dd.daleoNLOdelta+dd.daleoNLOplus+dd.daleoNLOrest
          #ax.plot(dd.qT,tot/dd.idisNLO,'%s-'%c)
  
          #ax.plot(dd.qT,tot/dd.idisNLO,'%s:'%c,label='')
      #ax.axvline(d.Q.values[0])
      #ax.semilogy()
      ax.set_ylim(1e-4,10)
      ax.set_xlim(0,50)
      ax.set_xticklabels([])
  
      ax.tick_params(axis='both', which='major', labelsize=30, direction='in')
      #ax.set_xticks([2,4,6])
      ax.set_yticks([2,6,10,14,18])
      if cnt!=1:  ax.set_yticklabels([])


      ax.set_title(r'$x_{\rm bj}=%0.2f~~Q^2=%0.1f{~\rm GeV}^2$'%(d.x.mean(),d.Q2.mean()),size=13)

      ax.set_ylim(0,20);
      ax.xaxis.set_tick_params(labelsize=15)
      ax.yaxis.set_tick_params(labelsize=15)
      if cnt==1: ax.set_ylabel(r'$\rm data/theory(LO)$',size=20)
      #ax.text(0.05,0.90,r'$\rm \left<Q\right>=%0.1f{\rm GeV}$'%(0.5*(Qmin+Qmax))\
      #       ,transform=ax.transAxes, size=20)
      #if cnt==2: ax.legend(loc=2,frameon=0,fontsize=20,bbox_to_anchor=(-0.01,0.9))
 
      ax.plot([d.Q.values[0]**2,8**2],[1,1],c='k',lw=1,alpha=0.5)
      #if cnt==3: 
      #    #LO1b = mpatches.Rectangle((0,0), 0.05, 0.1, ec="none",color='k',alpha=0.1)
      #    #LO1l=py.Line2D([0,2], [0,0], color='k',ls=':')
      #    #LO2=py.Line2D([0], [0], color='k',ls='-')
      #    qTrange = mpatches.Rectangle((0,0), 0, 0, ec="none",color='y',alpha=0.5)
      #    #qTrange=py.Line2D([0], [0], color='y',ls='-',lw=10,alpha=0.5)
      #    ax.legend([qTrange],[r'$q_{\rm T}>Q$']\
      #              ,bbox_to_anchor=[-0.35, 0.3]\
      #              ,loc=4,fontsize=13,frameon=0)
      #if cnt==1: ax.text(0.8,0.8,r'$\rm a.1)$',transform=ax.transAxes,size=30)
      #if cnt==2: ax.text(0.8,0.8,r'$\rm b.1)$',transform=ax.transAxes,size=30)
      #if cnt==3: ax.text(0.8,0.8,r'$\rm c.1)$',transform=ax.transAxes,size=30)
      #if cnt==4: ax.text(0.8,0.8,r'$\rm d.1)$',transform=ax.transAxes,size=30)
      #    msg=r'${\rm COMPASS~17}~h^+$'
      #    ax.text(-2,2.8,msg,transform=ax.transAxes,size=50)
      #    msg =r'$\frac{d\sigma}{dx_{\rm bj}dQ^2dzdP^2_{\rm T}}/\frac{d\sigma}{dx_{\rm bj}dQ^2}'
      #    msg+=r'({\tiny {\rm GeV}^{-2}})'
      #    msg+=r'~{\rm vs.}~q_{\rm T}{~\rm (GeV)}$'
      #    ax.text(-2,2.2,msg,transform=ax.transAxes,size=50)
      if cnt==2: 
          ax.legend(loc=4,fontsize=11,frameon=0\
                   ,markerscale=1,handletextpad=0.1,bbox_to_anchor=(0.75,0.46))
      #ax.text(0.6,0.8,'(%d,%d)'%(ir,ic),transform=ax.transAxes,size=20)

  for k in [(2,6),(1,6),(0,7)]:
      cnt+=1
      ir,ic=k
      ax = py.subplot(nrows,ncols,cnt)
      d=data.query('%s and  had=="h+"'%bins[k])
      for i in range(len(zbins)):
          dd=d.query('%f<z and z<%f'%(zbins[i][0],zbins[i][1]))

          thy=(dd.DDSLO+dd.DDSNLO)/dd.idisNLO
          ax.errorbar(dd.qT**2,dd.value/thy,dd.alpha/thy,\
                      fmt='.',label=r'$%0.2f<z<%0.2f$'%(zbins[i][0],zbins[i][1]))
  
          #ax.plot(dd.qT,tot/dd.idisNLO,'%s:'%c,label='')
      #ax.axvline(d.Q.values[0])
      #ax.semilogy()
      ax.set_ylim(1e-4,10)
      ax.set_xlim(0,50)
      #ax.set_yticklabels([])
  
      ax.tick_params(axis='both', which='major', labelsize=30, direction='in')
      #ax.set_xticks([2,4,6])
      ax.set_yticks([2,6,10,14,18])
      if cnt!=4:  ax.set_yticklabels([])


      ax.set_ylim(0,20);
      ax.xaxis.set_tick_params(labelsize=15)
      ax.yaxis.set_tick_params(labelsize=15)
      if cnt==4: ax.set_ylabel(r'$\rm data/theory(NLO)$',size=20)
      ax.set_xlabel(r'$q_{\rm T}^2({\rm GeV}^2)$',size=20)
      #ax.text(0.05,0.90,r'$\rm \left<Q\right>=%0.1f{\rm GeV}$'%(0.5*(Qmin+Qmax))\
      #       ,transform=ax.transAxes, size=20)
      #if cnt==2: ax.legend(loc=2,frameon=0,fontsize=20,bbox_to_anchor=(-0.01,0.9))

 
      ax.plot([d.Q.values[0]**2,8**2],[1,1],c='k',lw=1,alpha=0.5)
      if cnt==6: 
          #LO1b = mpatches.Rectangle((0,0), 0.05, 0.1, ec="none",color='k',alpha=0.1)
          qTrange=py.Line2D([0,2], [0,0], color='k',ls='-',alpha=0.5)
          #LO2=py.Line2D([0], [0], color='k',ls='-')
          #qTrange = mpatches.Rectangle((0,0), 0, 0, ec="none",color='y',alpha=0.5)
          #qTrange=py.Line2D([0], [0], color='y',ls='-',lw=10,alpha=0.5)
          ax.legend([qTrange],[r'$q_{\rm T}>Q$']\
                    ,bbox_to_anchor=[-0.22, 0.3]\
                    ,loc=4,fontsize=13,frameon=0)
      if cnt==5: 
          ax.legend(loc=4,fontsize=13,frameon=0\
                   ,markerscale=1,handletextpad=0.1,bbox_to_anchor=(1.0,0.46))
      #ax.text(0.6,0.8,'(%d,%d)'%(ir,ic),transform=ax.transAxes,size=20)
  
  py.tight_layout()    
  py.savefig('gallery/compass17-rat.pdf')

def plot_hermes():

  data=pd.read_excel('expdata/hermes.xlsx')

  bins=[]
  bins.append("Q2<1.5 and x<0.04 ")
  bins.append("1.5<Q2 and Q2<1.8 and 0.04<x and x<0.065 ")
  bins.append("1.8<Q2 and Q2<2.8 and 0.065<x and x<0.1 ")
  bins.append("2.8<Q2 and Q2<5.0 and 0.1<x and x<0.2 ")
  bins.append("5.0<Q2 and Q2<9.0 and 0.2<x and x<0.3 ")
  bins.append("9.0<Q2 and 0.3<x ")

  zbins=[]
  zbins.append("z<0.20")
  zbins.append("0.20<z and z<0.27")
  #zbins.append("0.27<z and z<0.32")
  zbins.append("0.32<z and z<0.40")
  #zbins.append("0.40<z and z<0.50")
  zbins.append("0.50<z and z<0.60")
  #zbins.append("0.60<z and z<0.75")
  zbins.append("0.75<z")

  nrows,ncols=2,4
  fig = py.figure(figsize=(ncols*3,nrows*3))
  py.subplots_adjust(left=0.11, bottom=None, right=0.99, top=0.9, wspace=0.1, hspace=0.2)
  #gs = gridspec.GridSpec(5,8)
  #gs.update(wspace=0.,hspace=0,left=0.12, right=0.96,bottom=0.13,top=0.96)
  AX={}
  cnt=0
  for k in range(len(bins)):
      cnt+=1
      if cnt==4: cnt+=1
      ax = py.subplot(nrows,ncols,cnt)
      d=data.query('%s and  had=="pi+"'%bins[k])
      for i in range(len(zbins)):
          dd=d.query('%s'%(zbins[i]))
          if dd.index.size==0: continue
        
          #if   i==0: c='r'
          #elif i==1: c='b'
          #elif i==2: c='g'
          #elif i==3: c='m'

          e=ax.errorbar(dd.qT,dd.value,dd.alpha,fmt='.',label=r'$<z>=%0.1f$'%(dd.z.mean()))
          c=e[0].get_color()
          ax.fill_between(dd.qT,dd.LOmin/dd.idisNLO,dd.LOmax/dd.idisNLO,color=c,alpha=0.1)

          ax.plot(dd.qT,dd.DDSLO/dd.idisNLO,color=c,ls=':')

          tot=dd.DDSLO+dd.DDSNLOdelta+dd.DDSNLOplus+dd.DDSNLOrest
          #tot=tot.values
          ##for itot in range(tot.size):
          ##  if np.isnan(tot[itot]): 
          ##    tot[itot]=tot[itot+1]*5.5
          #tot=smooth(dd.qT.values,tot)
          ax.plot(dd.qT,tot/dd.idisNLO,color=c,ls='-',label='')
          #ax.plot(dd.qT,dd.LO/dd.idisNLO,'%s--'%c)
          ##ax.plot(dd.qT,tot/dd.idisNLO,'%s-'%c)
  
          #ax.plot(dd.qT,tot/dd.idisNLO,'%s:'%c,label='')
      #ax.axvline(d.Q.values[0])
      ax.semilogy(nonposy='clip')
      ax.set_ylim(4e-3,3e1)
      ax.set_xlim(0,7)

      if cnt!=1 and cnt!=5: ax.set_yticklabels([])
      if cnt<5: ax.set_xticklabels([])
      if cnt>4: 
        ax.set_xlabel(r'$q_{\rm T}(\rm {GeV})$',size=20)
        ax.set_xticks([0,2,4,6])

      if cnt==1: 
        msg =r'$\frac{d\sigma}{dx_{\rm bj}dQ^2dzdP^2_{\rm T}}/\frac{d\sigma}{dx_{\rm bj}dQ^2}$'
        ax.set_ylabel(msg,size=30)
        ax.yaxis.set_label_coords(-0.2, -0.02)
      if cnt==7:
        ax.text(1.3,0,r'${\rm HERMES}~\pi^+$',size=20,transform=ax.transAxes)

      if cnt==3: 
          ax.legend(bbox_to_anchor=[1.5, 0.1], loc='center',fontsize=15,frameon=0\
                   ,markerscale=1,handletextpad=0.1)

      ax.set_title(r'$x_{\rm bj}=%0.2f~~Q^2=%0.1f~({\rm GeV}^2)$'%(dd.x.mean(),dd.Q2.mean()),size=13)

  
      ax.tick_params(axis='both', which='both', labelsize=15, direction='in')
      #if all(k!=_ for _ in [(4,0),(4,1),(4,2),(4,3),(4,4),(3,5),(2,6),(1,7)]): 
      #    ax.set_xticklabels([])
      #else:
      #    ax.set_xticks([2,4,6])
      #    #ax.set_xlabel(r'$q_{\rm T}$',size=35)
      #    #ax.xaxis.set_label_coords(0.95, -0.02)
      #if all(k!=_ for _ in [(4,0),(3,0),(2,2),(1,4),(0,6)]): 
      #    ax.set_yticklabels([])
      #    ax.set_yticks([1e-3,1e-2,1e-1,1])
      #else:
      #    ax.set_yticks([1e-3,1e-2,1e-1,1])
      
      #if k==(4,0):
      #    ax.annotate('', xy=(-0.35, 5.2), xycoords='axes fraction', xytext=(-0.35, -0.1),
      #                arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))
      #    ax.annotate('', xy=(8.2,-0.3), xycoords='axes fraction', xytext=(-0.1, -0.3),
      #                arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))        
      #    ax.annotate(r'$Q^2~({\rm GeV}^2)$', xy=(-1,3),xycoords='axes fraction',size=50,rotation=90)
      #    ax.annotate(r'$x_{\rm bj}$', xy=(3.9,-0.7),xycoords='axes fraction',size=50)
      #            
      #    for i in range(8):
      #        if xb[i]<2e-2: msg=r'$%0.3f$'%xb[i]
      #        else:msg=r'$%0.2f$'%xb[i]
      #        ax.text(0.5+i,-0.5,msg,transform=ax.transAxes,size=30,ha="center")
      #        ax.annotate('',xy=(i,-0.35),xycoords='axes fraction',xytext=(i+1, -0.35), \
      #                    arrowprops=dict(arrowstyle="<->", color='k'))
  
      #    for i in range(5):
      #        ax.text(-0.6,0.5+i,r'$%0.1f$'%Q2b[i],transform=ax.transAxes,size=30,rotation=90,va="center")
      #        ax.annotate('',xy=(-0.4,i),xycoords='axes fraction',xytext=(-0.4,i+1), \
      #                    arrowprops=dict(arrowstyle="<->", color='k'))
  
      ax.plot([d.Q.values[0],7],[7e-3,7e-3],c='y',lw=10,alpha=0.5)

      if cnt==7: 
          LO1b = mpatches.Rectangle((0,0), 0.05, 0.1, ec="none",color='k',alpha=0.1)
          LO1l=py.Line2D([0,2], [0,0], color='k',ls=':')
          LO2=py.Line2D([0], [0], color='k',ls='-')
          qTrange = mpatches.Rectangle((0,0), 0, 0, ec="none",color='y',alpha=0.5)
      #    #qTrange=py.Line2D([0], [0], color='y',ls='-',lw=10,alpha=0.5)
          ax.legend([(LO1b,LO1l),LO2,qTrange],[r'$\rm DDS~(LO)$',r'$\rm DDS~(NLO)$',r'$q_{\rm T}>Q$']\
                    ,bbox_to_anchor=[1.6, 0.5]\
                    ,loc='center',fontsize=15,frameon=0)
      #    msg=r'${\rm COMPASS~17}~h^+$'
      #    ax.text(-2,2.8,msg,transform=ax.transAxes,size=50)
      #    msg =r'$\frac{d\sigma}{dx_{\rm bj}dQ^2dzdP^2_{\rm T}}/\frac{d\sigma}{dx_{\rm bj}dQ^2}'
      #    msg+=r'({\tiny {\rm GeV}^{-2}})'
      #    msg+=r'~{\rm vs.}~q_{\rm T}{~\rm (GeV)}$'
      #    ax.text(-2,2.2,msg,transform=ax.transAxes,size=50)
      #ax.text(0.6,0.8,'(%d,%d)'%(ir,ic),transform=ax.transAxes,size=20)
  
  #py.tight_layout()        
  py.savefig('gallery/hermes.pdf')

def plot_hermes_rat():

  data=pd.read_excel('expdata/hermes.xlsx')

  bins=[]
  #bins.append("Q2<1.5 and x<0.04 ")
  #bins.append("1.5<Q2 and Q2<1.8 and 0.04<x and x<0.065 ")
  #bins.append("1.8<Q2 and Q2<2.8 and 0.065<x and x<0.1 ")
  bins.append("2.8<Q2 and Q2<5.0 and 0.1<x and x<0.2 ")
  bins.append("5.0<Q2 and Q2<9.0 and 0.2<x and x<0.3 ")
  bins.append("9.0<Q2 and 0.3<x ")

  zbins=[]
  zbins.append("z<0.20")
  zbins.append("0.20<z and z<0.27")
  #zbins.append("0.27<z and z<0.32")
  zbins.append("0.32<z and z<0.40")
  #zbins.append("0.40<z and z<0.50")
  zbins.append("0.50<z and z<0.60")
  #zbins.append("0.60<z and z<0.75")
  zbins.append("0.75<z")

  nrows,ncols=2,3
  fig = py.figure(figsize=(ncols*3,nrows*3))
  #py.subplots_adjust(left=0.11, bottom=None, right=0.99, top=0.9, wspace=0.1, hspace=0.2)

  cnt=0
  for k in range(len(bins)):
      cnt+=1
      if cnt==4: cnt+=1
      ax = py.subplot(nrows,ncols,cnt)
      d=data.query('%s and  had=="pi+"'%bins[k])
      for i in range(len(zbins)):
          dd=d.query('%s'%(zbins[i]))
          if dd.index.size==0: continue
          tot=dd.DDSLO/dd.idisNLO#+dd.DDSNLOdelta+dd.DDSNLOplus+dd.DDSNLOrest
          ax.errorbar(dd.qT**2,dd.value/tot,dd.alpha/tot,fmt='.',label=r'$<z>=%0.1f$'%(dd.z.mean()))
      ax.set_ylim(0,40)
      ax.set_xlim(0,7**2)
      if cnt==1: ax.set_ylabel(r'$\rm data/theory(LO)$',size=20)
      ax.set_yticks([10,20,30])
      ax.set_xticklabels([])
      if cnt>1: ax.set_yticklabels([])
      ax.set_title(r'$x_{\rm bj}=%0.2f~~Q^2=%0.1f~{\rm GeV}^2$'%(dd.x.mean(),dd.Q2.mean()),size=13)

      #if cnt!=1 and cnt!=5: ax.set_yticklabels([])
      #if cnt<5: ax.set_xticklabels([])
      #if cnt>4: 
      #  ax.set_xlabel(r'$q_{\rm T}(\rm {GeV})$',size=20)
      #  ax.set_xticks([0,2,4,6])

      #if cnt==1: 
      #  msg =r'$\frac{d\sigma}{dx_{\rm bj}dQ^2dzdP^2_{\rm T}}/\frac{d\sigma}{dx_{\rm bj}dQ^2}$'
      #  ax.set_ylabel(msg,size=30)
      #  ax.yaxis.set_label_coords(-0.2, -0.02)
      #if cnt==7:
      #  ax.text(1.3,0,r'${\rm HERMES}~\pi^+$',size=20,transform=ax.transAxes)

      #if cnt==3: 
      #    ax.legend(bbox_to_anchor=[1.5, 0.1], loc='center',fontsize=15,frameon=0\
      #             ,markerscale=1,handletextpad=0.1)

      #ax.set_title(r'$x_{\rm bj}=%0.2f~~Q^2=%0.1f~({\rm GeV}^2)$'%(dd.x.mean(),dd.Q2.mean()),size=13)

  
      ax.tick_params(axis='both', which='both', labelsize=15, direction='in')
      ax.plot([d.Q.values[0]**2,7**2],[1,1],c='k',lw=1,alpha=0.5)

      #if cnt==7: 
      #    LO1b = mpatches.Rectangle((0,0), 0.05, 0.1, ec="none",color='k',alpha=0.1)
      #    LO1l=py.Line2D([0,2], [0,0], color='k',ls=':')
      #    LO2=py.Line2D([0], [0], color='k',ls='-')
      #    qTrange = mpatches.Rectangle((0,0), 0, 0, ec="none",color='y',alpha=0.5)
      #    ax.legend([(LO1b,LO1l),LO2,qTrange],\
      #              [r'$\rm DDS~(LO)$',r'$\rm DDS~(NLO)$',r'$q_{\rm T}>Q$']\
      #              ,bbox_to_anchor=[1.6, 0.5]\
      #              ,loc='center',fontsize=15,frameon=0)


  for k in range(len(bins)):
      cnt+=1
      ax = py.subplot(nrows,ncols,cnt)
      d=data.query('%s and  had=="pi+"'%bins[k])
      for i in range(len(zbins)):
          dd=d.query('%s'%(zbins[i]))
          if dd.index.size==0: continue
          tot=(dd.DDSLO+dd.DDSNLOdelta+dd.DDSNLOplus+dd.DDSNLOrest)/dd.idisNLO
          ax.errorbar(dd.qT**2,dd.value/tot,dd.alpha/tot,fmt='.',label=r'$<z>=%0.1f$'%(dd.z.mean()))
      if cnt==4: ax.set_ylabel(r'$\rm data/theory(NLO)$',size=20)
      ax.set_ylim(0,40)
      ax.set_xlim(0,7**2)
      ax.set_yticks([10,20,30])
      if cnt>4: ax.set_yticklabels([])

      #if cnt!=1 and cnt!=5: ax.set_yticklabels([])
      #if cnt<5: ax.set_xticklabels([])
      #if cnt>4: 
      #  ax.set_xlabel(r'$q_{\rm T}(\rm {GeV})$',size=20)
      #  ax.set_xticks([0,2,4,6])

      #if cnt==1: 
      #  msg =r'$\frac{d\sigma}{dx_{\rm bj}dQ^2dzdP^2_{\rm T}}/\frac{d\sigma}{dx_{\rm bj}dQ^2}$'
      #  ax.set_ylabel(msg,size=30)
      #  ax.yaxis.set_label_coords(-0.2, -0.02)
      #if cnt==7:
      #  ax.text(1.3,0,r'${\rm HERMES}~\pi^+$',size=20,transform=ax.transAxes)

      if cnt==4: 
          ax.legend(loc=2,fontsize=15,frameon=0\
                   ,markerscale=1,handletextpad=0.1)

  
      ax.tick_params(axis='both', which='both', labelsize=15, direction='in')
      ax.plot([d.Q.values[0]**2,7**2],[1,1],c='k',lw=1,alpha=0.5)
      ax.set_xlabel(r'$q_{\rm T}^2({\rm GeV}^2)$',size=20)

      if cnt==5: 
          qTrange=py.Line2D([0,2], [0,0], color='k',ls='-',alpha=0.5)
          ax.legend([qTrange],[r'$q_{\rm T}>Q$']\
                    ,bbox_to_anchor=[0.5, 0.8]\
                    ,loc='center',fontsize=15,frameon=0)

  py.tight_layout()
  py.savefig('gallery/hermes-rat.pdf')

def plot_compass17_nlo():

  data=pd.read_excel('expdata/compass17.xlsx')


  bins={}
  bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  bins[(0,7)]="17<Q2 and 0.2<x "
  bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  bins[(1,5)]="6<Q2 and Q2<17 and 0.06<x  and x<0.09 "
  bins[(1,6)]="6<Q2 and Q2<17 and 0.09<x  and x<0.2 "
  bins[(1,7)]="6<Q2 and Q2<17 and 0.2<x "
  bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  bins[(2,3)]="3<Q2 and Q2<6 and 0.02<x  and x<0.03 "
  bins[(2,4)]="3<Q2 and Q2<6 and 0.03<x  and x<0.06 "
  bins[(2,5)]="3<Q2 and Q2<6 and 0.06<x  and x<0.09 "
  bins[(2,6)]="3<Q2 and Q2<6 and 0.09<x  and x<0.2 "
  bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  bins[(3,1)]="1.5<Q2 and Q2<3.0 and 0.009<x and x<0.012 "
  bins[(3,2)]="1.5<Q2 and Q2<3.0 and 0.012<x and x<0.02 "
  bins[(3,3)]="1.5<Q2 and Q2<3.0 and 0.02<x  and x<0.03 "
  bins[(3,4)]="1.5<Q2 and Q2<3.0 and 0.03<x  and x<0.06 "
  bins[(3,5)]="1.5<Q2 and Q2<3.0 and 0.06<x  and x<0.09 "
  bins[(4,0)]="Q2<1.5 and x<0.009 "
  bins[(4,1)]="Q2<1.5 and 0.009<x and x<0.012 "
  bins[(4,2)]="Q2<1.5 and 0.012<x and x<0.02 "
  bins[(4,3)]="Q2<1.5 and 0.02<x  and x<0.03 "
  bins[(4,4)]="Q2<1.5 and 0.03<x  and x<0.06 "

  Q2bins={}
  Q2bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  Q2bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  Q2bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  Q2bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  Q2bins[(4,0)]="Q2<1.5 and x<0.009 "
  Q2b=[]
  for k in Q2bins:
      d=data.query('%s and  had=="h+"'%Q2bins[k])
      Q2b.append(d.Q2.mean())
  Q2b=np.sort(np.unique(Q2b))

  xbins={}
  xbins[(4,0)]="x<0.009 "
  xbins[(4,1)]="0.009<x and x<0.012 "
  xbins[(4,2)]="0.012<x and x<0.02 "
  xbins[(4,3)]="0.02<x  and x<0.03 "
  xbins[(4,4)]="0.03<x  and x<0.06 "
  xbins[(2,3)]="0.02<x  and x<0.03 "
  xbins[(2,4)]="0.03<x  and x<0.06 "
  xbins[(2,5)]="0.06<x  and x<0.09 "
  xbins[(2,6)]="0.09<x  and x<0.2 "
  xbins[(1,7)]="0.2<x "
  
  xb=[]
  for k in xbins:
      d=data.query('%s and  had=="h+"'%xbins[k])
      xb.append(d.x.mean())
  xb=np.sort(np.unique(xb))


  zbins=[[0.24,0.3],[0.3,0.4],[0.4,0.5],[0.65,0.7]]

  nrows,ncols=5,8
  fig = py.figure(figsize=(ncols*3,nrows*3))
  gs = gridspec.GridSpec(5,8)
  gs.update(wspace=0.,hspace=0,left=0.12, right=0.96,bottom=0.13,top=0.96)
  AX={}
  for k in bins:
      ir,ic=k
      ax = py.subplot(gs[ir,ic])
      d=data.query('%s and  had=="h+"'%bins[k])
      for i in range(len(zbins)):
          dd=d.query('%f<z and z<%f'%(zbins[i][0],zbins[i][1]))
          if dd.index.size==0: continue
        
          #if   i==0: c='r'
          #elif i==1: c='b'
          #elif i==2: c='g'
          #elif i==3: c='m'

          #e=ax.errorbar(dd.qT,dd.value,dd.alpha,fmt='.',label=r'$%0.2f<z<%0.2f$'%(zbins[i][0],zbins[i][1]))
          #c=e[0].get_color()

          #ax.plot(dd.qT,dd.LO/dd.idisNLO,color='k',ls='--')
          #ax.plot(dd.qT,dd.DDSLO/dd.idisNLO,color=c,ls=':')

 
          tot=dd.LO+dd.NLOdelta+dd.NLOplus+dd.NLOregular
          tot=tot.values
          #tot=smooth(dd.qT.values,tot)
          ax.plot(dd.qT,tot/dd.idisNLO,ls='-',label='')


          tot=dd.DDSLO+dd.DDSNLOdelta+dd.DDSNLOplus+dd.DDSNLOrest
          tot=tot.values
          #tot=smooth(dd.qT.values,tot)
          #ax.plot(dd.qT,tot/dd.idisNLO,color=c,ls='-',label='')
          ax.plot(dd.qT,tot/dd.idisNLO,color='k',ls='--',label='')
 
          #ax.plot(dd.qT,tot/dd.idisNLO,'%s:'%c,label='')
      #ax.axvline(d.Q.values[0])
      ax.semilogy()
      ax.set_ylim(1e-4,10)
      ax.set_xlim(0,8)
      #ax.set_yticklabels([])
  
      ax.tick_params(axis='both', which='major', labelsize=30, direction='in')
      if all(k!=_ for _ in [(4,0),(4,1),(4,2),(4,3),(4,4),(3,5),(2,6),(1,7)]): 
          ax.set_xticklabels([])
      else:
          ax.set_xticks([2,4,6])
          #ax.set_xlabel(r'$q_{\rm T}$',size=35)
          #ax.xaxis.set_label_coords(0.95, -0.02)
      if all(k!=_ for _ in [(4,0),(3,0),(2,2),(1,4),(0,6)]): 
          ax.set_yticklabels([])
          ax.set_yticks([1e-3,1e-2,1e-1,1])
      else:
          ax.set_yticks([1e-3,1e-2,1e-1,1])
      
      if k==(4,0):
          ax.annotate('', xy=(-0.35, 5.2), xycoords='axes fraction', xytext=(-0.35, -0.1),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))
          ax.annotate('', xy=(8.2,-0.3), xycoords='axes fraction', xytext=(-0.1, -0.3),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))        
          ax.annotate(r'$Q^2~({\rm GeV}^2)$', xy=(-1,3),xycoords='axes fraction',size=50,rotation=90)
          ax.annotate(r'$x_{\rm bj}$', xy=(3.9,-0.7),xycoords='axes fraction',size=50)
                  
          for i in range(8):
              if xb[i]<2e-2: msg=r'$%0.3f$'%xb[i]
              else:msg=r'$%0.2f$'%xb[i]
              ax.text(0.5+i,-0.5,msg,transform=ax.transAxes,size=30,ha="center")
              ax.annotate('',xy=(i,-0.35),xycoords='axes fraction',xytext=(i+1, -0.35), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
          for i in range(5):
              ax.text(-0.6,0.5+i,r'$%0.1f$'%Q2b[i],transform=ax.transAxes,size=30,rotation=90,va="center")
              ax.annotate('',xy=(-0.4,i),xycoords='axes fraction',xytext=(-0.4,i+1), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
      ax.plot([d.Q.values[0],8],[2e-4,2e-4],c='y',lw=10,alpha=0.5)
      if k==(2,2): 
          #LO1b = mpatches.Rectangle((0,0), 0.05, 0.1, ec="none",color='k',alpha=0.1)
          LO1=py.Line2D([0,2], [0,0], color='k',ls='--')
          LO2=py.Line2D([0], [0], color='k',ls='-')
          qTrange = mpatches.Rectangle((0,0), 0, 0, ec="none",color='y',alpha=0.5)
          #qTrange=py.Line2D([0], [0], color='y',ls='-',lw=10,alpha=0.5)
          ax.legend([LO1,LO2,qTrange],[r'$\rm DDS~(NLO)$',r'$\rm ours~(NLO)$',r'$q_{\rm T}>Q$']\
                    ,bbox_to_anchor=[-1.2, 1.]\
                    ,loc='center',fontsize=40,frameon=0)
          msg=r'${\rm COMPASS~17}~h^+$'
          ax.text(-2,2.8,msg,transform=ax.transAxes,size=50)
          msg =r'$\frac{d\sigma}{dx_{\rm bj}dQ^2dzdP^2_{\rm T}}/\frac{d\sigma}{dx_{\rm bj}dQ^2}'
          msg+=r'({\tiny {\rm GeV}^{-2}})'
          msg+=r'~{\rm vs.}~q_{\rm T}{~\rm (GeV)}$'
          ax.text(-2,2.2,msg,transform=ax.transAxes,size=50)
      if k==(1,4): 
          ax.legend(bbox_to_anchor=[3, -2.5], loc='center',fontsize=40,frameon=0\
                   ,markerscale=2,handletextpad=0.1)
      #ax.text(0.6,0.8,'(%d,%d)'%(ir,ic),transform=ax.transAxes,size=20)
  
          
  py.savefig('gallery/compass17-nlo.pdf')

def plot_compass17_lo_diff():

  data=pd.read_excel('expdata/compass17B.xlsx')


  bins={}
  bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  bins[(0,7)]="17<Q2 and 0.2<x "
  bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  bins[(1,5)]="6<Q2 and Q2<17 and 0.06<x  and x<0.09 "
  bins[(1,6)]="6<Q2 and Q2<17 and 0.09<x  and x<0.2 "
  bins[(1,7)]="6<Q2 and Q2<17 and 0.2<x "
  bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  bins[(2,3)]="3<Q2 and Q2<6 and 0.02<x  and x<0.03 "
  bins[(2,4)]="3<Q2 and Q2<6 and 0.03<x  and x<0.06 "
  bins[(2,5)]="3<Q2 and Q2<6 and 0.06<x  and x<0.09 "
  bins[(2,6)]="3<Q2 and Q2<6 and 0.09<x  and x<0.2 "
  bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  bins[(3,1)]="1.5<Q2 and Q2<3.0 and 0.009<x and x<0.012 "
  bins[(3,2)]="1.5<Q2 and Q2<3.0 and 0.012<x and x<0.02 "
  bins[(3,3)]="1.5<Q2 and Q2<3.0 and 0.02<x  and x<0.03 "
  bins[(3,4)]="1.5<Q2 and Q2<3.0 and 0.03<x  and x<0.06 "
  bins[(3,5)]="1.5<Q2 and Q2<3.0 and 0.06<x  and x<0.09 "
  bins[(4,0)]="Q2<1.5 and x<0.009 "
  bins[(4,1)]="Q2<1.5 and 0.009<x and x<0.012 "
  bins[(4,2)]="Q2<1.5 and 0.012<x and x<0.02 "
  bins[(4,3)]="Q2<1.5 and 0.02<x  and x<0.03 "
  bins[(4,4)]="Q2<1.5 and 0.03<x  and x<0.06 "

  Q2bins={}
  Q2bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  Q2bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  Q2bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  Q2bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  Q2bins[(4,0)]="Q2<1.5 and x<0.009 "
  Q2b=[]
  for k in Q2bins:
      d=data.query('%s and  had=="h+"'%Q2bins[k])
      Q2b.append(d.Q2.mean())
  Q2b=np.sort(np.unique(Q2b))

  xbins={}
  xbins[(4,0)]="x<0.009 "
  xbins[(4,1)]="0.009<x and x<0.012 "
  xbins[(4,2)]="0.012<x and x<0.02 "
  xbins[(4,3)]="0.02<x  and x<0.03 "
  xbins[(4,4)]="0.03<x  and x<0.06 "
  xbins[(2,3)]="0.02<x  and x<0.03 "
  xbins[(2,4)]="0.03<x  and x<0.06 "
  xbins[(2,5)]="0.06<x  and x<0.09 "
  xbins[(2,6)]="0.09<x  and x<0.2 "
  xbins[(1,7)]="0.2<x "
  
  xb=[]
  for k in xbins:
      d=data.query('%s and  had=="h+"'%xbins[k])
      xb.append(d.x.mean())
  xb=np.sort(np.unique(xb))


  zbins=[[0.24,0.3],[0.3,0.4],[0.4,0.5],[0.65,0.7]]

  nrows,ncols=5,8
  fig = py.figure(figsize=(ncols*3,nrows*3))
  gs = gridspec.GridSpec(5,8)
  gs.update(wspace=0.,hspace=0,left=0.12, right=0.96,bottom=0.13,top=0.96)
  AX={}
  for k in bins:
      ir,ic=k
      ax = py.subplot(gs[ir,ic])
      d=data.query('%s and  had=="h+"'%bins[k])
      for i in range(len(zbins)):
          dd=d.query('%f<z and z<%f'%(zbins[i][0],zbins[i][1]))
          if dd.index.size==0: continue
        

 
          tot1=dd.LO#+dd.NLOdelta+dd.NLOplus+dd.NLOregular
          tot1=tot1.values

          tot2=dd.DDSLO#+dd.DDSNLOdelta+dd.DDSNLOplus+dd.DDSNLOrest
          tot2=tot2.values
          ax.plot(np.abs((tot1-tot2)/tot1*100),label=r'$%0.2f<z<%0.2f$'%(zbins[i][0],zbins[i][1]))
 
          #ax.plot(dd.qT,tot/dd.idisNLO,'%s:'%c,label='')
      #ax.axvline(d.Q.values[0])
      #ax.semilogy()
      #ax.set_ylim(1e-4,10)
      ax.set_xlim(0,8)
      ax.set_ylim(0,50)
      #ax.set_yticklabels([])
  
      ax.tick_params(axis='both', which='major', labelsize=30, direction='in')
      if all(k!=_ for _ in [(4,0),(4,1),(4,2),(4,3),(4,4),(3,5),(2,6),(1,7)]): 
          ax.set_xticklabels([])
      else:
          ax.set_xticks([2,4,6])
          #ax.set_xlabel(r'$q_{\rm T}$',size=35)
          #ax.xaxis.set_label_coords(0.95, -0.02)
      if all(k!=_ for _ in [(4,0),(3,0),(2,2),(1,4),(0,6)]): 
          ax.set_yticklabels([])
          ax.set_yticks([10,20,30,40])
      else:
          ax.set_yticks([10,20,30,40])
      if k==(4,0):
          ax.annotate('', xy=(-0.35, 5.2), xycoords='axes fraction', xytext=(-0.35, -0.1),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))
          ax.annotate('', xy=(8.2,-0.3), xycoords='axes fraction', xytext=(-0.1, -0.3),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))        
          ax.annotate(r'$Q^2~({\rm GeV}^2)$', xy=(-1,3),xycoords='axes fraction',size=50,rotation=90)
          ax.annotate(r'$x_{\rm bj}$', xy=(3.9,-0.7),xycoords='axes fraction',size=50)
                  
          for i in range(8):
              if xb[i]<2e-2: msg=r'$%0.3f$'%xb[i]
              else:msg=r'$%0.2f$'%xb[i]
              ax.text(0.5+i,-0.5,msg,transform=ax.transAxes,size=30,ha="center")
              ax.annotate('',xy=(i,-0.35),xycoords='axes fraction',xytext=(i+1, -0.35), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
          for i in range(5):
              ax.text(-0.6,0.5+i,r'$%0.1f$'%Q2b[i],transform=ax.transAxes,\
                       size=30,rotation=90,va="center")
              ax.annotate('',xy=(-0.4,i),xycoords='axes fraction',xytext=(-0.4,i+1), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
      if k==(2,2): 
          msg=r'${\rm COMPASS~17}~h^+$'
          ax.text(-2,2.8,msg,transform=ax.transAxes,size=50)
          msg =r'$|{\rm (ours-DDS)/ours}|\times 100~{\rm vs.}~q_{\rm T}{~\rm (GeV)}$'
          ax.text(-2,2.2,msg,transform=ax.transAxes,size=50)
          msg =r'$\rm at~LO$'
          ax.text(-2,1.8,msg,transform=ax.transAxes,size=50)
      if k==(1,4): 
          ax.legend(bbox_to_anchor=[3, -2.5], loc='center',fontsize=40,frameon=0\
                   ,markerscale=2,handletextpad=0.1)
  
          
  py.savefig('gallery/compass17-lo-diff.pdf')

def plot_compass17_nlo_diff():

  data=pd.read_excel('expdata/compass17B.xlsx')


  bins={}
  bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  bins[(0,7)]="17<Q2 and 0.2<x "
  bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  bins[(1,5)]="6<Q2 and Q2<17 and 0.06<x  and x<0.09 "
  bins[(1,6)]="6<Q2 and Q2<17 and 0.09<x  and x<0.2 "
  bins[(1,7)]="6<Q2 and Q2<17 and 0.2<x "
  bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  bins[(2,3)]="3<Q2 and Q2<6 and 0.02<x  and x<0.03 "
  bins[(2,4)]="3<Q2 and Q2<6 and 0.03<x  and x<0.06 "
  bins[(2,5)]="3<Q2 and Q2<6 and 0.06<x  and x<0.09 "
  bins[(2,6)]="3<Q2 and Q2<6 and 0.09<x  and x<0.2 "
  bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  bins[(3,1)]="1.5<Q2 and Q2<3.0 and 0.009<x and x<0.012 "
  bins[(3,2)]="1.5<Q2 and Q2<3.0 and 0.012<x and x<0.02 "
  bins[(3,3)]="1.5<Q2 and Q2<3.0 and 0.02<x  and x<0.03 "
  bins[(3,4)]="1.5<Q2 and Q2<3.0 and 0.03<x  and x<0.06 "
  bins[(3,5)]="1.5<Q2 and Q2<3.0 and 0.06<x  and x<0.09 "
  bins[(4,0)]="Q2<1.5 and x<0.009 "
  bins[(4,1)]="Q2<1.5 and 0.009<x and x<0.012 "
  bins[(4,2)]="Q2<1.5 and 0.012<x and x<0.02 "
  bins[(4,3)]="Q2<1.5 and 0.02<x  and x<0.03 "
  bins[(4,4)]="Q2<1.5 and 0.03<x  and x<0.06 "

  Q2bins={}
  Q2bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  Q2bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  Q2bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  Q2bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  Q2bins[(4,0)]="Q2<1.5 and x<0.009 "
  Q2b=[]
  for k in Q2bins:
      d=data.query('%s and  had=="h+"'%Q2bins[k])
      Q2b.append(d.Q2.mean())
  Q2b=np.sort(np.unique(Q2b))

  xbins={}
  xbins[(4,0)]="x<0.009 "
  xbins[(4,1)]="0.009<x and x<0.012 "
  xbins[(4,2)]="0.012<x and x<0.02 "
  xbins[(4,3)]="0.02<x  and x<0.03 "
  xbins[(4,4)]="0.03<x  and x<0.06 "
  xbins[(2,3)]="0.02<x  and x<0.03 "
  xbins[(2,4)]="0.03<x  and x<0.06 "
  xbins[(2,5)]="0.06<x  and x<0.09 "
  xbins[(2,6)]="0.09<x  and x<0.2 "
  xbins[(1,7)]="0.2<x "
  
  xb=[]
  for k in xbins:
      d=data.query('%s and  had=="h+"'%xbins[k])
      xb.append(d.x.mean())
  xb=np.sort(np.unique(xb))


  zbins=[[0.24,0.3],[0.3,0.4],[0.4,0.5],[0.65,0.7]]

  nrows,ncols=5,8
  fig = py.figure(figsize=(ncols*3,nrows*3))
  gs = gridspec.GridSpec(5,8)
  gs.update(wspace=0.,hspace=0,left=0.12, right=0.96,bottom=0.13,top=0.96)
  AX={}
  for k in bins:
      ir,ic=k
      ax = py.subplot(gs[ir,ic])
      d=data.query('%s and  had=="h+"'%bins[k])
      for i in range(len(zbins)):
          dd=d.query('%f<z and z<%f'%(zbins[i][0],zbins[i][1]))
          if dd.index.size==0: continue
        

 
          tot1=dd.LO+dd.NLOdelta+dd.NLOplus+dd.NLOregular
          tot1=tot1.values

          tot2=dd.DDSLO+dd.DDSNLOdelta+dd.DDSNLOplus+dd.DDSNLOrest
          tot2=tot2.values
          ax.plot(np.abs((tot1-tot2)/tot1*100),label=r'$%0.2f<z<%0.2f$'%(zbins[i][0],zbins[i][1]))
 
          #ax.plot(dd.qT,tot/dd.idisNLO,'%s:'%c,label='')
      #ax.axvline(d.Q.values[0])
      #ax.semilogy()
      #ax.set_ylim(1e-4,10)
      ax.set_xlim(0,8)
      ax.set_ylim(0,50)
      #ax.set_yticklabels([])
  
      ax.tick_params(axis='both', which='major', labelsize=30, direction='in')
      if all(k!=_ for _ in [(4,0),(4,1),(4,2),(4,3),(4,4),(3,5),(2,6),(1,7)]): 
          ax.set_xticklabels([])
      else:
          ax.set_xticks([2,4,6])
          #ax.set_xlabel(r'$q_{\rm T}$',size=35)
          #ax.xaxis.set_label_coords(0.95, -0.02)
      if all(k!=_ for _ in [(4,0),(3,0),(2,2),(1,4),(0,6)]): 
          ax.set_yticklabels([])
          ax.set_yticks([10,20,30,40])
      else:
          ax.set_yticks([10,20,30,40])
      if k==(4,0):
          ax.annotate('', xy=(-0.35, 5.2), xycoords='axes fraction', xytext=(-0.35, -0.1),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))
          ax.annotate('', xy=(8.2,-0.3), xycoords='axes fraction', xytext=(-0.1, -0.3),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))        
          ax.annotate(r'$Q^2~({\rm GeV}^2)$', xy=(-1,3),xycoords='axes fraction',size=50,rotation=90)
          ax.annotate(r'$x_{\rm bj}$', xy=(3.9,-0.7),xycoords='axes fraction',size=50)
                  
          for i in range(8):
              if xb[i]<2e-2: msg=r'$%0.3f$'%xb[i]
              else:msg=r'$%0.2f$'%xb[i]
              ax.text(0.5+i,-0.5,msg,transform=ax.transAxes,size=30,ha="center")
              ax.annotate('',xy=(i,-0.35),xycoords='axes fraction',xytext=(i+1, -0.35), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
          for i in range(5):
              ax.text(-0.6,0.5+i,r'$%0.1f$'%Q2b[i],transform=ax.transAxes,\
                       size=30,rotation=90,va="center")
              ax.annotate('',xy=(-0.4,i),xycoords='axes fraction',xytext=(-0.4,i+1), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
      if k==(2,2): 
          msg=r'${\rm COMPASS~17}~h^+$'
          ax.text(-2,2.8,msg,transform=ax.transAxes,size=50)
          msg =r'$|{\rm (ours-DDS)/ours}|\times 100~{\rm vs.}~q_{\rm T}{~\rm (GeV)}$'
          ax.text(-2,2.2,msg,transform=ax.transAxes,size=50)
          msg =r'$\rm at~NLO$'
          ax.text(-2,1.8,msg,transform=ax.transAxes,size=50)
      if k==(1,4): 
          ax.legend(bbox_to_anchor=[3, -2.5], loc='center',fontsize=40,frameon=0\
                   ,markerscale=2,handletextpad=0.1)
  
          
  py.savefig('gallery/compass17-nlo-diff.pdf')

def plot_compass17_relative_ours():

  data=pd.read_excel('expdata/compass17-channels.xlsx')


  bins={}
  bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  bins[(0,7)]="17<Q2 and 0.2<x "
  bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  bins[(1,5)]="6<Q2 and Q2<17 and 0.06<x  and x<0.09 "
  bins[(1,6)]="6<Q2 and Q2<17 and 0.09<x  and x<0.2 "
  bins[(1,7)]="6<Q2 and Q2<17 and 0.2<x "
  bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  bins[(2,3)]="3<Q2 and Q2<6 and 0.02<x  and x<0.03 "
  bins[(2,4)]="3<Q2 and Q2<6 and 0.03<x  and x<0.06 "
  bins[(2,5)]="3<Q2 and Q2<6 and 0.06<x  and x<0.09 "
  bins[(2,6)]="3<Q2 and Q2<6 and 0.09<x  and x<0.2 "
  bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  bins[(3,1)]="1.5<Q2 and Q2<3.0 and 0.009<x and x<0.012 "
  bins[(3,2)]="1.5<Q2 and Q2<3.0 and 0.012<x and x<0.02 "
  bins[(3,3)]="1.5<Q2 and Q2<3.0 and 0.02<x  and x<0.03 "
  bins[(3,4)]="1.5<Q2 and Q2<3.0 and 0.03<x  and x<0.06 "
  bins[(3,5)]="1.5<Q2 and Q2<3.0 and 0.06<x  and x<0.09 "
  bins[(4,0)]="Q2<1.5 and x<0.009 "
  bins[(4,1)]="Q2<1.5 and 0.009<x and x<0.012 "
  bins[(4,2)]="Q2<1.5 and 0.012<x and x<0.02 "
  bins[(4,3)]="Q2<1.5 and 0.02<x  and x<0.03 "
  bins[(4,4)]="Q2<1.5 and 0.03<x  and x<0.06 "

  Q2bins={}
  Q2bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  Q2bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  Q2bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  Q2bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  Q2bins[(4,0)]="Q2<1.5 and x<0.009 "
  Q2b=[]
  for k in Q2bins:
      d=data.query('%s and  had=="h+"'%Q2bins[k])
      Q2b.append(d.Q2.mean())
  Q2b=np.sort(np.unique(Q2b))

  xbins={}
  xbins[(4,0)]="x<0.009 "
  xbins[(4,1)]="0.009<x and x<0.012 "
  xbins[(4,2)]="0.012<x and x<0.02 "
  xbins[(4,3)]="0.02<x  and x<0.03 "
  xbins[(4,4)]="0.03<x  and x<0.06 "
  xbins[(2,3)]="0.02<x  and x<0.03 "
  xbins[(2,4)]="0.03<x  and x<0.06 "
  xbins[(2,5)]="0.06<x  and x<0.09 "
  xbins[(2,6)]="0.09<x  and x<0.2 "
  xbins[(1,7)]="0.2<x "
  
  xb=[]
  for k in xbins:
      d=data.query('%s and  had=="h+"'%xbins[k])
      xb.append(d.x.mean())
  xb=np.sort(np.unique(xb))


  zbins=[[0.24,0.3],[0.3,0.4],[0.4,0.5],[0.65,0.7]]

  nrows,ncols=5,8
  fig = py.figure(figsize=(ncols*3,nrows*3))
  gs = gridspec.GridSpec(5,8)
  gs.update(wspace=0.,hspace=0,left=0.12, right=0.96,bottom=0.13,top=0.96)
  AX={}



  for k in bins:
      ir,ic=k
      ax = py.subplot(gs[ir,ic])
      d=data.query('%s and  had=="h+"'%bins[k])
      for i in [2]:#range(len(zbins)):
          dd=d.query('%f<z and z<%f'%(zbins[i][0],zbins[i][1]))
          if dd.index.size==0: continue

          #if   i==0: c='r'
          #elif i==1: c='b'
          #elif i==2: c='g'
          #elif i==3: c='m'

       
          # OURS
          tot1=dd.LO+dd.NLOdelta+dd.NLOplus+dd.NLOregular
 
          tot2=dd.LO1+dd.NLOdelta1+dd.NLOplus1+dd.NLOregular1
          ax.plot(tot2.values/tot1.values,label=r'$\rm channel=1$')

          tot2=dd.LO2+dd.NLOdelta2+dd.NLOplus2+dd.NLOregular2
          ax.plot(tot2.values/tot1.values,label=r'$\rm channel=2$')

          tot2=dd.LO3+dd.NLOdelta3+dd.NLOplus3+dd.NLOregular3
          ax.plot(tot2.values/tot1.values,label=r'$\rm channel=3$')

          tot2=dd.LO4+dd.NLOdelta4+dd.NLOplus4+dd.NLOregular4
          ax.plot(tot2.values/tot1.values,label=r'$\rm channel=4$')

          tot2=dd.LO5+dd.NLOdelta5+dd.NLOplus5+dd.NLOregular5
          ax.plot(tot2.values/tot1.values,label=r'$\rm channel=5$')

          tot2=dd.LO6+dd.NLOdelta6+dd.NLOplus6+dd.NLOregular6
          ax.plot(tot2.values/tot1.values,label=r'$\rm channel=6$')

          #tot2=dd.LO3+dd.NLOdelta3+dd.NLOplus3+dd.NLOregular3
          #tot2=dd.LO4+dd.NLOdelta4+dd.NLOplus4+dd.NLOregular4
          #tot2=dd.LO5+dd.NLOdelta5+dd.NLOplus5+dd.NLOregular5
          #tot2=dd.LO6+dd.NLOdelta6+dd.NLOplus6+dd.NLOregular6

          #ax.plot(tot1,label='')
          #ax.plot(tot2,label='')
 
          #ax.plot(dd.qT,tot/dd.idisNLO,'%s:'%c,label='')
      #ax.axvline(d.Q.values[0])
      #ax.semilogy()
      #ax.set_ylim(1e-4,10)
      ax.set_xlim(0,8)
      ax.set_ylim(-0.1,1.0)
      #ax.set_yticklabels([])
  
      ax.tick_params(axis='both', which='major', labelsize=30, direction='in')
      if all(k!=_ for _ in [(4,0),(4,1),(4,2),(4,3),(4,4),(3,5),(2,6),(1,7)]): 
          ax.set_xticklabels([])
      else:
          ax.set_xticks([2,4,6])
          #ax.set_xlabel(r'$q_{\rm T}$',size=35)
          #ax.xaxis.set_label_coords(0.95, -0.02)
      if all(k!=_ for _ in [(4,0),(3,0),(2,2),(1,4),(0,6)]): 
          ax.set_yticklabels([])
          ax.set_yticks([0,0.2,0.4,0.6,0.8])
      else:
          ax.set_yticks([0,0.2,0.4,0.6,0.8])
      if k==(4,0):
          ax.annotate('', xy=(-0.35, 5.2), xycoords='axes fraction', xytext=(-0.35, -0.1),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))
          ax.annotate('', xy=(8.2,-0.3), xycoords='axes fraction', xytext=(-0.1, -0.3),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))        
          ax.annotate(r'$Q^2~({\rm GeV}^2)$', xy=(-1,3),xycoords='axes fraction',size=50,rotation=90)
          ax.annotate(r'$x_{\rm bj}$', xy=(3.9,-0.7),xycoords='axes fraction',size=50)
                  
          for i in range(8):
              if xb[i]<2e-2: msg=r'$%0.3f$'%xb[i]
              else:msg=r'$%0.2f$'%xb[i]
              ax.text(0.5+i,-0.5,msg,transform=ax.transAxes,size=30,ha="center")
              ax.annotate('',xy=(i,-0.35),xycoords='axes fraction',xytext=(i+1, -0.35), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
          for i in range(5):
              ax.text(-0.6,0.5+i,r'$%0.1f$'%Q2b[i],transform=ax.transAxes,size=30,rotation=90,va="center")
              ax.annotate('',xy=(-0.4,i),xycoords='axes fraction',xytext=(-0.4,i+1), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
      #ax.plot([d.Q.values[0],8],[2e-4,2e-4],c='y',lw=10,alpha=0.5)
      if k==(2,2): 
          msg=r'${\rm COMPASS~17}~h^+$'
          ax.text(-2,2.8,msg,transform=ax.transAxes,size=50)
          msg =r'$\sigma_{\rm ch}/\sigma_{\rm tot}(\rm ours)~{\rm vs.}~q_{\rm T}{~\rm (GeV)}$'
          ax.text(-2,2.2,msg,transform=ax.transAxes,size=50)
          zmean=dd.z.mean()
          msg =r'$z=%0.2f$'%zmean
          ax.text(-2,1.6,msg,transform=ax.transAxes,size=50)
      if k==(1,4): 
          ax.legend(bbox_to_anchor=[3, -2.2], loc='center',fontsize=40,frameon=0\
                   ,markerscale=2,handletextpad=0.1)
      #ax.text(0.6,0.8,'(%d,%d)'%(ir,ic),transform=ax.transAxes,size=20)
  
          
  py.savefig('gallery/compass17-relative-ours.pdf')

def plot_compass17_relative_dds():

  data=pd.read_excel('expdata/compass17-channels.xlsx')


  bins={}
  bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  bins[(0,7)]="17<Q2 and 0.2<x "
  bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  bins[(1,5)]="6<Q2 and Q2<17 and 0.06<x  and x<0.09 "
  bins[(1,6)]="6<Q2 and Q2<17 and 0.09<x  and x<0.2 "
  bins[(1,7)]="6<Q2 and Q2<17 and 0.2<x "
  bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  bins[(2,3)]="3<Q2 and Q2<6 and 0.02<x  and x<0.03 "
  bins[(2,4)]="3<Q2 and Q2<6 and 0.03<x  and x<0.06 "
  bins[(2,5)]="3<Q2 and Q2<6 and 0.06<x  and x<0.09 "
  bins[(2,6)]="3<Q2 and Q2<6 and 0.09<x  and x<0.2 "
  bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  bins[(3,1)]="1.5<Q2 and Q2<3.0 and 0.009<x and x<0.012 "
  bins[(3,2)]="1.5<Q2 and Q2<3.0 and 0.012<x and x<0.02 "
  bins[(3,3)]="1.5<Q2 and Q2<3.0 and 0.02<x  and x<0.03 "
  bins[(3,4)]="1.5<Q2 and Q2<3.0 and 0.03<x  and x<0.06 "
  bins[(3,5)]="1.5<Q2 and Q2<3.0 and 0.06<x  and x<0.09 "
  bins[(4,0)]="Q2<1.5 and x<0.009 "
  bins[(4,1)]="Q2<1.5 and 0.009<x and x<0.012 "
  bins[(4,2)]="Q2<1.5 and 0.012<x and x<0.02 "
  bins[(4,3)]="Q2<1.5 and 0.02<x  and x<0.03 "
  bins[(4,4)]="Q2<1.5 and 0.03<x  and x<0.06 "

  Q2bins={}
  Q2bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  Q2bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  Q2bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  Q2bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  Q2bins[(4,0)]="Q2<1.5 and x<0.009 "
  Q2b=[]
  for k in Q2bins:
      d=data.query('%s and  had=="h+"'%Q2bins[k])
      Q2b.append(d.Q2.mean())
  Q2b=np.sort(np.unique(Q2b))

  xbins={}
  xbins[(4,0)]="x<0.009 "
  xbins[(4,1)]="0.009<x and x<0.012 "
  xbins[(4,2)]="0.012<x and x<0.02 "
  xbins[(4,3)]="0.02<x  and x<0.03 "
  xbins[(4,4)]="0.03<x  and x<0.06 "
  xbins[(2,3)]="0.02<x  and x<0.03 "
  xbins[(2,4)]="0.03<x  and x<0.06 "
  xbins[(2,5)]="0.06<x  and x<0.09 "
  xbins[(2,6)]="0.09<x  and x<0.2 "
  xbins[(1,7)]="0.2<x "
  
  xb=[]
  for k in xbins:
      d=data.query('%s and  had=="h+"'%xbins[k])
      xb.append(d.x.mean())
  xb=np.sort(np.unique(xb))


  zbins=[[0.24,0.3],[0.3,0.4],[0.4,0.5],[0.65,0.7]]

  nrows,ncols=5,8
  fig = py.figure(figsize=(ncols*3,nrows*3))
  gs = gridspec.GridSpec(5,8)
  gs.update(wspace=0.,hspace=0,left=0.12, right=0.96,bottom=0.13,top=0.96)
  AX={}



  for k in bins:
      ir,ic=k
      ax = py.subplot(gs[ir,ic])
      d=data.query('%s and  had=="h+"'%bins[k])
      for i in [2]:#range(len(zbins)):
          dd=d.query('%f<z and z<%f'%(zbins[i][0],zbins[i][1]))
          if dd.index.size==0: continue


          tot1 = dd.DDSLO+dd.DDSNLOdelta+dd.DDSNLOplus+dd.DDSNLOrest

          # channel 1
          tot2 = dd.DDSLO5+dd.DDSNLOdelta5+dd.DDSNLOplus5+dd.DDSNLOrest5
          ax.plot(tot2.values/tot1.values,label=r'$\rm channel=1$')

          # channel 2
          tot2 = dd.DDSLO1+dd.DDSNLOdelta1+dd.DDSNLOplus1+dd.DDSNLOrest1
          tot2+= dd.DDSLO2+dd.DDSNLOdelta2+dd.DDSNLOplus2+dd.DDSNLOrest2
          tot2+= dd.DDSLO3+dd.DDSNLOdelta3+dd.DDSNLOplus3+dd.DDSNLOrest3
          ax.plot(tot2.values/tot1.values,label=r'$\rm channel=2$')

          # channel 3  +
          tot2 = dd.DDSLO4+dd.DDSNLOdelta4+dd.DDSNLOplus4+dd.DDSNLOrest4
          ax.plot(tot2.values/tot1.values,label=r'$\rm channel=3$')

          # channel 4
          tot2 = dd.DDSLO6+dd.DDSNLOdelta6+dd.DDSNLOplus6+dd.DDSNLOrest6
          ax.plot(tot2.values/tot1.values,label=r'$\rm channel=4$')

          # channel 5:  ours is zero???
          tot2 = dd.DDSLO7+dd.DDSNLOdelta7+dd.DDSNLOplus7+dd.DDSNLOrest7
          ax.plot(tot2.values/tot1.values,label=r'$\rm channel=5$')

          # channel 6: 
          tot2 = dd.DDSLO8+dd.DDSNLOdelta8+dd.DDSNLOplus8+dd.DDSNLOrest8
          tot2+= dd.DDSLO9+dd.DDSNLOdelta9+dd.DDSNLOplus9+dd.DDSNLOrest9
          tot2+= dd.DDSLO10+dd.DDSNLOdelta10+dd.DDSNLOplus10+dd.DDSNLOrest10
          ax.plot(tot2.values/tot1.values,label=r'$\rm channel=6$')


 
          #ax.plot(dd.qT,tot/dd.idisNLO,'%s:'%c,label='')
      #ax.axvline(d.Q.values[0])
      #ax.semilogy()
      #ax.set_ylim(1e-4,10)
      ax.set_xlim(0,8)
      ax.set_ylim(-0.1,1.0)
      #ax.set_yticklabels([])
  
      ax.tick_params(axis='both', which='major', labelsize=30, direction='in')
      if all(k!=_ for _ in [(4,0),(4,1),(4,2),(4,3),(4,4),(3,5),(2,6),(1,7)]): 
          ax.set_xticklabels([])
      else:
          ax.set_xticks([2,4,6])
          #ax.set_xlabel(r'$q_{\rm T}$',size=35)
          #ax.xaxis.set_label_coords(0.95, -0.02)
      if all(k!=_ for _ in [(4,0),(3,0),(2,2),(1,4),(0,6)]): 
          ax.set_yticklabels([])
          ax.set_yticks([0,0.2,0.4,0.6,0.8])
      else:
          ax.set_yticks([0,0.2,0.4,0.6,0.8])
      if k==(4,0):
          ax.annotate('', xy=(-0.35, 5.2), xycoords='axes fraction', xytext=(-0.35, -0.1),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))
          ax.annotate('', xy=(8.2,-0.3), xycoords='axes fraction', xytext=(-0.1, -0.3),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))        
          ax.annotate(r'$Q^2~({\rm GeV}^2)$', xy=(-1,3),xycoords='axes fraction',size=50,rotation=90)
          ax.annotate(r'$x_{\rm bj}$', xy=(3.9,-0.7),xycoords='axes fraction',size=50)
                  
          for i in range(8):
              if xb[i]<2e-2: msg=r'$%0.3f$'%xb[i]
              else:msg=r'$%0.2f$'%xb[i]
              ax.text(0.5+i,-0.5,msg,transform=ax.transAxes,size=30,ha="center")
              ax.annotate('',xy=(i,-0.35),xycoords='axes fraction',xytext=(i+1, -0.35), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
          for i in range(5):
              ax.text(-0.6,0.5+i,r'$%0.1f$'%Q2b[i],transform=ax.transAxes,size=30,rotation=90,va="center")
              ax.annotate('',xy=(-0.4,i),xycoords='axes fraction',xytext=(-0.4,i+1), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
      #ax.plot([d.Q.values[0],8],[2e-4,2e-4],c='y',lw=10,alpha=0.5)
      if k==(2,2): 
          msg=r'${\rm COMPASS~17}~h^+$'
          ax.text(-2,2.8,msg,transform=ax.transAxes,size=50)
          msg =r'$\sigma_{\rm ch}/\sigma_{\rm tot}(\rm DDS)~{\rm vs.}~q_{\rm T}{~\rm (GeV)}$'
          ax.text(-2,2.2,msg,transform=ax.transAxes,size=50)
          zmean=dd.z.mean()
          msg =r'$z=%0.2f$'%zmean
          ax.text(-2,1.6,msg,transform=ax.transAxes,size=50)

      if k==(1,4): 
          ax.legend(bbox_to_anchor=[3, -2.2], loc='center',fontsize=40,frameon=0\
                   ,markerscale=2,handletextpad=0.1)
      #ax.text(0.6,0.8,'(%d,%d)'%(ir,ic),transform=ax.transAxes,size=20)
  
          
  py.savefig('gallery/compass17-relative-dds.pdf')

def plot_compass17_nlo_channels_diff(channel):

  data=pd.read_excel('expdata/compass17-channels.xlsx')


  bins={}
  bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  bins[(0,7)]="17<Q2 and 0.2<x "
  bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  bins[(1,5)]="6<Q2 and Q2<17 and 0.06<x  and x<0.09 "
  bins[(1,6)]="6<Q2 and Q2<17 and 0.09<x  and x<0.2 "
  bins[(1,7)]="6<Q2 and Q2<17 and 0.2<x "
  bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  bins[(2,3)]="3<Q2 and Q2<6 and 0.02<x  and x<0.03 "
  bins[(2,4)]="3<Q2 and Q2<6 and 0.03<x  and x<0.06 "
  bins[(2,5)]="3<Q2 and Q2<6 and 0.06<x  and x<0.09 "
  bins[(2,6)]="3<Q2 and Q2<6 and 0.09<x  and x<0.2 "
  bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  bins[(3,1)]="1.5<Q2 and Q2<3.0 and 0.009<x and x<0.012 "
  bins[(3,2)]="1.5<Q2 and Q2<3.0 and 0.012<x and x<0.02 "
  bins[(3,3)]="1.5<Q2 and Q2<3.0 and 0.02<x  and x<0.03 "
  bins[(3,4)]="1.5<Q2 and Q2<3.0 and 0.03<x  and x<0.06 "
  bins[(3,5)]="1.5<Q2 and Q2<3.0 and 0.06<x  and x<0.09 "
  bins[(4,0)]="Q2<1.5 and x<0.009 "
  bins[(4,1)]="Q2<1.5 and 0.009<x and x<0.012 "
  bins[(4,2)]="Q2<1.5 and 0.012<x and x<0.02 "
  bins[(4,3)]="Q2<1.5 and 0.02<x  and x<0.03 "
  bins[(4,4)]="Q2<1.5 and 0.03<x  and x<0.06 "

  Q2bins={}
  Q2bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  Q2bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  Q2bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  Q2bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  Q2bins[(4,0)]="Q2<1.5 and x<0.009 "
  Q2b=[]
  for k in Q2bins:
      d=data.query('%s and  had=="h+"'%Q2bins[k])
      Q2b.append(d.Q2.mean())
  Q2b=np.sort(np.unique(Q2b))

  xbins={}
  xbins[(4,0)]="x<0.009 "
  xbins[(4,1)]="0.009<x and x<0.012 "
  xbins[(4,2)]="0.012<x and x<0.02 "
  xbins[(4,3)]="0.02<x  and x<0.03 "
  xbins[(4,4)]="0.03<x  and x<0.06 "
  xbins[(2,3)]="0.02<x  and x<0.03 "
  xbins[(2,4)]="0.03<x  and x<0.06 "
  xbins[(2,5)]="0.06<x  and x<0.09 "
  xbins[(2,6)]="0.09<x  and x<0.2 "
  xbins[(1,7)]="0.2<x "
  
  xb=[]
  for k in xbins:
      d=data.query('%s and  had=="h+"'%xbins[k])
      xb.append(d.x.mean())
  xb=np.sort(np.unique(xb))


  zbins=[[0.24,0.3],[0.3,0.4],[0.4,0.5],[0.65,0.7]]

  nrows,ncols=5,8
  fig = py.figure(figsize=(ncols*3,nrows*3))
  gs = gridspec.GridSpec(5,8)
  gs.update(wspace=0.,hspace=0,left=0.12, right=0.96,bottom=0.13,top=0.96)
  AX={}



  for k in bins:
      ir,ic=k
      ax = py.subplot(gs[ir,ic])
      d=data.query('%s and  had=="h+"'%bins[k])
      for i in range(len(zbins)):
          dd=d.query('%f<z and z<%f'%(zbins[i][0],zbins[i][1]))
          if dd.index.size==0: continue
       
          # OURS
          #tot1=dd.LO+dd.NLOdelta+dd.NLOplus+dd.NLOregular
          #tot2 =dd.LO1+dd.NLOdelta1+dd.NLOplus1+dd.NLOregular1
          #tot2+=dd.LO2+dd.NLOdelta2+dd.NLOplus2+dd.NLOregular2
          #tot2+=dd.LO3+dd.NLOdelta3+dd.NLOplus3+dd.NLOregular3
          #tot2+=dd.LO4+dd.NLOdelta4+dd.NLOplus4+dd.NLOregular4
          #tot2+=dd.LO5+dd.NLOdelta5+dd.NLOplus5+dd.NLOregular5
          #tot2+=dd.LO6+dd.NLOdelta6+dd.NLOplus6+dd.NLOregular6

          # DDS
          #tot1 = dd.DDSLO+dd.DDSNLOdelta+dd.DDSNLOplus+dd.DDSNLOrest
          #tot2 = dd.DDSLO1+dd.DDSNLOdelta1+dd.DDSNLOplus1+dd.DDSNLOrest1
          #tot2+= dd.DDSLO2+dd.DDSNLOdelta2+dd.DDSNLOplus2+dd.DDSNLOrest2
          #tot2+= dd.DDSLO3+dd.DDSNLOdelta3+dd.DDSNLOplus3+dd.DDSNLOrest3
          #tot2+= dd.DDSLO4+dd.DDSNLOdelta4+dd.DDSNLOplus4+dd.DDSNLOrest4
          #tot2+= dd.DDSLO5+dd.DDSNLOdelta5+dd.DDSNLOplus5+dd.DDSNLOrest5
          #tot2+= dd.DDSLO6+dd.DDSNLOdelta6+dd.DDSNLOplus6+dd.DDSNLOrest6
          #tot2+= dd.DDSLO7+dd.DDSNLOdelta7+dd.DDSNLOplus7+dd.DDSNLOrest7
          #tot2+= dd.DDSLO8+dd.DDSNLOdelta8+dd.DDSNLOplus8+dd.DDSNLOrest8
          #tot2+= dd.DDSLO9+dd.DDSNLOdelta9+dd.DDSNLOplus9+dd.DDSNLOrest9
          #tot2+= dd.DDSLO10+dd.DDSNLOdelta10+dd.DDSNLOplus10+dd.DDSNLOrest10



          if channel== 1:
            tot1 = dd.LO1+dd.NLOdelta1+dd.NLOplus1+dd.NLOregular1
            tot2 = dd.DDSLO5+dd.DDSNLOdelta5+dd.DDSNLOplus5+dd.DDSNLOrest5

          if channel== 2:
            tot1 = dd.LO2+dd.NLOdelta2+dd.NLOplus2+dd.NLOregular2
            tot2 = dd.DDSLO1+dd.DDSNLOdelta1+dd.DDSNLOplus1+dd.DDSNLOrest1
            tot2+= dd.DDSLO2+dd.DDSNLOdelta2+dd.DDSNLOplus2+dd.DDSNLOrest2
            tot2+= dd.DDSLO3+dd.DDSNLOdelta3+dd.DDSNLOplus3+dd.DDSNLOrest3

          if channel== 3:
            tot1 = dd.LO3+dd.NLOdelta3+dd.NLOplus3+dd.NLOregular3
            tot2 = dd.DDSLO4+dd.DDSNLOdelta4+dd.DDSNLOplus4+dd.DDSNLOrest4

          if channel== 4:
            tot1 = dd.LO4+dd.NLOdelta4+dd.NLOplus4+dd.NLOregular4
            tot2 = dd.DDSLO6+dd.DDSNLOdelta6+dd.DDSNLOplus6+dd.DDSNLOrest6

          if channel== 5:
            tot1 = dd.LO5+dd.NLOdelta5+dd.NLOplus5+dd.NLOregular5
            tot2 = dd.DDSLO7+dd.DDSNLOdelta7+dd.DDSNLOplus7+dd.DDSNLOrest7

          if channel== 6:
            tot1 = dd.LO6+dd.NLOdelta6+dd.NLOplus6+dd.NLOregular6
            tot2 = dd.DDSLO8+dd.DDSNLOdelta8+dd.DDSNLOplus8+dd.DDSNLOrest8
            tot2+= dd.DDSLO9+dd.DDSNLOdelta9+dd.DDSNLOplus9+dd.DDSNLOrest9
            tot2+= dd.DDSLO10+dd.DDSNLOdelta10+dd.DDSNLOplus10+dd.DDSNLOrest10

          # totals (using combined channels)
          #tot1 = dd.LO+dd.NLOdelta+dd.NLOplus+dd.NLOregular
          #tot2 = dd.DDSLO+dd.DDSNLOdelta+dd.DDSNLOplus+dd.DDSNLOrest

          # totals (using summed channels)
          #tot1 =dd.LO1+dd.NLOdelta1+dd.NLOplus1+dd.NLOregular1
          #tot1+=dd.LO2+dd.NLOdelta2+dd.NLOplus2+dd.NLOregular2
          #tot1+=dd.LO3+dd.NLOdelta3+dd.NLOplus3+dd.NLOregular3
          #tot1+=dd.LO4+dd.NLOdelta4+dd.NLOplus4+dd.NLOregular4
          #tot1+=dd.LO5+dd.NLOdelta5+dd.NLOplus5+dd.NLOregular5
          #tot1+=dd.LO6+dd.NLOdelta6+dd.NLOplus6+dd.NLOregular6

          #tot2 = dd.DDSLO1+dd.DDSNLOdelta1+dd.DDSNLOplus1+dd.DDSNLOrest1
          #tot2+= dd.DDSLO2+dd.DDSNLOdelta2+dd.DDSNLOplus2+dd.DDSNLOrest2
          #tot2+= dd.DDSLO3+dd.DDSNLOdelta3+dd.DDSNLOplus3+dd.DDSNLOrest3
          #tot2+= dd.DDSLO4+dd.DDSNLOdelta4+dd.DDSNLOplus4+dd.DDSNLOrest4
          #tot2+= dd.DDSLO5+dd.DDSNLOdelta5+dd.DDSNLOplus5+dd.DDSNLOrest5
          #tot2+= dd.DDSLO6+dd.DDSNLOdelta6+dd.DDSNLOplus6+dd.DDSNLOrest6
          #tot2+= dd.DDSLO7+dd.DDSNLOdelta7+dd.DDSNLOplus7+dd.DDSNLOrest7
          #tot2+= dd.DDSLO8+dd.DDSNLOdelta8+dd.DDSNLOplus8+dd.DDSNLOrest8
          #tot2+= dd.DDSLO9+dd.DDSNLOdelta9+dd.DDSNLOplus9+dd.DDSNLOrest9
          #tot2+= dd.DDSLO10+dd.DDSNLOdelta10+dd.DDSNLOplus10+dd.DDSNLOrest10

          ax.plot(np.abs((tot1.values-tot2.values)/tot1.values*100),label=r'$%0.2f<z<%0.2f$'%(zbins[i][0],zbins[i][1]))
          #ax.plot(tot1,label='')
          #ax.plot(tot2,label='')
 
          #ax.plot(dd.qT,tot/dd.idisNLO,'%s:'%c,label='')
      #ax.axvline(d.Q.values[0])
      #ax.semilogy()
      #ax.set_ylim(1e-4,10)
      ax.set_xlim(0,8)
      ax.set_ylim(0,50)
      #ax.set_yticklabels([])
  
      ax.tick_params(axis='both', which='major', labelsize=30, direction='in')
      if all(k!=_ for _ in [(4,0),(4,1),(4,2),(4,3),(4,4),(3,5),(2,6),(1,7)]): 
          ax.set_xticklabels([])
      else:
          ax.set_xticks([2,4,6])
          #ax.set_xlabel(r'$q_{\rm T}$',size=35)
          #ax.xaxis.set_label_coords(0.95, -0.02)
      if all(k!=_ for _ in [(4,0),(3,0),(2,2),(1,4),(0,6)]): 
          ax.set_yticklabels([])
          ax.set_yticks([10,20,30,40])
      else:
          ax.set_yticks([10,20,30,40])
          pass
      if k==(4,0):
          ax.annotate('', xy=(-0.35, 5.2), xycoords='axes fraction', xytext=(-0.35, -0.1),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))
          ax.annotate('', xy=(8.2,-0.3), xycoords='axes fraction', xytext=(-0.1, -0.3),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))        
          ax.annotate(r'$Q^2~({\rm GeV}^2)$', xy=(-1,3),xycoords='axes fraction',size=50,rotation=90)
          ax.annotate(r'$x_{\rm bj}$', xy=(3.9,-0.7),xycoords='axes fraction',size=50)
                  
          for i in range(8):
              if xb[i]<2e-2: msg=r'$%0.3f$'%xb[i]
              else:msg=r'$%0.2f$'%xb[i]
              ax.text(0.5+i,-0.5,msg,transform=ax.transAxes,size=30,ha="center")
              ax.annotate('',xy=(i,-0.35),xycoords='axes fraction',xytext=(i+1, -0.35), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
          for i in range(5):
              ax.text(-0.6,0.5+i,r'$%0.1f$'%Q2b[i],transform=ax.transAxes,size=30,rotation=90,va="center")
              ax.annotate('',xy=(-0.4,i),xycoords='axes fraction',xytext=(-0.4,i+1), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
      #ax.plot([d.Q.values[0],8],[2e-4,2e-4],c='y',lw=10,alpha=0.5)
      if k==(2,2): 
          #LO1b = mpatches.Rectangle((0,0), 0.05, 0.1, ec="none",color='k',alpha=0.1)
          #LO1=py.Line2D([0,2], [0,0], color='k',ls='--')
          #LO2=py.Line2D([0], [0], color='k',ls='-')
          #qTrange = mpatches.Rectangle((0,0), 0, 0, ec="none",color='y',alpha=0.5)
          #qTrange=py.Line2D([0], [0], color='y',ls='-',lw=10,alpha=0.5)
          #ax.legend([LO1,LO2,qTrange],[r'$\rm DDS~(LO)$',r'$\rm ours~(LO)$',r'$q_{\rm T}>Q$']\
          #          ,bbox_to_anchor=[-1.2, 1.]\
          #          ,loc='center',fontsize=40,frameon=0)
          msg=r'${\rm COMPASS~17}~h^+$'
          ax.text(-2,2.8,msg,transform=ax.transAxes,size=50)
          #msg =r'$\frac{d\sigma}{dx_{\rm bj}dQ^2dzdP^2_{\rm T}}/\frac{d\sigma}{dx_{\rm bj}dQ^2}'
          #msg+=r'({\tiny {\rm GeV}^{-2}})'
          msg =r'$|{\rm (ours-daleo)/ours}|\times 100~{\rm vs.}~q_{\rm T}{~\rm (GeV)}$'
          ax.text(-2,2.2,msg,transform=ax.transAxes,size=50)
          msg =r'$\rm channel~%d$'%channel
          ax.text(-2,1.6,msg,transform=ax.transAxes,size=50)
      if k==(1,4): 
          ax.legend(bbox_to_anchor=[3, -2.5], loc='center',fontsize=40,frameon=0\
                   ,markerscale=2,handletextpad=0.1)
      #ax.text(0.6,0.8,'(%d,%d)'%(ir,ic),transform=ax.transAxes,size=20)
  
          
  py.savefig('gallery/compass17-nlo-channel-%d.pdf'%channel)

def plot_compass17_nlo_channels(channel):

  data=pd.read_excel('expdata/compass17-channels.xlsx')


  bins={}
  bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  bins[(0,7)]="17<Q2 and 0.2<x "
  bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  bins[(1,5)]="6<Q2 and Q2<17 and 0.06<x  and x<0.09 "
  bins[(1,6)]="6<Q2 and Q2<17 and 0.09<x  and x<0.2 "
  bins[(1,7)]="6<Q2 and Q2<17 and 0.2<x "
  bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  bins[(2,3)]="3<Q2 and Q2<6 and 0.02<x  and x<0.03 "
  bins[(2,4)]="3<Q2 and Q2<6 and 0.03<x  and x<0.06 "
  bins[(2,5)]="3<Q2 and Q2<6 and 0.06<x  and x<0.09 "
  bins[(2,6)]="3<Q2 and Q2<6 and 0.09<x  and x<0.2 "
  bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  bins[(3,1)]="1.5<Q2 and Q2<3.0 and 0.009<x and x<0.012 "
  bins[(3,2)]="1.5<Q2 and Q2<3.0 and 0.012<x and x<0.02 "
  bins[(3,3)]="1.5<Q2 and Q2<3.0 and 0.02<x  and x<0.03 "
  bins[(3,4)]="1.5<Q2 and Q2<3.0 and 0.03<x  and x<0.06 "
  bins[(3,5)]="1.5<Q2 and Q2<3.0 and 0.06<x  and x<0.09 "
  bins[(4,0)]="Q2<1.5 and x<0.009 "
  bins[(4,1)]="Q2<1.5 and 0.009<x and x<0.012 "
  bins[(4,2)]="Q2<1.5 and 0.012<x and x<0.02 "
  bins[(4,3)]="Q2<1.5 and 0.02<x  and x<0.03 "
  bins[(4,4)]="Q2<1.5 and 0.03<x  and x<0.06 "

  Q2bins={}
  Q2bins[(0,6)]="17<Q2 and 0.09<x  and x<0.2 "
  Q2bins[(1,4)]="6<Q2 and Q2<17 and 0.03<x  and x<0.06 "
  Q2bins[(2,2)]="3<Q2 and Q2<6 and 0.012<x and x<0.02 "
  Q2bins[(3,0)]="1.5<Q2 and Q2<3.0 and x<0.009 "
  Q2bins[(4,0)]="Q2<1.5 and x<0.009 "
  Q2b=[]
  for k in Q2bins:
      d=data.query('%s and  had=="h+"'%Q2bins[k])
      Q2b.append(d.Q2.mean())
  Q2b=np.sort(np.unique(Q2b))

  xbins={}
  xbins[(4,0)]="x<0.009 "
  xbins[(4,1)]="0.009<x and x<0.012 "
  xbins[(4,2)]="0.012<x and x<0.02 "
  xbins[(4,3)]="0.02<x  and x<0.03 "
  xbins[(4,4)]="0.03<x  and x<0.06 "
  xbins[(2,3)]="0.02<x  and x<0.03 "
  xbins[(2,4)]="0.03<x  and x<0.06 "
  xbins[(2,5)]="0.06<x  and x<0.09 "
  xbins[(2,6)]="0.09<x  and x<0.2 "
  xbins[(1,7)]="0.2<x "
  
  xb=[]
  for k in xbins:
      d=data.query('%s and  had=="h+"'%xbins[k])
      xb.append(d.x.mean())
  xb=np.sort(np.unique(xb))


  zbins=[[0.24,0.3],[0.3,0.4],[0.4,0.5],[0.65,0.7]]

  nrows,ncols=5,8
  fig = py.figure(figsize=(ncols*3,nrows*3))
  gs = gridspec.GridSpec(5,8)
  gs.update(wspace=0.,hspace=0,left=0.12, right=0.96,bottom=0.13,top=0.96)
  AX={}



  for k in bins:
      ir,ic=k
      ax = py.subplot(gs[ir,ic])
      d=data.query('%s and  had=="h+"'%bins[k])
      for i in range(len(zbins)):
          dd=d.query('%f<z and z<%f'%(zbins[i][0],zbins[i][1]))
          if dd.index.size==0: continue
       
          # OURS
          #tot1=dd.LO+dd.NLOdelta+dd.NLOplus+dd.NLOregular
          #tot2 =dd.LO1+dd.NLOdelta1+dd.NLOplus1+dd.NLOregular1
          #tot2+=dd.LO2+dd.NLOdelta2+dd.NLOplus2+dd.NLOregular2
          #tot2+=dd.LO3+dd.NLOdelta3+dd.NLOplus3+dd.NLOregular3
          #tot2+=dd.LO4+dd.NLOdelta4+dd.NLOplus4+dd.NLOregular4
          #tot2+=dd.LO5+dd.NLOdelta5+dd.NLOplus5+dd.NLOregular5
          #tot2+=dd.LO6+dd.NLOdelta6+dd.NLOplus6+dd.NLOregular6

          # DDS
          #tot1 = dd.DDSLO+dd.DDSNLOdelta+dd.DDSNLOplus+dd.DDSNLOrest
          #tot2 = dd.DDSLO1+dd.DDSNLOdelta1+dd.DDSNLOplus1+dd.DDSNLOrest1
          #tot2+= dd.DDSLO2+dd.DDSNLOdelta2+dd.DDSNLOplus2+dd.DDSNLOrest2
          #tot2+= dd.DDSLO3+dd.DDSNLOdelta3+dd.DDSNLOplus3+dd.DDSNLOrest3
          #tot2+= dd.DDSLO4+dd.DDSNLOdelta4+dd.DDSNLOplus4+dd.DDSNLOrest4
          #tot2+= dd.DDSLO5+dd.DDSNLOdelta5+dd.DDSNLOplus5+dd.DDSNLOrest5
          #tot2+= dd.DDSLO6+dd.DDSNLOdelta6+dd.DDSNLOplus6+dd.DDSNLOrest6
          #tot2+= dd.DDSLO7+dd.DDSNLOdelta7+dd.DDSNLOplus7+dd.DDSNLOrest7
          #tot2+= dd.DDSLO8+dd.DDSNLOdelta8+dd.DDSNLOplus8+dd.DDSNLOrest8
          #tot2+= dd.DDSLO9+dd.DDSNLOdelta9+dd.DDSNLOplus9+dd.DDSNLOrest9
          #tot2+= dd.DDSLO10+dd.DDSNLOdelta10+dd.DDSNLOplus10+dd.DDSNLOrest10



          if channel== 1:
            tot1 = dd.LO1+dd.NLOdelta1+dd.NLOplus1+dd.NLOregular1
            tot2 = dd.DDSLO5+dd.DDSNLOdelta5+dd.DDSNLOplus5+dd.DDSNLOrest5

          if channel== 2:
            tot1 = dd.LO2+dd.NLOdelta2+dd.NLOplus2+dd.NLOregular2
            tot2 = dd.DDSLO1+dd.DDSNLOdelta1+dd.DDSNLOplus1+dd.DDSNLOrest1
            tot2+= dd.DDSLO2+dd.DDSNLOdelta2+dd.DDSNLOplus2+dd.DDSNLOrest2
            tot2+= dd.DDSLO3+dd.DDSNLOdelta3+dd.DDSNLOplus3+dd.DDSNLOrest3

          if channel== 3:
            tot1 = dd.LO3+dd.NLOdelta3+dd.NLOplus3+dd.NLOregular3
            tot2 = dd.DDSLO4+dd.DDSNLOdelta4+dd.DDSNLOplus4+dd.DDSNLOrest4

          if channel== 4:
            tot1 = dd.LO4+dd.NLOdelta4+dd.NLOplus4+dd.NLOregular4
            tot2 = dd.DDSLO6+dd.DDSNLOdelta6+dd.DDSNLOplus6+dd.DDSNLOrest6

          if channel== 5:
            tot1 = dd.LO5+dd.NLOdelta5+dd.NLOplus5+dd.NLOregular5
            tot2 = dd.DDSLO7+dd.DDSNLOdelta7+dd.DDSNLOplus7+dd.DDSNLOrest7

          if channel== 6:
            tot1 = dd.LO6+dd.NLOdelta6+dd.NLOplus6+dd.NLOregular6
            tot2 = dd.DDSLO8+dd.DDSNLOdelta8+dd.DDSNLOplus8+dd.DDSNLOrest8
            tot2+= dd.DDSLO9+dd.DDSNLOdelta9+dd.DDSNLOplus9+dd.DDSNLOrest9
            tot2+= dd.DDSLO10+dd.DDSNLOdelta10+dd.DDSNLOplus10+dd.DDSNLOrest10


          e=ax.plot(tot1.values,label=r'$%0.2f<z<%0.2f$'%(zbins[i][0],zbins[i][1]))
          c=e[0].get_color()
          ax.plot(tot2.values,color=c,ls='--',label='')

          #ax.plot(tot1,label='')
          #ax.plot(tot2,label='')
 
          #ax.plot(dd.qT,tot/dd.idisNLO,'%s:'%c,label='')
      #ax.axvline(d.Q.values[0])
      #ax.semilogy()
      #ax.set_ylim(1e-4,10)
      ax.set_xlim(0,8)
      #ax.set_ylim(0,50)
      #ax.set_yticklabels([])
  
      ax.tick_params(axis='both', which='major', labelsize=30, direction='in')
      if all(k!=_ for _ in [(4,0),(4,1),(4,2),(4,3),(4,4),(3,5),(2,6),(1,7)]): 
          ax.set_xticklabels([])
      else:
          ax.set_xticks([2,4,6])
          #ax.set_xlabel(r'$q_{\rm T}$',size=35)
          #ax.xaxis.set_label_coords(0.95, -0.02)
      if all(k!=_ for _ in [(4,0),(3,0),(2,2),(1,4),(0,6)]): 
          ax.set_yticklabels([])
          #ax.set_yticks([10,20,30,40])
      else:
          #ax.set_yticks([10,20,30,40])
          pass
      if k==(4,0):
          ax.annotate('', xy=(-0.35, 5.2), xycoords='axes fraction', xytext=(-0.35, -0.1),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))
          ax.annotate('', xy=(8.2,-0.3), xycoords='axes fraction', xytext=(-0.1, -0.3),
                      arrowprops=dict(arrowstyle="-|>, head_width=1, head_length=2", color='k',lw=3))        
          ax.annotate(r'$Q^2~({\rm GeV}^2)$', xy=(-1,3),xycoords='axes fraction',size=50,rotation=90)
          ax.annotate(r'$x_{\rm bj}$', xy=(3.9,-0.7),xycoords='axes fraction',size=50)
                  
          for i in range(8):
              if xb[i]<2e-2: msg=r'$%0.3f$'%xb[i]
              else:msg=r'$%0.2f$'%xb[i]
              ax.text(0.5+i,-0.5,msg,transform=ax.transAxes,size=30,ha="center")
              ax.annotate('',xy=(i,-0.35),xycoords='axes fraction',xytext=(i+1, -0.35), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
          for i in range(5):
              ax.text(-0.6,0.5+i,r'$%0.1f$'%Q2b[i],transform=ax.transAxes,size=30,rotation=90,va="center")
              ax.annotate('',xy=(-0.4,i),xycoords='axes fraction',xytext=(-0.4,i+1), \
                          arrowprops=dict(arrowstyle="<->", color='k'))
  
      #ax.plot([d.Q.values[0],8],[2e-4,2e-4],c='y',lw=10,alpha=0.5)
      if k==(2,2): 
          LO1b = mpatches.Rectangle((0,0), 0.05, 0.1, ec="none",color='k',alpha=0.1)
          LO1=py.Line2D([0,2], [0,0], color='k',ls='--')
          LO2=py.Line2D([0], [0], color='k',ls='-')
          ax.legend([LO1,LO2],[r'$\rm DDS~(LO)$',r'$\rm ours~(LO)$']\
                    ,bbox_to_anchor=[-1.2, 1.5]\
                    ,loc='center',fontsize=40,frameon=0)
          msg=r'${\rm COMPASS~17}~h^+$'
          ax.text(-2,2.8,msg,transform=ax.transAxes,size=50)
          msg =r'$\frac{d\sigma}{dx_{\rm bj}dQ^2dzdP^2_{\rm T}}'
          msg+=r'({\tiny {\rm GeV}^{-2}})'
          msg+=r'~{\rm vs.}~q_{\rm T}{~\rm (GeV)}$'
          ax.text(-2,2.2,msg,transform=ax.transAxes,size=50)
      if k==(1,4): 
          ax.legend(bbox_to_anchor=[3, -2.5], loc='center',fontsize=40,frameon=0\
                   ,markerscale=2,handletextpad=0.1)
      #ax.text(0.6,0.8,'(%d,%d)'%(ir,ic),transform=ax.transAxes,size=20)
  
          
  py.savefig('gallery/compass17-nlo-channel-%d-abs.pdf'%channel)



if __name__=="__main__":


  #plot_compass17()
  #plot_compass17_K()
  #plot_compass17_rat()
  #plot_hermes()
  #plot_hermes_rat()


  #plot_compass17_lo_diff()
  #plot_compass17_nlo_diff()
  #plot_compass17_relative_ours()
  #plot_compass17_relative_dds()
  #for channel in range(1,7): plot_compass17_nlo_channels_diff(channel)
  plot_compass17_nlo_channels(5)
















