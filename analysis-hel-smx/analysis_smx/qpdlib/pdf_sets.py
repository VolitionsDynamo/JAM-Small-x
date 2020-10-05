#--setup canvas
nrows,ncols=3,2
fig = py.figure(figsize=(ncols*5,nrows*2.5))
AX=[py.subplot(nrows,ncols,cnt+1) for cnt in range(6)]


#--uv and dv    
ax=AX[0]
H,L={},{}
H['JAM']=plot_bandA(ax,'um','r',r'$\rm JAM$')
L['JAM']=r'$\rm\bf JAM19$'
for _ in ['MMHT', 'CJ', 'ABMP', 'NNPDF','CSKK']: 
    if _=='CJ': label=r'$\rm CJ15$'
    elif _=='NNPDF': label=r'$\rm NNPDF3.1$'
    elif _=='MMHT': label=r'$\rm MMHT14$'
    elif _=='ABMP': label=r'$\rm ABMP16$'
    else: label=r'$\rm %s$'%_
    L[_]=label
    H[_]=plot_bandB(ax,groups[_],'uv',color=groups[_]['color'],alpha=0.2,label=r'$\rm %s$'%_)
plot_bandA(ax,'dm','r',r'$\rm JAM$')
for _ in groups:
    plot_bandB(ax,groups[_],'dv',color=groups[_]['color'],alpha=0.2,label=r'$\rm %s$'%_)
ax.text(0.81,0.78,r'\boldmath$xu_v$',color='k',transform=ax.transAxes,size=28)
ax.text(0.60,0.08,r'\boldmath$xd_v$',color='k',transform=ax.transAxes,size=28)
ax.set_ylim(0,0.9)
ax.set_yticks([0,0.2,0.4,0.6,0.8])
ax.set_yticklabels([r'$0$',r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$'])
K=['JAM','CSKK','CJ']
ax.legend([H[_] for _ in K],[L[_] for _ in K],frameon=0,loc='best',bbox_to_anchor=[0.46,1.02],handletextpad=0.4,fontsize=19)


#--glue
ax=AX[1]
H,L={},{}
H['JAM']=plot_bandA(ax,'g','r',r'$\rm JAM19$',scale=1)
for _ in ['MMHT', 'CJ', 'ABMP', 'NNPDF','CSKK']: 
    if _=='CJ': label=r'$\rm CJ15$'
    elif _=='NNPDF': label=r'$\rm NNPDF3.1$'
    elif _=='MMHT': label=r'$\rm MMHT14$'
    elif _=='ABMP': label=r'$\rm ABMP16$'
    else: label=r'$\rm %s$'%_
    L[_]=label
    H[_]=plot_bandB(ax,groups[_],'g',color=groups[_]['color'],alpha=0.2,label=label,scale=1)
ax.text(0.8,0.8,r'\boldmath$xg$',color='k',transform=ax.transAxes,size=31)
ax.set_ylim(0,2.5)
ax.set_yticks([0,0.5,1,1.5,2])
ax.set_yticklabels([r'$0$',r'$0.5$',r'$1$',r'$1.5$',r'$2$'])
K=['ABMP','NNPDF','MMHT']
ax.legend([H[_] for _ in K],[L[_] for _ in K],frameon=0,loc='best',bbox_to_anchor=[0.52,0.6],handletextpad=0.4,fontsize=19)


#--db+ub
ax=AX[2]
H,L={},{}
H['JAM']=plot_bandA(ax,'db+ub','r',r'$\rm JAM$')
L['JAM']=r'$\rm JAM19$'
for _ in ['MMHT', 'CJ', 'ABMP', 'NNPDF','CSKK']: 
    if _=='CJ': label=r'$\rm CJ15$'
    elif _=='NNPDF': label=r'$\rm NNPDF3.1$'
    elif _=='MMHT': label=r'$\rm MMHT14$'
    elif _=='ABMP': label=r'$\rm ABMP16$'
    else: label=r'$\rm %s$'%_
    L[_]=label
    H[_]=plot_bandB(ax,groups[_],'db+ub',color=groups[_]['color'],alpha=0.2,label=r'$\rm %s$'%_)
ax.text(0.60,0.77,r'\boldmath$x(\bar{d}\!+\!\bar{u})$',color='k',transform=ax.transAxes,size=28)
ax.set_ylim(0,0.48)
ax.set_yticks([0,0.1,0.2,0.3,0.4])
ax.set_yticklabels([r'$0$',r'$0.1$',r'$0.2$',r'$0.3$',r'$0.4$'])


#--db-ub
ax=AX[3]
plot_bandA(ax,'db-ub','r',r'$\rm JAM$')
for _ in groups:
    plot_bandB(ax,groups[_],'db-ub',color=groups[_]['color'],alpha=0.2,label=r'$\rm %s$'%_)
ax.text(0.05,0.78,r'\boldmath$x(\bar{d}\!-\!\bar{u})$',color='k',transform=ax.transAxes,size=28)
ax.axhline(y=0.0,c='k',ls='-',lw=1,alpha=0.5)
ax.set_ylim(-0.07,0.10)
ax.set_yticks([-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.10])
ax.set_yticklabels([r'',r'-$0.04$',r'',r'$0$',r'',r'$0.04$',r'',r'$0.08$',r''])


#--s+sb
ax=AX[4]
plot_bandA(ax,'s+sb','r',r'$\rm JAM$')
for _ in groups:
    plot_bandB(ax,groups[_],'s+sb',color=groups[_]['color'],alpha=0.2,label=r'$\rm %s$'%_)
ax.text(0.62,0.77,r'\boldmath$x(s\!+\!\bar{s})$',color='k',transform=ax.transAxes,size=28)
ax.set_ylim(0,0.35)
ax.set_yticks([0,0.1,0.2,0.3])
ax.set_yticklabels([r'$0$',r'$0.1$',r'$0.2$',r'$0.3$'])


#--Rs
ax=AX[5]
plot_bandA(ax,'s+sb/db+ub','r',r'$\rm JAM$')
for _ in groups:
    plot_bandB(ax,groups[_],'Rs',color=groups[_]['color'],alpha=0.2,label=r'$\rm %s$'%_)
ax.text(0.5,0.6,r'\boldmath$R_s$',color='k',transform=ax.transAxes,size=30)
ax.set_ylim(0,1.3)
ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2])
ax.set_yticklabels([r'$0$',r'',r'',r'',r'$0.4$',r'',r'',r'',r'$0.8$',r'',r'',r'',r'$1.2$'])


for ax in AX:
    ax.semilogx()
    ax.tick_params(axis='both',which='both',right=True,top=True,direction='in',labelsize=20)
    ax.xaxis.set_label_coords(0.7, -0.05)
    ax.set_xlim(8e-3,9e-1)
    ax.set_xticks([1e-2,1e-1])
    if ax==AX[4] or ax==AX[5]:
        ax.set_xlabel(r'\boldmath$x$',size=30)   
        ax.set_xticklabels([r'$0.01$',r'$0.1$'])
        ax.text(0.84,-0.096,r'$0.5$',color='k',transform=ax.transAxes,size=20)
    else:
        ax.set_xticklabels([r'',r''])
    
py.tight_layout()
py.subplots_adjust(left=0.05, bottom=0.06, right=0.99, top=0.99, wspace=0.17, hspace=0.04)
checkdir('gallery')
py.savefig('gallery/pdf_sets.pdf')
