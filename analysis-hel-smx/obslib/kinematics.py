






ifig = py.figure(figsize=(6,4))
ax1=py.subplot(121)
ax2=py.subplot(122)

data=load('%s/%s'%('data/cov-exp/calc','data_and_thy.dat'))
NA10=data['dy-pion']['tabs']['NA10']
E615=data['dy-pion']['tabs']['E615']
H1=data['ln']['tabs']['H1']
ZEUS=data['ln']['tabs']['ZEUS']

ax1.plot(E615['xF'],E615['Q2'],'bx',label=r'$\rm E615$')
ax1.plot(NA10['xF'],NA10['Q2'],'rx',label=r'$\rm NA10$')
ax1.plot(H1['xpi'],H1['Q2'],'gx',label=r'$\rm H1$')
ax1.plot(ZEUS['xpi'],ZEUS['Q2'],'mx',label=r'$\rm ZEUS$')
ax1.semilogy()
ax1.semilogx()
ax1.set_xlim(None,0.1)
ax1.spines['right'].set_visible(False)
ax1.legend(loc=2,fontsize=15,frameon=0)

ax2.plot(E615['xF'],E615['Q2'],'bx',label=r'$\rm E615$')
ax2.plot(NA10['xF'],NA10['Q2'],'rx',label=r'$\rm NA10$')
ax2.plot(H1['xpi'],H1['Q2'],'gx',label=r'$\rm H1$')
ax2.plot(ZEUS['xpi'],ZEUS['Q2'],'mx',label=r'$\rm ZEUS$')
ax2.semilogy()
ax2.set_xlim(0.1,1)
ax2.set_yticks([])
ax2.set_yticklabels([])
ax2.tick_params(axis='y',which='both',left='off')
ax2.spines['left'].set_visible(False)
py.subplots_adjust(wspace=0, hspace=0)

ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.set_yticks([10,100])
ax1.set_yticklabels([r'$10$',r'$100$'])
ax1.yaxis.set_label_coords(-0.14,0.85)
ax1.set_ylabel(r'$Q^2$',size=30,rotation=0)
ax1.set_xticks([0.001,0.01,0.1])
ax1.set_xticklabels([r'$0.001$',r'$0.01$',r'$0.1$'])


ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.xaxis.set_label_coords(0.95, -0.03)
ax2.set_xticks([0.3,0.5,0.7])
ax2.set_xticklabels([r'$0.3$',r'$0.5$',r'$0.7$',r'$0.9$'])
ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.set_xlabel(r'$x_{\pi}$',size=30)

