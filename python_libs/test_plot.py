import numpy as np 
import matplotlib.pyplot as plt 
import plot_settings as ps

fig, ax = plt.subplots()

C = ps.CURVEFAMILY( 10, ax )
C.set_individual_colors( color_map = "cmptugreen" )
C.set_individual_markerstyles()
C.set_individual_markeveries()
C.set_shared_kwargs( markersize=10 )

C_fit = ps.CURVEFAMILY( 10, ax )
C_fit.set_individual_colors( color_map = "cmptugreen" )
C_fit.set_shared_linestyle("--")
#C_fit.set_individual_markerstyles()
C_fit.set_individual_markeveries()

x = np.linspace(0,10,100)
for i in range(0,10):
    d = np.sin( x + i*0.1 )
    C.plot( x, d, label = str(i) )

x = np.linspace(0,10,100)
for i in range(0,10):
    d = np.cos( x + i*0.1 )
    C_fit.plot( x, d, label = str(i) )

plt.ylim(-2,2)
#plt.legend(loc='best')
plt.savefig('test.pdf')