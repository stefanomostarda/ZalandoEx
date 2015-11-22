#!/bin/bash
# -*- coding: utf-8 -*-

from numpy import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

## spree river in degrees
spree = [[52.529198,13.274099],[52.531835,13.29234],[52.522116,13.298541],[52.520569,13.317349],[52.524877,13.322434],[52.522788,13.329],[52.517056,13.332075],[52.522514,13.340743],[52.517239,13.356665],[52.523063,13.372158],[52.519198,13.379453],[52.522462,13.392328],[52.520921,13.399703],[52.515333,13.406054],[52.514863,13.416354],[52.506034,13.435923],[52.496473,13.461587],[52.487641,13.483216],[52.488739,13.491456],[52.464011,13.503386]]

S	= [52.464011, 13.274099]	## reference point
B	= [52.516288, 13.377689]	## brandeburg tor
Sa  = [[52.590117, 13.39915], [52.437385, 13.553989]] ## satellite start and arrival

## parameters distributions
Mean = 4.700
Mode = 3.877
s1 = 0.6825											## sigma normal distr from spree
s2	= sqrt(2.0/3.0*(log(Mean)-log(Mode)))	## sigma log-normal distr from brand-tor
m2	= (2*log(Mean)+log(Mode))/3.0				## mu 	log-normal distr from brand-tor
s3	= 0.6000											## sigma normal distr from satellite

## transform degrees to Km
def degrees_to_km( P, S ):
	P_x = (P[0] - S[0]) * cos(S[1] * pi / 180.0) * 111.323
	P_y = (P[1] - S[1]) * 111.323
	return P_x, P_y

## transform Km to degrees
def km_to_degrees( P, S ):
	P_x = S[0] + P[0]/( cos(S[1] * pi / 180.0) * 111.323)
	P_y = S[1] + P[1]/111.323
	return P_x, P_y

## distance between two points
def dist_p1_p2( p1, p2):
	return math.hypot( p2[0]-p1[0], p2[1]-p1[1] )

## gaussian distribution (didn't use the built-on one to check functionality)
def prob_gauss( x, s ):
	return 1.0/(s * sqrt(2*pi)) * exp (-(x**2)/(2*s**2))

## log-normal distribution (didn't use the built-on one to check functionality)
def prob_log_norm( x, s, m):
	return 1.0/( x*s*sqrt(2*pi)) * exp (-(log(x)-m)**2/(2*s**2))

## minimum distance from Spree river:
## This method is EXACT. Given a point in space X, it takes one segment of Spree s_i at the time,
## it calculates the projection of X on s_i and checks if it belongs to the segment.
## If it belongs, that is the local minimum; if it doesn't it takes the end of s_i which is closest to X.
## The global minimum will be the minimum after all Spree segments are taken into account.
def min_dist_from_spree( X, spree ):
	dist_spree = []                    ## record the distances from each segment
	for ii3 in range(len(spree)-1):
		spree_v = [spree[ii3+1][0]-spree[ii3][0], spree[ii3+1][1]-spree[ii3][1]]      ## spree vector i to i+1
		ii_P_vec = [ X[0]-spree[ii3][0], X[1]-spree[ii3][1]]                          ## vector from i to X
		scale = dot( ii_P_vec, spree_v) / dot(spree_v, spree_v)                       ## projection of ii_P_vec on spree_v
		proj_x_on_line = [ spree[ii3][0] + scale * spree_v[0] , spree[ii3][1] + scale * spree_v[1] ]  ## projection of X on spree_v 
		## is the projection belonging to the segment? If yes, take that one, otherwise take the closest vertix
		if dist_p1_p2( proj_x_on_line, spree[ii3]) + dist_p1_p2( proj_x_on_line, spree[ii3+1]) > 1.001*dist_p1_p2( spree[ii3], spree[ii3+1]):
			dist_spree.append( min(dist_p1_p2(X, spree[ii3]), dist_p1_p2(X, spree[ii3+1])) )
		else:
			dist_spree.append( dist_p1_p2(X,proj_x_on_line) )
	return min(dist_spree)

## transform degrees in km
spree   = [ degrees_to_km( ii, S ) for ii in spree  ]
B 		=   degrees_to_km( B, S  )									## becomes 5.66, 11.53
Sa		= [ degrees_to_km( ii, S ) for ii in Sa     ]

precision = 0.1    ## in Km

######## Probabilities Calculation are a grid of 30.0K x 30K

## probability from spree
p1 = array([[prob_gauss( min_dist_from_spree( [ii*precision , ii2*precision], spree ), s1  ) for ii2 in range(int(30.0/precision))] for ii in range(int(30.0/precision))])       
## probability from brandeburg tor
p2 = array([[prob_log_norm( dist_p1_p2( [ii*precision , ii2*precision], B ), s2, m2  ) for ii2 in range(int(30.0/precision))] for ii in range(int(30.0/precision))])
## probability from staellite
p3 = array([[prob_gauss( min_dist_from_spree( [ii*precision , ii2*precision], Sa ), s3 ) for ii2 in range(int(30.0/precision))] for ii in range(int(30.0/precision))])
## the total probability is the product of the 3
pt = array([[p1[ii][ii2]*p2[ii][ii2]*p3[ii][ii2] for ii2 in range(int(30.0/precision))] for ii in range(int(30.0/precision)) ]) ## max in 52.492622310971008, 13.49238273292132




## Plot results
dx, dy = precision, precision
#y, x = mgrid[slice(dy, 30 + dy, dy), slice(dx, 30 + dx, dx)]
x, y = mgrid[slice(dy, 30 + dy, dy), slice(dx, 30 + dx, dx)]

def preparation_plot( data, cmap ):
    data    = data[:-1, :-1]           ## remove the last row/column
    levels  = MaxNLocator(nbins=15).bin_boundaries(data.min(), data.max())
    norm    = BoundaryNorm(levels , ncolors=cmap.N, clip=True)
    return data, levels, norm
    
fig, ((ax1, ax2), (ax3, ax4))= plt.subplots(nrows = 2, ncols = 2)
cmap 		                 = plt.get_cmap('bwr')

pt, levelspt, normpt = preparation_plot( pt , cmap )
p1, levelsp1, normp1 = preparation_plot( p1 , cmap )
p2, levelsp2, normp2 = preparation_plot( p2 , cmap )
p3, levelsp3, normp3 = preparation_plot( p3 , cmap )



cf = ax1.contourf(y[:-1, :-1] + dy/2., x[:-1, :-1] + dx/2., pt, levels=levelspt,cmap=cmap)
fig.colorbar(cf, ax=ax1)
#ax1.set_xlim([0,8])
ax1.set_title('Total probability')
ax1.set_ylabel('Km')
ax1.set_xlabel('Km')
#ax1.set_ylim([20,28])

cf = ax2.contourf(y[:-1, :-1] + dy/2., x[:-1, :-1] + dx/2., p1, levels=levelsp1,cmap=cmap)
fig.colorbar(cf, ax=ax2)
#ax2.set_xlim([0,8])
ax2.set_title('Probability from Spree')
ax2.set_ylabel('Km')
ax2.set_xlabel('Km')
#ax2.set_ylim([20,28])

cf = ax3.contourf(y[:-1, :-1] + dy/2., x[:-1, :-1] + dx/2., p2, levels=levelsp2,cmap=cmap)
fig.colorbar(cf, ax=ax3)
#ax3.set_xlim([0,8])
ax3.set_title('Probability from Brandeburg Tor')
ax3.set_ylabel('Km')
ax3.set_xlabel('Km')
#ax3.set_ylim([20,28])

cf = ax4.contourf(y[:-1, :-1] + dy/2., x[:-1, :-1] + dx/2., p3, levels=levelsp3, cmap=cmap)
fig.colorbar(cf, ax=ax4)
#ax4.set_xlim([0,8])
ax4.set_title('Probability from satellite')
ax4.set_ylabel('Km')
ax4.set_xlabel('Km')
#ax4.set_ylim([20,28])

plt.show()
