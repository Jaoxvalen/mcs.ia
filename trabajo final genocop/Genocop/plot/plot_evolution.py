import numpy as np
import matplotlib.pyplot as plt

from benchmarks import *

plt.ion()

#_bmFunction = BMAckley(20,0.2, 2*np.pi ,2, BenchmarkFunction.PLOTTING_MODE_SURFACE )
#_bmFunction.setRange(0.0,25.0)
#file = open("t01_2_plot.txt", "r") 

#_bmFunction = BMSchwefel(2, BenchmarkFunction.PLOTTING_MODE_SURFACE )
#_bmFunction.setRange(0.0,2000.0)
#file = open("t02_2_plot.txt", "r") 

_bmFunction = BMShafferFcn6( 2, BenchmarkFunction.PLOTTING_MODE_SURFACE )
_bmFunction.setRange(0,1)
file = open("t03_2_plot.txt", "r") 


#read file and fill the matrix
_buff = []

_buff_gen = []

_optimal = []

for line in file: 
	spl = line.split()
	
	if spl[0] == "pop_generation":
		if len(_buff_gen)>0:
			_buff.append(_buff_gen)
			_buff_gen = [] 
		
	if spl[0] == "values":
		_buff_xy = []
		_buff_xy.append( float(spl[2]) )
		_buff_xy.append( float(spl[3]) )
		_buff_gen.append(_buff_xy)

	if spl[0] == "optimal":
		_buff_xy_o = []
		_buff_xy_o.append( float(spl[1]) )
		_buff_xy_o.append( float(spl[2]) )
		_optimal.append(_buff_xy_o)

_buff.append(_buff_gen)


_buff = np.asarray(_buff)
_optimal = np.asarray(_optimal)



for i in range(0,100):
	 _bmFunction.axes().cla()
	 _bmFunction.axesContour().cla()
	 _bmFunction.plotBase()
	 _bmFunction.plotSequence( np.reshape( _buff[i], ( len( _buff[i] ), _bmFunction.ndim() ) ) )
	 _bmFunction.plotSingle(_optimal, 'bo' , 'bo')

	 plt.figure(0);
	 plt.savefig('plot_t03/gen/'+str(i)+'.png')
	 plt.figure(1);
	 plt.savefig('plot_t03/con/'+str(i)+'.png')
	 plt.pause( 0.00001 )


'''
while True :

    _bmFunction.axes().cla()

    _bmFunction.plotBase()

    _p = np.random.uniform( _bmFunction.min(), _bmFunction.max(), ( 1, _bmFunction.ndim() ) )
    _buff.append( _p )

    _bmFunction.plotSequence( np.reshape( _buff, ( len( _buff ), _bmFunction.ndim() ) ) )

    plt.pause( 0.02 )
'''