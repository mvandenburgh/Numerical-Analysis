# from rp import *
import numpy as np
import sys
sys.path.append("..") # Adds higher directory to python modules path.
from t32 import uniform_random as random_float

def random_floats(n, a, b):
    return [random_float(a,b) for i in range(n)]

τ=2*np.pi           #Tau is better than pi. Fight me lol xD
def simulated_probability(needle_length):
	n=500               #Number of points on the quadrifolium curve's approximation
	θ=np.linspace(0,τ,n)#[0τ÷n, 1τ÷n, 2τ÷n, ... nτ÷n]
	r=np.sin(2*θ)
	x=r*np.cos(θ)
	y=r*np.sin(θ)
	quadrifolium=np.transpose([x,y])#A list of n x,y pairs

	# display_clear()
	# display_path(quadrifolium                         ,color='black')#Draw the quadrifolium
	# display_path([[-1,-1],[-1,1],[1,1],[1,-1],[-1,-1]],color='black')#Draw the bounding box

	class Needle:
		def __init__(self):
			self.angle =random_float (   0,τ)#A random angle
			self.center=random_floats(2,-1,1)#2 random floats between -1 and 1

			#Calculate start and end point coordinates for this needle
			#In my rp library, we call lists of 2d points "paths"
			r  =needle_length/2
			Δx =r * np.cos(self.angle)
			Δy =r * np.sin(self.angle)
			self.start =self.center+[Δx,Δy]
			self.end   =self.center-[Δx,Δy]
			self.path=[self.start,self.end]


		def intersects_quadrifolium(self)->bool:
			#Return True if this needle intersects the quadrifolium, False otherwise
			number_of_tests=100
			x=np.linspace(self.start[0],self.end[0],number_of_tests)
			y=np.linspace(self.start[1],self.end[1],number_of_tests)
			θ=np.arctan2(x,y)
			r=np.sqrt(x**2+y**2)
			in_quadrifolium=r<np.abs(np.sin(2*θ))
			return np.any(in_quadrifolium) and not np.all(in_quadrifolium)#Some points are inside, and others are outside

		def plot(self):

			if self.intersects_quadrifolium():
				color = 'red'
			else:
				color = 'blue'

			# display_path(self.path,color=color,block=None)

	needles=[]
	N=60000#Number of needles
	for _ in range(N):
		needle=Needle()
		needle.intersects_quadrifolium()#Do that calculation now, so we can see how long its taking...
		needles.append(needle)
		# Needle().plot() #Drawing too many needles can cause python to crash
		# if not _%500:
			# update_display()
			# print("Simulated",_,"of",N,"needles")
	# update_display()

	probability=np.mean(list(needle.intersects_quadrifolium() for needle in needles))
	return probability

for length in [.1,.25,.5,1]:
	print("For length",length,"probability is",simulated_probability(length))