import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import spdiags
import time
import matplotlib.animation as animation

def laplacian(N):
	"""
	Create a sparse matrix to calculate the laplacian
	    N = grid size
	"""
	e = np.ones(N**2)
	e2 = ([1]*(N-1)+[0])*N
	e3 = ([0]+[1]*(N-1))*N
	L = spdiags([-4*e,e2,e3,e,e],[0,-1,1,-N,N],N**2,N**2)
	return L

def initialize(N):
	"""
	Setting up the initial condition
	    N = grid size
	"""
	N2, r = np.int64(N/2), 16

	# We start with every grid cell that has a lot of chemical A
	# and no chemical B except at the small region in the center
	A = np.ones((N, N))
	B = np.zeros((N, N))

	B[N2-r:N2+r, N2-r:N2+r] = 1.0

	A = A.reshape((N*N))
	B = B.reshape((N*N))

	return A, B

def update(A, B, dT, Da, Db, F, K, L):
	"""
	Update the iteration via Euler method
		A = Concentration of chemical A
		B = Concentration of chemical B
		dT = Time step
		Da = Diffusion coefficient of chemical A
		Db = Diffusion coefficient of chemical B
		F = Feed rate
		K = Kill rate
		L = Sparse matrix for Laplacian
	"""
	A = A + (Da*L.dot(A) - A*B*B + F*(1-A)) * dT
	B = B + (Db*L.dot(B) + A*B*B - (F+K)*B) * dT

	return A, B

def drawPattern(A, ax=None):
	"""
	Draw the pattern based on the chemical concentrations
		A = Concentration of chemical A
		ax = axis object
	"""
	ax.pcolor(A.reshape((N,N)), cmap=plt.cm.RdBu)
	ax.set_axis_off()
	# ax.axis('tight')

if __name__ == '__main__':

	do_pattern_growth_sim = True
	do_full_pattern_sim = False

	"""
	For pattern growth simulation
	"""
	if do_pattern_growth_sim:
		N = 256       # Grid size
		T = 5000     # Total time
		dT = 1.0      # Time step
		n = int(T/dT) # Number of iterations
		Da, Db, F, K = 0.14, 0.06, 0.035, 0.065  # Coefficients

		fig, axes = plt.subplots(3,3, figsize=(8,8), dpi=400, facecolor='w', edgecolor='k')
		step_plot = n // 9

		# Simulation part
		t1 = time.time()

		A, B = initialize(N)
		L = laplacian(N)

		for i in range(n):
			# Iteration update
			A, B = update(A, B, dT, Da, Db, F, K, L)

			# Plot at 9 different times
			if i % step_plot == 0 and i < 9 * step_plot:
				ax = axes.flat[i // step_plot]
				drawPattern(A, ax=ax)
				ax.set_title(f'iter = {i}')
		    
		fig.suptitle(f'Cell Division Pattern',fontsize=20)
		t2 = time.time()

		print(f'Simulation time = {t2-t1} seconds')
		plt.show()

	"""
	For full pattern simulation
	"""
	if do_full_pattern_sim:
		N = 256       # Grid size
		T = 3.2e4     # Total time
		dT = 1.0      # Time step
		n = int(T/dT) # Number of iterations
		Da, Db, F, K = 0.14, 0.06, 0.035, 0.065  # Coefficients for Cell Division Pattern
		# Da, Db, F, K = 0.16, 0.08, 0.060, 0.062  # Coefficients for Coral Pattern
		# Da, Db, F, K = 0.12, 0.08, 0.020, 0.050  # Coefficients for Spiral Pattern
		# Da, Db, F, K = 0.16, 0.08, 0.035, 0.060  # Coefficients for Middle East Pattern

		fig, ax = plt.subplots(1, dpi=400, figsize=(4,4), facecolor='w', edgecolor='k')

		# Simulation part
		t1 = time.time()

		A, B = initialize(N)
		L = laplacian(N)

		for i in range(n):
			# Iteration update
			A, B = update(A, B, dT, Da, Db, F, K, L)

		t2 = time.time()

		print(f'Simulation time = {t2-t1} seconds')

		# Plot at the last iteration
		drawPattern(A, ax=ax)
		fig.suptitle(f'Cell Division Pattern',fontsize=15);
		fig.savefig('../../Figures/Cell_Division_Pattern.png',format='png',dpi=300)
		plt.show()