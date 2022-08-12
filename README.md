Batik Pattern Inspiration from the Gray-Scott Model
==============================================================================

[Batik](https://en.wikipedia.org/wiki/Batik) is a technique of wax-resist dyeing applied to the whole cloth. This technique originated from the island of Java, Indonesia. Batik is made either by drawing dots and lines of the resist with a spouted tool called canting, or by printing the resist with a copper stamp called a cap.

In this project, I aim to introduce a new way to create batik patterns from the famous
mathematical model called the Gray-Scott (GS) model. The idea is to numerically solve the GS
equations and track the pattern as the iteration increases. By varying the equation coefficients, some interesting patterns that are proper for batik motifs can be found at certain iterations. Batik designers can use the selected pattern as a building block that can be repeated to form a batik motif. They can also play around with the parameters to discover new patterns.

This code is adapted from the Reaction-Diffusion Tutorial by Karl Sims [1]. Instead of performing a 3x3 convolution operation to calculate the Laplacian (as in the tutorial), a 5-point stencil approximation is used in this code.

# Dependencies
NumPy, SciPy, Matplotlib

# Usage
********************************************************************************
```python
		N = 256       # Grid size
		T = 5000     # Total time
		dT = 1.0      # Time step
		n = int(T/dT) # Number of iterations
		Da, Db, F, K = 0.14, 0.06, 0.035, 0.065  # Coefficients for Cell Division Pattern
    # Da, Db, F, K = 0.16, 0.08, 0.060, 0.062  # Coefficients for Coral Pattern
		# Da, Db, F, K = 0.12, 0.08, 0.020, 0.050  # Coefficients for Spiral Pattern
		# Da, Db, F, K = 0.16, 0.08, 0.035, 0.060  # Coefficients for Middle East Pattern

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
```
# Results

### Original Image
<p float="center">
  <img src="pattern_growth.GIF" width="100%"/>
</p>

# Reference
[1] Karl Sims. Reaction-Diffusion Tutorial. (2013). [Link](http://www.karlsims.com/rd.html)
