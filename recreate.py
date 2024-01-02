import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

n = 10 #number of iterations
b = 0.0245
c = 0.4
a = np.linspace(0.1, 2, n)

x_e = np.sqrt(2*a*b*(2*a*c + 3 - np.sqrt(8*a*c + 9))) / (2*a)
y_e = np.sqrt(6)*(2*a*c + 3 - np.sqrt(8*a*c + 9)) / (np.sqrt(8*a*c + 9) - 3)
z_e = 2*np.sqrt(3*a*b*(2*a*c + 3 - np.sqrt(8*a*c + 9))) / (np.sqrt(8*a*c + 9) - 3)

a0 = c - b + a
a1 = (-1/3)*y_e**2 + (1/6)*y_e*np.sqrt(6) + 1 + x_e**2 + c*a - c*b + (1/3)*z_e**2 - b*a
a2 = (2/3)*y_e*z_e*x_e - (1/6)*np.sqrt(6)*x_e*z_e + (1/3)*b*y_e**2 - (1/6)*y_e*b*np.sqrt(6) - b + a*x_e**2 + (1/3)*c*z_e**2 - c*b
ones = np.ones(n)

eig = []
for i in range(n):
    eig.append(np.roots([ones[i], a2[i], a1[i], a0[i]]))
eig = np.array(eig)

real = []
imag = []
for i in range(n):
    real.append(eig[i].real)
    imag.append(eig[i].imag)
real = np.array(real)
imag = np.array(imag)


real1 = real[:, 2]
imaginary1 = imag[:, 2]

real2 = real[:, 1]
imaginary2 = np.flip(imag[:, 1])

real3 = real[:, 0]
imaginary3 = np.flip(imag[:, 0])
print(eig)

plt.plot(real1, imaginary1)
plt.plot(real2, imaginary2)
plt.plot(real3, imaginary3)
plt.grid()
plt.ylabel("Imag(eignevalue)")
plt.xlabel("Real(eignevalue)")
plt.show()





