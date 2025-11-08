import matplotlib.pyplot as plt

k_values = [1, 5, 10, 15, 30, 50, 75, 100, 150, 175, 180]
norms = [
    39045.483388, 
    20501.518212, 
    14919.064448, 
    12218.031020, 
    8456.577677, 
    6124.659909, 
    4567.945162, 
    3630.536462, 
    2469.341613, 
    2085.406675, 
    2018.163274
]

plt.figure(figsize=(10, 6))
plt.plot(k_values, norms, marker='o', linestyle='-', color='b')
plt.title('Image Reconstruction Error (Frobenius Norm) vs. k')
plt.xlabel('k (Number of Singular Values Used)')
plt.ylabel('Frobenius Norm of Error (||A - A_k||_F)')

plt.grid(True)

plt.show()
