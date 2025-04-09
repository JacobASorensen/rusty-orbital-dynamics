import rusty_orbital_dynamics as rod
import matplotlib.pyplot as plt
import numpy as np 

def Basic_Graphic(data,size):
    y_values = np.array(data[1])
    # Static Plot
    # plt.clf()
    plt.figure( figsize=(5,5) )

    current = int(len(y_values[0])/2)

    x_vals = []
    y_vals = []
    for x in range(size):
        x_val = y_values[:,current]
        x_vals.append(x_val)
        y_val = y_values[:,current+1]
        y_vals.append(y_val)
        plt.plot( x_val , y_val , label = 'orbit #' + str(x+1) )
        current += 3

    plt.legend()
    plt.grid()
    plt.axis('equal')
    plt.show()
    plt.savefig('test1.png')

def Formatter(vector):
    VN = []
    RN = []
    for i in range(len(vector)):
        if (i % 2 == 0):
            RN += vector[i]
        else:
            VN += vector[i]
    return VN + RN

# Mass Data
masses = [10**29,10**29,10**29,10**25,10**25]

# Positional and Velocity data; Of format Rx, Vx
vectors = [
    [-300000,0,0],
    [0,0,0],
    [0,0,0],
    [250,250,0],
    [300000,0,0],
    [0,0,0],
    [-100000,-300000,0],
    [300,300,0],
    [0,6000000,0],
    [50,0,0]
]

data = rod.rkf45(Formatter(vectors),masses,0,150000,1,6.67259*10**-20,10**-10,10,0.1)

# Generates a static plot of the data
Basic_Graphic(data,len(masses))