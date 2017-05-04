import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig = plt.figure()

ax = fig.add_axes([0.1, 0.1, .9, .9,])

c = [(0,0,0), (0,0,0), (0,0,0)]
x = [.1,.2,.3]
y = [.1,.1,.1]
p = ax.scatter(x,y,facecolors=c)

def update(i):
    index = i % 3
    if c[index][0] == 0: 
        c[index] = (1,.5,1)
    else:
        c[index] = (0,0,0)
    
    p.set_facecolors(c)

an = animation.FuncAnimation(fig, update, interval=100)
plt.show()


plt.show()
