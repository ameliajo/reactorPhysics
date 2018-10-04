import matplotlib.pyplot as plt

def plotRay(r,color):
    x = [r.x0,r.xMax]
    y = [r.y0,r.yMax]
    plt.plot(r.x0,r.y0,color,marker='bo')
    plt.plot(x,y,color)

def plotRaySegment(r,color,firstIntersection):
    x = [r.x0,firstIntersection["x-int"]]
    y = [r.y0,firstIntersection["y-int"]]
    #plt.plot(r.x0,r.y0,color,marker='o')
    plt.plot(x,y,color)


def plotBox(ax,box):
    plt.plot([box.L.x,box.L.x,box.R.x,box.R.x,box.L.x],
             [box.U.y,box.D.y,box.D.y,box.U.y,box.U.y])
    circles = box.C
    for circle in circles:
        c = plt.Circle((circle.x, circle.y), circle.r,fill=False)
        ax.add_artist(c)



