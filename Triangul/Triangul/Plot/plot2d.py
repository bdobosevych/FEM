import matplotlib.pyplot as plt
from decimal import Decimal
f = open("./../file/treangulation/vertices.txt", "r")
res_dist={}

for info in f.read().split('\n'):
    if info != "":
        r = info.split()
        res_dist[r[0]]=[float(r[1].replace(",", ".")),float(r[2].replace(",", "."))]
print(res_dist)

fig, ax = plt.subplots()
f1 = open("./../file/treangulation/triangles.txt", "r")
for info in f1.read().split('\n'):
    if info != "":
        r = info.split()
        x=[res_dist[r[1]][0],res_dist[r[2]][0],res_dist[r[3]][0],res_dist[r[1]][0]]
        y=[res_dist[r[1]][1],res_dist[r[2]][1],res_dist[r[3]][1],res_dist[r[1]][1]]
        ax.plot(x, y,color='black')
plt.show()