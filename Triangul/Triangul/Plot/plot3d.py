import matplotlib.pyplot as plt
from decimal import Decimal
f = open("./../file/treangulation/vertices.txt", "r")
res_dist={}

for info in f.read().split('\n'):
    if info != "":
        r = info.split()
        res_dist[r[0]]=[float(r[1].replace(",", ".")),float(r[2].replace(",", "."))]
print(res_dist)
f1 = open("./../file/CLAR/Un.txt", "r")
for info in f1.read().split('\n'):
    if info != "":
        r = info.split()
        for k in res_dist:
            if res_dist[k]==[float(r[0].replace(",", ".")),float(r[1].replace(",", "."))]:
                res_dist[k].append(float(r[2].replace(",", ".")))
print(res_dist)

fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(111, projection='3d')
f1 = open("./../file/treangulation/triangles.txt", "r")
for info in f1.read().split('\n'):
    if info != "":
        r = info.split()
        x=[res_dist[r[1]][0],res_dist[r[2]][0],res_dist[r[3]][0],res_dist[r[1]][0]]
        y=[res_dist[r[1]][1],res_dist[r[2]][1],res_dist[r[3]][1],res_dist[r[1]][1]]
        z=[res_dist[r[1]][2],res_dist[r[2]][2],res_dist[r[3]][2],res_dist[r[1]][2]]
        ax.plot(x, y, z,color='black')

plt.show()