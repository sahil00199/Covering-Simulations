import matplotlib.pyplot as plt

name = "right"

file = open(name, 'r')
data = ""
for chunk in file:
	data += chunk

file.close()
data = data.split("\n")
data = [x for x in data if x != ""]

x = []
y = []
for element in data:
	x.append(float(element.split()[0]))
	y.append(float(element.split()[1]))

plt.plot(x, y)
plt.axis([0, 3.8, 0, 1.1])
plt.grid(True)
plt.savefig("../plots/"+name+".png")
