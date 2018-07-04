import matplotlib
import matplotlib.pyplot as plt 

file = open("input", 'r')
data = ""
for chunk in file:
	data += chunk
file.close()
data = data.split("\n")
data = [x for x in data if x != ""]

ax = plt.subplot(aspect = 'equal')
for element in data:
	radius = float(element.split()[0])
	coods = element.split()[1][1:-1].split(",")
	x = float(coods[0])
	y = float(coods[1])
	circle = plt.Circle((x, y), radius, fill = False)
	ax.add_artist(circle)


plt.xlim((-0.5, 3.5))
plt.ylim((-0.5, 3.5))

plt.grid(True)

plt.show()
