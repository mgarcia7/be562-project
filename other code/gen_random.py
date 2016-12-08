import random

def main():
	b = "AGTC"
	g = []

	stop = random.randint(0, (5000 - 20))

	for i in range(stop):
		g.append(b[random.randint(0, 3)])

	g.append("AAAAAAAAAAAAAAAAAAAA")

	for i in range((5000 - 20) - stop):
		g.append(b[random.randint(0, 3)])

	out = "".join(g)

	print(out)

if __name__ == "__main__":
	main()