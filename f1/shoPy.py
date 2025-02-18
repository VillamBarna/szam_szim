import matplotlib.pyplot as plt
import time
import numpy as np
import sho

def feladat1():
	file_name = "keteres-ido.txt"
	sho.run_simulation(1, 0, 1, 3, 100, file_name, True)
	data = np.loadtxt(file_name).T
	plt.plot(data[0], data[1])
	plt.show()

def feladat2():
	file_name = "kieters-sebesseg.txt"
	sho.run_simulation(.5, 0, 1, 100, 100, file_name, True)
	data = np.loadtxt(file_name).T
	fig, ax = plt.subplots()
	ax.plot(data[1], data[2], '.-')
	ax.set_aspect('equal')
	plt.show()

def feladat3_ec():
	file_name = "energia-ec.txt"
	sho.run_simulation(.5, 0, 1, 100, 100, file_name, True)
	data = np.loadtxt(file_name).T
	fig, ax = plt.subplots()
	ax.plot(data[0], 0.5*(.5**2*np.power(data[1], 2)+np.power(data[2], 2)))
	plt.show()

def feladat3_e():
	file_name = "energia-e.txt"
	sho.run_simulation(.5, 0, 1, 100, 100, file_name, False)
	data = np.loadtxt(file_name).T
	fig, ax = plt.subplots()
	ax.semilogy(data[0], 0.5*(.5**2*np.power(data[1], 2)+np.power(data[2], 2)))
	plt.show()

def feladat4():
	file_name = "tmp.txt"
	number_of_periods_list = np.array([x for x in range(1000)])
	running_time = []
	print("starting")
	for number_of_periods in number_of_periods_list:
		print("Progress:", number_of_periods/len(number_of_periods_list)*100, "%", end="\r")
		start_time = time.time()
		sho.run_simulation(1, 0, 1, number_of_periods, 100, file_name, True)
		running_time.append(time.time()-start_time)
	plt.plot(number_of_periods_list*100, running_time)
	plt.show()
