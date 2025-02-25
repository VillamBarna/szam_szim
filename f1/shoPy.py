import matplotlib.pyplot as plt
import time
import numpy as np
import sho

plt.rcParams.update({'font.size':14})
def feladat1():
	sho.run_simulation(3, 0, 1, 6, 100, "keteres_ido_ec.txt", True)
	sho.run_simulation(3, 0, 1, 6, 100, "keteres_ido_e.txt", False)
	data = np.loadtxt( "keteres_ido_ec.txt").T
	data2 = np.loadtxt( "keteres_ido_e.txt").T
	plt.plot(data[0], data[1], label="Euler-Cromer")
	plt.plot(data2[0], data2[1], "--", label="Euler")
	plt.xlabel("t[s]")
	plt.ylabel("x[m]")
	plt.legend()
	plt.grid()
	plt.tight_layout()
	plt.savefig("kiteres_ido.png")

def feladat2(ec):
	file_name = "kiteres_sebesseg_ec.txt" if ec else "kiteres_sebesseg_e.txt"
	sho.run_simulation(.5, 0, 1, 1000, 100, file_name, ec)
	data = np.loadtxt(file_name).T
	fig, ax = plt.subplots()
	ax.plot(data[1], data[2])
	plt.grid()
	plt.xlabel("x[m]")
	plt.ylabel("y[m]")
	plt.tight_layout()
	plt.savefig(file_name[:-3]+"png")

def feladat3(ec):
	file_name = "energia_ec.txt" if ec else "energia_e.txt"
	sho.run_simulation(.5, 0, 1, 20, 100, file_name, ec)
	data = np.loadtxt(file_name).T
	fig, ax = plt.subplots()
	ax.semilogy(data[0], 0.5*(.5**2*np.power(data[1], 2)+np.power(data[2], 2)))
	plt.xlabel("t[s]")
	plt.ylabel("E[J]")
	plt.tight_layout()
	plt.savefig(file_name[:-3]+"png")

def feladat4():
	file_name = "tmp.txt"
	number_of_periods_list = np.array([x for x in range(0, 10000, 100)])
	running_time = []
	print("starting")
	for number_of_periods in number_of_periods_list:
		print("Progress:", number_of_periods/len(number_of_periods_list), "%", end="\r")
		start_time = time.time()
		sho.run_simulation(1, 0, 1, number_of_periods, 100, file_name, True)
		running_time.append(time.time()-start_time)
	plt.plot(number_of_periods_list*100, running_time)
	plt.xlabel("Kiértékelések száma[db]")
	plt.ylabel("Futási idő[s]")
	plt.grid()
	plt.tight_layout()
	plt.savefig("futasi_ido.png")

feladat4()
