import pendulum
import time
import itertools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import multiprocessing

plt.rcParams.update({'font.size': 14})

def energy(L, theta, omega):
	return 1/2*L**2*omega**2 + 9.8*(1-np.cos(theta))

def matematikai():
    L = 1
    q = 0
    Omega_D = 0
    F_D = 0
    nonlinear = False
    theta = 3/180*np.pi
    omega = 0
    tmax = 40
    pendulum.run_pendulum("RK45_adaptive", nonlinear, L, q,
                            Omega_D, F_D, theta, omega, tmax)

    data = np.loadtxt('pendulum.data').T
    fig, axes = plt.subplots(4, 1, figsize=(10, 12), gridspec_kw={
                             'height_ratios': [1, 1, 1, 3]})
    axes[0].plot(data[0], data[1])
    axes[0].set_xlabel("t[s]")
    axes[0].set_ylabel("$\\theta[rad]$")
    axes[0].grid()

    axes[1].plot(data[0], data[2])
    axes[1].set_xlabel("t[s]")
    axes[1].set_ylabel("$\\omega[\\frac{1}{s}]$")
    axes[1].grid()

    axes[2].plot(data[0], energy(L, data[1], data[2]))
    axes[2].set_xlabel("t[s]")
    axes[2].set_ylabel("E[J]")
    axes[2].grid()

    axes[3].plot(data[1], data[2])
    axes[3].set_xlabel("$\\theta[rad]$")
    axes[3].set_ylabel("$\\omega[\\frac{1}{s}]$")
    axes[3].grid()
    plt.tight_layout()
    plt.show()


def csillapitott():
    L = 1
    q = 0.5
    Omega_D = 0
    F_D = 0
    nonlinear = False
    theta = 3/180*np.pi
    omega = 0
    tmax = 20
    pendulum.run_pendulum("RK45_adaptive", nonlinear, L, q,
                            Omega_D, F_D, theta, omega, tmax)

    data = np.loadtxt('pendulum.data').T
    fig, axes = plt.subplots(4, 1, figsize=(10, 12), gridspec_kw={
                             'height_ratios': [1, 1, 1, 3]})
    axes[0].plot(data[0], data[1])
    axes[0].set_xlabel("t[s]")
    axes[0].set_ylabel("$\\theta[rad]$")
    axes[0].grid()

    axes[1].plot(data[0], data[2])
    axes[1].set_xlabel("t[s]")
    axes[1].set_ylabel("$\\omega[\\frac{1}{s}]$")
    axes[1].grid()

    axes[2].plot(data[0], energy(L, data[1], data[2]))
    axes[2].set_xlabel("t[s]")
    axes[2].set_ylabel("E[J]")
    axes[2].grid()

    axes[3].plot(data[1], data[2])
    axes[3].set_xlabel("$\\theta[rad]$")
    axes[3].set_ylabel("$\\omega[\\frac{1}{s}]$")
    axes[3].grid()
    plt.tight_layout()
    plt.show()


def gerjesztett():
    L = 1
    q = .5
    Omega_D = 1
    F_D = 0.1
    nonlinear = False
    theta = 0
    omega = 6/180*np.pi
    tmax = 40
    pendulum.run_pendulum("RK45_adaptive", nonlinear, L, q,
                            Omega_D, F_D, theta, omega, tmax)

    data = np.loadtxt('pendulum.data').T
    print(data)
    fig, axes = plt.subplots(4, 1, figsize=(10, 12), gridspec_kw={
                             'height_ratios': [1, 1, 1, 3]})
    axes[0].plot(data[0], data[1])
    axes[0].set_xlabel("t[s]")
    axes[0].set_ylabel("$\\theta[rad]$")
    axes[0].grid()

    axes[1].plot(data[0], data[2])
    axes[1].set_xlabel("t[s]")
    axes[1].set_ylabel("$\\omega[\\frac{1}{s}]$")
    axes[1].grid()

    axes[2].plot(data[0], energy(L, data[1], data[2]))
    axes[2].set_xlabel("t[s]")
    axes[2].set_ylabel("E[J]")
    axes[2].grid()

    axes[3].plot(data[1], data[2])
    axes[3].set_xlabel("$\\theta[rad]$")
    axes[3].set_ylabel("$\\omega[\\frac{1}{s}]$")
    axes[3].grid()
    print("gerj")
    plt.tight_layout()
    plt.show()


def fizikai():
    L = 1
    q = 0
    Omega_D = 0
    F_D = 0
    nonlinear = True
    theta = 100/180*np.pi
    omega = 0
    tmax = 40
    pendulum.run_pendulum("RK45_adaptive", nonlinear, L, q,
                            Omega_D, F_D, theta, omega, tmax)

    data = np.loadtxt('pendulum.data').T
    fig, axes = plt.subplots(4, 1, figsize=(10, 12), gridspec_kw={
                             'height_ratios': [1, 1, 1, 3]})
    axes[0].plot(data[0], data[1])
    axes[0].set_xlabel("t[s]")
    axes[0].set_ylabel("$\\theta[rad]$")
    axes[0].grid()

    axes[1].plot(data[0], data[2])

    axes[1].set_xlabel("t[s]")
    axes[1].set_ylabel("$\\omega[\\frac{1}{s}]$")
    axes[1].grid()

    axes[2].plot(data[0], energy(L, data[1], data[2]))
    axes[2].set_xlabel("t[s]")
    axes[2].set_ylabel("E[J]")
    axes[2].grid()

    axes[3].plot(data[1], data[2])
    axes[3].set_xlabel("$\\theta[rad]$")
    axes[3].set_ylabel("$\\omega[\\frac{1}{s}]$")
    axes[3].grid()
    plt.tight_layout()
    plt.show()
    # plt.savefig("fizikai.png")


def dupla():
    L = 1
    q = 0
    Omega_D = 0
    F_D = 0
    nonlinear = True
    theta = np.pi/2
    theta2 = np.pi/2
    omega = 0
    omega2 = 0
    tmax = 100
    pendulum.run_double_pendulum("RK45_adaptive", nonlinear, L, q,
                            Omega_D, F_D, theta, theta2, omega, omega2, tmax)

    data = np.loadtxt('pendulum.data').T
    plt.plot(data[0], data[3])
    plt.show()


L = 1.0
res = (100, 100)
theta = np.linspace(-np.pi, np.pi, res[0])
theta2 = np.linspace(-np.pi, np.pi, res[1])
omega = 0.0
omega2 = 0.0

def f(t1, t2):
	print(f"{(t1*len(theta)+t2)/(len(theta)*len(theta2))*100:.2f}%", end="\r")
	start = time.time()
	tmp = pendulum.flipover(L, theta[t1], omega, theta2[t2], omega2, 10000*np.sqrt(1/9.8))
	a = time.time() - start
	return tmp, a

def flipover():
	pairs = itertools.product(range(res[0]), range(res[1]))
	num_workers = multiprocessing.cpu_count()

	start = time.time()
	with multiprocessing.Pool(processes=num_workers) as pool:
		results = [pool.apply_async(f, args=(arg1, arg2)) for arg1, arg2 in pairs]
		results = [r.get() for r in results]  # Collect results dynamically
	print(time.time() - start)

	results = np.array(results)
	times = results[:, 0]
	a = results[:, 1]
	times[times==-1] = np.nan
	times = np.reshape(times, res)
	print(np.average(a))
	print(np.std(a))
	np.save("times", times)
	# tmp = pendulum.flipover(L, theta[i], omega, theta2[j], omega2, 1000.0)
	
	norm = LogNorm(vmin=np.nanmin(times)+1e-6, vmax=np.nanmax(times))
	plt.imshow(times.T, origin="lower", norm=norm)
	plt.savefig("times.png", dpi=300)
flipover()
