from sys import *
from math import *
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt 

if len(argv) < 2 :
	print("Il n'y a pas assez d'arguments")
	exit(1)

dossier = "../image/Python/Resultats/"

for i in range(1, len(argv)):
	fichier = argv[i]
	A = np.loadtxt(fichier)

	nom = fichier[27:len(fichier) - 4]
	# print(nom[19:])
	N = int(nom[19:])

	tabX = np.reshape(A[:, 0], (N, N))
	tabY = np.reshape(A[:, 1], (N, N))
	U1 = np.reshape(A[:, 2], (N, N))
	sol1 = np.reshape(A[:, 3], (N, N))
	resJ1 = np.reshape(A[:, 4], (N, N))
	resGS1 = np.reshape(A[:, 5], (N, N))
	resS1 = np.reshape(A[:, 6], (N, N))
	resMJ1 = np.reshape(A[:, 7], (N, N))
	resMGS1 = np.reshape(A[:, 8], (N, N))
	U2 = np.reshape(A[:, 9], (N, N))
	sol2 = np.reshape(A[:, 10], (N, N))
	resJ2 = np.reshape(A[:, 11], (N, N))
	resGS2 = np.reshape(A[:, 12], (N, N))
	resS2 = np.reshape(A[:, 13], (N, N))
	resMJ2 = np.reshape(A[:, 14], (N, N))
	resMGS2 = np.reshape(A[:, 15], (N, N))
	
	axe = plt.axes(projection ='3d')
	axe.plot_surface(tabX, tabY, U1)
	plt.title("Vraie solution pour F1 avec N = " + str(N))
	plt.savefig(dossier + "Vraie solution pour F1 avec N = " + str(N) + ".png")
	# plt.show()

	axe = plt.axes(projection ='3d')
	axe.plot_surface(tabX, tabY, sol1)
	plt.title("Solution directe pour F1 avec N = " + str(N))
	plt.savefig(dossier + "Solution directe pour F1 avec N = " + str(N) + ".png")
	# plt.show()

	axe = plt.axes(projection ='3d')
	axe.plot_surface(tabX, tabY, resJ1)
	plt.title("Solution Jacobi pour F1 avec N = " + str(N))
	plt.savefig(dossier + "Solution Jacobi pour F1 avec N = " + str(N) + ".png")
	# plt.show()

	axe = plt.axes(projection ='3d')
	axe.plot_surface(tabX, tabY, resGS1)
	plt.title("Solution Gauss-Seidel pour F1 avec N = " + str(N))
	plt.savefig(dossier + "Solution Gauss-Seidel pour F1 avec N = " + str(N) + ".png")
	# plt.show()

	axe = plt.axes(projection ='3d')
	axe.plot_surface(tabX, tabY, resS1)
	plt.title("Solution SOR (avec w = 0.5) pour F1 avec N = " + str(N))
	plt.savefig(dossier + "Solution SOR (avec w = 0.5) pour F1 avec N = " + str(N) + ".png")
	# plt.show()

	axe = plt.axes(projection ='3d')
	axe.plot_surface(tabX, tabY, resMJ1)
	plt.title("Solution multigrid avec Jacobi pour F1 avec N = " + str(N))
	plt.savefig(dossier + "Solution multigrid avec Jacobi pour F1 avec N = " + str(N) + ".png")
	# plt.show()

	axe = plt.axes(projection ='3d')
	axe.plot_surface(tabX, tabY, resMGS1)
	plt.title("Solution multigrid avec Gauss-Seidel pour F1 avec N = " + str(N))
	plt.savefig(dossier + "Solution multigrid avec Gauss-Seidel pour F1 avec N = " + str(N) + ".png")
	# plt.show()

	axe = plt.axes(projection ='3d')
	axe.plot_surface(tabX, tabY, U2)
	plt.title("Vraie solution pour F2 avec N = " + str(N))
	plt.savefig(dossier + "Vraie solution pour F2 avec N = " + str(N) + ".png")
	# plt.show()

	axe = plt.axes(projection ='3d')
	axe.plot_surface(tabX, tabY, sol2)
	plt.title("Solution directe pour F2 avec N = " + str(N))
	plt.savefig(dossier + "Solution directe pour F2 avec N = " + str(N) + ".png")
	# plt.show()

	axe = plt.axes(projection ='3d')
	axe.plot_surface(tabX, tabY, resJ2)
	plt.title("Solution Jacobi pour F2 avec N = " + str(N))
	plt.savefig(dossier + "Solution Jacobi pour F2 avec N = " + str(N) + ".png")
	# plt.show()

	axe = plt.axes(projection ='3d')
	axe.plot_surface(tabX, tabY, resGS2)
	plt.title("Solution Gauss-Seidel pour F2 avec N = " + str(N))
	plt.savefig(dossier + "Solution Gauss-Seidel pour F2 avec N = " + str(N) + ".png")
	# plt.show()

	axe = plt.axes(projection ='3d')
	axe.plot_surface(tabX, tabY, resS2)
	plt.title("Solution SOR (avec w = 0.5) pour F2 avec N = " + str(N))
	plt.savefig(dossier + "Solution SOR (avec w = 0.5) pour F2 avec N = " + str(N) + ".png")
	# plt.show()

	axe = plt.axes(projection ='3d')
	axe.plot_surface(tabX, tabY, resMJ2)
	plt.title("Solution multigrid avec Jacobi pour F2 avec N = " + str(N))
	plt.savefig(dossier + "Solution multigrid avec Jacobi pour F2 avec N = " + str(N) + ".png")
	# plt.show()

	axe = plt.axes(projection ='3d')
	axe.plot_surface(tabX, tabY, resMGS2)
	plt.title("Solution multigrid avec Gauss-Seidel pour F2 avec N = " + str(N))
	plt.savefig(dossier + "Solution multigrid avec Gauss-Seidel pour F2 avec N = " + str(N) + ".png")
	# plt.show()