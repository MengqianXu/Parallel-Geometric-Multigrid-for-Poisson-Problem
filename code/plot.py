from sys import *
from math import *
import numpy as np
import matplotlib.pyplot as plt 

if len(argv) != 2 :
	print("Il n'y a pas assez d'arguments")
	exit(1)

_, fichier = argv
A = np.loadtxt(fichier)

nom = fichier[17:len(fichier) - 4]
dossier = "../image/Python/"

if nom == "Temps_des_méthodes_itératives":
	tabN = []
	tabN2 = []
	tabN3 = []
	tabTD1 = []
	tabTJ1 = []
	tabTGS1 = []
	tabTS1 = []
	tabTD2 = []
	tabTJ2 = []
	tabTGS2 = []
	tabTS2 = []
	
	for i in range(len(A)):
		tabN.append(A[i, 0])
		tabN2.append(A[i, 1])
		tabN3.append(A[i, 2])
		tabTD1.append(A[i, 3])
		tabTJ1.append(A[i, 4])
		tabTGS1.append(A[i, 5])
		tabTS1.append(A[i, 6])
		tabTD2.append(A[i, 7])
		tabTJ2.append(A[i, 8])
		tabTGS2.append(A[i, 9])
		tabTS2.append(A[i, 10])

	plt.plot(tabN, tabN, label = "N")
	plt.plot(tabN, tabN2, label = "N^2")
	plt.plot(tabN, tabN3, label = "N^3")
	plt.plot(tabN, tabTD1, label = "Direct")
	plt.plot(tabN, tabTJ1, label = "Jacobi")
	plt.plot(tabN, tabTGS1, label = "Gauss Seidel")
	plt.plot(tabN, tabTS1, label = "SOR avec w = 0.5")
	plt.title("Temps d'exécution des méthodes itératives")
	plt.ylabel("Temps en secondes")
	plt.xlabel("N")
	plt.legend()
	plt.xscale("log")
	plt.yscale("log")
	plt.savefig(dossier + "Temps d'exécution des méthodes itératives.png")

elif nom == "Convergence_des_méthodes_itératives":
	tabN = []
	tabH = []
	tabH2 = []
	tabC = []
	tabAD1 = []
	tabAJ1 = []
	tabAGS1 = []
	tabAS1 = []
	tabRD1 = []
	tabRJ1 = []
	tabRGS1 = []
	tabRS1 = []
	tabAD2 = []
	tabAJ2 = []
	tabAGS2 = []
	tabAS2 = []
	tabRD2 = []
	tabRJ2 = []
	tabRGS2 = []
	tabRS2 = []
	
	for i in range(len(A)):
		tabN.append(A[i, 0])
		tabH.append(A[i, 1])
		tabH2.append(A[i, 2])
		tabC.append(A[i, 3])
		
		tabAD1.append(A[i, 4])
		tabAJ1.append(A[i, 5])
		tabAGS1.append(A[i, 6])
		tabAS1.append(A[i, 7])
		
		tabRD1.append(A[i, 8])
		tabRJ1.append(A[i, 9])
		tabRGS1.append(A[i, 10])
		tabRS1.append(A[i, 11])
		
		tabAD2.append(A[i, 12])
		tabAJ2.append(A[i, 13])
		tabAGS2.append(A[i, 14])
		tabAS2.append(A[i, 15])
		
		tabRD2.append(A[i, 16])
		tabRJ2.append(A[i, 17])
		tabRGS2.append(A[i, 18])
		tabRS2.append(A[i, 19])

	plt.plot(tabN, tabC)
	plt.ylabel("Conditionnement de la matrice A")
	plt.xlabel("N")
	plt.title("Conditionnement de la matrice A en fonction de N")
	plt.savefig(dossier + "Conditionnement de la matrice A en fonction de N.png")
	plt.show()

	plt.plot(tabN, tabH, label = "H")
	plt.plot(tabN, tabH2, label = "H^2")
	plt.plot(tabN, tabAD1, label = "Direct")
	plt.plot(tabN, tabAJ1, label = "Jacobi")
	plt.plot(tabN, tabAGS1, label = "Gauss Seidel")
	plt.plot(tabN, tabAS1, label = "SOR avec w = 0.5")
	plt.title("Erreur absolue des méthodes itératives pour F1")
	plt.ylabel("Temps en secondes")
	plt.xlabel("N")
	plt.legend()
	plt.xscale("log")
	plt.yscale("log")
	plt.savefig(dossier + "Erreur absolue des méthodes itératives pour F1.png")
	plt.show()
	
	plt.plot(tabN, tabH, label = "H")
	plt.plot(tabN, tabH2, label = "H^2")
	plt.plot(tabN, tabRD1, label = "Direct")
	plt.plot(tabN, tabRJ1, label = "Jacobi")
	plt.plot(tabN, tabRGS1, label = "Gauss Seidel")
	plt.plot(tabN, tabRS1, label = "SOR avec w = 0.5")
	plt.title("Erreur relative des méthodes itératives pour F1")
	plt.ylabel("Temps en secondes")
	plt.xlabel("N")
	plt.legend()
	plt.xscale("log")
	plt.yscale("log")
	plt.savefig(dossier + "Erreur relative des méthodes itératives pour F1.png")
	plt.show()
	
	plt.plot(tabN, tabH, label = "H")
	plt.plot(tabN, tabH2, label = "H^2")
	plt.plot(tabN, tabAD2, label = "Direct")
	plt.plot(tabN, tabAJ2, label = "Jacobi")
	plt.plot(tabN, tabAGS2, label = "Gauss Seidel")
	plt.plot(tabN, tabAS2, label = "SOR avec w = 0.5")
	plt.title("Erreur absolue des méthodes itératives pour F2")
	plt.ylabel("Temps en secondes")
	plt.xlabel("N")
	plt.legend()
	plt.xscale("log")
	plt.yscale("log")
	plt.savefig(dossier + "Erreur absolue des méthodes itératives pour F2.png")
	plt.show()
	
	plt.plot(tabN, tabH, label = "H")
	plt.plot(tabN, tabH2, label = "H^2")
	plt.plot(tabN, tabRD2, label = "Direct")
	plt.plot(tabN, tabRJ2, label = "Jacobi")
	plt.plot(tabN, tabRGS2, label = "Gauss Seidel")
	plt.plot(tabN, tabRS2, label = "SOR avec w = 0.5")
	plt.title("Erreur relative des méthodes itératives pour F2")
	plt.ylabel("Temps en secondes")
	plt.xlabel("N")
	plt.legend()
	plt.xscale("log")
	plt.yscale("log")
	plt.savefig(dossier + "Erreur relative des méthodes itératives pour F2.png")

else:
	if fichier[17] == "T":
		N = fichier[38 + 17:len(fichier) - 4]
		tabPre = []
		tabPost = []
		tabTMJ1 = []
		tabTMGS1 = []
		tabTMJ2 = []
		tabTMGS2 = []

		for i in range(len(A)):
			tabPre.append(A[i, 0])
			tabPost.append(A[i, 1])
			tabTMJ1.append(A[i, 2])
			tabTMGS1.append(A[i, 3])
			tabTMJ2.append(A[i, 4])
			tabTMGS2.append(A[i, 5])

		tabPre = np.reshape(tabPre, (10, 10))
		tabPost = np.reshape(tabPost, (10, 10))
		tabTMJ1 = np.reshape(tabTMJ1, (10, 10))
		tabTMGS1 = np.reshape(tabTMGS1, (10, 10))
		tabTMJ2 = np.reshape(tabTMJ2, (10, 10))
		tabTMGS2 = np.reshape(tabTMGS2, (10, 10))

		plt.contour(tabPost, tabPre, tabTMJ1)
		plt.title("Temps d'exécution de la méthode multigrid avec Jacobi et N = " + N + " pour F1")
		plt.savefig(dossier + "Temps d'exécution de la méthode multigrid avec Jacobi et N = " + N + " pour F1.png")
		plt.show()
		
		plt.contour(tabPost, tabPre, tabTMGS1)
		plt.title("Temps d'exécution de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F1")
		plt.savefig(dossier + "Temps d'exécution de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F1.png")
		plt.show()
		
		plt.contour(tabPost, tabPre, tabTMJ2)
		plt.title("Temps d'exécution de la méthode multigrid avec Jacobi et N = " + N + " pour F2")
		plt.savefig(dossier + "Temps d'exécution de la méthode multigrid avec Jacobi et N = " + N + " pour F2.png")
		plt.show()
		
		plt.contour(tabPost, tabPre, tabTMGS2)
		plt.title("Temps d'exécution de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F2")
		plt.savefig(dossier + "Temps d'exécution de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F2.png")
		plt.show()

	else:
		N = fichier[38 + 23:len(fichier) - 4]
		tabPre = []
		tabPost = []
		tabTAJ1 = []
		tabTAGS1 = []
		tabTRJ1 = []
		tabTRGS1 = []
		tabTAJ2 = []
		tabTAGS2 = []
		tabTRJ2 = []
		tabTRGS2 = []

		for i in range(len(A)):
			tabPre.append(A[i, 0])
			tabPost.append(A[i, 1])
			tabTAJ1.append(A[i, 2])
			tabTAGS1.append(A[i, 3])
			tabTRJ1.append(A[i, 4])
			tabTRGS1.append(A[i, 5])
			tabTAJ2.append(A[i, 6])
			tabTAGS2.append(A[i, 7])
			tabTRJ2.append(A[i, 8])
			tabTRGS2.append(A[i, 9])

		tabPre = np.reshape(tabPre, (10, 10))
		tabPost = np.reshape(tabPost, (10, 10))
		tabTAJ1 = np.reshape(tabTAJ1, (10, 10))
		tabTAGS1 = np.reshape(tabTAGS1, (10, 10))
		tabTRJ1 = np.reshape(tabTRJ1, (10, 10))
		tabTRGS1 = np.reshape(tabTRGS1, (10, 10))
		tabTAJ2 = np.reshape(tabTAJ2, (10, 10))
		tabTAGS2 = np.reshape(tabTAGS2, (10, 10))
		tabTRJ2 = np.reshape(tabTRJ2, (10, 10))
		tabTRGS2 = np.reshape(tabTRGS2, (10, 10))

		plt.contour(tabPost, tabPre, tabTAJ1)
		plt.title("Erreur absolue de la méthode multigrid avec Jacobi et N = " + N + " pour F1")
		plt.savefig(dossier + "Erreur absolue de la méthode multigrid avec Jacobi et N = " + N + " pour F1.png")
		plt.show()
		
		plt.contour(tabPost, tabPre, tabTAGS1)
		plt.title("Erreur absolue de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F1")
		plt.savefig(dossier + "Erreur absolue de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F1.png")
		plt.show()
		
		plt.contour(tabPost, tabPre, tabTRJ1)
		plt.title("Erreur relative de la méthode multigrid avec Jacobi et N = " + N + " pour F1")
		plt.savefig(dossier + "Erreur relative de la méthode multigrid avec Jacobi et N = " + N + " pour F1.png")
		plt.show()
		
		plt.contour(tabPost, tabPre, tabTRGS1)
		plt.title("Erreur relative de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F1")
		plt.savefig(dossier + "Erreur relative de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F1.png")
		plt.show()
		
		plt.contour(tabPost, tabPre, tabTAJ2)
		plt.title("Erreur absolue de la méthode multigrid avec Jacobi et N = " + N + " pour F2")
		plt.savefig(dossier + "Erreur absolue de la méthode multigrid avec Jacobi et N = " + N + " pour F2.png")
		plt.show()
		
		plt.contour(tabPost, tabPre, tabTAGS2)
		plt.title("Erreur absolue de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F2")
		plt.savefig(dossier + "Erreur absolue de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F2.png")
		plt.show()
		
		plt.contour(tabPost, tabPre, tabTRJ2)
		plt.title("Erreur relative de la méthode multigrid avec Jacobi et N = " + N + " pour F2")
		plt.savefig(dossier + "Erreur relative de la méthode multigrid avec Jacobi et N = " + N + " pour F2.png")
		plt.show()
		
		plt.contour(tabPost, tabPre, tabTRGS2)
		plt.title("Erreur relative de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F2")
		plt.savefig(dossier + "Erreur relative de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F2.png")
		plt.show()