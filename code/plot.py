from sys import *
from math import *
import numpy as np
import matplotlib.pyplot as plt 

if len(argv) != 2 :
	print("Il n'y a pas assez d'arguments")
	exit(1)

_, fichier = argv
A = np.loadtxt(fichier)

nom = fichier[30:len(fichier) - 4]
dossier = "../image/Python/Comparaisons/"

if nom == "Temps_des_méthodes_itératives":
	tabN = A[:, 0]
	tabN2 = A[:, 1]
	tabN3 = A[:, 2]
	tabTD1 = A[:, 3]
	tabTJ1 = A[:, 4]
	tabTGS1 = A[:, 5]
	tabTS1 = A[:, 6]
	tabTD2 = A[:, 7]
	tabTJ2 = A[:, 8]
	tabTGS2 = A[:, 9]
	tabTS2 = A[:, 10]
	
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
	tabN = A[:, 0]
	tabH = A[:, 1]
	tabN2 = A[:, 2]
	tabH2 = A[:, 3]
	tabC = A[:, 4]
	tabAD1 = A[:, 5]
	tabAJ1 = A[:, 6]
	tabAGS1 = A[:, 7]
	tabAS1 = A[:, 8]
	tabRD1 = A[:, 9]
	tabRJ1 = A[:, 10]
	tabRGS1 = A[:, 11]
	tabRS1 = A[:, 12]
	tabAD2 = A[:, 13]
	tabAJ2 = A[:, 14]
	tabAGS2 = A[:, 15]
	tabAS2 = A[:, 16]
	tabRD2 = A[:, 17]
	tabRJ2 = A[:, 18]
	tabRGS2 = A[:, 19]
	tabRS2 = A[:, 20]
	
	plt.plot(tabN, tabC, label = "Cond(A)")
	plt.plot(tabN, tabN, label = "N")
	plt.plot(tabN, tabN2, label = "N^2")
	plt.ylabel("Conditionnement de la matrice A")
	plt.xlabel("N")
	plt.title("Conditionnement de la matrice A en fonction de N")
	plt.legend()
	plt.xscale("log")
	plt.yscale("log")
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
	if nom[0] == "T":
		N = nom[38:]
		tabPre = A[:, 0]
		tabPost = A[:, 1]
		tabTMJ1 = A[:, 2]
		tabTMGS1 = A[:, 3]
		tabTMJ2 = A[:, 4]
		tabTMGS2 = A[:, 5]

		tabPre = np.reshape(tabPre, (10, 10))
		tabPost = np.reshape(tabPost, (10, 10))
		tabTMJ1 = np.reshape(tabTMJ1, (10, 10))
		tabTMGS1 = np.reshape(tabTMGS1, (10, 10))
		tabTMJ2 = np.reshape(tabTMJ2, (10, 10))
		tabTMGS2 = np.reshape(tabTMGS2, (10, 10))

		plt.contour(tabPre, tabPost, tabTMJ1)
		plt.xlabel("pre")
		plt.ylabel("post")
		plt.colorbar()
		plt.grid()
		plt.title("Temps d'exécution de la méthode multigrid avec Jacobi et N = " + N + " pour F1")
		plt.savefig(dossier + "Temps d'exécution de la méthode multigrid avec Jacobi et N = " + N + " pour F1.png")
		plt.show()
		
		plt.contour(tabPre, tabPost, tabTMGS1)
		plt.xlabel("pre")
		plt.ylabel("post")
		plt.colorbar()
		plt.grid()
		plt.title("Temps d'exécution de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F1")
		plt.savefig(dossier + "Temps d'exécution de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F1.png")
		plt.show()
		
		plt.contour(tabPre, tabPost, tabTMJ2)
		plt.xlabel("pre")
		plt.ylabel("post")
		plt.colorbar()
		plt.grid()
		plt.title("Temps d'exécution de la méthode multigrid avec Jacobi et N = " + N + " pour F2")
		plt.savefig(dossier + "Temps d'exécution de la méthode multigrid avec Jacobi et N = " + N + " pour F2.png")
		plt.show()
		
		plt.contour(tabPre, tabPost, tabTMGS2)
		plt.xlabel("pre")
		plt.ylabel("post")
		plt.colorbar()
		plt.grid()
		plt.title("Temps d'exécution de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F2")
		plt.savefig(dossier + "Temps d'exécution de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F2.png")
		plt.show()

	else:
		N = nom[44:]
		tabPre = A[:, 0]
		tabPost = A[:, 1]
		tabTAJ1 = A[:, 2]
		tabTAGS1 = A[:, 3]
		tabTRJ1 = A[:, 4]
		tabTRGS1 = A[:, 5]
		tabTAJ2 = A[:, 6]
		tabTAGS2 = A[:, 7]
		tabTRJ2 = A[:, 8]
		tabTRGS2 = A[:, 9]
		
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
		
		plt.contour(tabPre, tabPost, tabTAJ1)
		plt.xlabel("pre")
		plt.ylabel("post")
		plt.colorbar()
		plt.grid()
		plt.title("Erreur absolue de la méthode multigrid avec Jacobi et N = " + N + " pour F1")
		plt.savefig(dossier + "Erreur absolue de la méthode multigrid avec Jacobi et N = " + N + " pour F1.png")
		plt.show()
		
		plt.contour(tabPre, tabPost, tabTAGS1)
		plt.xlabel("pre")
		plt.ylabel("post")
		plt.colorbar()
		plt.grid()
		plt.title("Erreur absolue de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F1")
		plt.savefig(dossier + "Erreur absolue de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F1.png")
		plt.show()
		
		plt.contour(tabPre, tabPost, tabTRJ1)
		plt.xlabel("pre")
		plt.ylabel("post")
		plt.colorbar()
		plt.grid()
		plt.title("Erreur relative de la méthode multigrid avec Jacobi et N = " + N + " pour F1")
		plt.savefig(dossier + "Erreur relative de la méthode multigrid avec Jacobi et N = " + N + " pour F1.png")
		plt.show()
		
		plt.contour(tabPre, tabPost, tabTRGS1)
		plt.title("Erreur relative de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F1")
		plt.xlabel("pre")
		plt.ylabel("post")
		plt.colorbar()
		plt.grid()
		plt.savefig(dossier + "Erreur relative de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F1.png")
		plt.show()
		
		plt.contour(tabPre, tabPost, tabTAJ2)
		plt.xlabel("pre")
		plt.ylabel("post")
		plt.title("Erreur absolue de la méthode multigrid avec Jacobi et N = " + N + " pour F2")
		plt.colorbar()
		plt.grid()
		plt.savefig(dossier + "Erreur absolue de la méthode multigrid avec Jacobi et N = " + N + " pour F2.png")
		plt.show()
		
		plt.contour(tabPre, tabPost, tabTAGS2)
		plt.xlabel("pre")
		plt.ylabel("post")
		plt.colorbar()
		plt.grid()
		plt.title("Erreur absolue de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F2")
		plt.savefig(dossier + "Erreur absolue de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F2.png")
		plt.show()
		
		plt.contour(tabPre, tabPost, tabTRJ2)
		plt.xlabel("pre")
		plt.ylabel("post")
		plt.title("Erreur relative de la méthode multigrid avec Jacobi et N = " + N + " pour F2")
		plt.colorbar()
		plt.grid()
		plt.savefig(dossier + "Erreur relative de la méthode multigrid avec Jacobi et N = " + N + " pour F2.png")
		plt.show()
		
		plt.contour(tabPre, tabPost, tabTRGS2)
		plt.xlabel("pre")
		plt.ylabel("post")
		plt.title("Erreur relative de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F2")
		plt.colorbar()
		plt.grid()
		plt.savefig(dossier + "Erreur relative de la méthode multigrid avec Gauss Seidel et N = " + N + " pour F2.png")
		plt.show()