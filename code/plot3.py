from sys import *
from math import *
import numpy as np
import matplotlib.pyplot as plt 

if len(argv) != 2 :
	print("Il n'y a pas assez d'arguments")
	exit(1)

_, fichier = argv
A = np.loadtxt(fichier)
dossier = "../image/"

tabA = A[:, 0]
tabN = A[:, 1]
tabN2 = A[:, 2]
tabN4 = A[:, 3]
tabN16 = A[:, 4]
tabH = A[:, 5]
tabH2 = A[:, 6]

tabN3 = [n**3 for n in tabN]

tabED1 = A[:, 7]
tabEJ1 = A[:, 8]
tabEGS1 = A[:, 9]
tabES1 = A[:, 10]
tabEMJ1 = A[:, 11]
tabEMGS1 = A[:, 12]

tabED2 = A[:, 13]
tabEJ2 = A[:, 14]
tabEGS2 = A[:, 15]
tabES2 = A[:, 16]
tabEMJ2 = A[:, 17]
tabEMGS2 = A[:, 18]

tabRJ1 = A[:, 19]
tabRGS1 = A[:, 20]
tabRS1 = A[:, 21]
tabRMJ1 = A[:, 22]
tabRMGS1 = A[:, 23]

tabRJ2 = A[:, 24]
tabRGS2 = A[:, 25]
tabRS2 = A[:, 26]
tabRMJ2 = A[:, 27]
tabRMGS2 = A[:, 27]

tabIJ1 = A[:, 29]
tabIGS1 = A[:, 30]
tabIS1 = A[:, 31]

tabIJ2 = A[:, 32]
tabIGS2 = A[:, 33]
tabIS2 = A[:, 34]

tabTD1 = A[:, 35]
tabTJ1 = A[:, 36]
tabTGS1 = A[:, 37]
tabTS1 = A[:, 38]
tabTMJ1 = A[:, 39]
tabTMGS1 = A[:, 40]

tabTD2 = A[:, 41]
tabTJ2 = A[:, 42]
tabTGS2 = A[:, 43]
tabTS2 = A[:, 44]
tabTMJ2 = A[:, 45]
tabTMGS2 = A[:, 46]

# plt.plot(tabN4, tabN4, label = "N")
# plt.plot(tabN4, tabN16, label = "N^2")
# plt.plot(tabN4, tabA, label = "Cond(A)")
# plt.title("Conditionning number of A")
# plt.xlabel("N = matrix size")
# plt.ylabel("Conditionning number")
# plt.xscale("log", base = 2)
# plt.yscale("log", base = 2)
# plt.legend()
# plt.show()

# plt.plot(tabN2, tabN2, label = "N")
# plt.plot(tabN2, tabN4, label = "N^2")
# plt.plot(tabN2, tabA, label = "Cond(A)")
# plt.title("Conditionning number of A")
# plt.xlabel("N = domaine size")
# plt.ylabel("Conditionning number")
# plt.xscale("log", base = 2)
# plt.yscale("log", base = 2)
# plt.legend()
# plt.show()

# plt.plot(tabN, tabN, label = "N")
# plt.plot(tabN, tabN2, label = "N^2")
# plt.plot(tabN, tabA, label = "Cond(A)")
# plt.title("Conditionning number of A")
# plt.xlabel("N = discretisation size")
# plt.ylabel("Conditionning number")
# plt.xscale("log", base = 2)
# plt.yscale("log", base = 2)
# plt.legend()
# plt.show()

# p4 = 
plt.plot(tabN, tabH, label = "H")
plt.plot(tabN, tabH2, label = "H^2")
plt.plot(tabN, tabED1, label = "Direct")
plt.plot(tabN, tabEJ1, label = "Jacobi")
plt.plot(tabN, tabEGS1, label = "Gauss-Seidel")
plt.plot(tabN, tabES1, label = "SOR")
plt.plot(tabN, tabEMJ1, label = "Multigrid Jacobi")
plt.plot(tabN, tabEMGS1, label = "Multigrid Gauss-Seidel")
plt.title("Error for F1")
plt.xlabel("N")
plt.ylabel("Error")
plt.xscale("log", base = 2)
plt.yscale("log", base = 2)
plt.legend()
plt.show()

# p5 = 
plt.plot(tabN, tabN, label = "N")
plt.plot(tabN, tabH, label = "H")
plt.plot(tabN, tabH2, label = "H^2")
plt.plot(tabN, tabRJ1, label = "Jacobi")
plt.plot(tabN, tabRGS1, label = "Gauss-Seidel")
plt.plot(tabN, tabRS1, label = "SOR")
plt.plot(tabN, tabRMJ1, label = "Multigrid Jacobi")
plt.plot(tabN, tabRMGS1, label = "Multigrid Gauss-Seidel")
plt.title("Residual for F1")
plt.xlabel("N")
plt.ylabel("Residual")
plt.xscale("log", base = 2)
plt.yscale("log", base = 2)
plt.legend()
plt.show()

# p6 = 
plt.plot(tabN, tabH, label = "H")
plt.plot(tabN, tabH2, label = "H^2")
plt.plot(tabN, tabED2, label = "Direct")
plt.plot(tabN, tabEJ2, label = "Jacobi")
plt.plot(tabN, tabEGS2, label = "Gauss-Seidel")
plt.plot(tabN, tabES2, label = "SOR")
plt.plot(tabN, tabEMJ2, label = "Multigrid Jacobi")
plt.plot(tabN, tabEMGS2, label = "Multigrid Gauss-Seidel")
plt.title("Error for F2")
plt.xlabel("N")
plt.ylabel("Error")
plt.xscale("log", base = 2)
plt.yscale("log", base = 2)
plt.legend()
plt.show()

# p7 = 
plt.plot(tabN, tabN, label = "N")
plt.plot(tabN, tabH, label = "H")
plt.plot(tabN, tabH2, label = "H^2")
plt.plot(tabN, tabRJ2, label = "Jacobi")
plt.plot(tabN, tabRGS2, label = "Gauss-Seidel")
plt.plot(tabN, tabRS2, label = "SOR")
plt.plot(tabN, tabRMJ2, label = "Multigrid Jacobi")
plt.plot(tabN, tabRMGS2, label = "Multigrid Gauss-Seidel")
plt.title("Residual for F2")
plt.xlabel("N")
plt.ylabel("Residual")
plt.xscale("log", base = 2)
plt.yscale("log", base = 2)
plt.legend()
plt.show()


# # p8 = plt.plot(tabN, [tabN, tabN2, tabIJ1, tabIGS1, tabIS1], label = ["N" "N^2" "Jacobi" "Gauss-Seidel" "SOR"], title = "Number of iterations depending on N for F1", xlabel = "N", ylabel = "Number of iterations")
# plt.plot(tabN, tabN, label = "N")
# plt.plot(tabN, tabN2, label = "N^2")
# plt.plot(tabN, tabN3, label = "N^3")
# plt.plot(tabN, tabIJ1, label = "Jacobi")
# plt.plot(tabN, tabIGS1, label = "Gauss-Seidel")
# plt.plot(tabN, tabIS1, label = "SOR")
# plt.title("Number of iterations depending on N for F1")
# plt.xlabel("N")
# plt.ylabel("Number of iterations")
# plt.xscale("log", base = 2)
# plt.yscale("log", base = 2)
# plt.legend()
# plt.show()

# # p9 = plt.plot(tabN, [tabN, tabN2, tabIJ2, tabIGS2, tabIS2], label = ["N" "N^2" "Jacobi" "Gauss-Seidel" "SOR"], title = "Number of iterations depending on N for F2", xlabel = "N", ylabel = "Number of iterations")
# plt.plot(tabN, tabN, label = "N")
# plt.plot(tabN, tabN2, label = "N^2")
# plt.plot(tabN, tabN3, label = "N^3")
# plt.plot(tabN, tabIJ2, label = "Jacobi")
# plt.plot(tabN, tabIGS2, label = "Gauss-Seidel")
# plt.plot(tabN, tabIS2, label = "SOR")
# plt.title("Number of iterations depending on N for F2")
# plt.xlabel("N")
# plt.ylabel("Number of iterations")
# plt.xscale("log", base = 2)
# plt.yscale("log", base = 2)
# plt.legend()
# plt.show()


# # p10 = title = "Execution time depending on N for F1", xlabel = "N", ylabel = "Execution time (s)")
# plt.plot(tabN, tabN, label = "N")
# plt.plot(tabN, tabN2, label = "N^2")
# plt.plot(tabN, tabN3, label = "N^3")
# plt.plot(tabN, tabTD1, label = "Direct")
# plt.plot(tabN, tabTJ1, label = "Jacobi")
# plt.plot(tabN, tabTGS1, label = "Gauss-Seidel")
# plt.plot(tabN, tabTS1, label = "SOR")
# plt.plot(tabN, tabTMJ1, label = "Multigrid Jacobi")
# plt.plot(tabN, tabTMGS1, label = "Multigrid Gauss-Seidel")
# plt.title("Execution time depending on N for F1")
# plt.xlabel("N")
# plt.ylabel("Execution time (s)")
# plt.xscale("log", base = 2)
# plt.yscale("log", base = 2)
# plt.legend()
# plt.show()

# # p11 = plt.plot(tabN, [tabN, tabN2, tabTD2, tabTJ2, tabTGS2, tabTS2, tabTMJ2, tabTMGS2], label = ["N" "N^2" "Jacobi" "Gauss-Seidel" "SOR" "Multigrid_Jacobi" "Multigrid_Gauss-Seidel"], title = "Execution time depending on N for F2", xlabel = "N", ylabel = "Execution time (s)")
# plt.plot(tabN, tabN, label = "N")
# plt.plot(tabN, tabN2, label = "N^2")
# plt.plot(tabN, tabN3, label = "N^3")
# plt.plot(tabN, tabTD2, label = "Direct")
# plt.plot(tabN, tabTJ2, label = "Jacobi")
# plt.plot(tabN, tabTGS2, label = "Gauss-Seidel")
# plt.plot(tabN, tabTS2, label = "SOR")
# plt.plot(tabN, tabTMJ2, label = "Multigrid Jacobi")
# plt.plot(tabN, tabTMGS2, label = "Multigrid Gauss-Seidel")
# plt.title("Execution time depending on N for F2")
# plt.xlabel("N")
# plt.ylabel("Execution time (s)")
# plt.xscale("log", base = 2)
# plt.yscale("log", base = 2)
# plt.legend()
# plt.show()

# # savefig(p4, "Erreur pour F1.png")
# # savefig(p5, "Residu pour F1.png")
# # savefig(p6, "Erreur pour F2.png")
# # savefig(p7, "Residu pour F2.png")
# # savefig(p8, "Iterations pour F1.png")
# # savefig(p9, "Iterations pour F2.png")
# # savefig(p10, "Temps pour F1.png")
# # savefig(p11, "Temps pour F2.png")