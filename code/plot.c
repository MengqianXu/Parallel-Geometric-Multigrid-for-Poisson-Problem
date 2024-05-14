#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/time.h>

int main(int argc, char **argv){
	char fichier[1000000] = "../data/Version3/Comparaisons/", commande[1000000];
	DIR* dossier = opendir(fichier);
	// if (dossier){
	// 	struct dirent* dir;
	// 	while((dir = readdir(dossier))){
	// 		char* nom = dir->d_name;
	// 		if((strcmp(".", nom)) && (strcmp("..", nom)) && (strcmp("Conclusions.txt", nom))){
	// 			// printf("%s\n", nom);
	// 			strcpy(fichier, "../data/Version3/Comparaisons/");
	// 			strcpy(commande, "python plot.py ");
	// 			strcat(fichier, nom);
	// 			strcat(commande, fichier);
	// 			system(commande);
	// 		}
	// 	}
	// 	closedir(dossier);
	// }
	strcpy(fichier, "../data/Version3/Resultats/");
	strcpy(commande, "python plot2.py ");
	dossier = opendir(fichier);
	if (dossier){
		struct dirent* dir;
		while((dir = readdir(dossier))){
			char* nom = dir->d_name;
			if((strcmp(".", nom)) && (strcmp("..", nom))){
				printf("%s\n", nom);
				strcpy(fichier, "../data/Version3/Resultats/");
				strcat(fichier, nom);
				strcat(fichier, " ");
				strcat(commande, fichier);
			}
		}
		closedir(dossier);
		system(commande);
	}
	else{
		printf("Il y a un probleme quelque part !\n");
	}
	return 0;
}