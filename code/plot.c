#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/time.h>

int main(int argc, char **argv){
	char fichier[1000000] = "../data/Version3/", commande[1000000];
	DIR* dossier = opendir(fichier);
	if (dossier){
		struct dirent* dir;
		// printf("Les fichiers presents sont : \n");
		while((dir = readdir(dossier))){
			char* nom = dir->d_name;
			if((strcmp(".", nom)) && (strcmp("..", nom))){
				strcpy(fichier, "../data/Version3/");
				strcpy(commande, "python3 plot.py ");
				// printf("%s\n", nom);
				strcat(fichier, nom);
				strcat(commande, fichier);
				// strcat(commande, "'");
				// printf("%s\n", commande);
				system(commande);
			}
		}
		closedir(dossier);
	}
	else{
		printf("Il y a un probleme quelque part !\n");
	}
	return 0;
}