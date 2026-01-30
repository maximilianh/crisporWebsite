# changes in CRISPOR - assistant mode

## 17/12/25 ;

ajout du mode "assistant"

## 22/12/25 :

### corrections dans getGeneSeq():
 
cette fonction retourne deux listes contenant la séquence des exons de tous les transcrits d'un gène, dans les fichiers .gp disponibles.
- une liste ('all.exons')
- une liste ('first.exons') correspondant aux exons du premiers tiers du CDS

#### changements depuis le 17/12/25 :
 
- retourne la liste des exons dans le bon ordre pour les gènes du brin "-".
- les exons de plus de 2300bp (MAXSEQLEN) sont scidés en plusieurs séquences.
- retrait des 5' et 3'UTRs.
- Un exon trop long pour entrer dans le premier tiers du CDS mais pouvant être tronqué l'est désormais.
- si le CDS est de < 100bp, tous les exons sont considérés (seuil à augmenter ?)

#### problèmes / à faire :

- pas d'annotation des transcripts ID dans le fichier .gp -> les trancrits sont pour l'instant numérotés selon leur ordre d'appariation.
- prise en compte des autres fichier d'annotation ? (genes.tsv, .bb, segments.bed) ? 

### Ajout du menu déroulant de la liste des gènes disponible en fonction de l'organisme sélectionné

via AJAX avec select2
chargement un peu lent (~2s), à voir avec d'autres génomes.

## 23/12/25

## ajouts divers 

- ajout du logo Celphedia
- possibilité de sélectionner un gene ID depuis le menu principal.
- possibilité de sélectionner le n°de l'exon à cibler.
- possibilité de choisir un PAM / type d'enzyme / longueur de guide personnalisé.
- utilisation de select2 sur tous les menus déroulants.

## à faire

- ajouter la possibilité de relancer la recherche sur l'exon suivant si sélection depuis geneID
- choix du transcrit consensus -> utiliser l'intersection de tous les transcrits dispo (cf. CHOPCHOP). Attention à ce qu'un transcrit "rare" ne dicte pas le consensus
- vérifier que les custom PAMs ne provoquent pas de bugs. +- FAIT
- définir les longeurs de guides possibles (pour l'instant : 10-30bp)
- définir d'autres types d'enzymes ?

# 26/12/25

## classement des guides selon un score global

- calcul de l'énergie libre minimum des structures secondaires du sgRNA avec RNAfold (avertissment si celle-ci est < à 3 kcal/mol)
- transfert du calcul du GC content dans mergeGuideInfo()
- calcul du score global : normalisation min-max des CFD specificity et rs3 efficiency scores, puis somme de ces scores avec application d'un coef (arbitraire pour l'instant). (0.60*CFD + 0.40*rs3)
- le score global est mis à une échelle de 0 à 100.
- une pénalité de 20 points est appliquée si l'énergie libre minimum est inférieure à 3kcal/mol (cf. Riesenberg et al 2025, attention unités à vérifier).
- une pénalité de 40 points est appliquée si le %GC inférieur à 25% ou supérieur à 75%. Dans CHOPCHOP, ces limites sont de 40-70%.
- ajout de la colonne "Global Score" dans le tableau des guides. Classement des guides avec ce score par défaut.

## à faire

- implémenter le EVA score de Riesenberg et al 2025, pour guides synthétiques ? 
- Vérifier / ajuster les paramètres de calcul de Global Score.
- ajouter des pénalités ?
- adapter le calcul du global score au mode cpf1 / saCas9.
- ajouter binary RNAfold dans /data/www/crispor/bin/Linux-aarch64 FAIT
- utiliser package python viennarna au lieu de RNAfold ? la performance ne semble pas être un problème 

# 29/12/25

## 1ère tentative d'implémentation du mode knock-in

- installation de protoSpaceJam dans /bin (doc https://czbiohub-sf.github.io/protoSpaceJAM/)
- mais multiples problèmes : protoSpaceJam fonctionne avec des PAMs pré-calculés sur un nombre restreint de génomes (Homme, souris, rat, zebrafish)
- de nombreuses fonctions sont redondantes avec CRISPOR (notamment pour le calcul des scores)
- solutions : soit réutiliser des fonctions (obtention et recodage de la séquence de l'ADN donneur, calcul des pénalités (séquences répétées, %GC, homogénéité...))
 mais ces fonction requièrent souvent beaucoup de paramètres obtenus ailleurs dans la programme.. pas sûr que ce soit possible à implémenter simplement.
- soit partir de 0 et implémenter les mêmes concepts (plus simple, mais le résultat sera probablement moins complet).

## implémentation "manuelle" du mode knock-in
- ajout du formulaire : l'utilisateur entre un gene ID et choisit la position (Nter ou Cter), ou alors entre une séquence avec '//' au niveau du site d'insertion.
- choix de la séquence d'insert : soit sélection parmi une liste de tags / linkers (+ choix de l'ordre), soit séquence custom.
- si geneID, la séquence du premier / dernier exon est extraite avec getGeneSeq() puis scindée au niveau du codon START / STOP
- si séquence, celle-ci est scindée aux '//' (un seul '//' autorisé)

## à faire

- étendre la séquence de 800bp (+ définir input pour la longueur des bras d'homologie) en 5' et 3' (geneID) FAIT
- séquence : mapper sur le génome puis étendre de 800bp FAIT
- merger la séquence d'insert, puis calculer %GC et segments répétés sur des bins de (50bp ?) 

#30/12/25

- design de l'ADN donneur avec la fonction getDonorSeq()
- réinstallation complète de crispor à cause d'un bug causé par l'installation de protoSpaceJam (mauvaise version de numpy)
- calcul du %GC sur des bins de 50bp
- 

## à faire :

- commit depuis branche anton
- (pré)calcul nb match dans genome de taille donnée
- adapter nb mismatch à longueur spacer (=< 3)
- déplacer calcul énergie libre dans effScores
- modifier contraintes PAM custom (3-8 nt, 2 non Ns) FAIT
- créer nouveau jobType pour recherche batch : multiseq / multiPAMs. (modifier JobType, readBatchParams(), runQueueWorker(), getOddTargets() (ou nouvelle fonction getOffTargets() ) (mysql au lieu de sqlite pour gérer plusieurs recherches en mm temps)
- knock-in :
    - vérifier que le gène code pour une protéine si insertion Nter / Cter.
    - recherche des guides 60 bp de part et 'autre du site d'insertion.
    - affichage résultats : toutes enzymes dispo (+ checkbox enzymes)
    - prise en compte de la distance au site d'insertion dans le classement des guides
    - guides chevauchant le site d'insertion en priorité
    - donneur sb ( + court, + simple à mettre en oeuvre, non adapté aux longs inserts ) + guides db.
    - design primers de part et d'autre des bras d'homologie ( attention à la longueur de l'amplicon)
    - retirer le calcul des scores d'efficacité et réduire de nb de mismatch pour rechercher plusieurs PAMs en mm temps.
    - si site clivage éloigné su site d'insertion -> recodage.
    - vérifier le cadre de lecture des bras d'homologie ?
    - indiquer le brin (codant / non codant) correspondant au donneur

# 07/01/26

## Ajouts divers

- ajout d'une pénalité pour les guides décrits par Graf et al, 2019.
- déplacement du de calcFreeEnergy() dans crisporEffScores.py 
- calcul du EVA score (pour guides synthétiques avec la fonction calcEvaScore() 
- calcul de l'énergie libre avec RNAstructure (lib à installer, utilisé par EVAscore)

## Mode Knock-in 

- extension des séquences (donneur / target) dans la limite des coordonnées chromosomiques
- definition chevauchement : <16nt

## à faire

- réfléchir à l'interface formulaire / résultats Knock-in (+ "général" que ProtoSpaceJam 
- installer RNAstructure / finir la fonction associée
- importer EVA score

# 08/01/26

## Ajouts divers 
- ajout de fonctions permettant d'obtenir une liste des transcrits et des exons (par n°) du geneID sélectionné, via AJAX
- (ne fonctionne pas encore parfaitement)

## Modification de getGeneSeq

- fonction trop longue et pas assez généraliste
- replacée par getGenePos : à partir d'un geneID, renvoie un dict des modèles de gènes correspondant à l'ID :
- contient : un dict par fichier .gp contenant un dict pour chaque colonne correspondnat à l'ID d'input

	- coordonnées des 5' et 3' UTRs (si le gène code pour une protéine)
	- coordonnées des exons (5' et 3' UTRs tronqués)
	- coordonnées des introns

## à faire :

- faire une nouvelle fonction getGeneSeq() prenant comme input l'outupt de getGenePos()
- optionellement, renvoie le premier tiers de la séquence codante.

## 09/01/26

## Modification de GetGeneSeq()
- Correction de getGenePos()
- ajout de getGeneSeq() : à partir des coordonnées issues de getGenePos(), retourne la liste des séquences correspondantes.
	- possibilité de sélectionner exons, intron, UTRs.
- ajout de getFirstThird(): retourne le premiers tiers de la séquence d'une lsite d'exons.
- ajout de seqSplit(): découpe une séquence en fragments de taille donnée

## à faire :
- faire un système d'envoi des résultats par mail en différé
- modifier le système de batch pour gérer plusieurs séquences / PAMs
- créer une fonction pour éxécuter getGenePos(), getGeneSeq().. et l'appeller dans assistantDispatcher()

## 12/01/26

## Modification de GetGeneSeq()

- ajout de fetchSeqsFromID() : à partir d'un geneID, retourne les séquences correspondant à la feature d'input 
(en éxécutant les fonctions ajoutées le 09/01)
	- ignore les séquences < 23 bp
	- fragmente les séquences > MAXSEQLEN
	- assigne un index correspondant au n°de la séquence avant traitement (ex, n° d'exon)
- retourne une liste de tuples [(id, seq)]
 
## gestion de multiples séquences / PAMs en batch 

- pramètre CGI : multiseq = liste de tuples (out fetchSeqsFromID() )

- newBatch() accepte le paramètre optionnel multiseq, qui est ensuite stocké dans le json
- dans crisprSearch() : ajout d'un batch avec paramètre "multiseq" -> ajout d'un job "multisearch" 
- dans runQueueWorker : si le jobType == "multisearch"
	- lance un sub-batch par séquence -> getOffTargets() puis parseOffTargets par seq # supprimé le 13/01/26

# 13/01/26

## gestion multi seq / PAM 

- fonction submitMultiSearch() : ajoute un job "multisearch" dans la queue
- dans runQueueWorker() : 
	- ajout de writeMultiFasta : écrit un fichier fasta correspondant aux quides de toutes les séquences en input
	- fasta header : >seqId.exonId.pamId (où seqId est l'ID correspondant à la séquence d'input et exonID correspond au n° de l'exon)
	- recherche des off-targets / calcul des scores dans processMultiSubmission()
	- pour l'instant : recherche des offtargets OK

## à faire

- adapter createBatchEffScoreTable() (simplification, retrait de scores)
- implémenter l'affichage des résultats
- faire la recherche pour plusieurs PAMs (et plusieurs geneIDs ?)

## à faire (divers)

- installer RNAstructure, ajouter MIT score à EVA score dans crispor.py
- optimiser getGenePos(), utiliser des coordonnées plutôt que des séquences
- ajouter le paramètre maxLen à getFirstThird()


# 14/01/26

# gestion multi seq / PAM

- recherche pour plusieurs pams : modification de processMultiSubmission()
	- écriture d'un fichier bed vide
	- dans une boucle : 
		- call de setupPamInfo()
		- writeMultiFasta() (recherche des guides pour chaque séquence en input)
		- findOffTargetsBwa() >> ajout des résultats dans le fichier bed
 	- index des guides / offtargets: PAM.seqId.exonId.pamId
- réduction du maximum de mismatchs à 3 lors de la recherche des offtargets (~5-10x plus rapide)

- call des effscores dans la boucle avec createMultiBatchEffScoreTable()

# bugs / à faire
- certains geneIDs retournent une liste vide
- Pour certains geneIDs, le job "multisearch" n'est jamais lancé.
	- (huître) NM_001308865.1 -> fonctionne
- finaliser le calcul des effscores

# 15/01/26

## gestion multiseq / PAM

- correction de processMultiSubmission et de calcMultiSaveEffScores
- le job "multisearch" fontionne pour tous les geneIDs (besoin d'effacer la base de données avec ./crispor.py --clear)
- homogéinisation du tableau des effscores (1 score / PAM + oof) ~ non terminé

## divers

- mise a jour de la branche anton ~ non terminé

# à faire

- headers dupliqués dans tableau effscores + 3 scores pour pam NGG
- keras ne peut pas importer le backend theano ?! -> impossible d'utiliser cpf1

# 16/21/26

## gestion multiseq / PAM

- mise au propre du tableau effscores
- dans crisprSearch -> parseMultisearchInfo

## réécriture de getGenePos()

- stocke uniquement les donneés relaives au geneID d'input
- sotckage des coordonnées plutôt que des séquences
- en 4 fonctions : 
	- getGenePos() retourne chrom, start, (exonStarts, exonEnds)
	- getFirstThird() : filtre les exons correspondant au premiers tiers de la seq. codante
	- formatExonPos() : conversion des coordonnées en posStr(chrom:start-end:strand), retrait des exons < PAMLEN et assignation du n° de l'exon)
        - getExonsFromId() : définit maxLen, PAMLEN -> call des trois fonctions précédentes. retourne [(exonIDs, exonPosStr)]

## divers 

- réinstallation du docker container + image (à jour avec master)
- réinitialisation de la branche anton -> push de la nouvelle branche sur le repo

## à faire 

- malgrès réinstallation comlpète, pb de version de keras toujours présent ??
- séparer multiseq / multiPAM OK ~~

- corriger l'écriture des effScores dans processMultiSeqSubmission()
- corrigder : CFD toujours entre 97 et 100

#19/01/25

## Gestion multiseq / multiPAM 

- séparation des jobs multiseq / multipam (multipam nécéssaire uniquement pour knock-in)
- jobType : "multiseq" : 
        - call de processMultiSeqSubmission() : 
                - par exon : obtention et extension de la séquence et calcul des effscores (dans un seul fichier) 
                - écriture d'un fichier fasta contenant tous les guides (pour tous les exons)
                - recherche des offtargets 
- lorsque le job est terminé :
        - call de parseMultiSearchInfo():
                - lecture des fichiers offtargets et effscores
                - par exon : call de mergeGuideInfo() et aggrégation des données dans allGuideData (liste) allGuideScores (dict) et allPamIdToSeq (dict) 
		- ajout de showExonAndPams (affichage de la séquence de / n° / longueur de chaque exon)
		- ou alternativement, affichage de la séquence des exons mis bout à bout et séparés par '//' (supprimé, à refaire demain)
 	
- dans colonne "position / strand" du tableau des guides : ajout de "in exon n"

## divers : 

- si sélection d'un geneID depuis le menu principal, affichage / sélection du nombre d'exons correspondant (via ajax)

## à faire / bugs

- corriger la position des pams (startDict) lorsque les exons sont affichés bout à bout
- pour certains gènes, le calcul des effscores retourne un dict vide
- si beaucoup de guides (> 25 ?) filtre des 25 meilleurs guides
- makePamLines retournes moins de PAMs que le nombre de guides possible.

# 20/01/26

## divers

- correction de getFirstThird (retournait la séquence complète si > maxlen)
- installation de RNAstructue : 
	- installation de make / g++
	- compilation des programmes puis de l'interface python
	- note : besoin de compiler l'interface python avec PYTHON=/data/www/crispor/venv/bin/python3

- ajout du EVA score (calcul EVA score sans MIT dans crisporEffScores.py
 puis ajout du MIT dans crispor.py avec calcEVAscore() / mergeGuideInfo

## multiseq mode : 
- modification de microhompage() ->  recherche du pamId en fonction de l'exon dans lequel se trouve le pam.
- modification de saveOutcomeData() (pour scores oof / lindel), n'efface plus les données précédentes si la fonction est exécutée dans une boucle
- interface KO : OK. (reste à afficher les exons bout à bout ?)

## à faire / bugs
- finaliser l'affichage du mode KO
- ajouter description pour EVA score et vérifier le calcul (notamment structurenumber=1)
- lorsque un batch multiseq est relancé, perte du paramètre "exonSeqs".
- ajouter / finir les descriptions des fonctions

# 21/01/26

## divers

- dans calcFreeEnergyRNAStructure() : calcul de la mfe structure uniquement
- comparaison calcul EVA scores avec données de Riesenberg et al. 2025 (OK, quelques différences (négligeables ?), probablement dues aux arrondis, à tester sur + d'exemples)
- comparaison énergie libre entre RNAstructure et RNAfold : OK

## KO / multiseq mode

- modification de primerDetailsPage() -> obtention de la séquence / exonId depuis les données du batch.
- ajout d'options pour le type de knock-out : frameshift dans la séquence codante ou excision du locus : 
	- si excision du locus -> sélection de la taille de la région target n bp upstream/downstream du TSS/TES
	
## à faire
- limiter la recherche upstream / downstream dans les limites des coordonnées chromosomiques.
- ajouter target promoteur ?
- ajouter liens ucsc exons
- internal server error sur lien téléchargement données guides.

# 22/01/26

## mode knock-out

- ajout de liens vers le browser, par exon.
- modification de parseMultiSearchInfo : retourne optionellement les données des guides (pour call dans downloadFile() )
- correction d'un bug dans getFirstThird() (allongement le dernier exon) 
- affichage du gene model : 
	- ajout de getGeneModel (call dans getExonsFromId() ) : retourne geneModel : liste de la longueur des exons / introns (uniquement entre CDS start / end)
	- ajout de geneModel dans params de newMultiBatch()
	- ajout de printGeneModel() affiche les exons sous forme de bloc, et introns sous forme de lignes
		- coloration des exons target
		- si dernier exon tronqué, coloration de la partie target
## à faire
(suggestions JP)
- option : masquer l'affichage seq / PAM des exons DONE
- cliquer sur un exon du gene model affiche la seq / pams de l'exon DONE
- rechercher guides dans seq exon étendue (exon=uppercase, intron=lowercase) : pam 6+ bp splice site prio ++ ~DONE 
- dans tableau effscores :
	- filtrer et trier par exon / prio ++ DONE
	- ajouter lindel / oof dans global score ? / prio -
	- supprimer doench 2016 et scores optionneles / prio ++ DONE
	- sélection du effscore pour calcul du score global (formulaire avec dropdown menu sous header score global) / prio + DONE
	- ajouter titres (indiquer expType / geneID, etc..) / prio ++ DONE
	- vérifier la page clonage (warning si pam =/= spCas9 -> tracr diff.) et ajouter champ input seq / prio -

# 23/01/26

## knock-out mode 

- Correction de l'affichage des exons target dans le gene model
- sur le gene model : cliquer sur un des exons target affiche  la séquence / pams de l'exon correspondant (+ bouton permettant d'afficher tous les exons)
- ajout d'un titre à la page de résultats correspondant au type d'expérience, avec un lien NCBI vers le transcript id.
- correction erreur lorsque la page est rechargée:
	- si le json est déjà présent, newMultiBatch() ne l'overwrite pas
	- si le fichier offtarget et effscores et déjà présent, processMultiSeqSubmission ne l'overwrite pas et retourne le nom des fichiers

## à faire
(visio JP / Max)
- ajouter option production du guide (transcrit / synthétique ..) -> effscore par défault / prio +
- recherche .gp par altname (symbol) afficher tous les transcrits correspondant (et la taille de la prot.)  DONE
- pour l'homme : surligner MANE vs Basic + lien vers browser gtex (fq utilisation des exons / tissu) / prio -
- ajouter mode CRISPRa : recherche autour du TxStart -> retrouver guides dans données Broad Institue (Pooled lib) (stockage dans database) / prio --
- dans processmultiseqsubmission() -> retoruner Fname si Fname existe, pour éviter crash lors du rechargement de la page. DONE

(bugs)
- la numérotation des exons ne prend pas en compte les exons situés entièrement en 5' et 3' UTR (OK)
- sur le formulaire knock-out, cliquer sur "submit" ne redirige pas vers une nouvelle page 

# 26/01/26

## knock-out mode 

- extension de la séquence des exons à GUIDELEN - 6bp (pour que le site de coupure soit à +6bp du site d'épissage)
- modification de findPams() afin d'éviter de rechercher des pams en dehors (dans l'extension) de l'exon. 
- lorsqu'un exon est sélectionné, la fenêtre affichant la séquence et les pams correspndant est agrandie.
- suppression des scores optionnels et changement de l'ordre (rs3 / EVA / Mor-Mateos)
- Possibilité de sélectionner le score d'efficacité utilisé pour le calcul du global score (modifications dans printTableHead() )

## à faire
- ajouter variants dans parseAndPrintMultiSearchInfo() DONE
- certains transcrits (ex ENST00000591702.1_1 ne sont pas présents sur ucsc / ENSEMBL database)

# 27/01/26

## knock-out mode

- affichage des variants : déplacement de l'obtention de l'information des variants dans une fontion dédiée (getVariants() )
- correction de printGeneModel() : les exons pour lesquels aucun guide n'a été trouvé ne sont pas sélectionables.

## à faire

- MANE : ucsc tools / tableBrowser -> gencode (track allGENCODE V49 (comprehensive + MANE) + ncbiRefSeq (all +  select) -> classer ordre (ou surligne / annoter)
- afficher "coding exons" au lieu de "exon"
- afficher dropdown variants (en dehors de la bouche exons)
- dans printKoForm() : empêcher la ré-impression du formulaire lorsque submit=submit

# 28/01/26

## divers

- les formulaires knock-out / knock-in ne sont plus réimprimés après soumission

## knock-out mode

- ajout d'un message affichant sur quel exon les données sont filtrés
- ajout d'un mouseover affichant le n° de l'exon si celui-ci est trop petit pour l'afficher directement 

## knock-in mode

- correction de getDonorSeq
- ajout de la variable globale pamLists (dict) : contient les listes de PAMs à utiliser pour multipam job. key = params["multipam"]
- création de processmultiPamSubmission() : dans une boucle, écrit les effscores / offtargets pour plusieurs PAMs (1 effscore / pam pour l'instant)
- ajout de newMultiPamBatch : sauvegarde les paramètres nécéssaires au job "multipam" dans json
- modification de crisprSearch(), readBatchParams() 
- ajout de parseAndPrintMultiPamInfo() : récupère données offtargets / effscores -> affichage (à compléter)

## à faire

- données OK -> corriger l'affichage (tableau + seq)
- ajouter la distance PAM / site d'insertion (sauvegarder insertpos dans params)

# 29/01/26

## divers

- correction d'un bug dans crisprSearch() (execution du bloc dédié au mode signle seq / pam après affichage des résultats)

## knock-in mode

- corrections de bugs dans processMultiPamSubmission() et parseAndPrintMultiPamInfo() : les varables globales n'étaient pas redéfinies correctement pour chaque PAM.

## bugs

- GUIDELEN = 20 pour tous les PAMs dans le fichier offtarget, ce qui casse annotateOffTargets() ensuite. Pourtant, GUIDELEN est normalement redéfinie juste avant findOffTargetsBwa().

# 30/21/26

## divers

- corrections d'erreurs de syntaxe (à terminer)
- correction de writePamFlank() -> pass de pamFullName dans flankSeqIter() pour redifinir GUIDELEN

## mode kock-in

- modification de createMultiBatchEfffScoresTable() et calcMultiSaveEffScores() : sauvegarde de tous les scores disponibles,
si un score ne correspond pas au PAM traité -> 0 (si le score est un string, impossible de trier le tableau)

## à faire

- implémenter staggered cut pour eSpOT-ON (pam NGG-22) 
- travailler sur l'affichage knock-in
- modifier le formulaire knock-in (input seq -> desired seq)

