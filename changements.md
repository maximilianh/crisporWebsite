# changes in CRISPOR - assistant mode

# 17/12/25 ;

ajout du mode "assistant"

# 22/12/25 :

## corrections dans getGeneSeq():
 
cette fonction retourne deux listes contenant la séquence des exons de tous les transcrits d'un gène, dans les fichiers .gp disponibles.
- une liste ('all.exons')
- une liste ('first.exons') correspondant aux exons du premiers tiers du CDS

## changements depuis le 17/12/25 :
 
- retourne la liste des exons dans le bon ordre pour les gènes du brin "-".
- les exons de plus de 2300bp (MAXSEQLEN) sont scidés en plusieurs séquences.
- retrait des 5' et 3'UTRs.
- Un exon trop long pour entrer dans le premier tiers du CDS mais pouvant être tronqué l'est désormais.
- si le CDS est de < 100bp, tous les exons sont considérés (seuil à augmenter ?)

## problèmes / à faire :

- pas d'annotation des transcripts ID dans le fichier .gp -> les trancrits sont pour l'instant numérotés selon leur ordre d'appariation.
- prise en compte des autres fichier d'annotation ? (genes.tsv, .bb, segments.bed) ? 

- Ajout du menu déroulant de la liste des gènes disponible en fonction de l'organisme sélectionné
via AJAX avec select2
chargement un peu lent (~2s), à voir avec d'autres génomes.

# 23/12/25

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

# 30/12/25

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

# 09/01/26

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

# 12/01/26

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

## gestion multi seq / PAM

- recherche pour plusieurs pams : modification de processMultiSubmission()
	- écriture d'un fichier bed vide
	- dans une boucle : 
		- call de setupPamInfo()
		- writeMultiFasta() (recherche des guides pour chaque séquence en input)
		- findOffTargetsBwa() >> ajout des résultats dans le fichier bed
 	- index des guides / offtargets: PAM.seqId.exonId.pamId
- réduction du maximum de mismatchs à 3 lors de la recherche des offtargets (~5-10x plus rapide)

- call des effscores dans la boucle avec createMultiBatchEffScoreTable()

## bugs / à faire
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

## à faire

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

# 19/01/25

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
- corriger geneModel

# 02/02/26

## mode knock-in

- modification de parseAndPrintMultiPamInfo() : aggrège les données offtargets / effscores dans une boucle -> showSeqAndPams / showGuideTable avec données aggrégées
- modification de showSeqAndPams() : mode "multipam" optionnel : crée les pamLines dans une boucle, pour chaque pam
- modification de flankSeqIter() : séparation en deux fonctions
    - getPamLines() : génère les motifs à afficher sur la ligne, pour un PAM (exécuté dans un boucle)
    - layoutPamLines : place des motifs sur la séquence
    - conservation de flankSeqIter(), qui exécute les deux fonctions ci-dessus
 
- modification de showGuideTable() : ajout du PAM correspondant au guide + description PAM / enzyme (en mouseover) dans la première colonne, changement des '0' en '--' pour les scores non calculés
- ajout de showDonor() : affiche la séquence de l'ADN donneur (non modifié) + bouton copy to clipboard
(format des résultats OK, reste à ajuster l'affichage)

## à faire

- améliorer l'affichage (titre de la page, etc..) ~ DONE
- établir la liste des scores à utiliser
- créer une fonction pour obtenir la distance site de coupure <-> site d'insertion et afficher ces distances dans le tableau (+ prendre en compte pour global score) DONE
- en mode knock-in, le choix de l'effscore pour la calcul du global score ne fonctionne plus

# 03/02/26

## mode knock-in

réécriture de getDonorSeq() : séparation en deux fonctions
    - getTargetSeq : retourne les coordonnées (si user input = geneID) ou séquence target + position d'insertion (exécution avant crisprSearch()
    - writeDonorSeq : retourne séquence + coordonnées et écrit la séquence de l'ADN donneur dans les batch params (éxécution par le worker)
- getSeq() et findPerfectMatch() sont désormais exécutés par le worker (dans writeDonorSeq() )
- écriture de la position du site d'insertion, de la longueur des bras d'homologie et de la séquence d'insert dans les batch params
- dans newMultiPamBatch() : batchId unique en fonction de la position du site d'insertion et de la longueur des bras d'homologie

- affichage du site d'insertion +- 23bp sur la fenêtre séquence / PAMs
- affichage de la distance entre le site d'insertion et le site de coupure dans le tableau des guides
- prise en compte de cette distance dans le calcul du score global (-40 si distance > longueur du guide + PAM)

# 04/02/26

## divers

- en mode knock-out, affichage par défaut = séquence de l'exon 1 / tous les guides dans le tableau
- les trois guides non chevauchants au plus haut global score sont surlignés en vert dans le tableau
(ajout de flagBestGuides, call dans showGuideTable() )
- correction description global score


## mode knock-in

- correction de l'affichage: ajout titre avec lien vers NCBI/ENSEMBL si geneID, + autres détails
- ajout de la liste de PAMs + message indiquant où trouver l'info sur l'enzyme correspondante
- tri des scores affichés : rs3 / EVA / crisprScan / deepCpf1 / Najm (saCas9)

## à faire : 

- suggestions JP:
- ko :
    - afficher domaines protéine sur le gene model ? si faisable
    - mettre en surbrillance les 3 meilleurs guides non chevauchants
- ki : 
    - input = séquence départ / voulue (pour permettre substitutions et remplacement)
    - éventuellement, prédiction de la structure de la protéine ? (pour insertion =/= extrémités)
    - ne pas modifier le global score
    - par défaut : classer les guides en fonction de la distance site de coupue / site d'insertion
    - par défaut : recherche avec NGG, afficher les guides à +-20 bp du site d'insertion. 
    - si pas de guide avec NGG -> rediriger vers une nouvelle recherche
    - faire un schéma ADN donneur / site d'insertion au lieu d'afficher la séquence du donneur
    - si pas de chevauchement entre les 15 bp de l'extrémité du guide et le site d'insertion -> flag recodage pour éviter coupure du donneur

## bugs

- le score SaCas9 n'est jamais calculé

# 05/02/26

## divers

- ajout d'un bandeau en haut du menu principal : sélection de trois modes
    - classic
    - knock-out
    - knock-in
- transfert de assistanDispatcher() dans printBody()
- réorganisation des formulaires knock-out et knock-in

## mode knock-in 

- tri du tableau des guides par la distance site de coupure / site d'insertion (par défaut) + correction calcul distance
- modification de calcInsertDistance() : retourne doRecoding = True si les 15bp à l'extrémité du guide ne chevauchent pas le site d'insertion
- ajout d'un avertissement si  doRecoding = True.

# 06/02/26

## divers

- correction de readAnnGenomes()
- amélioration de l'affichage du bandeau pour le choix du mode
- amélioration de l'affichage des formulaires knock-out / knock-in 
- correction du mode "ko by excision of the gene locus" : 
    - affiche downstream + upstream region
    - affiche taille de la délétion
    - faut-il séparer le tableau en deux parties (upstream / downstream) ?
- correction du surlignage des guides : 
    - déplacement de flagBestGuides() dans ShowGuideTable() (évite double itération)
    - surlignage de toutes les cases du tableau
    - définition de la région occuppée par la RNP : 10bp en amont du PAM, 3bp en aval du guide.

## à faire

- biblio "competition" entre guides (cf Anders et al. 2014)
- dans le formulaire knock-in : custom séquence et insert par défault
 
# 09/02/26

## divers

- ajout du paramètre optionnel "maxlen=True" dans getSeq() : retire longueur maximum si False (pour ADN donneur > 2300 bp)
- reformatage de crispor.py avec black

# 10/02/26

## divers

- corrections d'erreurs html

## mode knock-in 

- choix de la longeur des bras d'homologie (+ donneur db/sb) -> step 3
- steps 4 et 5:
    - soit séquence de départ / séquence éditée
    - soit sélection geneID -> protein tagging en Nter ou Cter

# 11/02/26

## mode knock-in 

- réorganisaiton de printBody() pour expType == "ki" :
    - ajout de processCustomInsertSeq() depuis séquence de départ / séquence d'arrivée, extrait le type de knock-in (insertion, remplacement, substituion),
        la position d'insertion et la séquence d'insert (à compléter)

- modification de showDonor() : représentation du donneur par des lignes (+ longueur des bras d'homologie / séquence d'insert) au lieu de la séquence 

# 12/02/26

## mode knock-in 
- modification de processCustomInsertSeq() : 
    - prise en compte des modifications en uppercase 
    - retourne type de knock-in, position d'insertion, séquence d'insertion 

- dans formulaire knock-in :
    - les modifications sont indiquées en uppercase par l'utilisateur
    - ajout de deux boutons permettant de convertir la sélection en uppercase / lowercase 
## à faire

- stratégie knock-out : retrait d'un exon (en dehord du cardre de lecture) (précalculer oof exons)
- adapter pénalités global score au choix de l'effscore choisi
- ajouter liens vers CasPedia (http://caspedia.org/) 

# 13/02/26

## mode knock-in

- changement range longueur bras d'homologie : 50 -> 2000bp
- ajout de messages d'avertissement pour les type d'édition non supportés (remplacement, insertions multiples..)
- support du mode substitution : modification de writeDonorSeq(), parseAndPrintMultiPamInfo() et showSeqAndPams() : affiche le type de substitution (base WT -> base Edit) dans le titre, affichage de la substitution sur la séquence target.
- dans mode protein Tagging : ajout du mode "qTAG" (cf https://doi.org/10.1038/s44318-024-00337-5) :
    - contruction de la cassette qTAG
- déplacement de la séquence des linkers, tags et markers dans une variable globale
- stockage du nom et de l'ordre des séquences de tagging dans les paramères de batch (pour affichage sur geneModel)
- création d'un batchId unique pour différentes séquences d'insert
## à faire

-  liste PAMS : NGG > dispo commercialement > addGene > cas naturelles / cas engeneered (activité + faible + offtargets)
- ajouter qTAGs dans la liste + ajouter interface pour construire la cassette complète
- si mode "protein tagging" : afficher geneModel avec représentation de l'insertion (nom tag + couleur)
- dans le mode substitution -> ajouter l'acide aminé correspondant à la substitution dans le gene model

# bugs

- si un geneID n'est pas sélectionné, une page vide s'ouvre (devrait rediriger vers le formulaire)
- bug sélection EF1alpha

# 16/02/26

## mode knock-in

- réinitialisation de la sélection tag/linkers lors du changement de menu
- simplification de calcInsertDistance() + modification de la fenêtre de recodage si substitution
- si subtitution, recodage uniquement si la séquence du pam est altérée

## bugs

- double impression du formulaire ki 

# 14/02/26

## général

- retrait de l'autocorrection des éléments textarea

## mode knock-in

- correction de l'affichage qTAG
- modification de getInsertSeq() : pas de recodage si la séquence du pam est modifiée en mode insertion
- si insertion dans la séquence du pam mais que l'insert recrée le pam -> recodage (à adapter pour tous les types de pams)
- ajout de donorDesignPage() : lien depuis le tableau des guides, redirige vers un formulaire pour le design de l'ADN donneur

## notes Schubert et al (design ssODN)

- taille min bras d'homologie : 20-30nt

- amélioration efficacité HDR : 
    - donneur asymétrique
    - ajout phosphorothioate
    - double nicking : distance optimale ~40-68nt (dépend de l'enzyme)  
        - permet de réduire l'impact de la distance site de coupure / site d'insertion sur l'efficacité
        - si site d'insertion à mi-distance du site de coupure des deux guides, efficacité > à WT Cas9

- choix design : 
    - 200nt max (insert 160nt max)
    - proposer template sur le brin target ou non target (préférence varie en fonction lu locus / type cell)
    - 

- choix du guide:
    - efficacité > distance pour guides dont le site de coupe est ~< 10-14bp du site d'insertion.

- recodage:
    - pas de préférence transition / transversion
    - recodage PAM efficace, en particulier si insertion vers extrémité 3'
    - recodage guide en 3' efficace
    - 2 mutations suffisantes (3/4 n'améliorent pas l'efficacité)
    - si trop de mutations -> perte d'efficacité HDR

## bugs

- prendre en compte de global score pour mutated guides (dans page cloning)

# 18/02/26

## global

- ajout d'un message d'avertissement sur la page cloning / primers : tracr varie en fonction de l'enzyme

## knock-in mode

- correction du bug de la page clonig / primers en mode knock-in
- affichage du geneModel (si applicable) dans la page du design ADN donneur + affichage de l'insertion sur le geneModel
- ajout du mode délétion : détection de la délétion dans la séquence d'input et modification de l'affichage des résultats
- changement de la fenêtre de recodage en mode délétion
- correction bug : décalage de l'affichage de la base substituée si subtitution proche des extrémités de la séquence

## notes

- design donneur
    - proposer la séquence le donneur sous forme de plasmide ?
        - insérer target sgRNA de part et d'autre du donneur
    - ou ajout de biotine, amino-dT ou carbon spacers en 5'
        - évite multimérisation -> augmente efficacité (https://doi.org/10.1016/j.omtn.2024.102344) 

## à faire

- intersecion fenêtre recodage / fenêtre délétion mal calculée

# 19/02/26

- docker reset sur config root au démarrage! -> dockerd-rootless-setuptool.sh

## knock-in mode

- correction : en mode délétion, les guides dont le site de coupure est situé dans la délétion ont une insertDistance de 0
- correction : le changement de l'effscore pour le calcul du global score fonctionne en mode knock-in
- correction du chevauchement fenêtre de recodage du guide avec intervalle de la délétion

## notes

- structure prot -> préférence insertion Nter / Cter
- recodage : effets sur l'épissage
- input sq protéique (pour tagging)
- alphaGenome / VEP (ensembl)
- recodage entre site d'insertion / site de coupure
- ko par délétion du promoteur

# 20/02/26

## divers

- retrait de l'input "longueur des bras d'homologie" dans le formulaire knock-in

## knock-in mode

- ajout du formulaire pour design donor DNA
    - sélection ssODN / dsODN
    - si ssODN : 
        - choix du brin template
        - réduction de la longueur max des bras d'homologie (20 -> 200 - taille insert)
        - ajustement dynamique de la longueur des bras d'homologie afin que la taille max du ssODN de dépasse pas 200
    - sélection un nombre de design
    - sélection de la longueur des bras d'homologie (5' et 3')
    - options de recodage (PAM / guide / entre coupure-insertion / homopolymères / GC rich / repeat)

## à faire

- affichage gene model -> sélection longueur bras d'homologie sur gene model (éventuellement)
- enlever sélection PAM dans tableHead en mode knock-in
- option : recode only in coding seq
- "Recode between the cut site and insertion site" : reformuler
- symbole warning gif
- afficher lien ucsc correspondnat au bras d'homologie
- si pas de geneID -> annoter sq avec exons : ne pas recoder 5'UTR et splicing sites
- minVal à -564 pour si ssODN sélectionné ?!

# 23/02/26

# divers

- correction bug téléchargement des données en mode ko / ki
    - dans downloadFile() : parseAndPrintMulti(pam/seq)Info retourne les données nécéssaires au téléchargement
    - modification de xlsWrite(), iterGuideRows(), makeGuideHeaders() et intToExtPamId()
    - prise en compte des pamIds avec préfixe et ajout du global score / distance site de coupure - site d'insertion

## mode knock-in

- corrections dans writeDonorSeq() : 
    - prise en compte du mode délétion + recherche des bras d'homologie 5' et 3' indépendamment
    - choix du brin target / non-target
- modification de showDonor() : affiche la séquence du donneur en FASTA (bras d'homologie en lowercase, insert en uppercase)

# 24/02/26

## divers

- déplacement des binaires RNAstructure dans le dossier approprié (bin/Linux-x86_64)
- utilisation de RNAfold pour le calcul du EVA score (temporaire ? / pour permettre à crisporTest de fonctionner)

## mode knock-in

- ajout de tooltips dans donorDesignPage()
- changement de la sélection du brin modèle pour ssODN, selon Paix et al., 2017
- affichage de la distance / direction site de coupure / site d'insertion 

## bugs

- internal server error si downloadFile() est éxécuté en mode classic

# 25/02/26

## divers

- correction (définitive ?) de downloadFile()
- ajout de showSecondaryStructure() : exécute RNAfold et affiche une représentation graphique de la structure secondaire (lien depuis "Guide sequence + PAM")
- ajout du lien vers le genome browser en mode knock-in

## à faire
- lien -> vérification épissage
- masquer barcode pour substitutions
- ajout option donorName
- ajouter remplacement : définir taille remplacement -> 1 ou 2 guides (cf. double nicking)
- dans input "edited sequence" -> option pour afficher les 3 cadres de lecture, puis sélection 
- dans showDonor() -> ajouter séquence du guide (+ guides compatibles) + récap but de l'expérience
- ajouter bouton téléchargement fasta

# 26/02/26

## divers

- correction de la page off-target PCR (prise en compte du pamId complet en mode multipam)
- ajout sélection de la longueur du primer (15-30 bases)

## knock-in mode

- dans donorDesignPage() : 
    - ajout input fasta header pour ADN donneur
    - masquage des options "trim the donor to facilitate its synthesis" pour ssODN
    - masquage de l'option "barcode" en mode substitution
    - si distance site d'insertion / site de coupure < 10bp : sélection du brin target pour template ssOND

- dans showDonor() :
    - ajout de la séquence du guide + résumé de l'expérience
    - ajout d'un lien de retour vers la page de design du donneur (modification de printBackLink())
- correction d'un bug dans writeDonorSeq() : polarité assignée au donneur double brin + correcion coordonnées substitution

## à faire

- dans processMultiPamSubmission retourne une erreur si la séquence n'est pas dans le génomes -> afficher message à la place
- si rechargement de la page off-target primers -> "primers not found at this Tm"

# 02/03/26

## divers

- les en-têtes du tableau sont figées lors du défilement de la page (suggestion d'Axel Benchetrit)
- dans showSecondaryStructure : ajout d'un formulaire pour définir le température + extension de la séquence du guide en 3'
- correction du bug dans otPrimerPage() : ampLen, Tm et primerLen = None dans printHiddenFields()

## knock-in mode

- dans showDonor() : ajout d'un bouton pour le téléchargement des séquences du guide + ADN donneur en fasta

## à faire (avant lancement de la phase de test)

### global

- vérifier le workflow
- vérifier tous les textes (y compris tooltips!)
- vérifier qu'aucune information importante ne manque en mode knock-out / knock-in :
    - affichage des variants
    - affichage des coordonnées ucsc DONE
- rendre tableau scrollable au lieu de figer les en-têtes DONE
- corriger surlignage du PAM sur la séquence si sélection depuis le tableau
- "faux" off targets en fonction de l'assemblage
- ajouter message pour type d'exp sélectionné (surbrillance + explication)
- description au dessus -> sélection

- ajouter génomes lignées cellulaires

### knock-out

- ajouter "knock-out par délétion du promoteur" DONE
- ajouter surlignage du MANE transcript (cf. notes 27/01/26)
- ajouter affichage getGeneModels()
- alerter si un gène essentiel est ciblé (humain / souris)
- surligner ATG in frame dans showSeqAndPams() DONE

### knock-in

- créer les listes de PAMs OK + lien CasPedia
    - commercially available recombinant nucleases
        - spCas9 (NGG)
        - Cpf1 
        - Cas12f
        - saCas9
    - expression plasmid available from addGene
        (to be completed)
    - pam less / engeneered pam variants with low specificity
- + custom list (max 5) 

- obtenir les séquences des tags / linkers listés dans taggingSeqs
- ajouter sq tags 3FLAGSBP / SBP3FLAG OK
- annoter la région de l'ADN donneur avec séquences codantes / 5'UTR OK
    - soit depuis fichier genePred
    - soit en convertissant genePred en bigBed -> getGeneModels() OK
- implémenter recodage OK
- implémenter le mode "remplacement" (à préciser)
- ajouter custom track pour les bras d'homologie

# 03/03/26

## global

- Global score :
    - les pénalités varient en fonction du type de delivery (transcription in vivo / vitro / synthèse)
    - Description du calcul du global score plus détaillée (ajouter refs)

## bugs

- Tm 64 par défaut sur on target PCR 
- primerLen 22 par défaut 

## à faire

- trackHub : customisation couleur / épaisseur
- ou customtrack
- lignées cell : cellosaurus
- automatisation submission genome
- formulaire : tableau paramètres -> stock db -> exécution scripts

# 04/03/26

## global

- correction de otPrimerPage(), primerDetailsPage(), printValidationPcrSection() et designOffTargetPrimers()
- merge des changements dans la classe queue + adaptation de runQueueWorker()
- dans crisporTest : séparation du répertoire temp / crisporJobs.db avec la version publique

## knock-in mode 

- call de getExonInfo() dans donorDesignPage() pour annoter l'ADN donneur

## bugs 

- en mode classic, "guideScores" referenced before assignment dans crispor.py (en + du message "not found in genome")

## à faire

- lors de la sélection d'un PAM, center la rangée correspondante (sinon le header la masque)

# 05/03/26

## knock-in mode

- ajout de getArmCoords() : pour chaque bras d'homologie, obtention des coordonnées
    - du PAM
    - des bases PAM-proximales du spacer
    - de la région entre le sitre d'édition et le site de coupure
    - des exons
- ajout de recodeDonor(): 
    - crée une liste des bases à muter, puis une liste des codons recoupant ces positions
    - pour chaque codon, crée une liste des codons synonymes
    - pour chaque base de chaque codon synonyme, vérifie si cette celle-ci correspond à une base à muter
    - si oui, le codon synonyme est conservé, et est remplacé sur le bras d'homologie (base mutée en uppercase)

## à faire

- en mode knock-in, adapter la liste des effscores à la liste de pams sélectionnée

## bugs

- en mode knock-in, "show all" dans colonne off-target ne fonctionne pas OK

# 06/03/26

## knock-in mode 

- dans recodeDonor() : dans les coordonnées du PAM, retrait des positions correspondant aux bases "N" du motif

- ajout de codonFrequency.py dans /tools/usrLocalBin :
    - pour chaque génome dans le répertoire /genomes, calcule l'occurence de chaque codon à paritr du fichier genePred
    - mais :
        - pour l'instant : chaque transcrit est pris en compte, donc biais pour gènes à épissage alternatif ++ ? 
        - si prise en compte de tous les exons codants par gène -> risque de décalage du cardre de lecture ?
        - si prise en compte d'un seul transcrit, lequel choisir ?
        - script lent
    - -> utilisation d'un outil dédié ?

# 09/03/26

## knock-in mode

- ajout des listes de pams :
    - "commercially available nucleases",
    - "engineered pam variants with low specificity",
    - "expression plasmid available from addGene"

- codonFrequency.py : possiblité d'utiliser EMBOSS cusp (+ rapide)
    - mais requiert input CDS
    - donc, soit décompte des codons pour tous les transcrits, soit pour le transcrit le + long (ou MANE ou humain)

# 10/03/26

- récupération/réécriture des données après bug git dans VScode (git reset au commit 06/03/26)
- tentative de centrer la rangée correspondant au guide sélectionné depuis le sequence viewer

## knock-in mode 

- finalisation de codonFrequency.py : calcul de la fréquence des codons pour tous les transcrits
    - valeurs hg19 comparables aux données en lignes +- 0.05 -> OK pour classement des codons selon la fréquence ?

# 11/03/26

## knock-in mode

- correction de writeDonorSeq et getArmCoords : 
    - prise en compte de la polarité pour les ssODN
    - simplification de la conversion des coordonnées en fonction du brin

- dans recodeDonor : prise en compte de la fréquence d'utilisation des codon (lecture JSON) : 
    - recodage du donneur jusqu'au codon ayant la fréquence d'utilisation la plus élevée

- surlignage des guides chevauchant le site de coupure : option "pamIdToRecode" dans mergeGuideInfo(), showSeqAndPams() et makePamLines()

## à faire

    - ajouter checkbox affichage enzymes
    - couleur si freeEnergy > -5 kcal/mol
    - warning si la région sélectionnée n'a pas pu être recodée OK
    - pour ssODN : prise en compte des coordonnées en dehors des bras d'homologie! OK

# 12/03/26

## knock-in mode

- dans le tableau des guides : 
    - retrait de l'affichage des enzymes
    - retrait des checkbox filtrage selon séquence
    - rertait du message "secondary structure" -> ajout d'un /!\ ou d'une gomette verte à côté du lien "predicted secondary structure"

- correction bug : les coordonnées en dehors des bras d'homologie ne sont plus prise en compte lors du recodage des ssODN
- recodeDonor fonctionne sans fichier de fréquence des codons
- ajout de  printMutEventsTable() : affichage d'un tableau des mutations apportées à l'ADN donneur
- surlignage de la séquence d'insert en jaune dans showDonor()

## à faire

- refaire le sticky header : rendre le tableau scrollable ?
- afficher option "introduce a mutation to check for homozygous editing" après affichage donneur OK
- indiquer lorsque le recodage n'est pas possible OK
- ajouter un lien pour Gibson assembly / golden gate
- dans codonFrequency, sélectionner le transcrit le plus long OK
- explication détaillées dans donorDesignPage()
- recodage dans régions non codantes (à part 5 bases en amont des exons) OK
- ne pas recoder si création d'un site d'épissage
- dans processCustomInsertSeq(), gérer newlines et espaces (+ format fasta) dans textarea OK

# 13/03/26

## knock-in mode

- correction recodage pour ssODN à polarité inverse : revComp() après écriture et recodage du donneur
- modification de getArmCoords() : écriture des coordonnées des codons, 5'UTR, 3'UTR et sites d'épissage uniquement en amont / aval du site d'édition 
- dans recodeDonor() les coordonnées 5'UTR et sites d'épissage ne sont pas recodées + avertissement dans printMutEventsTable
- ajout du recodage dans régions non codantes : introduction d'une transition toutes les 3pb
    - mais : à quelle fréquence recoder et quelles mutation introduire ??
- dans formulaire knock-in : retrait des espaces, tabs et newlines dans input séquence

# 17/03/26

## knock-in mode

- modification de codonFrequency.py : getExonPos ne retourne que le transcrit le plus long (pour un gene Symbol donné)
- correction des coordonnées des mutations dans printMutEventsTable() pour ssODN à polarité inverse
- surlignage des régions de 20bp à +80%GC et des homopolymères dans showDonor() : modification de getHighlightedRow
- par défaut, seuls les guides provoquant une DSB à moins de 10bp du site d'édition sont affichés. ajout d'un curseur pour modifier cette valeur

## à faire

- autoriser recodage dans 5'UTR sauf au niveau de la séquence consensus kozak (ATG-5bp) OK
- ajouter un dropdown pour afficher les autres motifs de PAMs sur le sequence viewer (puis ajouter lien vers nouvelle recherche avec ce PAM) OK

# 18/03/26

- ajout des coordonnées de la séquence consensus kozak dans getArmCoords()
- dans recodeDonor(), retrait des coordonnées 5' et 3'UTR des sites d'épissage
- correction d'un bug lors du calcul des coordonnées codon START + affichage d'un message sur codon START dans tableau
- ajout de l'option "show the position of other PAMs on the sequence" dans la page de résultats Ki :
    - modification de showSeqAndPams() et makePamLines()
- dans downloadDonor() : si recodage, séquence du donneur non recodé + recodé

## à faire / bugs 

- surligner séquence guide  + PAM dans showDonor() OK
- les coordonnées 5'UTR ne sont pas toujours annotées si l'input est une séquence
- suggestion de Geneviève Tavernier : sur la page off-target primers, n'afficher que les off-targets correspondants si un filtre (exon / chr) a été coché OK

# 19/03/26

## knock-in mode

- dans showDonor() : surlignage de la séquence du PAM + spacer. Si guide /pam chevauchant le site d'édition, surlignage de la séquence de part et d'autre
- dans printTableHead: si filtre onlyExons / onlyChrom cochés, écriture des params correspondants dans le lien "off target primers"
- dans otPrimerPage() : filtre des primers en fonction de onlyExons / onlyChrom + affichage d'un message en fonction du filtre appliqué
- correction bug : en mode ko / ki, "show all" dans colonne off-targets du tableau refonctionne (pas de "." en CSS)

# 20/03/26

## knock-out mode 

- ajout des méthodes "ko par délétion du promoteur" et "ko par interférence avec l'épissage" dans le formulaire knock-out

## notes JP :

- phase si exons DONE
- chemical synthesis par défaut DONE
- liste !!
- PAM highlight blue pas dispo sur serveur OK normalement
- scroll tableau vertical DONE
- liste addGene : Cas naturelles / artificielles
- reset -> PAM sélectioné only DONE
- cliquer sur la séquence du PAM ajouté lance une recherche avec l'enzyme correspondante 
- ne pas afficher homopolymers / GC rich pour ssODN DONE
- copy seq to clipboard ne fonctionne pas sur le serveur OK via https
- insert seq surlignée -> substitued Nt DONE
- dowload au format txt seq + récap + hyperlien ucsc + infos supp

## notes Max :

- inverser ordre phases sur brin inverse DONE (exonFrame relatif à la séquence)
- mouseover pam en noir DONE
- options sous seq DONE
- annoter repeats avec repeatmasker DONE
- underline / italique DONE
- dictionnaire annotation (download)
- ajouter label pour PAMs non calculés DONE
- tableau résolution : fit width = taille minimale DONE
- si div scroll -> parents avec width: 100%
- MPC dans firefox -> ajout agent IA

# 23/03/26

## général

- correction du bug Keras / deepCpf1 : réinstallation de keras 2.12.0 / tensorflow 2.12.0, puis réinstallation de numpy 1.26.4
- dans calcGlobScore : pour guides Cas12a, seule l'efficacité est utilisée pour le calcul (aucun score de spécificité disponible)

## knock-in mode

- déplacement des options d'affichage / filtre des PAMs sur la séquence en dessous de la séquence + formatage
- ajout label + mouseover pour PAMs supplémentaires
- annotation des régions répétées dans writeDonorSeq() et surlignage dans showDonor()
- adaptation du calcul des cordonnées du PAM + guide pour pour le recodage avec guides Cas12a

## à faire 

- corriger bug score saCas9 -> OK, le score n'est délibérément pas calculé

# 24/03/26

## général

- correction de l'affichage en mode classic pour guides Cas12a
- entrer un transcript ID en mode classic refonctionne
- modification des pénalités pour le calcul du global score ( + mise à jour description / ajout références):
    - Pénalité énergie libre augmentée pour guides synthétiques
    - Ajout des pénalités pour guides Cas12a
- en-tête du tableau des guides dans un tableau séparé à position fixe
- encapsulation des deux tableaux dans un div scrollable à largeur fixe

## knock-out mode

- affichage des effScores appropriés pour Cas12a / saCas9

## knock-in mode

- dans showDonor(), surlignage de la séquence d'insert
- ajsutement de l'affichage du filtre des PAMs et modification du titre du tableau en fonction

# 25/03/26

## général

- dans makeExonLines() : prise en compte du décalage entre phase de l'exon entier / exon dans la séquence

## knock-out mode

- ajouts dans le mode "KO par perturbation de l'épissage" :
    - recherche 10bp de part et d'autre des jonction des exons
    - sélection d'un exon ou recherche sur tous les exons

## knock-in mode 

- correction de bugs dans getArmCoords :
    - les exons terminant dans une délétion ou un remplacement ne sont plus pris en compte
    - les coordonnées des exons situés entièrement dans le bras d'homologie 3' sont correctement calculées
    - mauvaise assignation du codon START

- ajout du mode "remplacement" :
    - prise en compte des remplacement de < 10pb à une seule position dans processcustomInsertSeq()
    - affichage du remplacement dans showSeqAndPams()
    - prise en compte des coordonnées du remplacement dans getArmCoords(), recodeDonor() et showDonor()
- correction du bug lorsque la séquence d'input n'est pas trouvée dans le génome + affichage d'un message dans la page de résultats si tel est le cas

## à faire

- le menu "Options to modify the display of PAMs on the sequence viewer" a parfois une taille réduite ?
- lorsque KO par une paire de guides -> séparer le tableau en 2 (guide 1 / guide 2)
- en mode "KO par perturbation de l'épissage" -> indiquer ne n° des exons, le type de site (donneur / accepteur) et retirer 5'/3'UTR

# 26/03/26

## global

- correction des coordonnées des exons dans makeExonLines() : offset par 3 - exonFrame
- prise en compte du changement de phase dans trimExonAndFlip() si la séquence débute après l'exon (après flip des coordonnées si brin - !!)

## Knock-out mode

- ajout de la séquence traduite correspondant au transcriptID sélectioné

- mode "splicing" :
    - surlignage des sites donneurs / accepteurs d'épissage sur la séquence
    - ajustement de printGeneModel : affichage des séquences donneur / accepteur si un exon est sélectionné

## à faire

- corriger le formulaire de sélection du transcrit dans donorDesignPage()
- filtre exons / tableau de fonctionne plus

# 27/03/26

## global

- prise en compte du MANE transcript 

## knock-out mode

- dans printGeneSelection / dbsearchGene : recherche par gene symbol, puis dropdown avec optgroup des transcrits correspondants
    - problème : pour génomes avec gene symbol annotés / LOCXXXX : recherche longue
    - arrêt de dbSearch si +30 gene symbols ont été trouvés

- réparation du filtre du tableau lorsqu'un exon est sélecitonné
- ajustement des noms des régions dans séquence / tableau en mode "promoteur" / "épissage"
- modification de sortGuideData() pour trier le tableau / exon (ne fonctionne pas encore)

## knock-in mode 

- correction du formulaire de sélection d'un transcriptId dans showDonor + retrait des params dupliqués dans KiResultsPage()

# notes

- KO mode : sélection gene symbol puis transcript ID (par défaut, MANE select ou exons communs) DONE
- mouseover celphedia DONE
- schéma étapes KI
- ajouter PAM à pamlist DONE

# 30/03/26

## global 

- ajout mouseover Celphedia (= texte https://celphedia.eu/en/)

## Knock-out mode

- correction d'un bug du calcul des coordonnées en mode "délétion du promoteur" sur le brin +
- finalisation du filtre de sortGuideData() par exon

## knock-in mode 

- ajout d'un récapitulatif des étapes de l'expérience
- ajout des séquences des tags / markers qTAG (https://www.addgene.org/Laurence_Pelletier/)
- affichage des séquences d'insert similaire à snapGene (Claude) + ajustement clipping du texte

# 31/03/26

## global

- ajout de l'option -g dans codonFrequency.py : restriction du calcul aux génomes indiqués (ou "all" pour tous les génomes)
- ajout d'une largeur minimale dans printAssistant / printKoForm et printKiForm, pour ne pas altérer l'affichage à basses résolutions (mais page zoomée dans ce cas)
- ajustement du bandeau pour la séléction du mode
- correction d'un bug dans crisprSearch() : message d'erreur si séquence non trouvée dans le génome en mode classic

## knock-out mode

- pour les délétions via une paire de guides (excision et délétion du promoteur), séparation des résultats en 2 pages
- ajout d'un bouton pour afficher les résultats pour la région en amont ou en aval de la délétion
- modification de la fonction Javascript de filtre des résultats en fonction de l'exon sélectionné pour prise en compte du mode "splicing",
    où deux séquences doivent être affichées simultanément

## knock-in mode 

- déplacement du récap des étapes du KI dans un fonction dédiée (printKiSteps) -> étapes avec tooltips + lien pour revenir à l'étape précédente
- dans showDonor(), ajout d'un message pour qTAG : affichage des éléments sélectionnés et d'un lien vers la librairie de plasmides

## à faire 

- inverser ordre sites loxP pout qTAG en C-ter! + vérifier ordre de tous les tags
- lister les changements dans /doc/changes.html

# 01/04/26

## knock-out mode

- correction de la logique de getExonsFromId() en mode "splicing" + écriture de l'exon sélectionné dans batch params
   + correction dénomination des sites donneurs / accepteurs
- correction de downloadFile() : tous les type des fichiers fonctionnent
- correction des index (exonIds) dans le tableau -> tous 0-based : plus de mismatch entre tableau / séquence en mode "splicing"

## knock-in mode 

- lors du surlignage du PAM + spacer dans showDonor(), inversion des coodonnées pour ssODN à polarité inverse
- ajout de tooltips dans printKiSteps et donorDesignPage() + vérification des textes

# 02/04/26

## knock-in mode

- correction d'un bug dans getPosAndSeq / runQueueWorker : crash du job (attente infinie) pour gènes non-codants en mode "protein tagging".
    - dans ce cas, affichage d'un message d'avertissement sur la page des résultats.
- dans showDonor : ajout d'un formulaire pour modifier la séquence d'insert en fonction de l'expérience :
    - ajout manuel d'une nouvelle séquence
    - choix tags / linkers
    - choix d'une autre substitution
- par défaut, les reapeats / homopolymères / GC rich ne sont pas surlignées : ajout de checkbox pour les afficher

# 07/04/26

## knock-out mode

- ajout de la sélection des exons communs à tous les transcrits :
    - ajout du paramètre commonExons à dbsearchGene() : retourne symbol
    - dans getGenePos : intersection de tous les transcrits + cds et tss start/end

## knock-in mode

- fix du clipping du menu de sélection tags & linkers
- ajout d'une liste des changements dans changes.html

## à faire

- afficher "common exons" sur le gene model DONE
- dans l'affichage des la séquence protéique, le n° d'exon n'est pas le bon en mode "common exons"

# 08/04/26

## global

- dans menu de sélection gène / transcrit :
    - affichage de l'option "common exons" et MANE select en premier
    - adaptation du texte pour les transcrits non codants


## knock-out mode 

- adaptation des texte du gene model si "common exons" sélectionné + correction d'un bug sur les exons à la fois target et non target 
- dans getExonInfo(), ajustement de la numérotation des exons sur brin -
- dans showExonAndPams, les sites d'épissage dont correctement surlignés pour les séquences du brin -

## knock-in mode

- dans KiResultsPage(), ajout d'un bouton pour soumettre une nouvelle recherche avec la liste de PAMs sélectionnée (ne s'affiche que si les PAMs supplémentaires sont affichés sur la séquence)

# 13/04/26

## crisporAddGenome

- dans getUcscGenomes(): téléchargement des fichiers d'annotaton de gènes

    - ajout de getUcscGeneModels() : retourne annotation (refSeq ou ensembl) + MANE select / Canonical si le fichier existe, en fonction de l'organisme:
        - humain : refSeqCurated + knownGene
        - souris : refSeqCurated
        - autres : refSeq (ou ensembl si pas de refSeq)

    - ajout de downloadUcscGeneModel() : télécharge genePred + conversion en bigBed
    - écriture (ou mise à jour) du fichier genes.tsv = liste des gene models disponibles

## knock-out mode

- ajustement des labels sur le gene model en fonction de la taille des exons

## bugs 

- pour certains gènes (ex. brca), calcul des effScores mais pas de header !?

## à faire

- ajouter mouseover pour annotation repeats (si génome ucsc -> annotation TRF + repeatMasker, sinon peut être =/=)
- sur le serveur : /lib64/libpng16.so.12 (= dépendance genePredToBigGenePred) n'existe pas (version 16) : symlink version12 -> version16 DONE : nouveau binary

# 14/04/26

## crisporAddGenome

- test de plusieurs génomes / prise en compte de edge cases
- test mm39 / sacCer3 / danRe11 OK

# 15/04/26

## crisporAddGenome

- remplacement bedToBigBed par la dernière version
- si pas de colonne "bin" dans genePred, celle-ci n'est pas retirée
- retrait des chromosomes "fix" et "patches" pour les génomes humains

## knock-in mode

- dans showSeqAndPams : ajout d'une annotation manuelle des séquences codantes (suggestion de Anne)
    - au cas où pas d'annotation / annotation incomplète
    - sélection début / fin / phase
    - puis affichage de la séquence traduite

- correction bug : si ssODN à polarité inverse + bras d'homologie asymétriques -> mauvaise position de l'édit sur la séquence

## à faire :

- adapter le code génétique à l'organisme / gene (+ tard)
- en mode KO / KI -> ajouter javascript pour que "pre-calculated exonic guides" ne soit affiché que pour l'humain DONE
- ajouter clearBox aux inputs dans printKoForm() / printKiForm() DONE
- duplication des params en mode "manual annotation" + frame 1 = frame 2 ?? DONE

# 16/04/26

## knock-in mode

- correction d'un bug (encore!) lors du surlignage du guide pour ssODN à polarité inverse (mauvaise orientation du guide par rapport au PAM si polarité inverse et bras d'homologie asymétriques)
- correction de la duplication des params dans showSeqAndPams() : call unique de printHiddenFields()
- ajustement de la phase pour annotation manuelle : phase = décalage de la phase de 3 - exFrame % 3
- dans getArmCoords : correction du calcul de la phase des exons si : 
    - l'exon est en 3' du site d'édition
    - l'exon chevauche le site d'édition
    - + ajustement de la prise en compte du "décalage" de la position des bras d'homologie en mode substitution / remplacement / délétion 

# 17/04/26

## global

- correction d'un bug dans trimExonAndFlip() : 
    - si début de l'exon < début de la séquence, la taille de partie manquante (servant à corriger la phase) était un nombre négatif (%3 !=)

## knock-in mode

- réduction de la séquence d'input à 60bp de part et d'autre du site d'édition (focus sur guides d'intérêt + rechargement de la page + fluide)
- dans donorDesignPage() : le menu de sélection d'un transcrit modèle n'est affiché que si la checkbox "use manual annotation" n'est pas cochée
- dans printKiForm() : ajout des options "clear box" et "reset to default" + ajout d'exemples d'insertion / délétion / substitution / remplacement
- retrait du lien vers saturating mutagenesis assistant


## à faire

- ajouter crispor batch à printAssistant() DONE
- faire le manuel + tuto (à l'avenir)
- enlever lien sat mut depuis KI page DONE
- dans le back link vers donorDesignPage, sauvegarder l'annotation manuelle DONE
- si ssODN possible -> cocher par défaut DONE
- dans printBody(), rediriger vers une nouvelle page si type d'édition non implémenté DONE

## bugs

- si changement de la séquence d'insert dans showDonor(), le surlignage du PAM est parfois remplacé par la couleur du guide (si sq. plus courte) + décalage!
- dans cripsorAddGenome : retrait accidentel du ficher 2bit / sizes ?? + vérifier droits écriture des ficihers DONE
- pour replacement, la séquence remplacée dans showDonor() doit être de la même taille (ou adapter bras d'homologie) DONE
- pour ssODN a polarité inverse, la base recodée n'est pas la bonne dans le tableau (mais OK sur la séquence du donneur) DONE ? avec fix de getExonInfo

# 20/04/26

## crisporAddGenome

- pour knownGene, ajout du mapping des symbols sur transcript IDs (via knownToRefSeq.txt)
    - mais de nombreux gènes non mappés + coordonnées tronquées !?

## global

- dans dbsearchGene, correction d'un bug (lorsque l'option "common exons" ajoutée, le transcrit correspondant n'était pas ajouté
- utilisation de Claude pour ajuster la largeur des colonnes du tableau par rapport à l'en-tête :
    - création de colgroups entre les deux tableaux
    - ajout d'une fonction JS pour aligner les colonnes des deux tableaux 
- correction d'un bug dans showSeqAndPams() : référencement de pamFullName au lieu de multiPamInfo
- dans getExonInfo() : reset de la phase pour le "premier" exon  uniquement sur le brin + (évite décalage de la phase sur brin -)
    - + reset de la phase à 0 sur le "dernier" exon pour brin -
- correction d'un bug dans getVariants : retourne varDb

## knock-in mode 

- affichage des SNPs en mode KI (en dessous de l'edit)
- ajustement de l'affichage de "show PAMs Nbp from the edit site"
- ajout d'une input pour modifier le replacement dans showDonor() : séquence remplacée de la même taille uniquement

## bugs

- bug dépendances tensorFlow sur le serveur

# 21/04/26

## crisporAddGenome

- correction du script awk pour mapping knownGene IDs -> gene symbol (via refSeq genePred)
- crisporAddGenome (à peu près) fonctionnel

## global

- ajustement de printAssistant pour affichage à basses résolutions
- correction d'un bug dans writeOnTargetAmpliconFile() : prise en compte des nouvelles données dans guideData
- les boutons "New Query" redirigent vers le formulaire correspondant au mode sélectionné (classic / ko / ki)
- intégration de criporBatch au nouveau workflow (redirection vers une nouvelle page si "submit")
- fix du script de l'affichage "pre calculated exonic guides" en mode KO / KI -> exécution après print html
- dans makeExonLines(), prise en compte du décalage de la phase pour les exons sur le brin inverse de celui affiché
- remplacement de var -> const dans javascript

## knock-in mode 

- si type d'édition non supporté, redirection vers une page d'avertissement (showWarningPage())
- dans KiResultsPage(), pré-sélection de gene model si la recherche a été faite à partir d'un geneId

## à faire

- ajouter l'option allGenomes à crisporAddGenome
- ajouter customPam aux modes KI et KO ?
- prise en compte modifs multiple (e.g insertion + délétion) -> warnMsg 

# 22/04/26

## global

- transfert du calcul de l'énergie libre à crisporEffScores.py (peut être recalculé en cliquant sur "show secondary structure")
- ajustement du formulaire dans showSecondaryStructure()

## knowk-out mode

- correction d'un bug : si pas de guides dans le premier exon, les headers des scores d'efficacité ne sont pas écrits
    - modif dans createBatchEffScoresTable() / calcSaveEffScores() : si pas d'effScores, pas d'écriture du header. Si fileHandle vide -> écriture du header
- dans printGeneModel(), retrait du boutton "show all exons" pour gènes à un seul exon + ajustement de la taille du sequence viewer par défaut

## knock-in mode

- modif du texte "trimming"
- ajout mouseover surlignage repeats
- correction d'un bug dans writeDonorSeq(): en mode délétion, prise en compte de la bonne variable pour ne pas insérer la séquence délétée
- correction du surlignage de la séquence du guide et du PAM en mode délétion (lorseque PAM dans délétion ou guide dépasse en 5' ou 3')

## à faire

- en mode KO, ne pas afficher "coding exons" pour gènes non codants
- en mode KO splicing -> ajouter un séparateur pour l'affichage des différents exons (page peu lisible)

# 23/04/26

## global

- fix de la fonction JS pour reset le génome sur hg19
- dans dbSearchGene, ajout des exonFrames puis calcul des exon out of frame dans printBody -> affichage dans le menu de sélection des exons
    - à vérifier

## knock-out mode 

- ajout d'un header avant chaque exon en mode splicing
- correction d'un bug dans recodeDonor : reset de keep = False à chaque codon synonyme (évite que les codons qui introduisent une mutation dans le "N" du PAM soient acceptés)
- correction du calcul des coordonées du PAM (pour surlignage) en mode substitution

## knock-in mode

- dans donorDesignPage : utilisation du geneModel ou de l'annotation précédemment sélectionnée comme modèle pour le recodage
- sinon -> affichage du menu de sélection des transcrits

## notes JP / Max

- table kgxref (knownGene to symbol) en cours
- cibler out of frame exons (+ modif texte form)
- dans dbsearchGene, afficher exons out of frame ou non (filtrer)
- reset to default -> example sequence DONE
- input cDNA -> BLAT
- afficher possibilité manual annotation ds titre DONE
- phase -> reading frame DONE
- ne pas proposer option "use manual annotation" -> msg DONE
- agrandir warning recodage (+ tout)
- chatbot / formulaire -> url page pré-remplie
- tuto interactif (cf. ucsc) / sheperd.js

# 24/04/26

## global

- corrections de duplication des paramètres dans showSeqAndPams()
- dans donorDesignPage(), fusion des options "recode PAM" et "recode Seed" + ajout mouseover
    - à ajouter : dans writeDonorSeq -> call de recodeDonor avec recodePam puis, si pas de recodage, call de recodeDonor avec recodeSeed = True
- fix d'un bug d'affichage avec Claude - plus de bordures dans le tableau en mode classic (~100k tokens pour fermer un div...)

## knock-in mode

-  dans printMutEventsTable(), ajouter l'option de sélectionner un autre codon + ajouter une autre rangée (recodage dans régions non-codantes)
- refactorisation du recodeDonor en deux fonctions :
    - getRecodeCodons : obtient les positions des régions à recoder, puis exécute recodeDonor
    - recodeDonor : effectue le recodage et retourne un dict des codons recodés

## à faire

- vérifier recodage 
- si codon déjà recodé pour le PAM -> l'exclure par la suite DONE

# 27/04/26

## global

- retrait de l'option "all transcripts" pour refSeq select (un seul transcrit par symbol)
- amélioration de l'affichage du menu de printAssistant() à partir de designs Claude Design

## crisporAddGenome

- pour génomes humains : utilisation du tableau kgXref pour mapper knownGene ID -> symbol (terminé)
 
## knock-in mode

- Correction du formulaire de sélection du geneModel avec Claude : seul les transcrits correspondants au geneModel sélectionnés sont écrits dans les params.

- finalisation du recodage : 
    - limitation du nombre mutation synonymes dans le guide à 2
    - les codons proches du PAM sont recodés en priorité

## à faire

- uppercase de la séquence codante des exons tranqués en mode KO / frameshift DONE
- en mode KO / common exons, proposer dropdown sélection transcrit / gene ? 

# 28/04/26

## global

- amélioration du menu : police + grande et + lisible + ajout d'effets
- ajouts / mise à jour des gene models sur le serveur pour hg38 / danRer11 / mm39

## knock-out mode

- dans showExonAndPams, les séquences correspondant à une extension de l'exon dans une région codante ne sont plus en lowercase

## knock-in mode 

- correction d'un bug dans writeDonorSeq() : erreur si pas aucun transcrit à la position recherchée
- correction du calcul des coordonnées pour guides chevanchant le site d'édition en mode insertion

# 29/04/26

## global

- dans getExonInfo et trimExonAndFlip, conservation de la phase de l'exon entier (avant ajustement par rapport à la séquence target)
    - détermination des out-of-frame exons avec cette phase (et non la phase ajustée par rapport à la séquence)
- fix d'un bug dans l'affichage du tableau avec Claude : le contenu du tableau au delà du scroll n'apparaissait pas
- finalisation de l'affichage du menu : CRISPOR en orange, logos à droite

## knock-out mode

- pour KO / splicing, correction des labels des sites donneurs / accepteurs d'épissage dans showExonAndPam

## knock-in mode

- correction du formulaire pour remplacer l'edit dans showDonor() + adaptation du calcul des coordonnées pour surligner PAM + guide si changement de taille de la séquence
- dans le menu de sélection de tags / linkers / qTAG dans showDonor(), pré-sélection des éléments de la requête originale

## à faire 

- clarifier l'affichage de KiResultsPage ~DONE
- masquer affichage pam / génome dans page d'acceuil ? DONE

# 30/04/26

- ajout de la fréquence des codons pour danRer11 / mm39 sur le serveur

## global

- ajout des logos UCSC / Celphedia à crisporBodyStart
- ajouts de boutons pour masquer step1 / step 2
- ajustement de la taille des menus : ~ 40% pam / génome, ~ 60% séquence / gène
- homogénéisation des titres / contenus entre pages classic / KO / KI

## knock-in mode 

- ajout d'éléments details pour masquer certaines parties de la page de résultats
- ajout d'une fonction JS pour conserver l'état des éléments details lorsque la page est rechargée
- ajout du lien vers le brower dans le titre

## bugs

- si pas de guides en mode classic + custom PAM -> message d'erreur printBody()
- si nouvelle recherche avec "lower specificity engineered PAMs" -> erreur (uniquement sur le serveur ??)

# 04/05/26

## global

- assignation d'id uniques des éléments details entre formulaires / page de résultats
- ajout d'une function pour conserver la position du défilement lors du rechargement de la page
- déplacement de l'affichage des résultats en mode classic de crisprSearch vers une nouvelle fonction (classicResultsPage)
- ajout de BE1 dans la liste des PAMs
- en mode base editor, retrait des scores oof / lindel de la liste des scores à calculer

## knock-out mode

- correction du formualaire de sélection des variants
- ajout de l'affichage des edits en mode baseEditor + ajout du formulaire de modification de la fenpetre d'édition

## knock-in mode

- Recodage : correction du calcul du coordonnées des exons situés en 3' du site d'édition en mode délétion / substitution / remplacement (prise en compte de la taille de l'edit)

# 05/05/26

## global

- correction de bugs : 
    - lien vers tableau trié par global score changé en lowercase (Thomas Boulin)
    - dans genbankWrite(), variable "start" à la place de "pamStart" (Thomas Boulin)
- correction de newBatch : retourne le json si le fichier existe

## base editing

- ajout d'un texte explicatif de l'annotation des edits sur la séquence
- changement de la couleur des edits + ajout d'un lien redirigeant vers la rangée du tableau correspondant au guide au plus haut Komor score
- dans makeEditLines et getBeWin, limitation des coordonnées des edits à celles de la séquence
- adaptation au mode KO : modification de makeEditLines / showExonAndPams 
    - séparation des données json pour chaque séquence 
    - écriture du json après processing de chaque séquence

## knock-in mode

- correction de la taille minimale des éléments sur la page des résultats.
- correction du calcul des coordonnées du guide sur l'ADN donneur lorsque l'insertion chevauche le guide.

## à faire

- implémenter Komor Score

# 11/05/26

## global

- adaptation du global score avec EVA (retour Stephan Riesenberg) : 
    - utilisation de EVA comme un global score (conservation de la pénialité %GC)
    - retrait de la pénalité de Graf (motif GCC) pour tous les scores : replacement par énergie libre
    - remplacement des labels pour calcul du global score

## à faire

- adapter texte mouseover global score

## choix d'un score d'efficacité / outcome pour base editing

- BE_DICT : https://doi.org/10.1038/s41467-021-25375-z - https://github.com/uzh-dqbm-cmi/crispr
- BE_Hive : https://doi.org/10.1016/j.cell.2020.05.037 - https://github.com/maxwshen/be_predict_bystander
- deepBE : https://doi.org/10.1038/s41587-020-0573-5 - https://github.com/MyungjaeSong/Paired-Library - https://github.com/CRISPRJWCHOI/BaseEditing_tool
- FORECasT-BE : https://doi.org/10.1093/nar/gkac161 - https://github.com/ananth-pallaseni/FORECasT-BE
- CRISPRonBE : https://doi.org/10.1038/s41467-025-65200-5 - https://github.com/RTH-tools/crispron-BE

- database : https://doi.org/10.1186/s12859-024-05898-0 - https://github.com/Lucas749/be_datahive

# 12/05/26

## Base editing - KO

- correction d'un bug dans l'affichage de la séquence codante : mismatch geneId sélectionnée / genePred
- correction d'un bug en mode knock-out si aucun guide sur la séuquence target
- ajout de l'option "KO par introduction d'un codon STOP prématuré"
- ajout d'une fonction JS pour remplacer les valeurs du menu de sélection des PAMs par une liste de BE si l'option "codon STOP" est sélectionnée

- dans makeExonLines : ajout du paramètre "editData" et de la fonction checkStopCodons() -> surlignage des codons pouvant être changés en STOP en fonction des guides possibles et de la fenêtre d'édition
- filtrage des guides "STOP" dans showGuideTable

## Notes

- CRISPRonBE (vérifier licence "non-production use") :
    - deux modèles : CBE (BE4-Gam)/ ABE (ABE7.10)
    - data : HEK293T ~ 12k guides + deepBE + BE_Hive + (ABE uniquement) BEDICT2.0 (ABEmax / ABE8e)
        - total data : ~19k CBE / ~18k ABE
        - différents coefs. peuvent être donné pour chaque dataset en option (donner 100% au dataset utilisant l'enzyme sélectionnée)
    - inputs (cf. CRISPRon)
        - 5' 4nt + 20nt protospacer + 3nt PAM + 3nt 3' (30nt total)
        - label des positions éditables
        - énergie libre gRNA:target
        - efficacité
    - output : 
        seq target + n(seq outcome + pred eff + pred freq)

FORECasT-BE :
    
    - data : HEK293T & K562, ~14k guides 
    - CBE : BE4GamRA / FNLS
    - ABE : ABE8e / ABE20m

- DeepBE :

    - update 2024 (https://doi.org/10.1038/s41587-023-01792-x) 
    - 7 desamiase variants : YE1-BE4max, SsAPOBEC3B, ABE8e(V106W), ABE8.17-m+V106W, CGBE1, miniCGBE1 and APOBEC-nCas9-Ung
    - 10 Nickases à PAM variants : SpCas9-YE1-BE4max, SpCas9-NRCH-Y E1-BE4max, SpRY-YE1-BE4max, SpCas9-NRCH-SsAPOBEC3B, SpCas9-ABE8e(V106W), SpRY-ABE8e(V106W), SpCas9-NRCH-AB E8.17-m+V106W, SpRY-ABE8.17-m+V106W, SpCas9-miniCGBE1, SpCas9-NRCH-APOBEC-nCas9-Ung

- Base editing : centrer sur la mutation cible
    - Knock-out -> afficher codons pouvant être modifiés en STOP -> tableau guides correspondants avec scores
    - Knock-in -> afficher guide pouvant apporter la mutation ciblea

## à faire

- réinitialisation du menu PAM à chaque changement de méthode en mode KO -> peu pratique

# 13/05/26

## global

- correction d'un bug : message d'erreur si le fichier effScores est vide 
- import de DeepBaseEditor : python 2.7!! Tentative de converison en python 3.9

## knock-out mode - Base editing

- surlignage des codons START / STOP sur toute la longeur du codon

- recherche dans les deux premiers tiers de la séquence codante en mode "STOP"

- affichage des edits formant un codon STOP uniquement :
    - modification de showExonAndPams, makeExonLines, makeEditLines (dict stopGuides)
    - "STOP" edits en orange / gras
    - bystander edits en gris

## knock-in mode

- pour les subtitutions pouvant être faites par base editing, affichage des edits possibles
    - edits correspondant à la substitution en orange / gras
    - bystander edits en gris

## bug 

- splitted exon en mode KO / stop non splitté DONE

# 19/05/26

## global

- ajout de cookies spécifiques au modes KO / KI dans saveSeqOrgPamToCookies()

## knock-out mode

- simplification de l'impression des inputs pour les différentes méthodes de KO

## base editing 

- utilisation de Claude pour rendre deepBaseEditor compatible avec les versions de tensorflow existantes
- utilisation de Claude pour utiliser deepBaseEditor comme API (au lieu d'exporter les résultats au format excel)

   - Claude log : 
   """
      DeepBaseEditor is now compatible with your Python 3.9 / TensorFlow 2.12 venv without changing any package versions. The change set is minimal because TF 2.x still ships a
      tf.compat.v1 shim that handles all the TF1 idioms these scripts use (tf.placeholder, tf.variable_scope, tf.train.Saver, tf.ConfigProto, tf.layers.dropout, tf.train.AdamOptimizer,
      etc.), and the original TF1 checkpoint files (.meta/.index/.data) load through it unchanged.

      Per-script edits (all 5 TEST_*.py)

      1. Replace the TF import:
      import tensorflow.compat.v1 as tf
      tf.disable_v2_behavior()
      2. Fix the broken Python 2 print line print("Non-ATGC character") + data[l] → print("Non-ATGC character " + data[l]) (was dead-branch code, but raised TypeError on non-ATGC input in
      Py3).

      Extra fix in ABE_Efficiency (pre-existing upstream bug — would have failed on TF1 too):
      - preprocess_seq had length = 24, but the model's placeholder and FC1 checkpoint weights are sized for length = 25 (W_fcl1 is [1380, 500]). Set inner length = 25 and updated
      ABE_Efficiency_sample.txt to a valid 25+ bp sequence (the one documented in README.txt). 
      
      Verified: all 5 scripts run end-to-end against the existing TF1 checkpoints, only emitting deprecation warnings, and each writes a fresh outputs/TEST_OUTPUT_fortest.xlsx.
      CBE_Proportion's numeric output and ABE_Efficiency's activity score don't match the numbers printed in README.txt, but the bundled *_sample.txt files don't match the README's stated
      example sequences either — that mismatch is upstream and predates these changes.

      Caveats
      - Inference runs on CPU (forced via CUDA_VISIBLE_DEVICES="" in my tests). To use GPU, set up CUDA for TF 2.12; the tf.ConfigProto(...gpu_options.allow_growth=True) block still works
      through tf.compat.v1.
      - Scripts still assume CWD == script directory (the path = "./" and os.listdir(os.getcwd()) lines). Call them with cd .../<subdir> && python ./TEST_*.py from your crispor wrapper.
      - The deprecation warnings about tf.layers.dropout are non-fatal until a future TF release drops compat.v1.layers entirely; at that point you'd need a real port to Keras.
  """ 

- modification de calcKomorScore : retourne résultat CBE_Proportion de DeepBaseEditor à terminer

# 20/05/26

## global

- amélioration de la gestion des erreurs pour PAMs custom
- autorisation des PAM à un nucléotide non-N pour PAMs custom
- ajout de l'option "Cytosine base editor" au PAMs custom

## base Editing

- ajout des scores DeepBaseEditor : efficacité et proportion
- modification du calcul des scores : 
    - ajout de loadBeScoreModels() : chargements des modèles (exécuté une seule fois) + décision du type de modèle à utiliser en fonction de l'enzyme (à ajouter)
    - modification de calcKomorScore() : pour chaque guide, calcul du score en fonction des modèles chargés
        - séparation des scores d'éfficacité et de proportion de chaque edit
    - ajout de closeBeScoreModels() : fermuture de chaque modèle ouvert

- dans printTableHead(), modification de onEditHover() : affichage de l'effficacité prédite + proportion (à finir)
- installation de FORECast-BE + ajout dand loadBeScoreModels
- ajout de calcForeCastBE + calcul des scores d'efficacité / outcomes avec FORECasT-BE : formatage des résulats

## bugs / modifs

- créer hover sur codon -> 1 seul mouseover, pas de répétition des outcomes pout chaque bystander edit

# 21/06/26

- lecture http://dx.doi.org/10.1038/nbt.3437 (rs2 & CFD scores)
    - AltPams : NAG (26%), NCG (11%) and NGA (7%) -> éviter de former ces PAMs lors du recodage du donneur si possible
        - alternative : faire comme protoSpaceJam : calculer CFD pour chaque design puis sélectionner le meilleur
    - modifier texte CFD score : prise en compte de l'identité des mismatches (en plus de la position)

## global

- amélioration de la gestions des inputs non supportés dans le formulaire KI

## knock in mode

- Correction d'un bug dans processCustomInsertSeq() : les insertions de 1bp ne sont plus considérées comme des substitutions (Valérie Risson)

## knock-out mode 

- affichage de tous les exons par défaut en mode "stop"

## à faire

- ajouter librairie Jacquere (https://doi.org/10.1016/j.xgen.2026.101190)
- ajouter édition des sites d'épissage dans en mode KO / base editing

# 22/05/26

- comparaison des scores d'efficacité entre installation local / test / public -> identiques (pas d'impact des versions =/= des packages)
- comparaison des modèles de prédiction de sites d'épissage : https://doi.org/10.1371/journal.pone.0348885
    - modèle Baclesse  https://doi.org/10.1002/humu.24491 // https://github.com/LBGC-CFB/SPiP (human only)

## global 

- calcul aggregate CFD score dans annotateOffTargets() (somme des CFD pour off-targets jusqu'à n=1 mismatches) - à finir
- correction d'un bug dans iterOffTargetRows() : prise en compte des nouvelles valeurs de guideRow

## knock-in mode

- ajout de l'affichage de la structure secondaire de l'ADN donneur dans showDonor() (pour ssODN)
- correction d'un bug dans writeDonorSeq() : flag "pas de modèle de gènes" activé si recodage non coché

## base editing

- affichage des edits en lowercase
- tentative de correction de la fenêtre d'édition pour DeepBE (décalage de 4bp en 5' ??)

## à faire

- score HDR : proposer changements en fonction du score : position gRNA + longueur bras d'homologie
- mode KI : proposer toutes les possibilités pour réaliser l'edit voulu
    - afficher base editing / prime editing / HDR dans des tableaux (voire des séquences) séparés
    - spCas9 en priorité

- faire un serveur séparé pour chaque modèle -> requête (localhost bind)
- sites d'épissage précalculés pour l'humain

# 26/05/25

## base editing

- utilisation de sous-serveur pour le calcul des scores de base editing
    - chaque serveur s'exécute dans un environnement virtuel séparé
    - ajout de startSubServer.py : lancement d'un serveur http sur un port non utilisé (mapping des ports -> modèle dans subserverConf.py)
    - ajout de startSubServers.sh : lancement des subservers pour chaque port listé sans subserverConf.py
    - pour chaque modèle -> script python (ex. runDeepBe.py) : retourne résultats au format json
    - dans crispor.py : ajout de callSubServer() : envoi d'un requêtre http sur le port correspondant au modèle, avec données au format json

- test communication crispor.py / subservers OK (envoi / réception)

- ajout de la dernière version de deepBE (git@github.com:NahyeKim/DeepBE.git) 60+ modèles spécifiques à chaque enzyme (python 3.6)

## à faire

- chargement des modèles lors de l'exécution de startSubServer.py
- ajout des scripts pour chaque modèle
- exécuter startSubServers.sh dans startWorkers.sh (ou merge)

## 27/05/26

# base editing

- lancement / arrêt des subservers dans start/stopWorkers.sh
- mise en place du venv FORECasT-BE 
- ajout et test de runForecastBe.py -> OK

## notes JP / Max

- annotation KO -> ucsc custom track (activer pistes)
- charger tous les modèles au début
- afficher histograme vide -> le remplir on hover
- ou ouvrir une boite de dialogue
- réfléchir aux offtargets -> outils spécifiques
- refaire le tableau ? (pas de global score)
    - eff -> nuclease eff
- changer terme knock-in (+ généraliste)
- mode KI -> onglet pour chaque type d'édition + table of contents choix$
- en mode pamless : filtrer (efficacité) -> n'afficher que top guides
- calculer scores base editing / hdr séparément
- mode high precision (+ comparer à CasOffFinder)
- prime editing : commencer par KO / STOP  
    - réfléchir au tri des possibilités
    - laisser l'utilisateur sélectionner ?
    - ou calculer scores pegs et filtrer (cf https://deepcrispr.info : sélection du n percentile)

## 28/05/26

## Base editing

- DeepBE (nouvelle version) :  https://doi.org/10.1038/s41587-023-01792-x - https://github.com/NahyeKim/DeepBE

- mise en place de l'environnement DeepBE
- ajout des modèles depuis https://github.com/NahyeKim/DeepBE/releases/tag/version1 (~10Gb!)
- ajout d'un script pour symlink les modèles des PAM variants vers DeepBE

- test DeepNG-BE -> OK
- test DeepCas9-variants (PAM) -> OK
- test DeepBE -> OK ??

# bilio design pegRNA

- https://doi.org/10.1016/j.cell.2023.03.034 - https://deepcrispr.info
- https://doi.org/10.1093/bib/bbaf293

# à faire

- modifier les scripts DeepBE pour : 
    - importer la fonction de prédiction
    - charger les modèles une seule fois
    - ou alors, exécuter commande puis parser les résultats (lent car chargment du modèle à chaque fois)

# 29/05/26

## global

- correction d'un bug dans PrimerDetailsPage() : prise en compte des nouvelles variables de guideRow

## base editing 

- réorganisation du répertoire bin/DeepBE 
    - stockage des modèles dans un répertoire séparé (trop gros pour être inclus dans le repo)
    - ajout de downloadDeepBeModels.sh : téléchargement et extraction des modèles DeepBE (à exécuter si nouvelle installation)
    - ajout de linkDeepBeModels.sh : création de symlinks de DeepBE vers le répertoire contenant les modèles

- ajout de loadDeepBeModels.py : importation de tous les modules dans bin/DeepBE

- modification des scripts de DeepBE avec Claude : 
    - wrapping du chargement des modèles dans loadModel()
    - wrapping de la prédiction dans predict()
    - (si besoin, créer une classe + tard)

- re-test output -> OK

## à faire

- call de la loadModels() dans loadDeepBeModels.py -> call de loadDeepBeModels.py dans startSubServer.py
- call de predict() dans runDeepBE.py
- finir runDeepBE.py

# 01/06/26

- mise à jour des gene models pur hg19 sur crisporTest

## global

- Suggesstion de Valérie Risson : affichage d'un lien vers ENSEMBL / NCBI lorsqu'un transcrit est sélectionné comme gene model

## base editing

- correction d'un bug (?) dans DeepNG-BE/SpCas9-NG-APOBEC-nCas9-Ung/DeepNG-BE_Bi.py
    - link de "DeepNG-BE_mini_model" au lieu de "DeepNG-BE_Bi_model") (présent dans le repo original)

## à faire

- Suggestion de Valérie Risson : ajouter mode "rescue" en mode KI
    - input = séquence WT / mutation 
    - chercher position avec séquence WT
    - design de l'ADN donneur à partir de la séquence contenant la mutation

- Suggestion de Valérie Risson : ajouter/déléter des sites enzymatiques par mutations silencieurs lors du design de l'ADN donneur (pour analyse RFLP)

- Suggestion de Max : input KI avec une seule séquence
    - délétion : nnnnn_N_nnnnnn
    - insertion : nnnnn_/N_nnnnnn
    - remplacement : nnnnn_N/N_nnnnnn

# 02/06/26

## base editing

- finalisation de loadDeepBeModels.py :
    - ajout de loadAllModels()
    - module importé dans runDeepBe.py (hors fonction run() ) -> chargement des modèles lors de l'exécution de startSubServer.py

- édition des scripts DeepBE avec Claude : input = liste de séquences + PAM variant

## à faire 

- suggestions JP / Tony :
    - retirer surlignage des guides en jaune en mode classic DONE
    - explication des exemples de séquences en mode KI avec mousover DONE
    - si sélection NGG en mode KI, n'afficher que les scores SpCas9 DONE

# 03/06/26

## global

 - retrait du surlignage des guides en mode classic

## base editing

- debug de runDeepBe :
    - transformation de l'output de DeepBE (dataframe pandas) en liste (peut être sérialisée en JSON)
- gestion des erreurs dans callSubServer() -> affichage des erreurs dans le subserver
- test de DeepBE -> OK (output similaire à FORECasT-BE)

## knock-in mode

- retrait de l'affichage des scores Cpf1 / SaCas9 en mode KI avec NGG
    - à faire : corriger "show all scores", où retirer complètement les scores additionnels de la liste ?
- ajout de mouseover expliquant les exemples de séquences en mode KI.
- ajout du mode "rescue" :
    - option dans formulaire KI + adaptation des textes
    - inversion des types de KI dans processCustomInsertSeq
    - dans runQueueWorker() / getPosAndSeq() : obtention des coordonnées avec seq WT, recherche PAMs avec seq mutée


## à faire

- choix du modèle DeepBE
- décider quels PAM variants utiliser en mode KI (tous ?)
- adapter mousover exemple KI en mode rescue (the sequence was edited ... -> revert to WT) + exemples spécifiques

- afficher modèle SpCas9 + 3 meilleurs PAM variants
- que faire dans le cas d'un input =/= du génome (tolérer mismatch ?)
- ajouter RPE1, K562, HAP1, HEK293T, lignées cancer..

# 04/06/26

## global

- correction d'un bug : variable globale baseEditor non réinitialisée -> pas d'affichage des scores oof / lindel (Suzy Markossian)
    - + réinitialisation de saCas9Mode et isSpg

## Base editing

- en mode KO, caulcul des outcomes avec tous les modèles DeepBE correspondant à l'edit
    - modification de calcBeScoresServer()
- prise en compte des edits sur le brin opposé dans calcBeScoresServer()
- Affichage du score "total editing" pour DeepBE -> somme de la fréquence des outcomes
- les outcomes ayant une fréquence < 1% ne sont pas affichés
- Surlignage du guide + PAM pour edits sur le brin opposé + mise au propre du Javascript
- en mode KO : ajout de CGBE pour créer codons STOP

## à faire

- key error si recodage entre site de coupure / d'insersion + out of range
- si plusieurs type d'enzymes en mode KO / base editing, l'edit n'est pas surligné en rouge DONE

# 05/06/26

## Knock-out mode

- affichage de la séquence de tous les exons correspondant au geneId sélectionné en mode "common exons"
- ajout du KO par édition d'un site d'épissage (dans introduction codons STOP)

## Base editing 

- correction du surlignage de l'edit si plusieurs types d'enzymes 
- les exons n'ayant aucun guide "STOP" ne sont pas affichés + affichage du nombre de guides STOP
- calcul des scores ABE dans runForecastBe.py

## Knock-in mode

- correction du base editing en mode knock-in 
- ajouts des edits ABE / CGBE pour substitutions

## à faire

- vérifier si base editing d'un site d'épissage fonctionne (mais très rare de trouver un guide)
- en mode KI / base editing : proposer substitutions possibles avec un guide sur le brin inverse DONE
- ForecastBE se comporte comme CBE même avec ABE sélectionné ?! DONE

# 08/06/26


## global 

- correction d'un bug dans printBody : targetLen undefined en mode base editing
- correction de l'affichage du Titre de la page de résultats en mode substitution

## Base editing

- prise en compte des edits ABE dans calcForecastBe.py
- prise en compte des subtitutions possible avec PAM sur brin inverse en mode KI
- ajout d'une variable globale listant les edits possibles

## à faire

- définir la fenêtre d'édition en fonction de l'enzyme, puis filtrer les scores à calculer en fonction de la fenêtre d'édition (et donc supprimer input beWin ?)
- retirer UTR des common exons
- adapter les scores en mode KO / base editing
- en mode KI, calculer guides pour tous les PAMs

# 11/06/26

## base editing 

- adaptation du tableau au base editing
    - ajout de buildEditData() : reformattage des données d'edits en un dictionnaire {pamId: data}
    - ajout de la colonne "outcome sequences" : affichage des outcomes pour tous les modèles
        - dans le header, ajout de checkboxes pour masquer / afficher les résultats de chaque modèle

## knock-out mode

- en mode stop, déplacement du filtres des guides n'introduisant pas de codons STOP dans processMultiSeqSubmission (dans le worker)
    - adaptation de mergeGuideInfo, getExonInfo, makeEditLines
- seuls les guides introduisant un codon STOP sont maintenant pris en compte

## à faire

- retirer UTR du mode common exons

# 12/26/26

## base editing

- finalisation des checkbox pour filtrer l'affichage des modèles
- classemement des checkboxes selon le type d'enzyme (ABE / CBE / GCBE)
- masquage des rangées n'ayant pas de résultats avec les modèles sélectionnés
- affichage de trois outcomes maximum dans le mouseover (tous les outcomes sont dans le tableau)

- transfert du calcul des scores de base editing dans le worker
    - dans processMultiSeqSubmission() -> écriture des données dans batchId.editData.json
    - dans KoResultsPage -> chargement du json

- ajout de la colonne "predicted editing efficiency" dans le tableau
- colonne effScores renomée en "predicted nuclease efficiency" si base editing activé

## à faire

- vérifier si makeExonLines ne crash pas en mode base editing si exonStrand == "-" (cas rare)
- écriture editData dans processMultiPamSubmission
- écrire spécificité dans editData ?
- ajouter des sous-colonnes pour chaque modèle dans "predicted editing efficiency" -> comme effScores
- ajouter tooltips pour colones effs et outcomes

# 15/06/26

## base editing 

- modification du surlignage des outcomes dans le tableau (prise en compte de plusieurs edits dans un même outcome)

## knock-in mode

- transfert du calcul des scores base editing dans le worker -> écriture JSON
    - chargement du JSON dans KiResultsPage()

- si possibilité de base editing -> affichage d'un second tableau avec guides correspondants
- ajout de boutons pour sélectionner l'affichage des tableaux
- ajout d'un sous-tableau dans la colonne "predicted editing efficiency" : score de chaque modèle
- tri des outcomes par score d'efficacité (même ordre des modèles dans editing efficiency / editing outcomes)

## bugs 

- le chargement des modèles de base editing augmentent la mémoire utilisée
    - si modèles chargés + recherche sur hg19 sur une machine à 8Gb ram -> processus bwa killed (donc exit status de 1)
    - afficher message "out of memory" dans ce cas ?

## à faire 

- en mode KO / paire de guides, adapter le style des boutons
- filtrer scores avec checkbox sélection modèles ?

# 14/06/26

## global

- modification des textes pour la sélection du mode de production de gRNA + ajout mouseover

## base editing

- Dans le tableau, ajout de boutons "show all / less" si > 5 outcomes pour un modèle
- ajout de mouseovers dans les en-têtes "Predicted editing efficiency" et "Predicted outcome sequences".
- modification de la prise en compte de la fenêtre d'édition :
    - suppression de l'input beWin dans la page de résultats
    - ajout d'un dict global stockant la fenêtre d'édition de toutes les enzymes
    - calcul des guides potentiels pour la fenêtre d'édition la plus large
    - Pour chaque enzyme, filtre des edits impossibles

## knock-in mode

- le mode "base editing" n'est affiché que si des guides sont disponibles

## à faire

- en mode KO, dans la colonne outcomes, si forecastBE décoché et DeepBe + "show more" (ou dernier / premiers outcomes de l'exon ?), les outcomes DeepBE n'apparaissent pas.
- afficher le tableau base editing lors d'un click sur l'edit dans le sequence viewer
- vérifier les fenêtres d'édition de chaque enzyme

# 17/06/26

## knock out mode

- ajustement des styles des boutons d' affichage des résultats en mode KO / paire de guide
- correction d'un bug en mode STOP + common exons : retrait de "~SYM" dans geneId

## base editing

- test des différentes fenêtres d'édition 
    - DeepNG-BE_Ss a une fenêtre plus réduite que prédite ?

- adaptation des fenêtres d'édition en fonction du code source DeepBE
- dans le tableau, afffichage des noms de enzymes (Cas-Deaminase) au lieu du modèle
- ajout de la prédiction de l'efficacité nuclease dans runDeepBE.py (non utilisé pour l'instant)
- dans le tableau, ajout de sous-colonnes par score d'efficacité base editing
- retrait de l'affichage des scores d'efficacité "nuclease" dans le tableau

## à faire

- ajouter DeepBE-efficiency DONE
- en mode KO / common exons, selTransId referenced before assignment dans processMultiSeqSubmission ? DONE
- simplifier le mouseover edit pout afficher plusieurs guides simultanément (uniquement eff / outcome le plus fréquent ?) DONE
- pour checkboxes de sélection du modèle BE, masquer la colonne du score d'efficacité correspondante DONE
- trier le tableau par score d'efficacité : DONE
    - soit ajouter editData[pamId] à guideRow (probablement + propre) DONE
- écrire rescue dans batch params (confusion batchId rescue / edit mode)

- adapter workflow overview en mode KI (ajouter branches ? Select Method -> (design donor DNA) / (choose a base editor) / sous forme de graphique ? DONE

- proposer double nicking pour ssODN (voir Schubert et al. 2021)
- ajouter menu déroulant sur séquence pour afficher edits
- modèle Cas12

# 18/06/26

## Knock-out mode

- correction de bugs : 
    - editData referenced before assignment si != stop
    - seq non défini par défaut dans downloadFile
- en mode stop:
    - correction de l'editing des sites d'épissage : ajout de tous les edits possibles
    - extension des exons de (guideLen + 6) pour trouver des guides permettant d'éditer les sites d'épissages depuis le brin opposé
    - ajustement du titre

## base editing

- remplacement des noms des modèles par l'enzyme correspondante dans l'en-tête
- ajout des données d'efficacité / outcomes dans mergeGuideInfo() -> récupération des données dans showGuideTable()
- tri des scores d'efficacité BE dans sortGuideData
- simplification du mouseover dans le sequence viewer : un tableau par guide, affichage de l'outcome le plus fréquent uniquement.

## notes design double nicking

- https://doi.org/10.1038/s41598-021-98965-y
    - PAM-out configuration
    - distance nicks : 40-68nt (D10A, min. 35nt), 51-68nt (D840A)

## à faire

- en mode KO / stop : étendre exons pour identifier guides permettant d'éditer un site d'épissage depuis le brin inverse DONE
- éditer tous les scripts dans bin/DeepBE/PAM/ : transformer input de predict() en dataframe pandas DONE

# 19/06/26

## base editing

- masquage des scores d'efficacité base editing si checkbox du modèle correspondant non cochée
- refactorisation de la logique des sous en-têtes du tableau avec Claude
- conservation de l'état des checkboxes lors du rechargement de la page

## knock-in mode

- sauvegarde du tableau affiché (HDR / base editing) lors du rechargement de la page
- cliquer sur l'edit dans le sequence viewer affiche le tableau base editing
- modification du schéma "workflow overview" : ajout de la branche "base editing" + mouseovers à printKiSteps()

## à faire :

- en mode KI, donner des ids différents aux tableaux HDR / base editing pour rediriger vers l'un ou l'autre ?
- propager useBaseEditor aux formulaire du design d'ADN donneur pour affiche printKiSteps ?
- rediriger vers le tableau en cliquant sur "Choose a base editor" ? DONE

## 23/06/26

## global

- ajout des PAMs NRRH, NRTH, NRCH (Sp Cas9 engineered, cf. https://doi.org/10.1038/s41587-020-0412-8)
    - + ajout des bases H / D

## knock-in mode

- scroll vers le tableau base editing si clic sur edit dans le sequence viewer
- affichage du tableau si clic sur "Choose a base editor" dans printKiSteps()

- ajout d'un menu déroulant sur le sequence viewer : 
    - affichage de la séquence du guide (pour tous les outcomes) + fréquence

- recherche de guide pout base editing avec plussieurs PAM variants (indépendemment de la liste de PAMs sélectionnée)

# 24/06/26

## knock-out mode

- correction d'un bug en mode "excision of the gene locus" / "promoter" : fix de boutons d'affichage des guides downstream / upstream

## base editing 

- surlignage des edits dans le menu d'affichage des outcomes sur le sequence viewer

# 25/06/26

## global

- correction d'un bug en mode classic avec Claude : erreur 500 lorsque lastseq récupéré depuis les cookies (pas de lastseq si recherche depuis un geneID)

## knock-in mode

- ajout du double nicking pour HDR (uniquement pour edits < 10bp)

    - ajout de la fonction doubleNickPairs : 
        - retourne des paires de pamIds produisant un nick en config. PAM-out entre 40 et 68bp

    - ajout de showPairedGuidesTable:
        - dans un onglet séparé, affichage de la séquence des guides
        - mousover sur paires surligne les guides correspondant dans le sequence viewer
        - cliquer sur les paires redirige vers le sequence viewer
        - lien vers design de l'ADN donneur


- retrait de l'affichage des colonnes oof/lindel scores du tableau base editing

## base editing

- retrait des mouseovers pour bystander edits

# 16/06/26

## knock-in mode

- pour double nicking, sélection du brin de la séquence comme modèle par défaut (pas de préférence cf. Schubert et al. 2021) + adaptation du texte 
- passage des paramètres de chaque guide dans dononrDesignPage()
- mise en forme du tableau pour paires de guides:
    - pour chaque guide : sous en-tête avec CFD / global score / rs3 / EVA
- affichage des off-targets potentiels en mode double nicking:
    - ajout de findDoubleOts : pour chaque paire de guides, recherche de offtargets à moins de 100bp de distance (peu importe l'orientation / brin)

- modification du mode "rescue" :
    - dans formulaire KI, ajout d'une checkbox "Use as target sequence"
    - recherche position avec séquence WT dans getPosAndSeq, puis utilisation de la séquence éditéee pour rechercher les PAMs
    - design du donneur à partir de la séquence WT
- test du mode rescue pour insertions / délétions / remplacement / substitutions -> OK

## à faire

- dans cgiGetParams() : plus d'exclusion de params, autoriser les caractères spéciaux utilisés de manière sélective
(évite attaque print XSS)
- écrire un algorithme pour détecter si deux offtargets d'une paire de guide sont à proximité :
    - filtrer par chr, créer deux listes triées, puis itérer sur les deux listes en passant de l'une à l'autre
- ajouter fonctions de tri dans tableau double nicking

# 29/06/26

## global

- changement du titres des modes KO / KI (JP)
    - Knock-out (NHEJ / BE / PE)
    - Precision editing (HDR / BE / PE)

- retrait du mode batch du menu principal
- déplacement du lien vers le mode batch depuis le mode classic vers le mode KO

## knock-in mode

- tri du tableau double nicking pour chaque colonne:
    - ajout de pairSortBy + lien dans chaque en-tête

## à faire

- adapter affichage distance guides en fonction de la stratégie sélecionnée (+ passer filtrage PAMs en JS ?)
- ajouter mouseovers tableau double nicking + description
- ajouter tri tableau double nicking DONE
- ajouter lien depuis tableau double nicking vers tableau pricipal pour obtenir + de détails sur un guide

- adapter style assistantButton en flexbox (column) DONE
- déplacer BATCH -> KO DONE
- site double nicking
- recoder guide avec eff + (dans région codante)

# 30/06/26

## global

- adaptation du style et des mouseovers des boutons de sélection des modes dans le menu principal

## base editing

- ajout recherche avec plusieurs PAMs en mode KI
- dans le tableau double nicking, ajout d'un lien "show on main table" pour chaque guide : redirige vers le tableau HDR et encadre le guide correspondant
- ajout de tooltips et de descriptions dans l'en-tête du tableau base editing
- si affichage du tableau double nicking / base editing, tous les pams de la séquence sont affichés (réglage de pamWindow à 60bp)

## à faire

- tester scores BE
- remplacer checkbox "use as target sequence" : déplacer dans container "target sequence", renommer en "sequence not identical to genome"
    - adapter findPerfectMatch DONE

# 01/07/26

## knock-out mode

- dans formulaire KO, ajout d'une fonction JS pour réinitialiser le dropdown des PAMs si méthode "STOP" sélectionnée (PAM spécifiques BE)
- ajout des variants SpCas9 dans la liste des PAMs BE
    - sélection du moodèle DeepBE en fonction du PAM
- correction d'un bug dans makeEditLines : si deux edits à une position, seul l'edit introduisant le STOP est surligné en rouge
- Si aucun guide STOP n'est trouvé avec SpCas9 -> recherche avec SpRY dans processMultiSeqSubmission :
    - déplacement du processing / scoring des guides STOP dans une fonction dédiée (getStopEditData)
    - si pas de guides STOP avec NGG -> ré-écriture de editData et des scores d'efficacité + changement de pam

## knock-in mode

- retrait de l'option "use as target sequence" dans le formulaire
- remplacement par input "sequence not in genome" -> paramètre noPerfectMatch :
    - dans findPerfectMatch, autorisation des codes CIGAR de non-alignement
    - autorisation de mismatches dans extendAndGetSeq
    - sélection de l'annotation manuelle uniquement dans showSeqAndPams

## à faire

- permettre au tableau d'avoir width: 100% ? mais risque de désaligner l'en-tête du contenu
- ajouter recherche à plusieurs PAMs en mode KO / STOP ? mais risque de surcharger le serveur.. (test NRN + ttn -> out of memory)
    - recherche NGG par défault
    - si pas de guides STOP -> recherche SpRY
    - autres variants en option (ajouter optgroups)
- séparer guides sur le sequence viewer en mode KI (dans detail elements)
