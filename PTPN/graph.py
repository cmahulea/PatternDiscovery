# Graph ploting
# Le Moigne, 2023-11-01

"""
Plot a graph with graphviz. Adapted from my Saph-104 project
"""

import graphviz


# Définition des classes

class Graphe():
    def __init__(self, nom):
        self.nom = nom
        self.sommetsNommer = {}
        self.sommets = []
        self.arcs = []

    def __str__(self):
        return str(self.nom)

    def __repr__(self):
        return f"Graphe : {self}"

    def affichage(self):
        """Affichage d'un résumé textuel du graphe"""  # utilisée pour l'affichage console
        print("Contenu du graphe", self, ":")
        print(f" - {len(self.sommets)} sommets : ", end="")
        print(*self.sommets, sep=", ")
        print(f" - {len(self.arcs)} arcs:", end=" ")
        print(*self.arcs, sep=", ")

    def construireAvecDonnee(self, donneesSommets, donneesArcs):
        """construire le graphe à partir de données en listes/tuples"""
        for nomSommet in donneesSommets:
            self.ajouterSommet(nomSommet)
        for donneesArc in donneesArcs:
            nomDepart, nomArrive, poids = donneesArc
            self.ajouterArc(nomDepart, nomArrive, poids)

    def ajouterSommet(self, nom):
        """création d'un sommet dans ce graphe"""
        sommet = Sommet(nom)
        self.sommetsNommer[nom] = sommet
        self.sommets.append(sommet)

    def ajouterArc(self, nomDepart, nomArrive, poids):
        """création d'un arc dans ce graphe
        si un sommet n'est pas défini la commande est ignorée"""
        if nomDepart not in self.sommetsNommer:
            print("Vous essayez de créer un arc depuis un sommet inexistant !")
        elif nomArrive not in self.sommetsNommer:
            print("Vous essayez de créer un arc vers un sommet inexistant !")
        else:
            self.arcs.append(Arc(self.sommetsNommer[nomDepart], self.sommetsNommer[nomArrive], poids))

    def affichageConsole(self):
        """affiche dans la console un tableau avec les plus courts chemin"""

        # préambule : contenu du graphe
        self.affichage()
        # plus long des plus courts chemins
#        print(f"""
# Plus long des plus courts chemins :
#  Entre {self.departPlusLongPlusCourtsChemins} et {self.arrivePlusLongPlusCourtsChemins} : distance de {self.longueurPlusLongPlusCourtsChemins}
#    {self.chainePlusLongPlusCourtsChemins}""")

    def tracerGraphe(self):
        """trace le graphe au fomat PNG avec le module graphviz"""
        dot = graphviz.Digraph(comment=self.nom, format='png')
        dot.attr(rankdir='LR')  # sens horizontal
        dot.attr('edge', arrowhead='open', fontname="Arial")  # style de pointe et police des flèche
        dot.attr('node', fontname="Arial")  # police des noeuds
        for sommet in self.sommets:
            dot.node(sommet.nom)
        # dot.node(self.arrivePlusLongPlusCourtsChemins.nom, color = 'red', fontcolor = 'red') # coloration du noeud final
        for arc in self.arcs:
            # if arc in self.plusLongPlusCourtsChemins:
            #    dot.edge(arc.depart.nom, arc.arrive.nom, label=str(arc.poids), color = 'red', fontcolor = 'red')
            #    dot.node(arc.depart.nom, color = 'red', fontcolor = 'red') # coloration des noeuds de départ et de leurs arcs
            # else:
            dot.edge(arc.depart.nom, arc.arrive.nom, label=str(round(arc.poids*100, 2))+" %", penwidth=str(arc.poids*5))
        dot.render(filename=self.nom)   # création du fichier .PNG

    def affichageHTML(self):
        """affiche dans une page HTML un tableau avec les plus courts chemin"""
        self.tracerGraphe()       # génération de l'image avec graphviz
        # préambule : image du graphe
        codeHTML = f"""<!DOCTYPE html>
<!-- Page générée par un code Python de Le Moigne-->
<html>
    <head>
        <meta charset="utf-8" />
        <title>Continuous Pattern Discovery</title>
    </head>
    <body>
        <h1 style="text-decoration: underline;">Continuous Pattern Discover</h1>
        <h2>Graphe étudié : {self}</h2>
        <figure>
            <img src="{self}.png" alt="représentation du graphique" />
            <figcaption>Graph : {self} (weights are the percentage of total patients who follow the arc, not only from patients in the previous event)</figcaption>
        </figure>
    </body>
</html>"""

        html = open(f"{self}.html", "w", encoding="utf-8")  # efface le contenu de l'éventuelle page existante
        html.write(codeHTML)
        html.close()


class Sommet():
    def __init__(self, nom):
        self.nom = nom
        self.arcsSortant = []

    def __str__(self):
        return str(self.nom)

    def __repr__(self):
        return f"Sommet : {self}"

    def affichage(self):
        """Affichage du sommet"""  # pour le débuggage
        print(f"Sommet {self.nom} avec les arcs sortant {self.arcsSortant}")

    def ajouterArcSortant(self, arc):
        """ajoute un arc sortant de ce sommet"""
        self.arcsSortant.append(arc)


class Arc():
    def __init__(self, depart, arrive, poids):
        depart.ajouterArcSortant(self)
        self.depart = depart
        self.arrive = arrive
        self.poids = poids

    def __str__(self):
        return f"({self.depart},{self.arrive},{self.poids*100} %)"

    def __repr__(self):
        return f"Arc : {self}"

    def affichage(self):
        """Affichage de l'arc"""  # pour le débuggage
        print(f"Arc de {self.depart} vers {self.arrive} de poids {self.poids}")


# tests
# donnees=("Article's example 1", ['1', '2', '3', '4', '5', '6', '7', '8'],
#         [('1', '2', 1), ('2', '3', 3/4), ('7', '8', 3/4), ('3', '4', 1/2), ('4', '7', 1/2)])

# donnees=("Article's example 2", ['1', '2', '3', '4', '5', '6', '7', '8'],
#         [('1', '2', 1), ('4', '5', 1), ('2', '3', 4/5), ('6', '7', 1/2), ('7', '8', 1/2), ('3', '4', 2/5)])

# nomGraphe, donneesSommets, donneesArcs = donnees

# graphe = Graphe(nomGraphe)
# graphe.construireAvecDonnee(donneesSommets, donneesArcs)


# graphe.affichageConsole()
# graphe.tracerGraphe()
# graphe.affichageHTML()
