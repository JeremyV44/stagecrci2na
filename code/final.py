import os
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from sklearn.decomposition import PCA
from gprofiler import GProfiler

class DataPreprocessor:
    def __init__(self, mass_file_dir, rna_file_dir, design_file_name, result_dir):
        self.mass_file_dir = mass_file_dir
        self.rna_file_dir = rna_file_dir
        self.design_file_name = design_file_name
        self.result_dir = result_dir

    def load_data(self, file_dir):
        return pd.read_csv(file_dir)

    def preprocess_data(self, data):
        # Preprocess data (if needed)
        return data

    def save_result(self, data, file_name):
        # Save result to specified directory
        data.to_csv(os.path.join(self.result_dir, file_name), index=False)

class GeneAnalysis:
    def __init__(self, list_files, stepp, group, result_file, title, pvalue_cutoff=0.05, lf_cutoff=0.58):
        self.list_files = list_files
        self.stepp = stepp
        self.group = group
        self.result_file = result_file
        self.title = title
        self.pvalue_cutoff = pvalue_cutoff
        self.lf_cutoff = lf_cutoff

    def filter_data(self, pvalue_cutoff, lf_cutoff, filter_type):
        if filter_type == 'pvalue_and_lfc':
            return self.df[(self.df['padj'] < pvalue_cutoff) & (abs(self.df['log2FoldChange']) > lf_cutoff)]
        elif filter_type == 'upregulated':
            return self.df[(self.df['padj'] < pvalue_cutoff) & (self.df['log2FoldChange'] > lf_cutoff)]
        elif filter_type == 'downregulated':
            return self.df[(self.df['padj'] < pvalue_cutoff) & (self.df['log2FoldChange'] < -lf_cutoff)]

    def plot_histogram(self, title, drop_na=False):
        if drop_na:
            data = self.df.dropna(subset=['padj'])['pvalue']
        else:
            data = self.df['pvalue']
        plt.hist(data, color='gray')
        plt.title(title)
        plt.show()

    def get_venn_diagram(self):
        for group, files in self.group.items():
            lvenn = []
            for file in files:
                df = pd.read_csv(file)
                if self.stepp == "DE":
                    filtered_df = self.filter_data(self.pvalue_cutoff, self.lf_cutoff, 'pvalue_and_lfc')
                    lvenn.append(set(filtered_df.index))
                elif self.stepp == "enrichplot":
                    lvenn.append(set(df.index))

            names = [file.split('/')[-1] for file in files]
            venn = self._create_venn(lvenn, names)

            venn.plot_venn_diagram()
            plt.title(f"{self.title} {group}")
            plt.savefig(f"{self.result_file}{group}.png")

    def get_jacquard_index(self):
        jacquard_indices = {}
        for group, files in self.group.items():
            l_jcrd = []
            for file in files:
                df = pd.read_csv(file)
                if self.stepp == "DE":
                    filtered_df = self.filter_data(self.pvalue_cutoff, self.lf_cutoff, 'pvalue_and_lfc')
                    l_jcrd.append(set(filtered_df.index))
                elif self.stepp == "enrichplot":
                    l_jcrd.append(set(df.index))

            jacquard_index = self._calculate_jacquard_index(l_jcrd)
            jacquard_indices[group] = jacquard_index

        jacquard_df = pd.DataFrame.from_dict(jacquard_indices, orient='index', columns=['jacq_index'])
        jacquard_df.plot(kind='bar')
        plt.title(self.title)
        plt.xlabel('Comparison')
        plt.ylabel('Jacquard index value')
        plt.xticks(rotation=90)
        plt.savefig(self.result_file)

    def _create_venn(self, lvenn, names):
        # Créer le diagramme de Venn
        plt.figure()
        venn = venn2(subsets=lvenn, set_labels=names)
        return venn

    def _calculate_jacquard_index(self, sets):
        # Calculer l'indice de Jaccard
        intersection = set.intersection(*sets)
        union = set.union(*sets)
        jacquard_index = len(intersection) / len(union)
        return jacquard_index

        
class DownstreamAnalysis:
    def __init__(self, base_dir):
        self.base_dir = base_dir

    def plot_treeplot(self, enr, result_file):
        plt.figure(figsize=(12, 8))

        # Créer un graphe dirigé
        G = nx.DiGraph()

        # Ajouter les nœuds et les arêtes au graphe en utilisant les données de 'enr'
        for _, row in enr.iterrows():
            parent = row['parent']
            child = row['term']
            G.add_node(parent)
            G.add_node(child)
            G.add_edge(parent, child)

        # Positionner les nœuds du graphe
        pos = nx.spring_layout(G)

        # Dessiner les nœuds, les arêtes et les étiquettes
        nx.draw(G, pos, with_labels=True, node_size=1500, node_color="skyblue", font_size=10, font_weight='bold', arrows=False)

        # Titre et légende
        plt.title('Tree Plot')
        plt.legend()

        # Enregistrer l'image
        plt.savefig(f"{result_file}_treeplot.png")

        # Afficher le graphe
        plt.show()

    def perform_gprofiler_operation(self, res, result_file, title, regul, pvaluecutoff=0.05, lfcutoff=0.58, universe=None):
        if universe is None:
            universe = res.index.tolist()
        if not os.path.exists(result_file):
            os.makedirs(result_file)
        if regul == "both":
            res = self.filter_data(res, pvaluecutoff, lfcutoff, 'both')
        elif regul == "up":
            res = self.filter_data(res, pvaluecutoff, lfcutoff, 'upregulated')
        elif regul == "down":
            res = self.filter_data(res, pvaluecutoff, lfcutoff, 'downregulated')
        if isinstance(res, pd.DataFrame):
            res = res.index.tolist()
        gp = GProfiler(return_dataframe=True)
        gostres = gp.profile(organism='mmusculus', query=res, no_iea=True, sources=['GO:BP'], user_threshold=pvaluecutoff)
        gostres.to_csv(f"{result_file}.csv", index=False)
        gostres2 = gostres[gostres['source'] == "GO:BP"].copy()
        gostres2['GeneRatio'] = gostres2['intersection_size'].astype(str) + "/" + gostres2['query_size'].astype(str)
        gostres2['BgRatio'] = gostres2['term_size'].astype(str) + "/" + gostres2['effective_domain_size'].astype(str)
        gostres2 = gostres2.sort_values('p_value')
        gostres2.to_csv(f"{result_file}_simplify.csv", index=False)
        if len(gostres2) > 4:
            self.plot_treeplot(gostres2, result_file)            
            pass

    def enrich_with_GO(self, res, result_file, title, regul, pvaluecutoff=0.05, lfcutoff=0.58, universe=None):
        if universe is None:
            universe = res.index.tolist()
        if not os.path.exists(result_file):
            os.makedirs(result_file)
        if regul == "both":
            res = self.filter_data(res, pvaluecutoff, lfcutoff, 'both')
        elif regul == "up":
            res = self.filter_data(res, pvaluecutoff, lfcutoff, 'upregulated')
        elif regul == "down":
            res = self.filter_data(res, pvaluecutoff, lfcutoff, 'downregulated')
        if isinstance(res, pd.DataFrame):
            res = res.index.tolist()
        gp = GProfiler(return_dataframe=True)
        enr = gp.profile(organism='mmusculus', query=res, no_iea=True, sources=['GO:BP'], user_threshold=pvaluecutoff)
        enr.to_csv(f"{result_file}.csv", index=False)
        if len(enr) > 0:
            enr = enr[enr['source'] == "GO:BP"].copy()
            enr = enr.sort_values('p_value')
            enr.to_csv(f"{result_file}_simplify.csv", index=False)
            if len(enr) >= 10:
                self.plot_treeplot(enr, result_file)
                pass
                    
class Analysis:
    def __init__(self, preprocessed_mass_data, preprocessed_rna_data, design_data, result_dir):
        self.preprocessed_mass_data = preprocessed_mass_data
        self.preprocessed_rna_data = preprocessed_rna_data
        self.design_data = design_data
        self.result_dir = result_dir

    def run_analysis(self):
        # Perform PCA on preprocessed data
        pca = PCA(n_components=2)
        combined_data = pd.concat([self.preprocessed_mass_data, self.preprocessed_rna_data], axis=1)
        principal_components = pca.fit_transform(combined_data)

        # Plot PCA results
        self.plot_pca(principal_components, "pca_plot.png")

    def plot_pca(self, principal_components, file_name):
        plt.figure()
        plt.scatter(principal_components[:, 0], principal_components[:, 1])
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.title('PCA Plot')
        plt.savefig(os.path.join(self.result_dir, file_name))
        plt.close()

def main():
    # Répertoires et fichiers d'entrée
    mass_file_dir = "chemin/vers/les/données_massiques"
    rna_file_dir = "chemin/vers/les/données_rna"
    design_file_name = "chemin/vers/le/fichier_design"
    result_dir = "chemin/vers/le/répertoire_de_résultats"

    # Initialisation du préprocesseur de données
    preprocessor = DataPreprocessor(mass_file_dir, rna_file_dir, design_file_name, result_dir)

    # Chargement des données
    mass_data = preprocessor.load_data(os.path.join(mass_file_dir, "mass_data.csv"))
    rna_data = preprocessor.load_data(os.path.join(rna_file_dir, "rna_data.csv"))
    design_data = preprocessor.load_data(design_file_name)

    # Prétraitement des données
    preprocessed_mass_data = preprocessor.preprocess_data(mass_data)
    preprocessed_rna_data = preprocessor.preprocess_data(rna_data)

    # Initialisation de l'analyse des gènes
    gene_analysis = GeneAnalysis(list_files=[mass_data, rna_data], stepp="DE", group={"Group1": [mass_data], "Group2": [rna_data]},
                                 result_file=result_dir, title="Gene Analysis")

    # Exemple d'utilisation des méthodes de l'analyse des gènes
    gene_analysis.volcano_plot(pvalue_cutoff=0.05, lf_cutoff=0.58, title="Volcano Plot")
    gene_analysis.plot_pvalue_histogram(title="P-value Histogram")
    gene_analysis.plot_no_na_histogram(title="Histogram without NA values")

    # Exemple d'utilisation des méthodes pour la création de diagrammes de Venn et d'indices de Jaccard
    gene_analysis.get_venn_diagram()
    gene_analysis.get_jacquard_index()

    # Initialisation de l'analyse en aval
    downstream_analysis = DownstreamAnalysis(base_dir=result_dir)

    # Exemple d'utilisation des méthodes de l'analyse en aval
    # Assurez-vous d'avoir des données à utiliser pour ces méthodes
    # downstream_analysis.plot_treeplot(enr=dataframe_enrichment_results, result_file=result_dir)
    # downstream_analysis.g_profiler_result(res=your_data, result_file=result_dir, title="Your Title", regul="both")
    # downstream_analysis.enrich_with_GO(res=your_data, result_file=result_dir, title="Your Title", regul="both")

    # Initialisation de l'analyse PCA
    pca_analysis = Analysis(preprocessed_mass_data, preprocessed_rna_data, design_data, result_dir)

    # Exécution de l'analyse PCA
    pca_analysis.run_analysis()

if __name__ == "__main__":
    main()
