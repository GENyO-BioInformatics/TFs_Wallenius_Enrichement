
- The random selection of collectri TFs lists cannot include TFs that they or their targets are not in the annotation databases. Does it makes sense to say select random lists of TF of size 3 for test when actually are not available for testing? No it makes no sense. It is two different questions:
	- Ours: *Does specific annotations appear commonly significant in different databases by analysing random TFs lists?*
	- The here mentioned problem, _Does Collectri has all, more or less TFs/targets than the annotation databases? And from which universe (set of TFs) should we extract the random lists?_
	- **We need to work only with the common set of TFs between Collectri and annotation databases, that at least have one target.**
	- Seleccionar TFs que no estan en las bases de datos de anotacion funcional, es lo mismo que sacar una estadistica de enriquecimiento nula, ya que no apuntan a ningun término es como contar en indio...  Hacer un enriquecimiento con elementos que no existen en la db funcional no tiene sentido y mete ruido a las simulaciones de montecarlo. Hay que filtrar para lo comun por pares de dbs, es decir GO_BP-Collectri, GO_BP-Collectri, GO_BP-Collectri, GO_BP-Collectri tanto a nivel de TFs como a nivel de Targets... Coger 1000 listas aleatorias para cada conjunto
	- **Todas las listas aleatorias de TFs tienen que apuntar cada TF, al menos, a 1 target gene que esté en la db funcional**
	- Ask JuanAn about the proper statistical nomenclature to write about this.
	- 
