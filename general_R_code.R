by_species <- iris %>% group_by(Species)

# One function
by_species %>% summarise_all(n_distinct)
by_species %>% summarise_all(mean)

# Use the _at and _if variants for conditional mapping.
by_species %>% summarise_if(is.numeric, mean)

# summarise_at() can use select() helpers with the vars() function:
by_species %>% summarise_at(vars(Petal.Width), mean)
by_species %>% summarise_at(vars(matches("Width")), mean)

# You can also specify columns with column names or column positions:
by_species %>% summarise_at(c("Sepal.Width", "Petal.Width"), mean)
by_species %>% summarise_at(c(1, 3), mean)

# You can provide additional arguments. Those are evaluated only once:
by_species %>% summarise_all(mean, trim = 1)
by_species %>% summarise_at(vars(Petal.Width), mean, trim = 1)

# You can provide an expression or multiple functions with the funs() helper.
by_species %>% mutate_all(funs(. * 0.4))
by_species %>% summarise_all(funs(min, max))
# Note that output variable name must now include function name, in order to
# keep things distinct.

# Function names will be included if .funs has names or whenever multiple
# functions are used.
by_species %>% mutate_all(funs("in" = . / 2.54))
by_species %>% mutate_all(funs(rg = diff(range(.))))
by_species %>% summarise_all(funs(med = median))
by_species %>% summarise_all(funs(Q3 = quantile), probs = 0.75)
by_species %>% summarise_all(c("min", "max"))

# Two functions, continued
by_species %>% summarise_at(vars(Petal.Width, Sepal.Width), funs(min, max))
by_species %>% summarise_at(vars(matches("Width")), funs(min, max))


txt <- c("The", "licenses", "for", "most", "software", "are",
				 "designed", "to", "take", "away", "your", "freedom",
				 "to", "share", "and", "change", "it.",
				 "", "By", "contrast,", "the", "GNU", "General", "Public", "License",
				 "is", "intended", "to", "guarantee", "your", "freedom", "to",
				 "share", "and", "change", "free", "software", "--",
				 "to", "make", "sure", "the", "software", "is",
				 "free", "for", "all", "its", "users")
( i <- grep("[gu]", txt) ) # indices
stopifnot( txt[i] == grep("[gu]", txt, value = TRUE) )
(ot <- sub("[b-e]",".", txt))
(ot <- gsub("[b-e]",".", txt))
txt[ot != gsub("[b-e]",".", txt)]#- gsub does "global" substitution

txt[gsub("g","#", txt) !=
			gsub("g","#", txt, ignore.case = TRUE)] # the "G" words

regexpr("en", txt)

gregexpr("e", txt)

make.cepnames(c("Aa maderoi", "Poa sp.", "Cladina rangiferina",
								"Cladonia cornuta", "Cladonia cornuta var. groenlandica",
								"Cladonia rangiformis", "Bryoerythrophyllum"))
data(BCI)
colnames(BCI) <- make.cepnames(colnames(BCI))

make.names(c("a and b", "a-and-b"), unique = TRUE)
# "a.and.b"  "a.and.b.1"
make.names(c("a and b", "a_and_b"), unique = TRUE)
# "a.and.b"  "a_and_b"
make.names(c("a and b", "a_and_b"), unique = TRUE, allow_ = FALSE)
# "a.and.b"  "a.and.b.1"
make.names(c("", "X"), unique = TRUE)
# "X.1" "X" currently; R up to 3.0.2 gave "X" "X.1"

data(varespec)
## Print only more common species 
freq <- apply(varespec > 0, 2, sum)
vegemite(varespec, scale="Hult", sp.ind = freq > 10)
## Order by correspondence analysis, use Hill scaling and layout:
dca <- decorana(varespec)
vegemite(varespec, dca, "Hill", zero="-")
## Show one class from cluster analysis, but retain the ordering above
clus <- hclust(vegdist(varespec))
cl <- cutree(clus, 3)
sel <- vegemite(varespec, use=dca, select = cl == 3, scale="Br")
## Re-create previous
vegemite(varespec, sp=sel$sp, site=sel$site, scale="Hult")
## Re-order clusters by ordination
clus <- as.dendrogram(clus)
clus <- reorder(clus, scores(dca, choices=1, display="sites"), agglo.FUN = mean)
vegemite(varespec, clus, scale = "Hult")

## Abundance values have such a wide range that they must be rescaled
## or all abundances will not be shown in tabasco
tabasco(decostand(varespec, "log"), dca)

## Classification trees for species
data(dune, dune.taxon)
taxontree <- hclust(taxa2dist(dune.taxon))
plotree <- hclust(vegdist(dune), "average")
## Automatic reordering of clusters
tabasco(dune, plotree, sp.ind = taxontree)
## No reordering of taxonomy
tabasco(dune, plotree, sp.ind = taxontree, Colv = FALSE)
## Species cluster: most dissimilarity indices do a bad job when
## comparing rare and common species, but Raup-Crick makes sense
sptree <- hclust(vegdist(t(dune), "raup"), "average")
tabasco(dune, plotree, sptree)
