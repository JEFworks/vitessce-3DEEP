data <- read.table('data/20241209.HF.Coordinate.Celltype.Structure.transformed.txt', sep="\t", header = TRUE)
#data <- read.table('~/Downloads/20241201.HF.Coordinate.Celltype.Structure.txt', sep="\t", header = TRUE)

head(data)
dim(data)

hf <- table(data$HF) ## number of molecules per follicle 
hf <- sort(hf, decreasing=TRUE)
length(hf)
#hf.test <- 'HF245' ## focus on one follicle for now
set.seed(0)
#hf.test <- sample(names(hf), 10) ## works well
#hf.test <- sample(names(hf), 100) ## sample 100 hair follicles for now (gets a bit slow)
## sample at even intervals -> looks skewed, better do at log intervals
hist(hf)
n <- length(hf)
num_samples <- 50
log_indices <- round(10^(seq(0, log10(n), length.out = num_samples)))
log_indices <- unique(pmin(log_indices, n))
length(log_indices)
log_indices
sampled_elements <- hf[log_indices]
hf.test <- names(sampled_elements)

#vi <- data$HF == hf.test
vi <- data$HF %in% hf.test
table(vi)
data <- data[vi,]
head(data)

plot(data$transformedY, data$transformedZ, pch=".")
plot(data$x, data$y, pch=".")

########## need to order data in time else will overlap too much
## time actually doesn't seem as good of an ordering factor as size?
#time <- sapply(hf.test, function(x) {
#  mean(data$Time[data$HF == x])
#})
time <- table(data$HF)
hf.test <- names(sort(time, decreasing=TRUE))
hist(time[hf.test]) ## roughly uniform?

#data$HF <- factor(data$HF, levels=hf.test)

## need to rename HFs...
max_digits <- nchar(as.character(length(hf.test)))
str_nums <- sprintf(paste0("Hair-Follicle-%0", max_digits, "d"), 1:length(hf.test))
sorted_str_nums <- sort(str_nums)

data.reorder <- do.call(rbind, lapply(1:length(hf.test), function(i) {
  print(i)
  df <- data[data$HF == hf.test[i],]
  df$HF = sorted_str_nums[i]
  return(df)
}))
data <- data.reorder

########## fix missing structure annotation
table(data$Structure == "")
data$Structure[data$Structure == ""] <- "Unclear"

########## make cells.json
par(mfrow=c(1,1))
plot(data$transformedY, data$transformedZ, pch=".")
#plot(data$y, data$z, pch=".")
cells <- unique(data$cell)
cols <- sample(rainbow(length(cells)), length(cells)); names(cols) <- cells
final <- lapply(cells, function(cell.test) {
  
  data.cell <- data[data$cell == cell.test,]
  
  library(geometry)  # for convhulln
  #points <- as.matrix(data.cell[, c('x', 'y', 'z')])
  points <- as.matrix(data.cell[, c('transformedY', 'transformedZ')])
  if(nrow(points) > 3) { 
    hull <- convhulln(points, options = "FA")  # 'FA' ensures all attributes are returned
    hull.vertices <- points[hull$hull, , drop = FALSE]
    #if(nrow(hull.vertices) > 8) {
    #  set.seed(42)
    #  cluster <- kmeans(hull.vertices, centers = 5)
    #  hull.vertices <- cluster$centers
    #}
  } else {
    hull.vertices <- points
  }
  colnames(hull.vertices) <- c('x','y') ## for vitescce 
  cell.pos <- as.numeric(colMeans(points))
  
  angles <- atan2(hull.vertices[, 2] - cell.pos[2], hull.vertices[, 1] - cell.pos[1])
  ordered_indices <- order(angles)
  hull.vertices <- hull.vertices[ordered_indices, ]
  
  out <- list(
    genes = as.list(table(data.cell$Gene)),
    xy = cell.pos,
    factors = list(
      HF = as.character(data.cell$HF[1]),
      Time = as.character(data.cell$Time[1]), 
      CellType = as.character(data.cell$CellType[1]), 
      Structure = as.character(data.cell$Structure[1])
    ),
    poly = as.matrix(hull.vertices)
  )
  
  ## plot (double check)
  #points(points, pch=16, col=cols[cell.test])
  #points(out$poly, type="l")
  #points(t(out$xy), col='red')
  
  return(out)
})
names(final) <- cells
library(jsonlite)
exportJSON <- toJSON(final, auto_unbox = TRUE)
write(exportJSON, "data/soichiro_skin_cells.json")


########## make cell-sets.json
annot1 <- unique(data$CellType)
cellannotlist <- lapply(annot1, function(an) {
  set = data$cell[data$CellType == an]
  list("name" = an, "set" = set)
})
annot1.list <- list(
  "name" = "Cell Type",
  "children" = cellannotlist
)

annot2 <- unique(data$Structure)
strucannotlist <- lapply(annot2, function(an) {
  set = data$cell[data$Structure == an]
  list("name" = an, "set" = set)
})
annot2.list <- list(
  "name" = "Structure",
  "children" = strucannotlist
)
annot3 <- unique(data$HF)
hfannotlist <- lapply(annot3, function(an) {
  print(an)
  set = data$cell[data$HF == an]
  list("name" = an, "set" = set)
})
annot3.list <- list(
  "name" = "Hair Follicle",
  "children" = hfannotlist
)


cellset <- list("version" = "0.1.2",
                "datatype" = "cell",
                "tree" = list(
                  annot1.list,
                  annot2.list,
                  annot3.list)
)
exportJSON <- toJSON(cellset, auto_unbox = TRUE)
write(exportJSON, "data/soichiro_skin_cell-sets.json")


########## make clusters.json
data$Gene <- as.factor(data$Gene)
genes <- levels(data$Gene)
mat <- do.call(cbind, lapply(cells, function(cell.test) {
  data.cell <- data[data$cell == cell.test,]
  x <- as.numeric(table(data.cell$Gene))
  x/max(x) ## normalize for visualization only?
}))
rownames(mat) <- genes
colnames(mat) <- cells

info <- list(rows=genes, cols=cells, matrix=mat)
exportJSON <- toJSON(info)
write(exportJSON, "data/soichiro_skin_clusters.json")


########## make molecules.json
length(genes)
plot(data$transformedY, data$transformedZ, pch=".")
cols <- sample(rainbow(length(genes)), length(genes)); names(cols) <- genes
genes.pos <- lapply(genes, function(g) {
  #as.matrix(data[data$Gene == g, c('x', 'y', 'z')])
  #out <- as.matrix(data[data$Gene == g, c('x', 'y')])
  out <- as.matrix(data[data$Gene == g, c('transformedY', 'transformedZ')])
  colnames(out) <- c('x', 'y') ## for vitescce 
  
  # plot to check
  points(out, col=cols[g], pch=".")
  
  return(out)
})
names(genes.pos) <- genes

library(jsonlite)
exportJSON <- toJSON(genes.pos, auto_unbox = TRUE)
write(exportJSON, "data/soichiro_skin_molecules.json")


###################### make app
library(vitessceR)

## need to launch python server
## http-server ./ --cors -p 8000
#base_url <- "http://localhost:8000/"
base_url <- "https://raw.githubusercontent.com/JEFworks/vitessce-3DEEP/refs/heads/main/"

# Create Vitessce view config
vc <- VitessceConfig$new(schema_version = "1.0.16", name = "3DEEP of Skin v2")
dataset <- vc$add_dataset("3DEEP-Skin")$add_file(
  url = paste0(base_url, "data/soichiro_skin_cells.json"),
  file_type = FileType$CELLS_JSON
)$add_file(
  url = paste0(base_url, "data/soichiro_skin_molecules.json"), ## transcript positions
  file_type = FileType$MOLECULES_JSON
)$add_file(
  url = paste0(base_url, "data/soichiro_skin_clusters.json"), ## needed for gene names?
  file_type = FileType$CLUSTERS_JSON
)$add_file(
  url = paste0(base_url, "data/soichiro_skin_cell-sets.json"), ## cell-type annotations
  file_type = FileType$CELL_SETS_JSON
)

desc <- vc$add_view(dataset, Component$DESCRIPTION)
desc <- desc$set_props(description = "Soichiro et al, 3DEEP of Skin v3")
spatial <- vc$add_view(dataset, Component$SPATIAL)
spatial_layers <- vc$add_view(dataset, Component$LAYER_CONTROLLER)

gene_list <- vc$add_view(dataset, Component$FEATURE_LIST)
cell_sets <- vc$add_view(dataset, Component$OBS_SETS)
heatmap <- vc$add_view(dataset, Component$HEATMAP)$set_props(transpose = TRUE)

vc$layout(hconcat(
  vconcat(desc, spatial_layers),
  vconcat(spatial),
  vconcat(heatmap),
  vconcat(cell_sets, gene_list)
)) # need to change layout in exported contig
# {"dataset":"A"},"x":0,"y":0,"w":2,"h":6,"props":{"description":"Soichiro et al, 3DEEP of Skin v2"}},
# {"component":"spatial","coordinationScopes":{"dataset":"A"},"x":2,"y":0,"w":8,"h":12},
# {"component":"layerController","coordinationScopes":{"dataset":"A"},"x":0,"y":6,"w":2,"h":6},
# {"component":"featureList","coordinationScopes":{"dataset":"A"},"x":10,"y":0,"w":2,"h":6},
# {"component":"obsSets","coordinationScopes":{"dataset":"A"},"x":10,"y":6,"w":2,"h":6}],"initStrategy":"auto"}

# Render the Vitessce widget
vc$widget(theme = "light", width = "100%")

vc$export(with_config = TRUE, out_dir = "src/")


