data <- read.table('~/Downloads/20241201.HF.Coordinate.Celltype.Structure.txt', sep="\t", header = TRUE)
head(data)
dim(data)

## focus on one follicle for now
hf <- table(data$HF)
head(sort(hf, decreasing=TRUE))
hf.test <- 'HF182'

vi <- data$HF == hf.test
table(vi)
data <- data[vi,]
head(data)

plot(data$x, data$y, pch=".")

########## make cells.json
cells <- unique(data$cell)
final <- lapply(cells, function(cell.test) {
  
  data.cell <- data[data$cell == cell.test,]
  
  library(geometry)  # for convhulln
  #points <- as.matrix(data.cell[, c('x', 'y', 'z')])
  points <- as.matrix(data.cell[, c('x', 'y')])
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
  cell.pos <- colMeans(points)
  
  centroid <- colMeans(hull.vertices)
  angles <- atan2(hull.vertices[, 2] - centroid[2], hull.vertices[, 1] - centroid[1])
  ordered_indices <- order(angles)
  hull.vertices <- hull.vertices[ordered_indices, ]
  
  out <- list(
    genes = as.list(table(data.cell$Gene)),
    poly = as.matrix(hull.vertices),
    factors = list(
      HF = as.character(data.cell$HF[1]),
      Time = as.character(data.cell$Time[1]), 
      CellType = as.character(data.cell$CellType[1]), 
      Structure = as.character(data.cell$Structure[1])
    ),
    xy = as.numeric(cell.pos)
  )
  return(out)
})
names(final) <- cells
library(jsonlite)
exportJSON <- toJSON(final, auto_unbox = TRUE)
write(exportJSON, "~/Desktop/vitessce-3DEEP/data/soichiro_skin_cells.json")


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


cellset <- list("version" = "1.0.2",
                "datatype" = "cell",
                "tree" = list(
                  annot1.list,
                  annot2.list,
                  annot3.list)
                )
exportJSON <- toJSON(cellset, auto_unbox = TRUE)
write(exportJSON, "~/Desktop/vitessce-3DEEP/data/soichiro_skin_cell-sets.json")


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
write(exportJSON, "~/Desktop/vitessce-3DEEP/data/soichiro_skin_clusters.json")


########## make molecules.json
length(genes)
genes.pos <- lapply(genes, function(g) {
  #as.matrix(data[data$Gene == g, c('x', 'y', 'z')])
  as.matrix(data[data$Gene == g, c('x', 'y')])
})
names(genes.pos) <- genes

library(jsonlite)
exportJSON <- toJSON(genes.pos, auto_unbox = TRUE)
write(exportJSON, "~/Desktop/vitessce-3DEEP/data/soichiro_skin_molecules.json")



###################### make app
library(vitessceR)

## need to launch python server
## http-server ./ --cors -p 8000
base_url <- "http://localhost:8000/"

# Create Vitessce view config
vc <- VitessceConfig$new(schema_version = "1.0.2", name = "3DEEP of Skin")
dataset <- vc$add_dataset("3DEEP-Skin")$add_file(
  url = paste0(base_url, "soichiro_skin_cells.json"),
  file_type = FileType$CELLS_JSON
)$add_file(
  url = paste0(base_url, "soichiro_skin_molecules.json"), ## transcript positions
  file_type = FileType$MOLECULES_JSON
)$add_file(
  url = paste0(base_url, "soichiro_skin_clusters.json"), ## needed for gene names?
  file_type = FileType$CLUSTERS_JSON
)$add_file(
  url = paste0(base_url, "soichiro_skin_cell-sets.json"), ## cell-type annotations
  file_type = FileType$CELL_SETS_JSON
)

desc <- vc$add_view(dataset, Component$DESCRIPTION)
desc <- desc$set_props(description = "Soichiro et al, 3DEEP of Skin")
spatial <- vc$add_view(dataset, Component$SPATIAL)
spatial_layers <- vc$add_view(dataset, Component$LAYER_CONTROLLER)
gene_list <- vc$add_view(dataset, Component$FEATURE_LIST)
cell_sets <- vc$add_view(dataset, Component$OBS_SETS)

vc$layout(hconcat(
  vconcat(desc, spatial_layers),
  vconcat(spatial),
  vconcat(gene_list, cell_sets)
))

# Render the Vitessce widget
vc$widget(theme = "dark", width = "100%")

#vc$export(with_config = TRUE, out_dir = "~/Desktop/vitessce-test")


