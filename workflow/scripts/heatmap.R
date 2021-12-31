# load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))

#load data
dat <- read_csv(snakemake@input[[1]])

# create matrix for plotting
mat <- dat %>%
  select(-c(assembly)) %>%
  as.matrix()

# set matrix row names
rownames(mat) <- dat$assembly

# colors for the plot
col_fun = colorRamp2(c(0, 1), c("lightblue", "red"))

# function to determine figure size
calc_hm_size = function(heatmap) {
    pdf(NULL)
    hm = draw(heatmap)
    w = ComplexHeatmap:::width(hm)
    h = ComplexHeatmap:::height(hm)
    dev.off()

    c(w, h)
}

# create heatmap
hm <- mat %>%
  Heatmap(., col = col_fun(c(0,1)),
          column_title = "Defence system family",
          column_title_side = "bottom",
          cluster_rows = FALSE, cluster_columns = TRUE,
          show_row_names = TRUE, show_column_dend = FALSE,
          row_title = "Bacterial samples",
	  row_title_rot = 90, 
          row_names_side = "left",
          row_names_gp = gpar(cex = 0.75),
          border = TRUE,
          heatmap_legend_param = list(title = "Defence\nsystem\nfamily\npresent",
                                      labels = c("Yes", "No")),
	  width = ncol(mat)*unit(4, "mm"), 
          height = nrow(mat)*unit(4, "mm"))

size <- calc_hm_size(hm)

png(snakemake@output[[1]], width = size[1], height = size[2], res = 600, units = "mm")
draw(hm)
dev.off()

