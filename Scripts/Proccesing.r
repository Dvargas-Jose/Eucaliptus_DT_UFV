# ============================================================
# PCA con vegan + ggplot2 + Red de rasgos con MultiTraits
# Datos anatómicos foliares y de pecíolo de Eucalipto - UFV
# ============================================================

library(readxl)
library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(MultiTraits)
library(ggraph)
library(tidygraph)
library(igraph)
library(ggcorrplot)

# Crear carpeta de salida si no existe
if (!dir.exists("Output")) dir.create("Output")

# -----------------------------------------------------------
# 1. Cargar datos
# -----------------------------------------------------------
df <- read_excel("Data/Dados_anat_folh_peciol_euca_LIMPO.xlsx")

# Columnas de metadatos/agrupamiento
meta <- df[, c("Arvore", "Gen", "rep", "Phen", "Sp")]

# Rasgos numéricos: excluir columnas de identificación y grupo
traits_raw <- df %>%
  select(-Arvore, -Gen, -rep, -Phen, -Sp) %>%
  as.data.frame()

# Eliminar columnas con varianza cero
traits_raw <- traits_raw[, apply(traits_raw, 2, function(x) var(x, na.rm = TRUE) > 0)]

cat("Rasgos incluidos en el análisis:", ncol(traits_raw), "\n")
cat("Observaciones:", nrow(traits_raw), "\n\n")

# -----------------------------------------------------------
# 2. Estandarizar con vegan::decostand (z-score)
# -----------------------------------------------------------
traits_std <- decostand(traits_raw, method = "standardize")

# -----------------------------------------------------------
# 3. PCA con vegan::rda (distancias euclidianas implícitas)
# -----------------------------------------------------------
pca <- rda(traits_std)

# Varianza explicada por cada PC
eig     <- eigenvals(pca)
var_exp <- round(eig / sum(eig) * 100, 1)

cat("=== Varianza explicada ===\n")
cat("PC1:", var_exp[1], "%\n")
cat("PC2:", var_exp[2], "%\n")
cat("PC3:", var_exp[3], "%\n\n")

# Distancias euclidianas entre individuos en el espacio de los 2 primeros PCs
sites_2d <- scores(pca, display = "sites", choices = 1:2)
dist_euc <- dist(sites_2d, method = "euclidean")

cat("=== Distancias euclidianas (submatriz 4x4) ===\n")
print(round(as.matrix(dist_euc)[1:4, 1:4], 3))
cat("\n")

# -----------------------------------------------------------
# 4. Extraer scores para ggplot
# -----------------------------------------------------------
site_df <- as.data.frame(scores(pca, display = "sites",   choices = 1:2))
sp_df   <- as.data.frame(scores(pca, display = "species", choices = 1:2))

site_df$Phen   <- meta$Phen
site_df$Sp     <- meta$Sp
site_df$Arvore <- as.character(meta$Arvore)

sp_df$rasgo <- rownames(sp_df)

# Escalar flechas para que queden dentro del espacio de los sitios
escala <- 0.75 * max(abs(site_df[, 1:2])) / max(abs(sp_df[, 1:2]))
sp_df[, 1:2] <- sp_df[, 1:2] * escala

# -----------------------------------------------------------
# 5. ggplot2: Biplot PCA con individuos y variables
# -----------------------------------------------------------
pal_phen <- c(
  "tol"       = "#1565C0",   # azul oscuro
  "medio tol" = "#42A5F5",   # azul claro
  "mod sus"   = "#EF6C00",   # naranja
  "sus"       = "#B71C1C"    # rojo oscuro
)

shapes_sp <- c("EUR" = 16, "GRUR" = 17, "EGR" = 15)

p_pca <- ggplot() +
  # Líneas de referencia
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  # Flechas de variables
  geom_segment(
    data = sp_df,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow     = arrow(length = unit(0.18, "cm"), type = "closed"),
    color     = "grey50",
    linewidth = 0.35,
    alpha     = 0.7
  ) +
  # Etiquetas de variables
  geom_text_repel(
    data        = sp_df,
    aes(x = PC1, y = PC2, label = rasgo),
    size        = 2.2,
    color       = "grey30",
    max.overlaps = 20,
    segment.size = 0.2,
    box.padding  = 0.2
  ) +
  # Puntos de individuos
  geom_point(
    data  = site_df,
    aes(x = PC1, y = PC2, color = Phen, shape = Sp),
    size  = 3.5,
    alpha = 0.9
  ) +
  # Etiquetas de individuos
  geom_text_repel(
    data         = site_df,
    aes(x = PC1, y = PC2, label = Arvore, color = Phen),
    size         = 2.8,
    show.legend  = FALSE,
    max.overlaps = 20,
    segment.size = 0.2
  ) +
  scale_color_manual(values = pal_phen, name = "Fenotipo") +
  scale_shape_manual(values = shapes_sp, name = "Especie") +
  labs(
    title    = "PCA - Rasgos anatómicos y morfológicos de Eucalipto",
    subtitle = paste0("PC1: ", var_exp[1], "%  |  PC2: ", var_exp[2], "% de varianza explicada"),
    x        = paste0("PC1 (", var_exp[1], "%)"),
    y        = paste0("PC2 (", var_exp[2], "%)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey40", size = 10),
    panel.grid.minor = element_blank()
  )

ggsave("Output/PCA_biplot_eucalipto.png", p_pca,
       width = 11, height = 7.5, dpi = 300, bg = "white")
cat("Biplot PCA guardado en: Output/PCA_biplot_eucalipto.png\n\n")

# -----------------------------------------------------------
# 6. Screeplot con ggplot2
# -----------------------------------------------------------
scree_df <- data.frame(
  PC    = paste0("PC", seq_along(var_exp)),
  Var   = as.numeric(var_exp),
  Acum  = cumsum(as.numeric(var_exp))
)[1:15, ]

p_scree <- ggplot(scree_df, aes(x = factor(PC, levels = PC))) +
  geom_col(aes(y = Var), fill = "#1565C0", alpha = 0.8, width = 0.6) +
  geom_line(aes(y = Acum, group = 1), color = "#B71C1C", linewidth = 0.8) +
  geom_point(aes(y = Acum), color = "#B71C1C", size = 2) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "grey50") +
  labs(
    title = "Screeplot - Varianza explicada por PC",
    x     = "Componente Principal",
    y     = "% Varianza"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Output/PCA_screeplot_eucalipto.png", p_scree,
       width = 8, height = 5, dpi = 300, bg = "white")
cat("Screeplot guardado en: Output/PCA_screeplot_eucalipto.png\n\n")

# -----------------------------------------------------------
# 7. Red de rasgos con MultiTraits::PTN
# -----------------------------------------------------------
# Matriz: filas = individuos, columnas = rasgos (estandarizados)
traits_net           <- as.data.frame(traits_std)
rownames(traits_net) <- paste0("Arv", meta$Arvore)

# Construir la red de correlaciones entre rasgos
# rThres: umbral mínimo de |r| para incluir un enlace
# pThres: umbral de p ajustado (FDR)
net_ptn <- PTN(
  traits_matrix    = traits_net,
  method           = "spearman",
  rThres           = 0.6,
  pThres           = 0.05,
  phylo_correction = FALSE
)

cat("=== Métricas de la Red de Rasgos (PTN) ===\n")
metricas <- PTN_metrics(net_ptn)
print(metricas)
cat("\n")

# Copiar el vector de métricas al portapapeles para pegar en Excel o similar
library(clipr)
write_clip(metricas)



# Guardar el plot de la red
png("Output/Red_PTN_eucalipto.png", width = 1600, height = 1300, res = 150)
PTN_plot(
  graph            = net_ptn,
  style            = 1,         
  vertex.size      = 14,
  vertex.label.cex = 0.55
)
dev.off()
cat("Red de rasgos guardada en: Output/Red_PTN_eucalipto.png\n")

# -----------------------------------------------------------
# 8. Visualizaciones alternativas de la red de rasgos
# -----------------------------------------------------------
# Extraer atributos de aristas del objeto igraph
edges_df  <- igraph::as_data_frame(net_ptn, what = "edges")
# El peso de correlación puede llamarse "weight", "r" u otro; detectar automáticamente
cor_col   <- intersect(c("r", "weight", "correlation"), names(edges_df))[1]

tg <- as_tbl_graph(net_ptn)

p_ggraph <- ggraph(tg, layout = "fr") +
  geom_edge_link(
    aes(color = .data[[cor_col]], width = abs(.data[[cor_col]])),
    alpha = 0.85
  ) +
  scale_edge_color_gradient2(
    low      = "#B71C1C",
    mid      = "grey88",
    high     = "#1565C0",
    midpoint = 0,
    name     = "Correlación (r)"
  ) +
  scale_edge_width(range = c(0.4, 2.8), guide = "none") +
  geom_node_point(size = 9, color = "#2E7D32", alpha = 0.85) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.2, color = "grey15") +
  labs(
    title    = "Red de Rasgos (PTN) — ggraph",
    subtitle = "Spearman | |r| ≥ 0.6 | FDR < 0.05",
    caption  = "Azul = correlación positiva  |  Rojo = correlación negativa"
  ) +
  theme_graph(base_family = "sans") +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey40", size = 10),
    legend.position = "right"
  )

ggsave("Output/Red_PTN_ggraph_eucalipto.png", p_ggraph,
       width = 11, height = 9, dpi = 300, bg = "white")
cat("Red ggraph guardada en: Output/Red_PTN_ggraph_eucalipto.png\n")
