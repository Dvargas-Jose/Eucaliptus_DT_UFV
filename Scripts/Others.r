# ============================================================
# Modelos estadísticos y ecofisiológicos adicionales
# Datos anatómicos foliares y de pecíolo de Eucalipto - UFV
# ============================================================
# Modelos incluidos:
#   1. PERMANOVA (adonis2)         -- diferencias multivariadas Phen/Sp
#   2. MANOVA + ANOVAs univariados -- medias por grupo
#   3. LDA + validacion cruzada    -- discriminacion de fenotipos
#   4. Random Forest               -- importancia de variables (Phen)
#   5. Cluster jerarquico          -- dendrograma de individuos
#   6. NMDS                        -- ordenacion no-metrica
#   7. Conductividad hidraulica    -- modelo Hagen-Poiseuille (vasos)
#   8. Indices funcionales foliares-- ratios anatomicos ecofisiologicos
#   9. Alometria del peciolo       -- SMA log-log xilema vs area total
#  10. Particion de varianza       -- efecto Sp vs Phen (varpart)
#  11. Coordinacion estoma-venacion-- stomatal economics
# ============================================================

# ── Paquetes ──────────────────────────────────────────────
pkgs <- c("readxl", "vegan", "MASS", "randomForest", "ggplot2",
          "ggrepel", "dplyr", "tidyr", "patchwork", "dendextend",
          "cluster", "smatr", "RColorBrewer")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

if (!dir.exists("Output")) dir.create("Output")

# ── Paletas ────────────────────────────────────────────────
pal_phen <- c(
  "tol"       = "#1565C0",
  "medio tol" = "#42A5F5",
  "mod sus"   = "#EF6C00",
  "sus"       = "#B71C1C"
)
pal_sp <- c("EUR" = "#4CAF50", "GRUR" = "#9C27B0", "EGR" = "#FF9800")

# ============================================================
# 1. CARGAR Y PREPARAR DATOS
# ============================================================
df <- read_excel("Data/Dados_anat_folh_peciol_euca_LIMPO.xlsx")

id_cols  <- c("Arvore", "Gen", "rep", "Phen", "Sp", "SPE", "FDP", "FDE", "FDO")
lum_cols <- grep("^lum_", names(df), value = TRUE)

traits_raw <- df %>%
  select(-all_of(id_cols)) %>%
  as.data.frame()

traits_raw <- traits_raw[, apply(traits_raw, 2, function(x) var(x, na.rm = TRUE) > 0)]

phen_lvls <- c("tol", "medio tol", "mod sus", "sus")
df$Phen   <- factor(df$Phen, levels = phen_lvls)
df$Sp     <- factor(df$Sp)

traits_std <- decostand(traits_raw, method = "standardize")

cat("Datos cargados:", nrow(df), "arboles x", ncol(traits_raw), "rasgos\n")
cat("Phen:", paste(table(df$Phen), collapse = " / "),
    "(", paste(names(table(df$Phen)), collapse = " / "), ")\n")
cat("Sp  :", paste(table(df$Sp),   collapse = " / "),
    "(", paste(names(table(df$Sp)),   collapse = " / "), ")\n\n")

# ============================================================
# 2. PERMANOVA (adonis2)
# ============================================================
cat("=== 2. PERMANOVA (adonis2) ===\n")

set.seed(42)
perm_phen <- adonis2(traits_std ~ Phen, data = df,
                     method = "euclidean", permutations = 9999)
cat("-- Por Fenotipo (Phen) --\n"); print(perm_phen)

set.seed(42)
perm_sp <- adonis2(traits_std ~ Sp, data = df,
                   method = "euclidean", permutations = 9999)
cat("\n-- Por Especie (Sp) --\n"); print(perm_sp)

set.seed(42)
perm_both <- adonis2(traits_std ~ Sp + Phen, data = df,
                     method = "euclidean", permutations = 9999)
cat("\n-- Modelo conjunto Sp + Phen --\n"); print(perm_both)

dist_mat  <- dist(traits_std, method = "euclidean")
beta_phen <- betadisper(dist_mat, df$Phen)
cat("\nHomogeneidad de dispersion (Phen):\n")
print(permutest(beta_phen, permutations = 999))

beta_sp <- betadisper(dist_mat, df$Sp)
cat("\nHomogeneidad de dispersion (Sp):\n")
print(permutest(beta_sp, permutations = 999))
cat("\n")

# ============================================================
# 3. MANOVA + ANOVAs univariados
# ============================================================
cat("=== 3. MANOVA (rasgos anatomicos ~ Phen) ===\n")

top_anat <- c("anat_Esp_Limbo", "anat_Esp_Mesofilo", "anat_Esp_Epiderme_Sup",
              "anat_Esp_Par_Palicadico", "anat_Esp_Par_Lacunoso",
              "pec_pct_Tec_Vascular", "pec_pct_Xilema",
              "NET_mm2", "Media_ET", "DVE")
top_anat <- intersect(top_anat, names(df))

Y_man   <- as.matrix(df[, top_anat])
man_mod <- manova(Y_man ~ df$Phen)
cat("Pillai trace:\n"); print(summary(man_mod, test = "Pillai"))

anova_res <- lapply(top_anat, function(var) {
  s <- summary(aov(df[[var]] ~ df$Phen))[[1]]
  p <- s[["Pr(>F)"]][1]
  data.frame(variable = var,
             F_val = round(s[["F value"]][1], 2),
             p_val = round(p, 4),
             sig   = ifelse(p < 0.05, "*", ""))
})
cat("\nANOVAs univariados:\n")
print(do.call(rbind, anova_res), row.names = FALSE)
cat("\n")

# ============================================================
# 4. LDA -- Analisis Discriminante Lineal (Phen)
# ============================================================
cat("=== 4. LDA ===\n")

lda_vars <- intersect(
  c("anat_Esp_Limbo", "anat_Esp_Mesofilo", "anat_Esp_Par_Palicadico",
    "anat_Esp_Par_Lacunoso", "pec_pct_Xilema", "pec_pct_Tec_Vascular",
    "NET_mm2", "Media_ET"),
  names(df)
)

lda_df  <- df[, c(lda_vars, "Phen")]
lda_mod <- lda(Phen ~ ., data = lda_df)
lda_cv  <- lda(Phen ~ ., data = lda_df, CV = TRUE)

cat("Proporcion varianza por LD:\n")
print(round(lda_mod$svd^2 / sum(lda_mod$svd^2) * 100, 1))

conf_mat <- table(Obs = df$Phen, Pred = lda_cv$class)
cat("\nMatriz de confusion (LOO-CV):\n"); print(conf_mat)
acc_loo  <- sum(diag(conf_mat)) / sum(conf_mat)
cat(sprintf("Precision LOO: %.1f%%\n\n", acc_loo * 100))

lda_scores       <- as.data.frame(predict(lda_mod)$x)
lda_scores$Phen  <- df$Phen
lda_scores$Sp    <- df$Sp
lda_scores$Arv   <- as.character(df$Arvore)

# Centroides por grupo
ld_cols <- intersect(c("LD1", "LD2"), names(lda_scores))
centroids <- lda_scores %>%
  group_by(Phen) %>%
  summarise(across(all_of(ld_cols), mean), .groups = "drop")

p_lda <- ggplot(lda_scores, aes(x = LD1, y = LD2, color = Phen, shape = Sp)) +
  stat_ellipse(aes(group = Phen, fill = Phen), geom = "polygon",
               alpha = 0.08, level = 0.68, show.legend = FALSE) +
  geom_point(size = 4, alpha = 0.9) +
  geom_point(data = centroids, aes(x = LD1, y = LD2, fill = Phen),
             shape = 23, size = 5, color = "black", show.legend = FALSE) +
  geom_text_repel(aes(label = Arv), size = 2.8, show.legend = FALSE,
                  segment.size = 0.2, max.overlaps = 20) +
  scale_color_manual(values = pal_phen, name = "Fenotipo") +
  scale_fill_manual(values  = pal_phen) +
  scale_shape_manual(values = c(EUR = 16, GRUR = 17, EGR = 15), name = "Especie") +
  labs(title    = "Analisis Discriminante Lineal (LDA)",
       subtitle = sprintf("Validacion cruzada LOO: %.1f%% de precision", acc_loo * 100),
       x = "LD1", y = "LD2") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13))

ggsave("Output/LDA_fenotipos_eucalipto.png", p_lda,
       width = 10, height = 7, dpi = 300, bg = "white")
cat("LDA guardado: Output/LDA_fenotipos_eucalipto.png\n\n")

# ============================================================
# 5. RANDOM FOREST -- importancia de variables para Phen
# ============================================================
cat("=== 5. Random Forest ===\n")

rf_df      <- traits_raw
rf_df$Phen <- df$Phen
set.seed(42)
rf_mod <- randomForest(Phen ~ ., data = rf_df,
                       ntree = 1000, importance = TRUE, na.action = na.omit)

cat(sprintf("OOB error: %.1f%%\n", rf_mod$err.rate[1000, "OOB"] * 100))
cat("Matriz de confusion (OOB):\n"); print(rf_mod$confusion); cat("\n")

imp_df <- as.data.frame(importance(rf_mod))
imp_df$Variable <- rownames(imp_df)
imp_df <- imp_df %>% arrange(desc(MeanDecreaseAccuracy)) %>% slice_head(n = 15)

p_rf <- ggplot(imp_df, aes(x = reorder(Variable, MeanDecreaseAccuracy),
                            y = MeanDecreaseAccuracy)) +
  geom_col(fill = "#1565C0", alpha = 0.85, width = 0.7) +
  geom_col(data = imp_df[1:5, ],
           aes(x = reorder(Variable, MeanDecreaseAccuracy),
               y = MeanDecreaseAccuracy),
           fill = "#B71C1C", alpha = 0.85, width = 0.7) +
  coord_flip() +
  labs(title    = "Random Forest -- Importancia de variables",
       subtitle = "Rojo: top-5 | Azul: top 6-15  |  Criterio: MDA",
       x = NULL, y = "Mean Decrease Accuracy") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("Output/RF_importancia_variables_eucalipto.png", p_rf,
       width = 9, height = 6, dpi = 300, bg = "white")
cat("RF guardado: Output/RF_importancia_variables_eucalipto.png\n\n")

# ============================================================
# 6. CLUSTER JERARQUICO + DENDROGRAMA
# ============================================================
cat("=== 6. Cluster Jerarquico ===\n")

dist_hclust <- dist(traits_std, method = "euclidean")
hc           <- hclust(dist_hclust, method = "ward.D2")
hc$labels    <- paste0(df$Arvore, " (", df$Phen, ")")

dend <- as.dendrogram(hc)
leaf_order  <- order.dendrogram(dend)
leaf_colors <- pal_phen[as.character(df$Phen)][leaf_order]

dend <- dend %>%
  set("labels_col", leaf_colors) %>%
  set("labels_cex", 0.65) %>%
  set("branches_lwd", 1.3)

# Silueta para elegir k optimo
sil_w <- sapply(2:6, function(k) {
  km <- kmeans(traits_std, centers = k, nstart = 25, iter.max = 100)
  mean(silhouette(km$cluster, dist_hclust)[, 3])
})
k_opt <- which.max(sil_w) + 1
cat("Ancho de silueta por k:\n")
for (i in seq_along(sil_w)) cat(sprintf("  k=%d: %.3f\n", i + 1, sil_w[i]))
cat(sprintf("k optimo: %d\n\n", k_opt))

png("Output/Cluster_dendrograma_eucalipto.png",
    width = 2400, height = 1600, res = 180)
par(mar = c(9, 4, 4, 2))
plot(dend,
     main = "Cluster Jerarquico (Ward D2) -- Rasgos anatomicos Eucalipto",
     ylab = "Distancia euclidiana (estandarizada)")
legend("topright", legend = phen_lvls,
       fill = pal_phen[phen_lvls], bty = "n", cex = 0.75, title = "Fenotipo")
dev.off()
cat("Dendrograma guardado: Output/Cluster_dendrograma_eucalipto.png\n\n")

# ============================================================
# 7. NMDS -- Ordenacion no-metrica multidimensional
# ============================================================
cat("=== 7. NMDS ===\n")

set.seed(42)
nmds <- metaMDS(traits_std, distance = "euclidean", k = 2,
                trymax = 100, trace = FALSE)
cat(sprintf("Stress NMDS: %.4f  (< 0.10 excelente; < 0.20 aceptable)\n\n",
            nmds$stress))

nmds_df <- as.data.frame(scores(nmds, display = "sites"))
nmds_df$Phen <- df$Phen
nmds_df$Sp   <- df$Sp
nmds_df$Arv  <- as.character(df$Arvore)

p_nmds <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Phen, shape = Sp)) +
  stat_ellipse(aes(group = Phen, fill = Phen), geom = "polygon",
               alpha = 0.08, level = 0.68, show.legend = FALSE) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = Arv), size = 2.8, show.legend = FALSE,
                  segment.size = 0.2, max.overlaps = 20) +
  annotate("text",
           x = max(nmds_df$NMDS1) * 0.8, y = min(nmds_df$NMDS2) * 0.92,
           label = sprintf("Stress = %.3f", nmds$stress),
           size = 3.5, color = "grey40") +
  scale_color_manual(values = pal_phen, name = "Fenotipo") +
  scale_fill_manual(values  = pal_phen) +
  scale_shape_manual(values = c(EUR = 16, GRUR = 17, EGR = 15), name = "Especie") +
  labs(title    = "NMDS -- Ordenacion no-metrica multidimensional",
       subtitle = sprintf("Distancia euclidiana | Stress = %.4f", nmds$stress),
       x = "NMDS1", y = "NMDS2") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave("Output/NMDS_eucalipto.png", p_nmds,
       width = 10, height = 7, dpi = 300, bg = "white")
cat("NMDS guardado: Output/NMDS_eucalipto.png\n\n")

# ============================================================
# 8. MODELO HIDRAULICO DE VASOS XILEMÁTICOS
#    Hagen-Poiseuille: Kh ∝ Σ(n_i × D_i^4)
# ============================================================
cat("=== 8. Analisis Hidraulico de Vasos ===\n")

lum_midpoints <- c(
  lum_lt2   =  1.0, lum_2a4  =  3.0, lum_4a6   =  5.0,
  lum_6a8   =  7.0, lum_8a10 =  9.0, lum_10a12 = 11.0,
  lum_12a14 = 13.0, lum_14a16 = 15.0, lum_16a18 = 17.0,
  lum_18a20 = 19.0, lum_gt20  = 22.0
)
lum_present <- names(lum_midpoints)[names(lum_midpoints) %in% names(df)]

if (length(lum_present) >= 5) {
  lum_matrix <- as.matrix(df[, lum_present])
  D_vec      <- lum_midpoints[lum_present]
  tot        <- rowSums(lum_matrix) + 1e-9

  # Diametro medio ponderado
  D_mean_wt <- rowSums(sweep(lum_matrix, 2, D_vec,   "*")) / tot

  # Diametro hidraulico: D_h = [Σ(n_i*D_i^5) / Σ(n_i*D_i)]^(1/4)
  D_h <- (rowSums(sweep(lum_matrix, 2, D_vec^5, "*")) /
          (rowSums(sweep(lum_matrix, 2, D_vec, "*")) + 1e-9))^(1/4)

  # Conductividad hidraulica relativa (proporcional a Σ n_i * r_i^4)
  Kh_rel <- rowSums(sweep(lum_matrix, 2, (D_vec / 2)^4, "*"))

  # Indice de vulnerabilidad hidraulica: proporcion vasos >= 10 um
  large_cols <- lum_present[lum_midpoints[lum_present] >= 10]
  pct_large  <- rowSums(lum_matrix[, large_cols, drop = FALSE]) / tot * 100

  hyd_df <- data.frame(
    Arvore             = df$Arvore,
    Phen               = df$Phen,
    Sp                 = df$Sp,
    D_media_pond       = round(D_mean_wt, 2),
    D_hidraulico       = round(D_h, 2),
    Kh_relativa        = round(Kh_rel, 4),
    pct_vasos_grandes  = round(pct_large, 1)
  )

  cat("Resumen por fenotipo:\n")
  print(
    hyd_df %>% group_by(Phen) %>%
      summarise(D_media = round(mean(D_media_pond, na.rm = TRUE), 2),
                D_h     = round(mean(D_hidraulico, na.rm = TRUE), 2),
                Kh_rel  = round(mean(Kh_relativa,  na.rm = TRUE), 4),
                pct_lg  = round(mean(pct_vasos_grandes, na.rm = TRUE), 1),
                .groups = "drop")
  )

  kw_Dh <- kruskal.test(D_hidraulico ~ Phen, data = hyd_df)
  kw_Kh <- kruskal.test(Kh_relativa  ~ Phen, data = hyd_df)
  cat(sprintf("\nKruskal-Wallis D_hidraulico ~ Phen: chi2=%.2f, p=%.4f\n",
              kw_Dh$statistic, kw_Dh$p.value))
  cat(sprintf("Kruskal-Wallis Kh_relativa   ~ Phen: chi2=%.2f, p=%.4f\n\n",
              kw_Kh$statistic, kw_Kh$p.value))

  # Box-plot diametro hidraulico
  p_Dh <- ggplot(hyd_df, aes(x = Phen, y = D_hidraulico, fill = Phen)) +
    geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.55) +
    geom_jitter(aes(color = Sp), width = 0.15, size = 2.8, alpha = 0.85) +
    scale_fill_manual(values  = pal_phen, guide = "none") +
    scale_color_manual(values = pal_sp, name = "Especie") +
    labs(title    = "Diametro hidraulico (D_h) por fenotipo",
         subtitle = sprintf("Kruskal-Wallis: p = %.3f", kw_Dh$p.value),
         x = "Fenotipo", y = "D_h (um)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  # Distribucion de lumenes apilada
  lum_long <- hyd_df %>%
    select(Arvore, Phen) %>%
    bind_cols(as.data.frame(lum_matrix)) %>%
    pivot_longer(cols = all_of(lum_present), names_to = "Clase",
                 values_to = "Frec") %>%
    mutate(D_mid = lum_midpoints[Clase]) %>%
    group_by(Arvore) %>%
    mutate(Frec_rel = Frec / (sum(Frec) + 1e-9) * 100) %>%
    ungroup()

  lum_sum <- lum_long %>%
    group_by(Phen, D_mid) %>%
    summarise(Frec_mean = mean(Frec_rel, na.rm = TRUE), .groups = "drop")

  p_lum <- ggplot(lum_sum, aes(x = D_mid, y = Frec_mean, fill = Phen)) +
    geom_col(position = "dodge", alpha = 0.85, width = 1.6) +
    scale_fill_manual(values = pal_phen, name = "Fenotipo") +
    scale_x_continuous(breaks = unname(lum_midpoints[lum_present])) +
    labs(title    = "Distribucion de diametros de lumenes xilemáticos",
         subtitle = "Frecuencia relativa media (%) por fenotipo",
         x = "Diametro medio de clase (um)", y = "Frecuencia relativa (%)") +
    theme_bw(base_size = 12) +
    theme(plot.title  = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))

  p_hyd <- (p_Dh / p_lum) +
    plot_annotation(
      title = "Analisis hidraulico de vasos xilemáticos -- Eucalipto",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )

  ggsave("Output/Hidraulica_vasos_eucalipto.png", p_hyd,
         width = 11, height = 11, dpi = 300, bg = "white")
  cat("Analisis hidraulico guardado: Output/Hidraulica_vasos_eucalipto.png\n\n")
} else {
  cat("No se encontraron columnas de lumenes suficientes. Omitiendo seccion 8.\n\n")
}

# ============================================================
# 9. INDICES FUNCIONALES FOLIARES Y DE PECIOLO
# ============================================================
cat("=== 9. Indices Funcionales Foliares ===\n")

func_df <- df %>%
  transmute(
    Arvore = Arvore, Phen = Phen, Sp = Sp,
    # Relacion parenquima palisadico / lacunoso (mas alto = mas eficiencia fotosintetica)
    ratio_pal_lac   = anat_Esp_Par_Palicadico / (anat_Esp_Par_Lacunoso + 1e-6),
    # Mesofilo como fraccion del limbo total
    pct_mesofilo    = anat_Esp_Mesofilo / (anat_Esp_Limbo + 1e-6) * 100,
    # Epidermis superior como fraccion del limbo (costo de proteccion)
    pct_epid_sup    = anat_Esp_Epiderme_Sup / (anat_Esp_Limbo + 1e-6) * 100,
    # Fraccion fotosintetica neta del limbo
    idx_fotosint    = (anat_Esp_Par_Palicadico + anat_Esp_Par_Lacunoso) /
                      (anat_Esp_Limbo + 1e-6) * 100,
    # Inversion en xilema del peciolo
    ratio_xil_vasc  = pec_pct_Xilema / (pec_pct_Tec_Vascular + 1e-6),
    # Relacion tejido vascular / fundamental (eficiencia transporte)
    ratio_vasc_fund = pec_pct_Tec_Vascular / (pec_pct_Tec_Fundamental + 1e-6),
    NET_mm2         = NET_mm2,
    Media_ET        = Media_ET
  )

cat("Indices por fenotipo:\n")
print(func_df %>% group_by(Phen) %>%
        summarise(across(ratio_pal_lac:Media_ET,
                         ~ round(mean(.x, na.rm = TRUE), 3)),
                  .groups = "drop"))
cat("\n")

idx_vars <- c("ratio_pal_lac", "pct_mesofilo", "idx_fotosint",
              "ratio_vasc_fund", "NET_mm2")
kw_func <- lapply(idx_vars, function(v) {
  kt <- kruskal.test(func_df[[v]] ~ func_df$Phen)
  data.frame(Indice = v, chi2 = round(kt$statistic, 2),
             p_val = round(kt$p.value, 4),
             sig   = ifelse(kt$p.value < 0.05, "*", ""))
})
cat("Kruskal-Wallis indices funcionales:\n")
print(do.call(rbind, kw_func), row.names = FALSE); cat("\n")

labels_idx <- c(
  ratio_pal_lac   = "Palisadico/Lacunoso",
  pct_mesofilo    = "% Mesofilo/Limbo",
  idx_fotosint    = "% Tejido fotosintetico",
  ratio_vasc_fund = "Vasc./Fund. (peciolo)",
  NET_mm2         = "Densidad venacion (mm-2)"
)

func_long <- func_df %>%
  select(Arvore, Phen, Sp, all_of(idx_vars)) %>%
  pivot_longer(cols = all_of(idx_vars),
               names_to = "Indice", values_to = "Valor") %>%
  mutate(Indice = labels_idx[Indice])

p_func <- ggplot(func_long, aes(x = Phen, y = Valor, fill = Phen)) +
  geom_boxplot(alpha = 0.72, outlier.shape = NA, width = 0.55) +
  geom_jitter(aes(color = Sp), width = 0.15, size = 2.2, alpha = 0.8) +
  facet_wrap(~ Indice, scales = "free_y", ncol = 2) +
  scale_fill_manual(values  = pal_phen, guide = "none") +
  scale_color_manual(values = pal_sp,   name  = "Especie") +
  labs(title = "Indices funcionales foliares y de peciolo por fenotipo",
       x = "Fenotipo de tolerancia", y = "Valor del indice") +
  theme_bw(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        axis.text.x   = element_text(angle = 35, hjust = 1),
        strip.text    = element_text(face  = "bold"),
        legend.position = "right")

ggsave("Output/Indices_funcionales_eucalipto.png", p_func,
       width = 12, height = 10, dpi = 300, bg = "white")
cat("Indices funcionales guardados: Output/Indices_funcionales_eucalipto.png\n\n")

# ============================================================
# 10. ALOMETRIA DEL PECIOLO -- SMA log-log
#     xilema ~ area total del peciolo (alometria hidraulica)
# ============================================================
cat("=== 10. Alometria del Peciolo (SMA) ===\n")

allom_df <- df %>%
  select(Arvore, Phen, Sp,
         pec_Area_Total, pec_Xilema, pec_Tec_Vascular, pec_Floema_Cambio) %>%
  filter(pec_Area_Total > 0, pec_Xilema > 0)

sma_xil  <- sma(log10(pec_Xilema)       ~ log10(pec_Area_Total), data = allom_df)
sma_vasc <- sma(log10(pec_Tec_Vascular) ~ log10(pec_Area_Total), data = allom_df)

slope_xil  <- sma_xil$groupsummary$Slope
slope_vasc <- sma_vasc$groupsummary$Slope
r2_xil     <- sma_xil$groupsummary$r2
r2_vasc    <- sma_vasc$groupsummary$r2

cat(sprintf("SMA Xilema ~ Area_Total:        pendiente = %.3f, R2 = %.3f\n", slope_xil,  r2_xil))
cat(sprintf("SMA Tec.Vasc ~ Area_Total:      pendiente = %.3f, R2 = %.3f\n", slope_vasc, r2_vasc))
cat("(pendiente = 1 -> isometria; > 1 -> alometria positiva)\n\n")

# SMA por especie
sma_xil_sp <- sma(log10(pec_Xilema) ~ log10(pec_Area_Total) * Sp, data = allom_df)
cat("SMA por especie (Xilema ~ Area):\n")
print(sma_xil_sp$groupsummary[, c("group", "Slope", "Int", "r2", "pval")])
cat("\n")

p_allom <- ggplot(allom_df,
                  aes(x = log10(pec_Area_Total), y = log10(pec_Xilema),
                      color = Phen, shape = Sp)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey45",
              linetype = "dashed", linewidth = 0.8, alpha = 0.12) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_text_repel(aes(label = Arvore), size = 2.8, show.legend = FALSE,
                  segment.size = 0.2, max.overlaps = 20) +
  annotate("text",
           x = min(log10(allom_df$pec_Area_Total)) + 0.04,
           y = max(log10(allom_df$pec_Xilema))     - 0.04,
           label = sprintf("SMA: b = %.2f | R2 = %.2f", slope_xil, r2_xil),
           hjust = 0, size = 3.5, color = "grey25") +
  scale_color_manual(values = pal_phen, name = "Fenotipo") +
  scale_shape_manual(values = c(EUR = 16, GRUR = 17, EGR = 15), name = "Especie") +
  labs(title    = "Alometria del peciolo -- SMA log10",
       subtitle = "Area xilema ~ Area total (escala logaritmica)",
       x = "log10(Area total peciolo, um2)",
       y = "log10(Area xilema, um2)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")

ggsave("Output/Alometria_peciolo_eucalipto.png", p_allom,
       width = 10, height = 7, dpi = 300, bg = "white")
cat("Alometria guardada: Output/Alometria_peciolo_eucalipto.png\n\n")

# ============================================================
# 11. PARTICION DE VARIANZA (varpart) -- Sp vs Phen
# ============================================================
cat("=== 11. Particion de Varianza (varpart) ===\n")

X_sp   <- model.matrix(~ Sp   - 1, data = df)
X_phen <- model.matrix(~ Phen - 1, data = df)

vp <- varpart(traits_std, X_sp, X_phen)
fracs <- vp$part$fract$R.square
cat(sprintf("Fraccion [a] Sp sola    : %.3f\n", fracs[1]))
cat(sprintf("Fraccion [b] compartida : %.3f\n", fracs[2]))
cat(sprintf("Fraccion [c] Phen sola  : %.3f\n", fracs[3]))
cat(sprintf("Fraccion residual       : %.3f\n\n", 1 - fracs[4]))

set.seed(42)
rda_sp   <- rda(traits_std ~ Sp   + Condition(Phen), data = df)
rda_phen <- rda(traits_std ~ Phen + Condition(Sp),   data = df)
pval_sp   <- anova(rda_sp,   permutations = 9999)$`Pr(>F)`[1]
pval_phen <- anova(rda_phen, permutations = 9999)$`Pr(>F)`[1]
cat(sprintf("Fraccion pura Sp   (RDA parcial): p = %.4f\n", pval_sp))
cat(sprintf("Fraccion pura Phen (RDA parcial): p = %.4f\n\n", pval_phen))

png("Output/Particion_varianza_eucalipto.png",
    width = 1200, height = 900, res = 150)
plot(vp,
     Xnames = c("Especie\n(Sp)", "Fenotipo\n(Phen)"),
     bg     = c("#9C27B0", "#1565C0"),
     alpha  = 0.55,
     digits = 3)
title(main = "Particion de varianza -- Efecto Especie vs Fenotipo",
      cex.main = 1.1, font.main = 2)
dev.off()
cat("Particion de varianza guardada: Output/Particion_varianza_eucalipto.png\n\n")

# ============================================================
# 12. COORDINACION ESTOMA-VENACION (stomatal economics)
# ============================================================
cat("=== 12. Coordinacion estomática y vascular ===\n")

if (all(c("NET_mm2", "Media_ET") %in% names(df))) {
  stom_df <- df %>%
    select(Arvore, Phen, Sp, NET_mm2, Media_ET) %>%
    filter(!is.na(NET_mm2), !is.na(Media_ET), NET_mm2 > 0, Media_ET > 0)

  sma_stom <- sma(log10(Media_ET) ~ log10(NET_mm2), data = stom_df)
  sl_st    <- sma_stom$groupsummary$Slope
  r2_st    <- sma_stom$groupsummary$r2
  cat(sprintf("SMA: log(Media_ET) ~ log(NET_mm2): b = %.3f, R2 = %.3f\n\n",
              sl_st, r2_st))

  # Correlacion de Spearman
  sp_cor <- cor.test(stom_df$NET_mm2, stom_df$Media_ET, method = "spearman")
  cat(sprintf("Spearman: rho = %.3f, p = %.4f\n\n",
              sp_cor$estimate, sp_cor$p.value))

  p_stom <- ggplot(stom_df,
                   aes(x = log10(NET_mm2), y = log10(Media_ET),
                       color = Phen, shape = Sp)) +
    geom_smooth(method = "lm", se = TRUE, color = "grey55",
                linetype = "dashed", linewidth = 0.8, alpha = 0.12) +
    geom_point(size = 4, alpha = 0.9) +
    geom_text_repel(aes(label = Arvore), size = 2.8, show.legend = FALSE,
                    segment.size = 0.2, max.overlaps = 20) +
    annotate("text",
             x = min(log10(stom_df$NET_mm2)) + 0.02,
             y = max(log10(stom_df$Media_ET)) - 0.02,
             label = sprintf("SMA: b = %.2f | R2 = %.2f", sl_st, r2_st),
             hjust = 0, size = 3.5, color = "grey25") +
    scale_color_manual(values = pal_phen, name = "Fenotipo") +
    scale_shape_manual(values = c(EUR = 16, GRUR = 17, EGR = 15), name = "Especie") +
    labs(title    = "Coordinacion venacion-estomas (stomatal economics)",
         subtitle = sprintf("Spearman: rho = %.2f | p = %.3f", sp_cor$estimate, sp_cor$p.value),
         x = "log10(Densidad de venacion, mm-2)",
         y = "log10(Tamano estoma / Media_ET)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")

  ggsave("Output/Coordinacion_estoma_venacion_eucalipto.png", p_stom,
         width = 10, height = 7, dpi = 300, bg = "white")
  cat("Grafico guardado: Output/Coordinacion_estoma_venacion_eucalipto.png\n\n")
}

# ============================================================
# RESUMEN FINAL
# ============================================================
cat("=======================================================\n")
cat("ARCHIVOS GENERADOS EN Output/\n")
cat("=======================================================\n")
outputs <- c(
  "LDA_fenotipos_eucalipto.png",
  "RF_importancia_variables_eucalipto.png",
  "Cluster_dendrograma_eucalipto.png",
  "NMDS_eucalipto.png",
  "Hidraulica_vasos_eucalipto.png",
  "Indices_funcionales_eucalipto.png",
  "Alometria_peciolo_eucalipto.png",
  "Particion_varianza_eucalipto.png",
  "Coordinacion_estoma_venacion_eucalipto.png"
)
for (f in outputs) {
  fp <- file.path("Output", f)
  cat(sprintf("  %-55s [%s]\n", f, ifelse(file.exists(fp), "OK", "no generado")))
}
cat("\n")
