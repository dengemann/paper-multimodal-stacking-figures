#'---
#'title: "2.1 MRI, fMRI and MEG"
#'author: "Denis A. Engemann"
#'date: "8/6/2019"
#'output:
#'    html_document:
#'        code_folding:
#'            hide
#'    md_document:
#'        variant:
#'            markdown_github
#'---

#'In this series of visualizations we will explore how MRI, fMRI and #'MEG interact.
#'We shall answer the following questions:
#'
#'1. How does MEG improve anatomy based age-prediction compared to #'fMRI?
#'
#'no, they are about equal
#'
#'2. Can age prediction be improved by comabining all modalities?
#'
#'yes,  significantly
#'
#'3. How do MEG and fMRI interact, or do they?
#'
#'
#'4. Can the modalities be combined without cost in the absence of #'complete data?

#+ config
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)

# imports `color_cats` and `get_label_from_marker`
source('./utils.r')
source('./config.r')

PREDICTIONS <- './input_data/age_stacked_predictions_megglobal.csv'
PREDICTIONS2 <- './input_data/age_stacked_predictions__na_coded.csv'
PREDICTIONS3 <- './input_data/age_stacked_predictions_meglocal.csv'

SCORES <- './input_data/age_stacked_scores_megglobal.csv'
IMPORTANCES <- "./outputs/age_stacked_importance_8.csv"
data_pred_wide <- read.csv(PREDICTIONS)
data_scores_wide <- read.csv(SCORES)

#'Preprocess data.

#'r data
stacked_keys <- c(
  "MEG_handcrafted",
  "MEG_powers",
  "MEG_powers_cross_powers",
  "MEG_powers_cross_powers_handrafted",
  "MEG_cat_powers_cross_powers_correlation",
  "MEG_cat_powers_cross_powers_correlation_handcrafted",
  "MEG_cross_powers_correlation",
  "MEG_powers_cross_powers_correlation",
  "MEG_all",
  "ALL",
  "ALL_no_fMRI",
  "MRI",
  "ALL_MRI",
  "fMRI"
)

names(data_scores_wide) <- fix_column_names(names(data_scores_wide))
reshape_sel <- names(data_scores_wide)[!names(data_scores_wide) %in% c(
  'repeat_', 'repeat_idx', 'age')]

data_scores <- reshape(data = data_scores_wide,
                       direction = 'long',
                       varying = reshape_sel,
                       v.names = 'MAE',
                       timevar = 'marker',
                       times = reshape_sel)

data_scores['modality'] <- 'MEG'
data_scores[grepl('MEG_', data_scores$marker),]['modality'] <- 'MEG'
data_scores[grepl('Connectivity_Matrix', 
                  data_scores$marker),]['modality'] <- 'fMRI'
data_scores[grepl('ALL', data_scores$marker),]['modality'] <- 'Multimodal'

data_scores$marker <- sub("stacked_", "", data_scores$marker)
data_scores$prediction <- factor(ifelse(
  data_scores$marker %in% stacked_keys, 'stacked', 'linear'))
data_scores$marker <- factor(data_scores$marker)

#'Select stacked data.

#+r stacked_data
data_stacked <- subset(data_scores, prediction == 'stacked')

stacked_selection <- c(
  "ALL",
  "ALL_MRI",
  "ALL_no_fMRI",
  "MRI"
)

data_stacked_sel <- within(
    data_stacked[data_stacked$marker %in% stacked_selection,],
    {
      family <- rep('Multimodal', length(marker))
      family[marker == 'ALL_MRI'] <- 'MRI & fMRI'
      family[marker == 'ALL_no_fMRI'] <- 'MRI & MEG'
      family[marker == 'MRI'] <- 'MRI'
      family <- factor(family)
    }
)

#'Plot error distibution.
#+ fig2a
colors <- setNames(
  with(color_cats, c(black, orange, `blueish green`, blue)),
  c('Multimodal', 'MRI & fMRI', 'MRI & MEG', 'MRI'))

sort_idx <- order(aggregate(MAE ~ marker,
                            data = data_stacked_sel, FUN = mean)$MAE)
data_stacked_sel$cv_idx <- rep(rep(1:10, times = 10), length(colors))

fig2a <- ggplot(
  data = data_stacked_sel,
  mapping = aes(y = MAE, x = reorder(marker, MAE, mean),
                color = family, fill = family)) +
  coord_flip() +
  geom_boxplot(show.legend = T, outlier.shape = NA, alpha = 0.5, size = 0.7) +
  stat_summary(geom = 'text',
               mapping = aes(label  = sprintf("%1.1f", ..y..)),
               fun.y= mean, size = 3.2, show.legend = FALSE,
               position = position_nudge(x=-0.49)) +
  geom_beeswarm(alpha=0.3, show.legend = F, size = 3) +
  guides(
    color = guide_legend(nrow = 3, title.position = "top")) +
  ylab("MAE (years)") +
  xlab("Multimodal stacking") +
  scale_color_manual(values = colors,
                     labels = names(colors),
                     breaks = names(colors), name = "stacking") +
  scale_fill_manual(values = colors,
                    labels = names(colors),
                    breaks = names(colors), name = "stacking") +
  theme(axis.text.y = element_blank(),
        legend.position = 'top',
        legend.justification = 'left',
        legend.text.align = 0) +
  scale_y_continuous(breaks = seq(3, 20, 1)) +
  guides(color = guide_legend(nrow = 1, title.position = 'top'),
         fill = guide_legend(nrow = 1))
print(fig2a)

fname <- "./figures/elements_fig2_stacking_mri_meg."
ggsave(paste0(fname, "pdf"),
       width = save_width, height = save_height, useDingbats = F)
embedFonts(file = paste0(fname, "pdf"), outfile = paste0(fname, "pdf"))
ggsave(paste0(fname, "png"),
       width = save_width, height = save_height, dpi = 300)
knitr::include_graphics(paste0(fname, "png"), dpi = 200)

#'Now let us look at the difference from baseline.

#+ fig2b
mri_mean <- aggregate(MAE ~ family, data_stacked_sel, mean)[1, 2]
# The power of R ... do multi-line assignment expressions inside data frame
# environment and return updated data.frame
data_diff <- within(data_stacked_sel,
    {
      MAE_diff <- c(MAE[family == 'Multimodal'] - MAE[family == 'MRI'],
                    MAE[family == 'MRI & MEG'] - MAE[family == 'MRI'],
                    mri_mean - MAE[family == 'MRI'],
                    MAE[family == 'MRI & fMRI'] - MAE[family == 'MRI'])
      family <- factor(gsub("Multimodal", "Fully multimodal", family))
    }
)

legend_name <- "Improvement over anatomical MRI"
colors_fig2b <- setNames(
  with(color_cats, c(black, `blueish green`, orange, blue)),
  c("MRI, fMRI, MEG", 'MRI, MEG', 'MRI, fMRI', 'MRI'))

sort_idx <- order(aggregate(MAE_diff ~ family, data_diff, mean)$MAE_diff)

levels(data_diff$family) <- c("MRI, fMRI, MEG", "MRI", "MRI, MEG", "MRI, fMRI")

fig2b <- ggplot(
  data = data_diff,
  mapping = aes(y = MAE_diff,
                x = reorder(family, MAE_diff, function(x) mean(x)),
                color = family, fill = family)) +
  coord_flip(ylim = c(-2.9, 1.7)) +
  stat_summary(geom = "boxplot", fun.data = my_quantiles,
               alpha = 0.5, size = 0.7, width = 0.8) +
  stat_summary(geom = "errorbar", fun.data = my_quantiles,
               alpha = 0.5, size = 0.7, width = 0.5) +
  stat_summary(geom = 'text',
               mapping = aes(label  = sprintf("%1.1f",
                                              ..y.. +
                                              mri_mean)),
               fun.y= mean, size = 3.2, show.legend = FALSE,
               position = position_nudge(x=-0.49)) +
  geom_beeswarm(alpha=0.3, show.legend = F, size = 3) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  guides(color = guide_legend(nrow = 3, title.position = "top")) +
  ylab("MAE difference (years)") +
  xlab("Multimodal stacking") +
  scale_color_manual(values = colors_fig2b,
                     labels = names(colors_fig2b),
                     breaks = names(colors_fig2b),
                     name = legend_name) +
  scale_fill_manual(values = colors_fig2b,
                    labels = names(colors_fig2b),
                    breaks = names(colors_fig2b),
                    name = legend_name) +
  guides(color = guide_legend(nrow = 1, title.position = 'top'),
         fill = guide_legend(nrow = 1)) +
  scale_y_continuous(breaks = seq(-3, 1.5, 0.5)) +
  theme(axis.text.y = element_blank(),
        legend.position = 'top',
        legend.justification = 'left',
        legend.text.align = 0)
print(fig2b)

fname <- "./figures/elements_fig2_stacking_mri_meg_diff."
ggsave(paste0(fname,  "pdf"),
       width = save_width, height = save_height, useDingbats = F)
embedFonts(file = paste0(fname, "pdf"), outfile = paste0(fname, "pdf"))
ggsave(paste0(fname,  "png"), width = save_width, height = save_height,
       dpi = 300)
knitr::include_graphics(paste0(fname, "png"), dpi = 200)

#'This is a clear big picture!
#'Let's investigate the interaction between the MEG and MRI.
#'For this we will have to extract the raw prediction errors.
#'
#'We will first prepare a subset of data that allows us to
#'compare MEG and fRMI.

#+ data_wide
stacked_selection2 <- c(stacked_selection, "fMRI", "MEG_all")
data_pred_stacked_sel <- preprocess_prediction_data(
  df_wide = data_pred_wide, stack_sel = stacked_selection2, drop_na = T)

#'Now we can package the data for plotting.

#+ data_diff
data_pred_stacked_comp <- rbind(
  data.frame(
    MAE_meg = subset(data_pred_stacked_sel, family == 'MEG')$MAE,
    MAE_fmri = subset(data_pred_stacked_sel, family == 'fMRI')$MAE,
    combined = F),
   data.frame(
    MAE_meg = subset(data_pred_stacked_sel, family == 'MRI & MEG')$MAE,
    MAE_fmri = subset(data_pred_stacked_sel, family == 'MRI & fMRI')$MAE,
    combined = T)
)

data_pred_stacked_comp <- cbind(
  data_pred_stacked_comp,
  rbind(subset(data_pred_stacked_sel, family == 'MEG', select = -c(MAE)),
        subset(data_pred_stacked_sel, family == 'MEG', select = -c(MAE))))

data_pred_stacked_comp$combined <- factor(
  ifelse(data_pred_stacked_comp$combined, "MRI[anat.]~added", "no~MRI[anat.]"),
  levels = c("no~MRI[anat.]", "MRI[anat.]~added"))

#'Now we can plot it.

#+ fig2c
fig2c <- ggplot(
  data = aggregate(
    cbind(MAE_fmri, MAE_meg, age) ~ X*combined, data = data_pred_stacked_comp, FUN = mean),
  mapping = aes(x = MAE_fmri, y =  MAE_meg, size = age, color = age)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(0.01, 3),
                        trans = 'sqrt') +
  ylab(expression(MAE[MEG] ~ (years))) +
  xlab(expression(MAE[fMRI] ~ (years))) +
  facet_wrap(~combined, labeller = label_parsed) +
  scale_color_viridis_c() +
  theme(legend.position = 'top',
        legend.justification = 'left',
        legend.text.align = 0) +
  guides(
        color = guide_legend(title.position = "left"),
        size = guide_legend(title.position = "left")) +
  coord_fixed(ylim = c(0, 30.5), xlim = c(0, 30.5)) +
  labs(color = "age", shape = "age")
print(fig2c)

fname <- "./figures/elements_fig2_supplement_mri_meg_scatter."
ggsave(paste0(fname, "pdf"), plot = fig2c,
       width = save_width, height = save_height, useDingbats = F)
embedFonts(file = paste0(fname, "pdf"), outfile = paste0(fname, "pdf"))
ggsave(paste0(fname, "png"), plot = fig2c,
        width = save_width, height = save_height,
       dpi = 300)
knitr::include_graphics(paste0(fname, "png"), dpi = 200)

#'One can see that a) the MEG and MRI errors are not
#'strongly related b) comabining each of them with MRI
#'Makes the error somewhat more similar, especially in old
#'people, yet leaves them rather  uncorrelated.
#'

#'Let us now plot the error by age group and modality.

#+ error_by_age
# make qualitative ordered age group
get_age_group <- function(age, step = 5){
  age_seq <- rev(seq(25, 90, step))
  age_group <- rep(age_seq[[1]], length(age))
  for (this_age in age_seq[-1])
    age_group[age < this_age] <- this_age
  return(age_group)
}

colors_fig2e <- setNames(
  with(color_cats, c(black, orange, `blueish green`, blue, violet, vermillon)),
  c("MRI, fMRI, MEG", 'MRI, fMRI', 'MRI, MEG', 'MRI', 'fMRI', 'MEG'))

data_pred_stacked_sel$age_group <- get_age_group(data_pred_stacked_sel$age)
data_pred_stacked_sel$family <- factor(
  data_pred_stacked_sel$family,
  levels = c("MRI", "fMRI", "MEG", "MRI & fMRI", "MRI & MEG", "Multimodal"),
  labels = c("MRI", "fMRI", "MEG", "MRI, fMRI", "MRI, MEG", "MRI, fMRI, MEG"))


data_pred_stacked_sel_agg <- aggregate(
    cbind(MAE, age) ~ X*family, data = data_pred_stacked_sel, FUN = mean)


fig2d <- ggplot(data = data_pred_stacked_sel_agg,
                mapping = aes(x = age, y = MAE,
                              size = MAE, 
                              color = family, fill = family)) +
  stat_smooth(size = 1.2, show.legend = F,
              method = loess, method.args = list(degree = 2),
              fill = NA, level = .9999) +
  geom_hex(mapping =  aes(x = age, y = MAE, alpha = ..density..),
           show.legend = F, size = 0.1, bins = 20) +
  # geom_point(alpha = 0.3, show.legend = F) +
  scale_size_continuous(range = c(0.5, 2.5), trans = 'sqrt') +
  scale_alpha_continuous(range = c(0.05, 1), trans = 'sqrt') +
  scale_color_manual(breaks = names(colors_fig2e),
                     labels = names(colors_fig2e),
                     values = colors_fig2e) +
  scale_fill_manual(breaks = names(colors_fig2e),
                    labels = names(colors_fig2e),
                    values = colors_fig2e) +
  scale_x_continuous(breaks = seq(20, 80, 10),
                     labels = seq(20, 80, 10)) +
  xlab("Age (year)") +
  ylab("MAE (years)") +
  facet_wrap(~family, ncol = 3)
print(fig2d)

fname <- "./figures/elements_fig2_error_by_age_group."
ggsave(paste0(fname, "pdf"), plot = fig2d,
       width = save_width, height = save_height, useDingbats = F)
ggsave(paste0(fname, "png"), plot = fig2d,
        width = save_width, height = save_height,
       dpi = 300)
knitr::include_graphics(paste0(fname, "png"), dpi = 200)

#' We see, what we knew before that, stacking helps improve the error
#' However, we see that this is probably due to at least two mechanisms.
#' 1) extreme error is unifomely reduced. 2) Error in young and old groups
#' is mitigated. However, we also learn that the best model still shows
#' considerable brain age bias, with young and old subjects systematically
#' suffering from more error.

#' Time to investigate the 2D dependence. To select the right variables,
#' we first need to consider the variable importance.

#+ importance_revisited
data_importance <- read.csv(IMPORTANCES)

stack_models <- c('ALL', 'ALL no fMRI', 'ALL MRI')
for (ii in seq_along(stack_models)){
  this_model <- stack_models[[ii]]
  this_result <- get_importances(
    subset(data_importance, stack_model == this_model &
           mod_type == 'permutation'))
  this_result <- this_result[rev(order(this_result$importance)),]
  # Some printing.
  cat(this_model, "\n")
  for (jj in 1:5)
    cat(with(this_result,
              sprintf("#%d - %s: %0.2f", jj, marker[jj], importance[jj])),
        "\n")
}

#' We can see from the printed outputs that MEG power always makes it under
#' the top markers. Thus, it makes sense to focus on power in a 2D dependence
#' analysis.

OUT_DEPENDENCE_2D <- './input_data/age_stacked_dependence_model-full-2d.csv'
data_dependence2d <- read.csv(OUT_DEPENDENCE_2D)
data_dependence2d$marker <- fix_marker_names(data_dependence2d$marker)
data_dependence2d$var_x <- fix_marker_names(data_dependence2d$var_x)
data_dependence2d$var_y <- fix_marker_names(data_dependence2d$var_y)
data_dependence2d$model <- gsub(" ", "_", data_dependence2d$model)

pdp2dmap <- list(
  "ALL" = list(
    cases = c(
      "Connectivity_Matrix,_MODL_256_tan--power_diag",
      "Cortical_Thickness--Connectivity_Matrix,_MODL_256_tan",
      "Cortical_Thickness--power_diag",
      "Cortical_Thickness--Subcortical_Volumes",
      "Subcortical_Volumes--Connectivity_Matrix,_MODL_256_tan",
      "Subcortical_Volumes--power_diag"),
    labels = c(
      "fMRI-P[cat]",
      "CrtT-fMRI",
      "CrtT-P[cat]",
      "CrtT-SbcV",
      "SbcV-fMRI",
      "SbcV-P[cat]")
  ),
  "ALL_no_fMRI" = list(
    cases = c(
      "Cortical_Thickness--mne_power_diag_beta_low",
      "Cortical_Thickness--power_diag",
      "Cortical_Thickness--Subcortical_Volumes",
      "mne_power_diag_beta_low--power_diag",
      "Subcortical_Volumes--mne_power_diag_beta_low",
      "Subcortical_Volumes--power_diag"),
    labels = c(
      "CrtT-P[beta~low]",
      "CrtT-P[cat]",
      "CrtT-SbcV",
      "P[beta~low]-P[cat]",
      "SbcV-P[beta~low]",
      "SbcV-P[cat]")
  ),
  "ALL_MRI" = list(
    cases = c(
      "Cortical_Surface_Area--Connectivity_Matrix,_MODL_256_tan",
      "Cortical_Surface_Area--Cortical_Thickness",
      "Cortical_Surface_Area--Subcortical_Volumes",
      "Cortical_Thickness--Connectivity_Matrix,_MODL_256_tan",
      "Cortical_Thickness--Subcortical_Volumes",
      "Subcortical_Volumes--Connectivity_Matrix,_MODL_256_tan"
    ),
    labels = c(
      "CrSA-fMRI",
      "CrSA-CrtT",
      "CrSA-SbcV",
      "CrtT-fMRI",
      "CrtT-SbcV",
      "SbcV-fMRI"
    )
  )
)

print(nrow(data_dependence2d))
print(unique(data_dependence2d$model))
models <- c("ALL", "ALL_MRI", "ALL_no_fMRI")
for (i_model in seq_along(models)){
  this_model <- models[[i_model]]
  print(this_model)
  this_data <- subset(data_dependence2d,
                      model == this_model & model_type == 'rf_msqrt')

  #  make nicer marker labels for titles.
  this_data$marker_label <- this_data$marker
  for (i_case in seq_along(pdp2dmap[[this_model]][['cases']])) {
    this_data$marker_label <- gsub(
      pdp2dmap[[this_model]][['cases']][i_case],
      pdp2dmap[[this_model]][['labels']][i_case],
      this_data$marker_label)
  }

  marker_split <- strsplit(this_data$marker_label, '-')
  this_data$marker_label_x <- sapply(marker_split, `[[`, 1)
  this_data$marker_label_y <- sapply(marker_split, `[[`, 2)

  if(this_model == "ALL_MRI"){
    this_breaks <- seq(30, 80, 4)
  }else{
    this_breaks <- seq(46, 62, 2)
  }

  fig2e <- ggplot(data = this_data,
                  mapping = aes(x = x, y = y, z = pred)) +
      geom_raster(aes(fill = pred), show.legend = T) +
      stat_contour(breaks = this_breaks,
                  color = "white", bins = 7, show.legend = F) +
      scale_fill_viridis_c(breaks = this_breaks,
                          name = expression(hat(y)),
                          guide = guide_colorbar(barheight = 10)) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      facet_wrap(marker_label_x ~ marker_label_y,
                scales = "free",
                labeller = function(x) label_parsed(x, multi_line = F)) +
      xlab("Input age 1") +
      ylab("Input age 2")
  
  fname <- paste0("./figures/elements_fig2e_meg_pdp_2d",
                  "_", this_model, ".")
  ggsave(paste0(fname, "pdf"), plot = fig2e,
        width = save_width, height = save_height, useDingbats = T)
  embedFonts(file = paste0(fname, "pdf"), outfile = paste0(fname, "pdf"))
  ggsave(paste0(fname, "png"), plot = fig2e,
          width = save_width, height = save_height,
        dpi = 300)
}
# read in main figure.
fname <- paste0("./figures/elements_fig2e_meg_pdp_2d",
                  "_", 'ALL', ".")
knitr::include_graphics(paste0(fname, "png"), dpi = 200)

fname_fig2 <- './figures/fig_2ad.png'
#'After combining all outputs with an illustrator software, we get
#'Figure 2 (including `r fname_fig2` it if exists).
#+ include_fig1_illu
if (file.exists(fname_fig2)) {
  knitr::include_graphics(fname_fig2, dpi = 200)
}

#'## Short caption
#'
#'A) Difference distributions showing improvments over anatomical MRI. 
#' MEG and fMRI are similarly beneficial, when everything is combine performance
#' is best 
#'B) Breakdown of different stacking 
#' models with regard to error across age. Stacking helps avoid extreme errors.
#'C) MEG vs fMRI errors compared. Errors are largely independent, extreme errors
#' occur in different age groups.
#'D) 2D partial dependence plots showing how
#' the stacker combines age inputs. Patterns suggest dominance of additive
#' effects.

#' ## Todos
#' 
#' ### Figure 2A
#'
#' - [ ] Consider carrying over the fMRI and MEG from Figure 2B
#'
#' ### Figure 2C
#'
#' - [x] make output square-like
#' - [x] improve ambgious "with/no fMRI" titles
#' - [ ] perhaps increase size of symbols, consider averaging over repeats  
#' 
#' ### global todo.
#'
#' - [x] make error by age group plot
#' - [x] make 2d partial dependency plot
#' - [x] missing value analysis 1: global common.
#' - [x] missing value analysis 2: local.
#' - [x] missing value analysis 3: multimodal in local index
#' - [ ] missing value analysis 4: the rest, why is it bad? Show that the 
#'       full model will simply perform at the level of remaining variables
#' - [ ] Make really consistent x/y labels, e.g. input, output age  
#' - [ ] decide which outputs become supplements
#' - [ ] compose, insert and caption the figure

#' ## Session info

#+ session_info
print(sessionInfo())